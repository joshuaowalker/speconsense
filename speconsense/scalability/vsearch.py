"""Vsearch-based candidate finding for scalable sequence comparison."""

import hashlib
import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Optional

from tqdm import tqdm


class VsearchCandidateFinder:
    """Vsearch-based candidate finding using usearch_global.

    This implementation uses vsearch to quickly find approximate sequence matches,
    which can then be refined with exact scoring. It is designed for large-scale
    datasets where O(n^2) pairwise comparisons become infeasible.

    The implementation uses SHA256-based deduplication to reduce the database size
    when many identical sequences are present.
    """

    def __init__(self,
                 batch_size: int = 1000,
                 num_threads: int = 1):
        """Initialize VsearchCandidateFinder.

        Args:
            batch_size: Number of sequences to query per batch
            num_threads: Number of threads for vsearch (default: 1 for backward compatibility)
        """
        self.batch_size = batch_size
        self.num_threads = num_threads
        self._db_path: Optional[str] = None
        self._hash_to_ids: Dict[str, List[str]] = {}
        self._cache_dir: Optional[str] = None

    @property
    def name(self) -> str:
        """Human-readable name of this backend."""
        return "vsearch"

    @property
    def is_available(self) -> bool:
        """Check if vsearch is installed and accessible."""
        try:
            result = subprocess.run(
                ['vsearch', '--version'],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def build_index(self,
                    sequences: Dict[str, str],
                    output_dir: str,
                    cache_id: Optional[str] = None) -> None:
        """Build vsearch database with SHA256-based deduplication.

        Args:
            sequences: Dict mapping sequence_id -> sequence_string
            output_dir: Directory for cache files
            cache_id: Unique identifier for this cache (e.g., sample name).
                      If not provided, uses process ID to avoid collisions.
        """
        # Use cache_id or PID to ensure parallel instances don't collide
        unique_id = cache_id if cache_id else str(os.getpid())
        self._cache_dir = os.path.join(output_dir, f".vsearch_cache_{unique_id}")
        os.makedirs(self._cache_dir, exist_ok=True)

        self._db_path = os.path.join(self._cache_dir, "sequences.fasta")

        # Deduplicate sequences using hash
        unique_seqs: Dict[str, tuple] = {}  # hash -> (list of ids, sequence)
        for seq_id, seq in sorted(sequences.items()):
            seq_hash = hashlib.sha256(seq.encode()).hexdigest()[:16]
            if seq_hash not in unique_seqs:
                unique_seqs[seq_hash] = ([], seq)
            unique_seqs[seq_hash][0].append(seq_id)

        # Write deduplicated FASTA
        with open(self._db_path, 'w') as f:
            for seq_hash, (ids, seq) in unique_seqs.items():
                f.write(f">{seq_hash}\n{seq}\n")

        # Store mapping for result lookup
        self._hash_to_ids = {
            seq_hash: ids for seq_hash, (ids, _) in unique_seqs.items()
        }

        logging.debug(f"Built vsearch index: {len(unique_seqs)} unique sequences "
                      f"(deduplicated from {len(sequences)} total)")

    def find_candidates(self,
                        query_ids: List[str],
                        sequences: Dict[str, str],
                        min_identity: float,
                        max_candidates: int) -> Dict[str, List[str]]:
        """Find candidate matches using vsearch usearch_global.

        Args:
            query_ids: List of sequence IDs to query
            sequences: Dict mapping sequence_id -> sequence_string
            min_identity: Minimum identity threshold (0.0-1.0)
            max_candidates: Maximum candidates to return per query

        Returns:
            Dict mapping query_id -> list of candidate target_ids
        """
        if not self._db_path or not os.path.exists(self._db_path):
            raise RuntimeError("Index not built. Call build_index() first.")

        all_results: Dict[str, List[str]] = defaultdict(list)

        # Process in batches with progress bar
        with tqdm(total=len(query_ids), desc="Finding candidates with vsearch") as pbar:
            for i in range(0, len(query_ids), self.batch_size):
                batch_ids = query_ids[i:i + self.batch_size]
                batch_results = self._run_batch(batch_ids, sequences, min_identity, max_candidates)

                for query_id, candidates in batch_results.items():
                    all_results[query_id].extend(candidates)

                pbar.update(len(batch_ids))

        # Validate results - detect likely vsearch failures
        total_candidates = sum(len(c) for c in all_results.values())
        seqs_with_candidates = sum(1 for c in all_results.values() if c)

        logging.debug(f"vsearch found {total_candidates} candidates for {len(query_ids)} sequences "
                      f"({seqs_with_candidates} sequences with ≥1 candidate)")

        # If zero candidates for a large dataset, vsearch likely failed
        if len(query_ids) > 100 and total_candidates == 0:
            raise RuntimeError(
                f"vsearch returned zero candidates for {len(query_ids)} sequences. "
                "This may indicate vsearch was killed due to resource contention. "
                "Try running with --threads 1 when using GNU parallel."
            )

        return dict(all_results)

    def _run_batch(self,
                   query_ids: List[str],
                   sequences: Dict[str, str],
                   min_identity: float,
                   max_candidates: int) -> Dict[str, List[str]]:
        """Run vsearch on a single batch of queries.

        Args:
            query_ids: List of sequence IDs to query
            sequences: Dict mapping sequence_id -> sequence_string
            min_identity: Minimum identity threshold
            max_candidates: Maximum candidates per query

        Returns:
            Dict mapping query_id -> list of candidate target_ids
        """
        # Create temporary query file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            for seq_id in query_ids:
                f.write(f">{seq_id}\n{sequences[seq_id]}\n")
            query_path = f.name

        try:
            cmd = [
                'vsearch',
                '--usearch_global', query_path,
                '--db', self._db_path,
                '--userout', '/dev/stdout',
                '--userfields', 'query+target+id',
                '--id', str(min_identity),
                '--maxaccepts', str(max_candidates),
                '--threads', str(self.num_threads),
                '--output_no_hits'
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse results
            results: Dict[str, List[str]] = defaultdict(list)
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 3:
                    continue

                query_id, target_hash, identity = parts

                # Map hash back to original IDs
                if target_hash in self._hash_to_ids:
                    for original_id in self._hash_to_ids[target_hash]:
                        if original_id != query_id:  # Skip self-matches
                            results[query_id].append(original_id)

            return dict(results)

        except FileNotFoundError:
            raise RuntimeError(
                "vsearch command not found. Please install vsearch:\n"
                "  conda install bioconda::vsearch\n"
                "or visit https://github.com/torognes/vsearch for installation instructions."
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"vsearch failed with return code {e.returncode}")
            logging.error(f"vsearch stderr: {e.stderr}")
            raise

        finally:
            # Clean up temporary query file
            if os.path.exists(query_path):
                os.unlink(query_path)

    def cleanup(self) -> None:
        """Clean up cache directory and temporary files."""
        if self._cache_dir and os.path.exists(self._cache_dir):
            import shutil
            shutil.rmtree(self._cache_dir)
        self._db_path = None
        self._hash_to_ids = {}
        self._cache_dir = None
