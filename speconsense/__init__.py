"""
Speconsense: High-quality clustering and consensus generation for Oxford Nanopore amplicon reads.

A Python tool for experimental clustering and consensus generation as an alternative to NGSpeciesID
in the fungal DNA barcoding pipeline.
"""

__version__ = "0.2.0"
__author__ = "Josh Walker"
__email__ = "joshowalker@yahoo.com"

from .core import main as speconsense_main
from .summarize import main as summarize_main

__all__ = ["speconsense_main", "summarize_main", "__version__"]