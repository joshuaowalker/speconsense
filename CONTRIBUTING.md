# Contributing to Speconsense

Thank you for your interest in contributing to Speconsense! This document provides guidelines for development, testing, and contributing to the project.

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- External dependencies:
  - [SPOA (SIMD POA)](https://github.com/rvaser/spoa) - Required (install via conda)
  - [MCL](https://micans.org/mcl/) - Optional but recommended (install via conda)

### Installation for Development

```bash
# Create and activate a virtual environment (recommended)
python -m venv speconsense-dev
source speconsense-dev/bin/activate  # On Windows: speconsense-dev\Scripts\activate

# Clone the repository
git clone https://github.com/joshuaowalker/speconsense.git
cd speconsense

# Install in editable mode
pip install -e .
```

This installs the package in editable mode, allowing you to modify the code and immediately test changes without reinstalling.

## Running Tests

The project includes pytest-based integration tests.

### Setting Up the Test Environment

```bash
# Set up test virtual environment
python -m venv test-venv
source test-venv/bin/activate  # On Windows: test-venv\Scripts\activate

# Install project and test dependencies
pip install -e .
pip install pytest pytest-cov
```

### Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_augment_input.py -v

# Run with coverage
python -m pytest tests/ --cov=speconsense --cov-report=html
```

### Available Tests

- `tests/test_augment_input.py`: Integration tests for --augment-input functionality
- `tests/test_orientation.py`: Tests for sequence orientation feature

## Code Style

The project follows standard Python conventions:

- Type annotations throughout (Union, Optional, List, Dict, NamedTuple)
- Snake_case for functions/variables, PascalCase for classes
- Comprehensive docstrings with parameter descriptions
- Context managers for file operations
- Max 100 character line length
- Explicit error handling with appropriate logging levels

## Making Contributions

### Before Starting

1. Check existing issues and pull requests to avoid duplicate work
2. For major changes, open an issue first to discuss your approach
3. Ensure you have the development environment set up correctly

### Pull Request Process

1. Create a new branch for your feature or bugfix
2. Write or update tests for your changes
3. Ensure all tests pass
4. Update documentation (README.md, docstrings) as needed
5. Follow the existing code style
6. Write clear commit messages describing your changes
7. Submit a pull request with a clear description of the changes

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in present tense (e.g., "Add feature", "Fix bug")
- Include context about why the change was made

## Testing Guidelines

- Write integration tests for new features
- Test edge cases and error conditions
- Ensure tests are reproducible
- Use meaningful test data that represents real-world scenarios
- Document what each test is validating

## Documentation

When adding features or making changes:

- Update relevant sections in README.md
- Add docstrings to new functions/classes
- Update help text for command-line arguments
- Include usage examples for new features
- Update CHANGELOG.md with your changes

## Questions?

If you have questions about contributing, please open an issue for discussion.

## License

By contributing to Speconsense, you agree that your contributions will be licensed under the BSD 3-Clause License.
