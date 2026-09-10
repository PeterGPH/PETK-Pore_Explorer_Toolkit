"""
Main entry point for running the SEM package as a module.

This allows the package to be run with:
    python -m sem config.json run
    python -m sem config.json preview_only
"""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())