"""Convenience entry point for peptide window and pool generation."""
import sys
from variant_peptide_tools import main

if __name__ == '__main__':
    raise SystemExit(main(['pools', *sys.argv[1:]]))
