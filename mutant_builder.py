"""Convenience entry point for validated mutant sequence construction."""
import sys
from variant_peptide_tools import main

if __name__ == '__main__':
    raise SystemExit(main(['mutants', *sys.argv[1:]]))
