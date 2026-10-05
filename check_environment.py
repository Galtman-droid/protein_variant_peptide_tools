"""Run a small synthetic end-to-end check without user data or external tools."""
import csv
import importlib.metadata
import platform
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    print(f'Python {platform.python_version()} on {platform.system()}')
    if sys.version_info < (3, 10):
        print('FAIL: Python 3.10 or newer is required')
        return 1
    try:
        version = importlib.metadata.version('openpyxl')
        print(f'openpyxl {version}')
        if version != '3.1.5':
            raise RuntimeError('Install the tested dependencies: python -m pip install -r requirements.txt')
        from openpyxl import load_workbook
        root = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix='variant_peptide_check_') as temp:
            work = Path(temp)
            ref, variants, windows = (work / name for name in ['reference.fasta', 'variants.txt', 'windows.csv'])
            ref.write_text('>artificial_alphabet\nACDEFGHIKLMNPQRSTVWY\n', encoding='utf-8')
            variants.write_text('CA2A;\nCA2A,KA9A;\nYA20F;\nWT;\n', encoding='utf-8')
            windows.write_text('core_id,window,window_len,native_window,start\nG1,1,5,ACDEF,1\nG1,2,5,FGHIK,5\n', encoding='utf-8')
            commands = [
                [sys.executable, str(root / 'mutant_builder.py'), '--sequence', str(ref), '--variants', str(variants), '--output', str(work / 'mutants')],
                [sys.executable, str(root / 'peptide_pool_generator.py'), '--sequence', str(ref), '--mutants', str(work / 'mutants/mutant_sequences.xlsx'), '--windows', str(windows), '--output', str(work / 'pools')],
            ]
            for command in commands:
                result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=30)
                if result.returncode:
                    raise RuntimeError(result.stderr or result.stdout)
            with (work / 'pools/synthesis_catalog.csv').open(newline='', encoding='utf-8') as handle:
                catalog = list(csv.DictReader(handle))
            if {p['sequence'] for p in catalog} != {'AADEF', 'FGHIA'}:
                raise RuntimeError('Unexpected catalog sequence results')
            with (work / 'pools/pool_summary.csv').open(newline='', encoding='utf-8') as handle:
                summaries = list(csv.DictReader(handle))
            if len(summaries) != 4 or sum(s['status'] == 'NO_AFFECTED_WINDOWS' for s in summaries) != 2:
                raise RuntimeError('Unexpected pool summary')
            wb = load_workbook(work / 'mutants/mutant_sequences.xlsx', rich_text=True)
            try:
                rich = wb['Mutant sequences']['E5'].value
                if str(rich) != 'AADEFGHIKLMNPQRSTVWY':
                    raise RuntimeError('Excel mutant sequence was not preserved')
                red = ''.join(block.text for block in rich if getattr(getattr(block, 'font', None), 'color', None) and block.font.color.rgb == 'FFFF0000')
                if red != 'A':
                    raise RuntimeError('Excel mutation highlighting failed')
            finally:
                wb.close()
        print('PASS: sequence validation, both script entry points, Excel handoff, red residue highlighting, pool mapping, and sequence deduplication')
        print('Only temporary artificial data were used; temporary files removed. No EvoEF2 run was required.')
        return 0
    except Exception as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
