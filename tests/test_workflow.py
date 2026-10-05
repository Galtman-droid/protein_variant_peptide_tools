import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import variant_peptide_tools as v

REF = 'ACDEFGHIKLMNPQRSTVWY'
ROOT = Path(__file__).resolve().parents[1]


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def file(self, name, text):
        path = self.work / name
        path.write_text(text, encoding='utf-8')
        return path

    def test_single_fasta_and_whitespace(self):
        self.assertEqual(v.read_sequence(self.file('a.fa', '>id\nacde f\ngh\n')), 'ACDEFGH')

    def test_multifasta_rejected(self):
        with self.assertRaises(v.InputError):
            v.read_sequence(self.file('a.fa', '>a\nACD\n>b\nEFG'))

    def test_unsupported_sequence_rejected(self):
        for seq in ['', 'ACDX', 'ACD*', 'AC1']:
            with self.subTest(seq=seq), self.assertRaises(v.InputError):
                v.read_sequence(self.file('a.fa', seq))

    def test_mutation_semicolon_combination(self):
        mutations = v.parse_mutations('KA9A,CA2A;', REF)
        self.assertEqual(v.apply_mutations(REF, mutations), 'AADEFGHIALMNPQRSTVWY')

    def test_invalid_tokens(self):
        for token in ['CA2A,', ',', 'CA2A;;', 'XA2A', 'CA2X', 'CA2A;KA9A']:
            with self.subTest(token=token), self.assertRaises(v.InputError):
                v.parse_mutations(token, REF)

    def test_chain_case_and_other_chain(self):
        for token in ['CB2A', 'Ca2A']:
            with self.subTest(token=token), self.assertRaises(v.InputError):
                v.parse_mutations(token, REF)
        self.assertEqual(v.parse_mutations('Ca2A', REF, 'a')[0].position, 2)

    def test_number_offset(self):
        mutation = v.parse_mutations('CB102A', REF, 'B', 100)[0]
        self.assertEqual(mutation.position, 2)

    def test_wrong_wt_out_of_range_noop_duplicates(self):
        for text in ['AA2G', 'AA0G', 'AA21G', 'CA2C', 'CA2A,CA2G', 'CA2A,CA2A']:
            with self.subTest(text=text), self.assertRaises(v.InputError):
                v.parse_mutations(text, REF)

    def test_all_standard_substitutions_roundtrip(self):
        for index, wt in enumerate(REF, 1):
            for target in sorted(v.AA - {wt}):
                mutations = v.parse_mutations(f'{wt}A{index}{target}', REF)
                expected = REF[:index-1] + target + REF[index:]
                self.assertEqual(v.apply_mutations(REF, mutations), expected)

    def test_no_fixed_six_mutation_cap(self):
        text = ','.join(f'{aa}A{i}W' for i, aa in enumerate(REF[:7], 1))
        path = self.file('variants.txt', text)
        self.assertEqual(v.read_variants(path, REF, 'A', 0)[0]['n_mutations'], 7)
        with self.assertRaises(v.InputError):
            v.read_variants(path, REF, 'A', 0, 6)

    def test_wt_and_comments(self):
        variants = v.read_variants(self.file('variants.txt', '# comment\nWT;\nCA2A; # note'), REF, 'A', 0)
        self.assertEqual(variants[0]['mutant_sequence'], REF)
        self.assertEqual(variants[1]['variant_id'], 'VAR002')

    def test_overlapping_repeat_rejected(self):
        with self.assertRaises(v.InputError):
            v.find_start('AAA', 'AAAA')

    def test_repeat_explicit_start(self):
        path = self.file('windows.csv', 'core_id,window,window_len,native_window,start\nG,1,3,AAA,2\n')
        self.assertEqual(v.read_windows(path, 'AAAA')[0]['start'], 2)

    def test_bad_windows_fail_entire_file(self):
        prefix = 'core_id,window,window_len,native_window,start\nG,1,3,ACD,1\n'
        for bad in ['G,2,3,ACD,2', 'G,1,3,ACD,1', 'G,2,4,ACD,1', 'G,0,3,ACD,1', ',2,3,ACD,1']:
            with self.subTest(row=bad), self.assertRaises(v.InputError):
                v.read_windows(self.file('windows.csv', prefix + bad), REF)

    def test_inconsistent_annotations(self):
        path = self.file('windows.csv', 'core_id,window,window_len,native_window,core_seq\nG,1,3,ACD,ACD\nG,2,3,EFG,EFG\n')
        with self.assertRaises(v.InputError):
            v.read_windows(path, REF)

    def test_pool_selection_and_global_deduplication(self):
        variants = v.read_variants(self.file('variants.txt', 'CA2A\nCA2A,KA9A\nWT'), REF, 'A', 0)
        windows = v.read_windows(self.file('windows.csv', 'core_id,window,window_len,native_window\nG,1,5,ACDEF\nG,2,5,FGHIK\nH,1,5,ACDEF\n'), REF)
        rows, summary, catalog = v.generate_pools(variants, windows)
        self.assertEqual(len(rows), 6)
        self.assertEqual(sum(r['status'] == 'UNCHANGED' for r in rows), 1)
        self.assertEqual({p['sequence'] for p in catalog}, {'AADEF', 'FGHIA'})
        peptide = next(p for p in catalog if p['sequence'] == 'AADEF')
        self.assertEqual(peptide['n_variants'], 2)
        self.assertEqual(peptide['window_sources'], 'G:1;H:1')
        self.assertEqual(summary[-1]['status'], 'NO_AFFECTED_WINDOWS')
        changed_rows, _, _ = v.generate_pools(variants, windows, 'changed-windows')
        self.assertEqual(len(changed_rows), 5)

    def test_tampered_builder_csv_rejected(self):
        rows = [{'variant_id':'V1', 'mutation_string':'CA2A', 'n_mutations':1, 'validation':'OK', 'mutant_sequence':REF}]
        path = self.work / 'mutants.csv'
        v.write_csv(path, v.VARIANT_FIELDS, rows)
        with self.assertRaises(v.InputError):
            v.read_built_variants(path, REF, 'A', 0)

    def test_invalid_and_duplicate_builder_rows_rejected(self):
        base = {'variant_id':'V1', 'mutation_string':'WT', 'n_mutations':0, 'validation':'OK', 'mutant_sequence':REF}
        for rows in [[base, base], [{**base, 'validation':'FAIL'}], [{**base, 'n_mutations':1}]]:
            path = self.work / 'mutants.csv'
            v.write_csv(path, v.VARIANT_FIELDS, rows)
            with self.assertRaises(v.InputError):
                v.read_built_variants(path, REF, 'A', 0)

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(ROOT / 'variant_peptide_tools.py'), *map(str, args)],
                              cwd=self.work, capture_output=True, text=True, timeout=30)

    def test_cli_csv_and_xlsx_handoff(self):
        mutants = self.work / 'mutants'
        result = self.run_cli('mutants', '--sequence', ROOT / 'examples/reference.fasta', '--variants', ROOT / 'examples/variants.txt', '--output', mutants)
        self.assertEqual(result.returncode, 0, result.stderr)
        for ext in ['csv', 'xlsx']:
            result = self.run_cli('pools', '--sequence', ROOT / 'examples/reference.fasta', '--mutants', mutants / f'mutant_sequences.{ext}', '--windows', ROOT / 'examples/windows.csv', '--output', self.work / ext)
            self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / 'csv/peptide_pools.csv').read_bytes(), (self.work / 'xlsx/peptide_pools.csv').read_bytes())
        manifest = json.loads((mutants / 'manifest.json').read_text())
        self.assertEqual(manifest['variant_count'], 4)
        self.assertTrue(all(len(i['sha256']) == 64 for i in manifest['inputs']))

    def test_output_overwrite_refused(self):
        result = self.run_cli('mutants', '--sequence', ROOT / 'examples/reference.fasta', '--variants', ROOT / 'examples/variants.txt', '--output', self.work)
        self.assertEqual(result.returncode, 2)
        self.assertIn('already exists', result.stderr)

    def test_invalid_input_does_not_create_output(self):
        result = self.run_cli('mutants', '--sequence', ROOT / 'examples/reference.fasta', '--variants', self.file('bad.txt', 'CB2A'), '--output', self.work / 'out')
        self.assertEqual(result.returncode, 2)
        self.assertFalse((self.work / 'out').exists())

    def test_empty_catalog_is_valid(self):
        variants = v.read_variants(self.file('variants.txt', 'WT\nYA20F'), REF, 'A', 0)
        windows = v.read_windows(ROOT / 'examples/windows.csv', REF)
        rows, summaries, catalog = v.generate_pools(variants, windows)
        self.assertEqual(rows, [])
        self.assertEqual(catalog, [])
        self.assertEqual(len(summaries), 2)

    def test_excel_labels_stored_as_text(self):
        from openpyxl import Workbook
        wb = Workbook()
        ws = v.table_sheet(wb, 'Data', ['label'], [{'label':'=1+1'}])
        self.assertEqual(ws['A2'].data_type, 's')


if __name__ == '__main__':
    unittest.main()
