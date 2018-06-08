import unittest
from samtools_lookup import *

class samtools_lookup_test(unittest.TestCase):

	def test_get_seq(self):
		self.assertEqual(str(get_seq("I", 1, 101, name = None, zero = False, species = "elegans").seq), "GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTA")

	def test_get_seq_zero_index(self):
		self.assertEqual(str(get_seq("I", 0, 101, name = None, zero = True, species = "elegans").seq), "GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTA")

if __name__ == '__main__':
    unittest.main()