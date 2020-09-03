

import unittest
import lisa

class AssertionTests(unittest.TestCase):

    def test_species_instantiation(self):
        self.assertRaises(AssertionError, lisa.LISA, 'hg39')

    def test_cores(self):
        self.assertRaises( AssertionError, lisa.LISA, 'hg38', -2)
        self.assertRaises( AssertionError, lisa.LISA, 'hg38', 'all')

    def test_isd_method(self):
        self.assertRaises(AssertionError, lisa.LISA, 'hg38',-1,'motifsorwhatever')

    def test_datasets(self):
        args = (lisa.LISA, 'hg38',-1,'chipseq')
        self.assertRaises(AssertionError, *args, 10, 200)
        self.assertRaises(AssertionError, *args, 200, -10)
        self.assertRaises(AssertionError, *args, 0, 0)
        self.assertRaises(AssertionError, *args, 1000, 200)


if __name__ == "__main___":
    unittest.main()