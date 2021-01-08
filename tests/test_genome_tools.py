import unittest
from lisa.core import genome_tools
import numpy as np

class TestRegionObject(unittest.TestCase):

    def test_invalid_region(self):
        with self.assertRaises(AssertionError):
            genome_tools.Region('chr1',10, 5)

    def test_overlap_different_chr(self):
        self.assertFalse(
            genome_tools.Region('chr1',5,10).overlaps(genome_tools.Region('chr2',5,10))
        )

    def test_overlap_negative(self):
        self.assertFalse(
            genome_tools.Region('chr1',5,10).overlaps(genome_tools.Region('chr1',20,25))
        )

    def test_overlap_negative_abutted(self):
        self.assertFalse(
            genome_tools.Region('chr1',5,10).overlaps(genome_tools.Region('chr1',10,20))
        )

    def test_overlap_any(self):
        self.assertTrue(
            genome_tools.Region('chr1',5,10).overlaps(genome_tools.Region('chr1',8,15))
        )

        self.assertFalse(
            genome_tools.Region('chr1',1,11).overlaps(genome_tools.Region('chr1',9,15), min_overlap_proportion=0.4)
        )

        self.assertTrue(
            genome_tools.Region('chr1',1,11).overlaps(genome_tools.Region('chr1',6,15), min_overlap_proportion=0.4)
        )

    def test_distance_function(self):
        self.assertEqual(
            genome_tools.Region('chr1',20, 30).get_genomic_distance(genome_tools.Region('chr1', 120, 130)), 100)

    def test_regions_equal(self):
        self.assertTrue(
            genome_tools.Region('chr1',10,15) == genome_tools.Region('chr1', 10, 15)
        )

class TestGenomeObject(unittest.TestCase):

    def setUp(self):
        self.genome = genome_tools.Genome(['chr1','chr2','chr3'], [450, 375, 600], window_size=100)
        
        self.unsorted_genome = genome_tools.Genome(['chr1','chr3','chr2'],[200, 300, 250], window_size= 100, _sort=False)
        self.sorted_genome = genome_tools.Genome(['chr1','chr3','chr2'],[200, 400, 250], window_size= 100, _sort=True)
        
        self.correct_mapping = np.array([
            (0,0),
            (1,1),
            (2,5), 
            (3,6),
            (4,7),
            (5,2),
            (6,3),
            (7,4),
        ])

    def test_genome_genome_mapping(self):

        self.assertTrue(
            np.all(self.unsorted_genome.map_genomes(self.sorted_genome) == self.correct_mapping)
        )

    def test_indptr(self):
        self.assertEqual(
            tuple(self.genome.indptr), (0,5,9,15)
        )
    
    def test_region_check(self):
        with self.assertRaises(genome_tools.BadRegionError):
            self.genome.check_region(genome_tools.Region('chr4',10,20))

        with self.assertRaises(genome_tools.BadRegionError):
            self.genome.check_region(genome_tools.Region('chr1',-1, 20))

        with self.assertRaises(genome_tools.BadRegionError):
            self.genome.check_region(genome_tools.Region('chr2', 300, 400))

    def test_chromlen(self):
        self.assertEqual(
            self.genome.get_chromlen('chr2'), 375
        )

    def test_get_num_windows(self):
        self.assertEqual(
            self.genome.get_num_windows(self.genome.get_chromlen('chr1'), self.genome.window_size), 5
        )

        self.assertEqual(
            self.genome.num_windows_in_genome(), 15
        )

    def test_get_window_from_position(self):
        region, window_idx = self.genome.get_window_from_position('chr2', 250)
        self.assertTrue(
            region == genome_tools.Region('chr2', 200, 300)
        )

        next_region, next_window_idx = self.genome.get_next_window(region)
        self.assertTrue(
            next_region == genome_tools.Region('chr2',300,375)
        )

        self.assertTrue(
            window_idx == (next_window_idx - 1)
        )

    def test_get_window_idx(self):
        self.assertTrue(
            self.genome.get_region(0)[0] == genome_tools.Region('chr1', 0, 100)
        )

        self.assertTrue(
            self.genome.get_region(4)[0] == genome_tools.Region('chr1',400,450)
        )

        self.assertTrue(
            self.genome.get_region(5)[0] == genome_tools.Region('chr2',0,100)
        )

        self.assertTrue(
            self.genome.get_region(14)[0] == genome_tools.Region('chr3',500,600)
        )

        with self.assertRaises(AssertionError):
            self.genome.get_region(15)


class TestRegionSet(unittest.TestCase):

    def setUp(self):
        self.genome = genome_tools.Genome(['chr1','chr2','chr3'], [450, 375, 600], window_size=50)
        self.scrambled_genome = genome_tools.Genome(['chr1','chr3','chr2'], [450, 600, 375], window_size=50)

        self.regions1A = [
            genome_tools.Region('chr1',20,40),
            genome_tools.Region('chr1',30,60),
            genome_tools.Region('chr1',210,230),
            genome_tools.Region('chr2',100,150),
            genome_tools.Region('chr2',220,233),
            genome_tools.Region('chr3',430,450),
        ]
        self.regions1B = [
            genome_tools.Region('chr1',20,40),
            genome_tools.Region('chr1',30,60),
            genome_tools.Region('chr1',210,230),
            genome_tools.Region('chr2',100,150),
            genome_tools.Region('chr2',220,233),
            genome_tools.Region('chr3',430,450),
        ]

        self.auto_distancing_truth = np.array(
            [[0,15,0,0,0,0],
             [15,0,175,0,0,0],
             [0,175,0,0,0,0],
             [0,0,0,0,101,0],
             [0,0,0,101,0,0],
             [0,0,0,0,0,0]
            ]
        )

        self.m2m_map_truth = np.array([
            (0,0),
            (0,1),
            (1,1),
            (4,2),
            (11,3),
            (13,4),
            (25,5)
        ])

    def test_auto_distancing(self):
        distance_matrix = genome_tools.RegionSet(self.regions1A, self.genome).map_intersects(genome_tools.RegionSet(self.regions1B, self.genome), 
            lambda x,y : x.get_genomic_distance(y), slop_distance=75)

        self.assertTrue(
            np.all(np.array(distance_matrix.todense()).astype(int) == self.auto_distancing_truth)
        )

    def test_auto_distancing_scrambled(self):

        distance_matrix = genome_tools.RegionSet(self.regions1A, self.scrambled_genome)\
            .map_intersects(genome_tools.RegionSet(self.regions1B, self.scrambled_genome), 
            lambda x,y : x.get_genomic_distance(y), slop_distance=75)

        self.assertTrue(
            np.all(np.array(distance_matrix.todense()).astype(int) == self.auto_distancing_truth)
        )

    def test_genome_bin_mapping(self):

        m2m_map = genome_tools.RegionSet(self.regions1A, self.genome)\
            .map_genomic_windows(min_window_overlap_proportion=0.0, regions_to_bins=False)

        self.assertTrue(
            np.all(m2m_map == self.m2m_map_truth)
        )

if __name__ == '__main__':
    unittest.main()
