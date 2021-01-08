import unittest
from lisa.core import data_interface
import numpy as np
import tempfile
import os
from scipy import sparse

class TestProjections(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = data_interface.DataInterface('mm10', window_size=1000, path = './test.h5', make_new=True)

    def test_array_projection(self):

        proj = np.array([
            (2,0),
            (3,1),
            (5,3),
        ])

        arr = np.array([0,0,1,2,0,3])
        correct = correct = np.array([1,2,0,3])
    
        self.assertTrue(
            np.all(self.data.project_array(arr, proj, 4) == correct)
        )

    def test_indices_projection(self):

        ind = [0,4,3]

        proj = np.array([
            (0,1),
            (1,2),
            (3,4),
            (4,5)
        ])

        correct = [1,5,4]

        ans = self.data.project_indices(ind, proj)

        self.assertTrue(
            len(set(ans)) == len(set(correct))
        )

        self.assertTrue(
            len(set(correct).intersection(set(ans))) == 3
        )

    def test_matrix_projection(self):

        arr = sparse.csr_matrix(np.array([
            (0,1),
            (1,0),
            (2,3),
            (3,4),
            (6,8),
            ])
        )

        proj = np.array([
            (0,1),
            (1,0),
            (4,2),
            (2,3),
        ])

        ans = np.array(self.data.project_sparse_matrix(arr, proj, 5).todense())

        correct = np.array([
            (1,0),
            (0,1),
            (6,8),
            (2,3),
            (0,0),
        ])

        self.assertTrue(
            np.all(correct == ans)
        )

    @classmethod
    def tearDownClass(cls):
        os.remove('./test.h5')


class TestDataConsistency(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = data_interface.DataInterface('mm10', window_size=100000, path = './test.h5', make_new=True)

    @staticmethod
    def sparse_matrix_equal(m1, m2):

        m1, m2 = m1.tocsr(), m2.tocsr()

        return np.all(m1.indices == m2.indices) and np.all(m1.indptr == m2.indptr) and \
            np.all(np.isclose(m1.data, m2.data))

    def test_rp_map_consistency(self):

        rp_map_shape = self.data.get_rp_map_shape()

        #random_rp_map = sparse.random(*rp_map_shape, density = 0.01)
        num_entries = 1000
        i = np.random.choice(rp_map_shape[0], num_entries, replace = False)
        j = np.random.choice(rp_map_shape[1], num_entries, replace = False)
        data = np.random.randn(num_entries)

        random_rp_map = sparse.csr_matrix((data, (i,j)), shape = rp_map_shape)

        self.data.add_rp_map('test', random_rp_map)

        self.assertTrue(
            tuple(self.data.get_rp_maps()) == ('test',)
        )

        reloaded = self.data.get_rp_map('test')

        self.assertTrue(
            self.sparse_matrix_equal(random_rp_map, reloaded)
        )

    def test_factor_binding_consistency(self):

        num_bins = len(self.data.genome)

        simulated_hits = np.random.choice(num_bins, 1000, replace = False)

        self.data.add_binding_data('ChIP-seq', '1', simulated_hits, factor = 'STAT', tissue = 'intestine',
            cell_line = '1893A', cell_type = 'blah', qc  = 1)

        reloaded, meta = self.data.get_binding_dataset('ChIP-seq', '1')

        self.assertTrue(
            len(set(simulated_hits).symmetric_difference(set(reloaded))) == 0
        )

        reloaded_matrix, ids, meta = self.data.get_binding_data('ChIP-seq')

        sim_hits_matrix = sparse.csc_matrix((np.ones_like(simulated_hits), simulated_hits, [0, len(simulated_hits)]),
                shape = (len(self.data.genome), 1))

        self.assertTrue(
            self.sparse_matrix_equal(sim_hits_matrix, reloaded_matrix)
        )

        self.assertTrue(
            tuple(ids) == ('1',)
        )

    def test_profile_consistency(self):

        simulated_profile = np.abs(np.random.randn(len(self.data.genome)))[:, np.newaxis]

        num_genes, num_bins = self.data.get_rp_map_shape()

        diag_rp_map = sparse.diags(
            np.ones(num_genes), shape = (num_genes, num_bins)
        )

        self.data.add_profile_data('DNase', '1', simulated_profile, [diag_rp_map], ['diag'], tissue = 'intestine',
            cell_line = '1893A', cell_type = 'blah', qc  = 1, norm_depth=None)

        reloaded, meta = self.data.get_profile('DNase', '1')

        self.assertTrue(
            np.all(np.isclose(simulated_profile, reloaded, rtol=1e-2))
        )

        rp_matrix, dataset_ids = self.data.get_rp_matrix('DNase','diag')

        self.assertTrue(
            np.all(np.isclose(rp_matrix[:num_bins], simulated_profile[:num_bins], rtol=1e-2))
        )


    @classmethod
    def tearDownClass(cls):
        os.remove('./test.h5')


if __name__ == '__main__':
    unittest.main()