import unittest
import cellmesh

class CellMeshTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_cellmesh_query_basic(self):
        genes = ['HLA-DRB5',
                 'HLA-DQA1',
                 'HCG24',
                 'ARHGEF1',
                 'HLA-DPB1',
                 'SCN4A',
                 'TEAD2',
                 'HLA-DRB1',
                 'LST1',
                 'HLA-DQB1',
                 'AIM1L']
        results = cellmesh.hypergeometric_test(genes)
        print(results[0])
        self.assertEqual(results[0][1], 'Dendritic Cells')
        genes = ['MYH11',
                 'CNN1',
                 'ACTA2',
                 'MYL9',
                 'OPTC']
        results = cellmesh.hypergeometric_test(genes)
        print(results[0])
        self.assertEqual(results[0][1], 'Myocytes, Smooth Muscle')

    def test_cellmesh_query_species(self):
        print('mouse_results')
        genes = ['HLA-DRB5',
                 'HLA-DQA1',
                 'HCG24',
                 'ARHGEF1',
                 'HLA-DPB1',
                 'SCN4A',
                 'TEAD2',
                 'HLA-DRB1',
                 'LST1',
                 'HLA-DQB1',
                 'AIM1L']
        results = cellmesh.hypergeometric_test(genes, species='mouse')
        print(results[0])
        genes = ['MYH11',
                 'CNN1',
                 'ACTA2',
                 'MYL9',
                 'OPTC']
        results = cellmesh.hypergeometric_test(genes, species='mouse')
        print(results[0])
        self.assertEqual(results[0][1], 'Myocytes, Smooth Muscle')
        results = cellmesh.hypergeometric_test(genes, species='both')
        print(results[0])
        self.assertEqual(results[0][1], 'Myocytes, Smooth Muscle')

    def test_cellmesh_prob_method(self):
        print('test_cellmesh_prob_method')
        from cellmesh import prob_method
        genes = ['MYH11',
                 'CNN1',
                 'ACTA2']
        results = prob_method.prob_test(genes, species='mouse')
        print(results[:3])
        self.assertEqual(results[0][1], 'Myocytes, Smooth Muscle')
        results = prob_method.prob_test(genes, species='both')
        print(results[:3])
        self.assertEqual(results[0][1], 'Myocytes, Smooth Muscle')

    def test_anatomy_query(self):
        genes = ['Tnnt2',
                 'Ifit1',
                 'Ifit3',
                 'Ifit3b',
                 'Rsad2',
                 'Psd',
                 'Tbx20',
                 'Gata4',
                 'Iigp1',
                 'Lum',
                 'Cd34',
                 'Corin',
                 'Pgam2',
                 'Dcn']
        results = cellmesh.hypergeometric_test(genes, db_dir=cellmesh.ANATOMY_DB_DIR)
        print(results[0])

    def test_cellmesh_prob_method_anatomy(self):
        print('test_cellmesh_prob_method_anatomy')
        from cellmesh import prob_method
        genes = ['MYH11',
                 'CNN1',
                 'ACTA2']
        results = prob_method.prob_test(genes, species='human', db_dir=cellmesh.ANATOMY_DB_DIR)
        print(results[:3])
        results = prob_method.prob_test(genes, species='mouse', db_dir=cellmesh.ANATOMY_DB_DIR)
        print(results[:3])
        self.assertEqual(results[0][1], 'Muscle, Smooth')
        results = prob_method.prob_test(genes, species='both', db_dir=cellmesh.ANATOMY_DB_DIR)
        print(results[:3])
        self.assertEqual(results[0][1], 'Muscle, Smooth')


    def test_get_descendants(self):
        # brain
        cell_id_names = cellmesh.get_all_cell_id_names(db_dir=cellmesh.ANATOMY_DB_DIR)
        cell_id_names = {c[0]: c[1] for c in cell_id_names}
        descendants = cellmesh.get_descendants(('D001921',))
        self.assertTrue(len(descendants) > 1)
        print('Brain descendants')
        #for d in descendants:
        #    print(d, cell_id_names[d])
        descendants = cellmesh.get_descendants(('D001921', 'D008168'))
        self.assertTrue(len(descendants) > 2)
        print('Brain + Lung descendants')
        #for d in descendants:
        #    print(d, cell_id_names[d])

    def test_go_query(self):
        from cellmesh import go_query
        genes = ['Tnnt2',
                 'Ifit1',
                 'Ifit3',
                 'Ifit3b',
                 'Rsad2',
                 'Psd',
                 'Tbx20',
                 'Gata4',
                 'Iigp1',
                 'Lum',
                 'Cd34',
                 'Corin',
                 'Pgam2',
                 'Dcn']
        results = go_query.gene_set_query(genes, return_header=False)
        self.assertTrue(len(results) > 0)
        print('mouse results')
        print(results[0])
        results = go_query.gene_set_query(genes, return_header=False, species='human')
        self.assertTrue(len(results) > 0)
        print('human results')
        print(results[0])

    def test_tree(self):
        tree = cellmesh.get_cellmesh_anatomy_tree()
        self.assertTrue(len(tree) > 0)
        #print(tree)


if __name__ == '__main__':
    unittest.main()
