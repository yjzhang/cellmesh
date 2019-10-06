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

    def test_anatomy_query(self):
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
        results = cellmesh.hypergeometric_test(genes, db_dir=cellmesh.ANATOMY_DB_DIR)
        for r in results:
            print(r)

    def test_get_descendants(self):
        # brain
        cell_id_names = cellmesh.get_all_cell_id_names(db_dir=cellmesh.ANATOMY_DB_DIR)
        cell_id_names = {c[0]: c[1] for c in cell_id_names}
        descendants = cellmesh.get_descendants(('D001921',))
        print('Brain descendants')
        for d in descendants:
            print(d, cell_id_names[d])
        descendants = cellmesh.get_descendants(('D001921', 'D008168'))
        print('Brain + Lung descendants')
        for d in descendants:
            print(d, cell_id_names[d])




if __name__ == '__main__':
    unittest.main()
