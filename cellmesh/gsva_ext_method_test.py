import unittest

from cellmesh.gsva_ext_method import calc_cell_gene_ranks

#python -m unittest test_module.TestClass.test_method

class TestGsvaExt(unittest.TestCase):
  def test_calc_cell_gene_ranks(self):
    cell_gene_count = [("a",2.0), ("b",1.0)]
    all_genes = set(["a", "b", "c", "d"])
    cell_gene_rank_sorted = calc_cell_gene_ranks(
      cell_gene_count, all_genes)
    print(cell_gene_rank_sorted)
    return

if __name__ == "__main__":
  unittest.main()