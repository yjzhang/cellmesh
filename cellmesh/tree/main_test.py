from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import pdb

from cellmesh.tree.imesh import \
  iMeSH_D, \
  iMeSH_D_Tree

from cellmesh.tree.path import path_iMeSH_D_Tree

def test_iMeSH_D_Tree():
  iMT = iMeSH_D_Tree(path_iMeSH_D_Tree["json_iMeSH_D"])

  # show whole tree
  iMT.show_all_tree()

  # show subtree for two mids
  # modify "1", "2" etc to any tag strings you want
  #mid_list = [("D011901", "1"), ("D019599", "2")]

  #http://uncurl-app.yjzhang.com:8888/user/test_10x_400_new/view
  # issues:
  # - query genes: may be missing (synonym?) or not occurred in literature?
  # - ref cell: naming strategy CD8+ cytotoxic t vs CD8+/CD45RA+ naive cytotoxic
  # cluster 1
  mid_list = [("D009000", "ref"), ("D002999", "1"), ("D016176", "2"), ("D007694", "3")]
  # cluster 5
  #mid_list = [("D013602", "ref"), ("D004912", "1"), ("D020298", "2"), ("D007694", "3"), ("D003713", "4"), ("D015603", "5")]
  # cluster 6 CD4+/CD45RO+ memory t
  #mid_list = [("D016176", "ref"), ("D004912", "1"), ("D015603", "2"), ("D019169", "3")]
  # cluster 7 CD8+/CD45RA+ naive cytotoxic
  #mid_list = [("D013602", "ref"), ("D020298", "1"), ("D007694", "2"), ("D003713", "3")]
  st = iMT.tree_str(mid_list)
  print("\n\nsubtree for mid_list=%s"%str(mid_list))
  print(st)
  return

def main():
  test_iMeSH_D_Tree()
  return

if __name__ == '__main__':
  main()
