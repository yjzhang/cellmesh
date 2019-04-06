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
  mid_list = [("D011901", "1"), ("D019599", "2")]
  st = iMT.tree_str(mid_list)
  print("\n\nsubtree for mid_list=%s"%str(mid_list))
  print(st)
  return

def main():
  test_iMeSH_D_Tree()
  return

if __name__ == '__main__':
  main()
