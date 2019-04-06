from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import pdb
import json
import re
import numpy as np

from cellmesh.tree.util import Timer, load_json
from cellmesh.tree.tree import TreeNode, Tree

class iMeSH_D():
  '''
  [i]ndexed [MeSH_D] object
  '''
  def __init__(self):
    '''
    members:
      self.dic = {
        DescriptorUI:{
          DescriptorName:str,
          TreeNumberList:[TreeNumber],
          ScopeNote:str
        }
      }    
    '''
    self.dic = {}
    return

  def get_mesh_ids(self):
    return self.dic.keys()

  def get_mesh_name(self, mesh_id):
    '''
    return mesh term name wrt mesh_id
      in case mesh_id not found, return empty string
    '''
    info = self.dic.get(mesh_id, {})
    if 'DescriptorName' in info:
      return info['DescriptorName']
    else:
      return ''

  def get_mesh_itm(self, mesh_id):
    return self.dic[mesh_id]

  def __str__(self, n=None):
    '''
    Input:
      n: number records to show
    Output:
      a json-formatted string with n or all(n=None) records
    '''
    if n is None:
      st = json.dumps(self.dic, indent=2)
    else:
      dic2 = {}
      cnt = 0
      for k, v in list(self.dic.items()):
        if cnt==n: break
        cnt += 1
        dic2[k] = v
      st = json.dumps(dic2, indent=2)
    return st
  
  def export_json(self, output_json):
    with open(output_json, 'w') as fo:
      timer = Timer()
      json.dump(self.dic, fo)

      used_time = timer.elapse()
      print('%s written (%d sec)'%(
        output_json, used_time))

  def load_json(self, input_json):
    with open(input_json, 'r') as fi:
      assert self.dic=={}
      timer = Timer()
      self.dic = json.load(fi)

      used_time = timer.elapse()
      print('%s loaded (%d sec)'%(
        input_json, used_time))

class iMeSH_D_Tree():
  def __init__(self, args):
    '''
    Build a MeSH tree, for retrieval analysis

    Input:
      json_iMeSH_D: a json file containing iMeSH_D
    Output:
      updated members:
        imesh: iMeSH_D object containing 
          {
            DescriptorUI:{
              DescriptorName:str,
              TreeNumberList:[TreeNumber],
              ScopeNote:str
            }
          }
    '''
    json_iMeSH_D = args

    imesh = iMeSH_D()
    imesh.load_json(json_iMeSH_D)

    treeNumber2name_dic = {} #key is treeNumber and val is cell name
    for mid in imesh.get_mesh_ids():
      itm = imesh.get_mesh_itm(mid)
      name = itm.get("DescriptorName", "")
      for treeNumber in itm.get("TreeNumberList", []):
        if len(treeNumber)>=len("A11") and treeNumber[:3]=="A11":
          treeNumber2name_dic[treeNumber] = mid + ":" + name

    treeNumber_list = list(treeNumber2name_dic.keys())

    tree = Tree()
    tree.build_tree(treeNumber_list, treeNumber2name_dic)

    self.imesh = imesh 
    self.tree = tree 

    return

  def show_all_tree(self):
    self.tree.traverse(self.tree.get_root(), show=True)

  def get_cellTreeNumbers(self, mid):
    '''
    Return a list of Cell treeNumbers given a MeSH term ID

    Input:
      mid: MeSH term ID
    Output:
      res: a list of Cell treeNumbers
    '''
    itm = self.imesh.get_mesh_itm(mid)
    res = []
    for treeNumber in itm.get("TreeNumberList", []):
      if len(treeNumber)>=len("A11") and treeNumber[:3]=="A11":
        res.append(treeNumber)
    return res

  def tree_str(self, mid_list):
    '''
    Input:
      mid_list: a list of (MeSH ids, tag)
        tag is arbitrary strings to describe additional information for MeSH ids
    Output:
      string of a minimum subtree containing these mesh ids in mid_list 
    '''
    treeNumber2tag_dic = {}

    for mid,tag in mid_list:
      itm = self.imesh.get_mesh_itm(mid)
      for treeNumber in itm.get("TreeNumberList", []):
        if len(treeNumber)>=len("A11") and treeNumber[:3]=="A11":
          treeNumber2tag_dic[treeNumber] = tag

    treeNumber_set = set(list(treeNumber2tag_dic.keys()))
    ancestor_treeNumber, candidates_to_print = self.tree.find_subtree(treeNumber_set)
    return self.tree.traverse(self.tree.get_node(ancestor_treeNumber),
      level=0, candidates_to_print=candidates_to_print,
      treeNumber2tag_dic=treeNumber2tag_dic)




