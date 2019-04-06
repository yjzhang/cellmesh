from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import pdb
import json
import re
import numpy as np

from cellmesh.tree.util import Timer, Log, iterCounter, write2json, load_json
from cellmesh.tree.tree import TreeNode, Tree
from cellmesh.tree.mesh import MeSH_D

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

def get_overlap(term, mesh_name):
  '''
  find list of overlap chunks between a term (e.g. cell name from cell ontology) and mesh term name

  Note:
    - term/mesh_name split by delimiters
    - augment term/mesh_name by synonyms
    - not consider stop words (e.g. cell/cells etc) as overlap
  '''

  delimiters = [' ', ',', '-','\(','\)']
  delimiter_st = '|'.join(delimiters)

  synonym_pairs = [['cancer', 'tumor']]
  stop_words = ['cell', 'cells', 'stem', 'like', 'derived', 'system',
    'associated', 'small', 'progenitor']  
  
  term_chunks = [plural2single(t) for t in re.split(delimiter_st, term) if t!='' and t not in stop_words]

  mesh_chunks = [plural2single(t) for t in re.split(delimiter_st, mesh_name) if t!='' and t not in stop_words]

  for s1, s2 in synonym_pairs:
    if s1 in term_chunks and s2 in mesh_chunks:
      term_chunks.append(s2)
    if s2 in term_chunks and s1 in mesh_chunks:
      term_chunks.append(s1)

  term_set = set(term_chunks)
  mesh_chunks = set(mesh_chunks) 
  return list(term_set.intersection(mesh_chunks))

def plural2single(t):
  '''
  convert a chunk from plural to single if necessary
  '''
  if t!='' and t[-1]=='s':
    if len(t)>=len('ies') and t[-3:]=='ies':
      return t.lower()[:-3]+'y'
    else:
      return t.lower()[:-1]
  else:
    return t.lower()

def get_plural_chunks(name):
  '''
  Input:
    name: a term name (e.g. Cell Ontology cell name, or MeSH cell name)
  Output:
    a list of name's lower-case chunks
      if chunk ends with s (plural form) and converted to single form
  Note:
    - term/mesh_name split by delimiters
    - '..ies' ==> '..y' (e.g. bodies to body)
  '''
  delimiters = [' ', ',', '-', '\(', '\)']
  delimiter_st = '|'.join(delimiters)
  #chunks = [t.lower() for t in re.split(delimiter_st, name) if t!='' and t[-1]=='s']
  chunks = []
  for t in re.split(delimiter_st, name):
    t = plural2single(t)
    if t!='':
      chunks.append(t)

  return chunks

def get_matched_MeSH_terms(args):
  '''
  Input:
    terms: a list of terms (e.g. cell names from Cell Ontology)
    json_iMeSH_D: a json file path to iMeSH_D indexing
  Output:
    res_dic: key as term in terms, val as {'MeSH term ID':['MeSH term name', list of overlap terms]}
      the MeSH term name shares overlap with original term
      it's possible res_dic[term]={}
    log: None or dump res_dic info for debug purpose
  '''
  terms, json_iMeSH_D, log = args
  imesh = iMeSH_D()
  imesh.load_json(json_iMeSH_D)

  logger = Log(log)
  '''
  debug - check chunks of a term ending at 's'
  '''
  if False:
    term_dic = {} #key as term, val as list of chunks ending at s
    for i in range(len(terms)):
      term = terms[i]
      plural_chunks = get_plural_chunks(term)
      if plural_chunks != []:
        term_dic[term] = plural_chunks
    if log is not None:
      logger.logging('----- term plural chunks -----\n')
      itms = sorted(list(term_dic.items()), key=lambda x: -len(x[1]))
      for itm in itms:
        logger.logging('%s\t%s\n'%(itm[0], str(itm[1])))
      logger.logging('\n')

  '''
  debug - check chunks of a MeSH term ending at 's'
  '''
  if False:
    term_dic = {} #key as MeSH term, val as list of chunks ending at s
    for mesh_id in imesh.get_mesh_ids():
      mesh_name = imesh.get_mesh_name(mesh_id)
      plural_chunks = get_plural_chunks(mesh_name)
      if plural_chunks != []:
        term_dic[mesh_name] = plural_chunks
    if log is not None:
      logger.logging('----- term plural chunks -----\n')
      itms = sorted(list(term_dic.items()), key=lambda x: -len(x[1]))
      for itm in itms:
        logger.logging('%s\t%s\n'%(itm[0], str(itm[1])))
      logger.logging('\n')

  res_dic = {}
  i = 0
  for term in terms:
    print('%d-th term: %s'%(i, term))
    i+=1
    for mesh_id in imesh.get_mesh_ids():
      mesh_name = imesh.get_mesh_name(mesh_id)
      if mesh_name == '': continue

      overlap_list = get_overlap(term, mesh_name)
      if overlap_list != []:
        if term not in res_dic: res_dic[term] = {}
        res_dic[term][mesh_id] = [mesh_name, overlap_list]
    if term not in res_dic:
      res_dic[term] = {}

  if log is not None:
    logger.logging('----- res dic -----\n')
    logger.logging(json.dumps(res_dic,indent=2)+'\n')
    logger.logging('\n')
  logger.close()
  return res_dic

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

  def add_dist(self, mid1, mid2, dist):
    '''
    Input:
      mid1: MeSH term ID 1
      mid2: MeSH term ID 2
      dist: calculated dist b/w mid1 and mid2
    Output:
      self.mesh_pairwise_dist updated
    '''
    if mid1==mid2:
      assert dist==0
      return
    elif mid1<mid2:
      k = "%s,%s"%(mid1,mid2)
      if k not in self.mesh_pairwise_dist:
        self.mesh_pairwise_dist[k] = dist
      return 
    else:
      k = "%s,%s"%(mid2,mid1)
      if k not in self.mesh_pairwise_dist:
        self.mesh_pairwise_dist[k] = dist
      return

  def get_dist(self, mid1, mid2):
    '''
    retrieve the distance b/w mid1 and mid2 from self.mesh_pairwise_dist
    if no dist available, return -1
    '''
    if mid1==mid2:
      return 0
    elif mid1<mid2:
      k = "%s,%s"%(mid1,mid2)
      return self.mesh_pairwise_dist.get(k, -1)
    else:
      k = "%s,%s"%(mid2,mid1)
      return self.mesh_pairwise_dist.get(k, -1)

  def calc_and_store_all_mesh_pairwise_dist(self):
    '''
    calc all (mid1, mid2) dist and update self.mesh_pairwise_dist

    Output:
      self.mesh_pairwise_dist
    '''
    mid_list = list(self.imesh.get_mesh_ids())
    L = len(mid_list)
    L_distinct_pairs = L*(L-1)/2

    #cnt = 0
    timer = Timer()
    ic = iterCounter(L_distinct_pairs,
      "calc_and_store_all_mesh_pairwise_dist: %d MeSH IDs, %d distinct pairs"%(
      L, L_distinct_pairs))

    for i in range(L):
      for j in range(i+1, L):
        ic.inc()

        mid1 = mid_list[i]
        mid2 = mid_list[j]
        d = self.calc_and_store_two_mid_dist(mid1, mid2)
    
    ic.finish()
    print("calc_and_store_all_mesh_pairwise_dist (%d sec)"%timer.elapse())

    return self.mesh_pairwise_dist

  def calc_and_store_two_mid_dist(self, mid1, mid2):
    '''
    Inputs:
      mid1: MeSH term ID 1
      mid2: MeSH term ID 2
    Output:
      distance between mid1 and mid2

    Note:
      - order of mid1, mid2 does not matter
      - will utilize self.mesh_pairwise_dist as a cache, also update it
      - Different from Tree.calc_and_store_two_mid_dist,
        each mesh id wrt several tree numbers, will calc
        min_t1,t2 dist(t1,t2) for t1 in mid1 and t2 in mid2
    '''
    dist = self.get_dist(mid1, mid2)

    if dist>=0:
      return dist

    # no cache value found
    mid1_tns = self.get_cellTreeNumbers(mid1)
    mid2_tns = self.get_cellTreeNumbers(mid2)
    assert len(mid1_tns)>0 and len(mid2_tns)>0

    min_dist = np.inf
    for tn1 in mid1_tns:
      for tn2 in mid2_tns:
          d = self.tree.calc_two_node_dist(tn1, tn2)
          if d<min_dist:
            min_dist = d

    self.add_dist(mid1, mid2, min_dist)

    return min_dist

  def export_all_mesh_pairwise_dist(self, json_file):
    '''
    Output:
      json_file stores self.mesh_pairwise_dist
    '''
    write2json(self.mesh_pairwise_dist, json_file)
    return

  def load_all_mesh_pairwise_dist(self, json_file):
    '''
    Input:
      json_file: containing the stored self.mesh_pairwise_dist info
    Output:
      self.mesh_pairwise_dist loaded from json file
    '''
    dic = load_json(json_file)
    for k in dic.keys():
      self.mesh_pairwise_dist[k] = int(dic[k])
    return

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




