from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

from cellmesh.tree.util import Timer, iterCounter

class TreeNode():
  def __init__(self, tn):
    '''
    Inputs:
      tn: tree number
    '''
    self.dic = {}
    self.dic["tn"] = tn
    self.dic["name"] = tn #name to describe the TreeNode (e.g. tn is MeSH TreeNumber and name is MeSH tern)
    self.dic["height"] = -1 #distance from TreeNode to a Tree's root node

    self.dic["parent"] = None # another TreeNode
    self.dic["children"] = [] # list of TreeNodes
    return

  def add_child(self, treeNode):
    self.dic["children"].append(treeNode)
    return

  def add_parent(self, treeNode):
    assert self.dic["parent"] is None
    self.dic["parent"] = treeNode
    return

  def no_parent(self):
    return self.dic["parent"] is None

  def get_treeNumber(self):
    return self.dic["tn"]

  def get_name(self, show_details=False, tag=""):
    st = self.dic["name"]
    if show_details:
      st += "(tn=%s,height=%d)"%(
        self.get_treeNumber(), self.get_height())
    if tag != "":
      st += " -->%s"%tag
    return st

  def set_name(self, name):
    self.dic["name"] = name

  def set_height(self, h):
    self.dic["height"] = h 

  def get_height(self):
    return self.dic["height"]

  def iter_children(self):
    for c in self.dic["children"]:
      yield c

class Tree():
  def __init__(self):
    '''
    members:
      "nodes": {MeSH tree number: TreeNode}
      "root": root TreeNode
        exists in "nodes"
      "pairwise_dist": {"tn1,tn2": dist}
        could be lazily updated, serving as a cache
        to save space, tn1<tn2
    '''
    self.dic = {
      "nodes": {},
      "root": None,
      "pairwise_dist": {}
    }
    return

  def add_node(self, tn):
    '''
    Input:
      tn: tree number
    '''
    self.dic["nodes"][tn] = TreeNode(tn)
    return

  def get_node(self, tn):
    return self.dic["nodes"].get(tn, None)

  def add_dist(self, tn1, tn2, dist):
    '''
    Input:
      tn1: treeNumber 1
      tn2: treeNumber 2
      dist: calculated dist b/w tn1 and tn2
    Output:
      self.dic["pairwise_dist"] updated
    '''
    if tn1==tn2:
      assert dist==0
      return
    elif tn1<tn2:
      k = "%s,%s"%(tn1,tn2)
      self.dic["pairwise_dist"][k] = dist
      return 
    else:
      k = "%s,%s"%(tn2,tn1)
      self.dic["pairwise_dist"][k] = dist
      return

  def get_dist(self, tn1, tn2):
    '''
    retrieve the distance b/w tn1 and tn2 from self.dic["pairwise_dist"]
    if no dist available, return -1
    '''
    if tn1==tn2:
      return 0
    elif tn1<tn2:
      k = "%s,%s"%(tn1,tn2)
      return self.dic["pairwise_dist"].get(k, -1)
    else:
      k = "%s,%s"%(tn2,tn1)
      return self.dic["pairwise_dist"].get(k, -1)
      

  def get_treeNumber_list(self):
    return list(self.dic["nodes"].keys())

  def get_treeNode_list(self):
    return self.dic["nodes"].values()

  def exists(self,tn):
    return tn in self.dic["nodes"]

  def connect(self, tn_p, tn_c):
    '''
    Input:
      tn_p: parent tree number
      tn_c: child tree number
    Output:
      tree_node_p (wrt tn_p) and tree_node_c (wrt tn_c) connected
    '''
    tree_node_p = self.get_node(tn_p)
    tree_node_c = self.get_node(tn_c)
    assert tree_node_p is not None and tree_node_c is not None
    tree_node_p.add_child(tree_node_c)
    tree_node_c.add_parent(tree_node_p)
    return

  def set_root(self, treeNode):
    self.dic["root"] = treeNode
    return

  def get_root(self):
    return self.dic["root"]

  def build_tree(self, tn_list, tn2name_dic=None):
    '''
    Input:
      tn_list: a list of tree numbers
        For example, tree number can be A11, A11.031 etc
        if tree_node_1 has tree_number tn1, and tree_node_2 has tn2,
          and tn1 + '.<num>' = tn2, then tree_node_2 is a child of tree_node_1
      tn2name_dic: None or a dic {treeNumber: treeName}
    '''

    # nodes
    for tn in tn_list:
      self.add_node(tn)

    # edges
    for tn2 in tn_list:
      tn1 = [t for t in tn2.split(".") if t!=""]
      if len(tn1)<=1: continue

      tn1 = ".".join(tn1[:-1])
      if self.exists(tn1):
        self.connect(tn1, tn2)

    # root
    sources = []
    for treeNode in self.get_treeNode_list():
      if treeNode.no_parent():
        sources.append(treeNode)
    assert len(sources)>0
    if len(sources)==1:
      self.set_root(sources[0])
    else:
      tn = "root"
      self.add_node(tn)
      for s in sources:
        self.connect(tn, s.get_treeNumber())
      self.set_root(self.get_node(tn))

    # heights and names
    self.traverse_to_update_tree_nodes(self.get_root(), 0, tn2name_dic)

    return

  def traverse_to_update_tree_nodes(self, treeNode, level, tn2name_dic=None):
    '''
    traverse the tree, starting from treeNode, to update
      each node's level/height and name (if tn2name_dic is not None)

    Inputs:
      treeNode: current "root" node to start the traverse
      level: or height from the true tree root node
      tn2name_dic: None or a dic {treeNumber: treeName}
    '''
    treeNode.set_height(level)
    
    tn = treeNode.get_treeNumber()
    if tn2name_dic is not None and tn in tn2name_dic: 
      treeNode.set_name(tn2name_dic[tn])

    for c in treeNode.iter_children():
      self.traverse_to_update_tree_nodes(c, level+1, tn2name_dic)

    return

  def traverse(self, treeNode, level=0, candidates_to_print=None, show=False, treeNumber2tag_dic=None):
    '''
    pre-order traverse the tree starting from treeNode.

    candidates_to_print (a set of treeNumbers)
      if None, only print treeNodes whose treeNumber is in candidates_to_print 

    show=True to print the traverse process

    treeNumber2tag_dic: {treeNumber: tag str}
      if None, the node's name will be augmented by tag str if the node has treeNumber in treeNumber2tag_dic
      would be useful to tag the retrieved cells in the tree-structured string
    '''
    st = ''
    if level>0:
      st += "|" + "--"*level

    tn = treeNode.get_treeNumber()

    if treeNumber2tag_dic is not None and tn in treeNumber2tag_dic:
      tag_st = treeNumber2tag_dic[tn]
      st += treeNode.get_name(tag=tag_st)
    else:
      st += treeNode.get_name()

    res_str = ''
    if candidates_to_print is None or tn in candidates_to_print:
      res_str += st + "\n"
      if show is True:
        print(st)
    for c in treeNode.iter_children():
      res_str += self.traverse(c, level+1, candidates_to_print, show, treeNumber2tag_dic)
    return res_str

  def find_subtree(self, treeNumber_set):
    '''
    find a treeNumber which is the lowest common ancestor of the treeNumber_set

    Output:
      tn: the treeNumber which is the lowest common ancestor of the treeNumber_set
      candidate_set: the treeNumbers that have themselves and their children overlapping with treeNumber_set
        traverse starting from tn, they shall limit the candidates to be printed

    Note:
      a typical usage is:
        ancestor_treeNumber, candidates_to_print = find_subtree(treeNumber_set)
        traverse(self.get_node(ancestor_treeNumber), level=0, candidates_to_print)
    '''

    #key: treeNumber tn, val: {"coverage": number of treeNumber_set covered by tn's subtree (including tn),
    #                          "level": level}
    status_dic = {}
    for tn in self.get_treeNumber_list():
      status_dic[tn] = {"coverage":0, "level":0}

    def update(treeNode, level, status_dic, treeNumber_set):
      '''
      recursively update the status_dic for treeNode's children and itself
      '''
      coverage = 0
      treeNumber = treeNode.get_treeNumber()

      for c in treeNode.iter_children():
        update(c, level+1, status_dic, treeNumber_set)
        child_cov = status_dic[c.get_treeNumber()]["coverage"]
        coverage += child_cov

        # early stop
        if child_cov == len(treeNumber_set):
          break

      if treeNumber in treeNumber_set:
        coverage += 1

      status_dic[treeNumber]["coverage"] = coverage
      status_dic[treeNumber]["level"] = level
      return

    update(self.get_root(), level=0, status_dic=status_dic, treeNumber_set=treeNumber_set)

    itms = sorted([y for y in list(status_dic.items()) if y[1]["coverage"]==len(treeNumber_set)], \
      key=lambda x: -x[1]["level"])

    candidate_set = set([y[0] for y in list(status_dic.items()) if y[1]["coverage"]>0])

    return itms[0][0], candidate_set

  def calc_two_node_dist(self, tn1, tn2):
    '''
    Inputs:
      tn1: treeNumber 1
      tn2: treeNumber 2
    Output:
      distance between tn1 and tn2

    Note:
      - order of tn1, tn2 does not matter
      - will utilize self.dic["pairwise_dist"] as a cache, also update it
      - Basic idea is to calculate LCA of tn1 and tn2, then
        distance = Level(tn1) + Level(tn2) - 2 Level(LCA)
      - another method is to use tree number structure
    '''
    dist = self.get_dist(tn1, tn2)

    if dist>=0:
      return dist

    # no cache value found
    '''
    tn_ann, _ = self.find_subtree(set([tn1, tn2]))
    '''
    # another method is to use tree number directly
    tn1_chunks = tn1.split("."); L1 = len(tn1_chunks)
    tn2_chunks = tn2.split("."); L2 = len(tn2_chunks)

    #print(tn1_chunks)
    #print(tn2_chunks)

    cnt = min(L1, L2) #need this init val to avoid bug
    for i in range(min(L1, L2)):
      if tn1_chunks[i]!=tn2_chunks[i]:
        cnt = i
        #print("i=%d break"%i)
        break
    tn_ann = ".".join(tn1_chunks[:cnt])

    d1 = self.get_node(tn1).get_height()
    d2 = self.get_node(tn2).get_height()
    #try:
    d_ann = self.get_node(tn_ann).get_height()
    #except:
    #  pdb.set_trace()
    d = d1 + d2 - 2*d_ann

    self.add_dist(tn1, tn2, d)

    return d

  def calc_pairwise_node_dist(self, max_cnt=-1):
    '''
    Input:
      max_cnt: -1 (disable) or >0 to limit cnts processed (for debug purpose)
    Output:
      self.dic["pairwise_dist"] = {"tn1,tn2":dist} (tn1<tn2)

    Note:
      Basic idea is to calculate LCA of tn1 and tn2, then
        distance = Level(tn1) + Level(tn2) - 2 Level(LCA)
    '''
    treeNumber_list = self.get_treeNumber_list()
    N = len(treeNumber_list)
    N_distinct_pairs = N*(N-1)/2

    cnt = 0
    timer = Timer()
    ic = iterCounter(N_distinct_pairs,
      "calc_pairwise_node_dist: %d treeNumber, %d distinct pairs"%(N, N_distinct_pairs))

    for i in range(len(treeNumber_list)):
      for j in range(i+1, len(treeNumber_list)):
        cnt += 1
        if max_cnt>0 and cnt>=max_cnt:
          break
        ic.inc()

        tn_i = treeNumber_list[i]
        tn_j = treeNumber_list[j]
        d_ij = self.calc_two_node_dist(tn_i, tn_j)

    ic.finish()

    return self.dic["pairwise_dist"]





