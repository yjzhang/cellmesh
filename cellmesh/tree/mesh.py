from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import pdb

from cellmesh.tree.util import Timer, Count
'''
dic with tag:[default_val, path_for_tag]
  path_for_tag is used to find content string of the tag in XML Element
'''
MESH_D_TARGETS = { 
    'DescriptorUI':['', 'DescriptorUI'],
    'DescriptorName':['', 'DescriptorName/String'],
    'TreeNumberList':[[], 'TreeNumber'],
    'ScopeNote':['', 'ConceptList/Concept/ScopeNote']
  } 

MESH_D_TARGETS_LIST = [
    'DescriptorUI',
    'DescriptorName',
    'TreeNumberList',
    'ScopeNote'
  ]

MESH_D_FORMATTED_TARGETS_LIST = [
    'DescriptorUI',
    'DescriptorName'
  ]

MESH_D_UID='DescriptorUI'
MESH_D_NAME='DescriptorName'

def extract_target(XML_Element, targets):
  '''
  Input:
    XML_Element: an XML Element as source
    targets: dic with key as tag,
      val as [default tag val, path_for_tag]
  Output:
    dic: key as tag and val as extracted from XML_Element
      based on targets
  '''

  dic = {}
  for k, v_p in list(targets.items()):
    v, p = v_p
    if isinstance(v, str):
      dic[k] = ''
      q = XML_Element.find(p)
      if q is not None:
        dic[k] = q.text
    elif isinstance(v, list):
      dic[k] = []
      q = XML_Element.find(k)
      if q is not None:
        for r in q.iter(p):
          dic[k].append(r.text)
  return dic

class MeSH_D():
  '''
  an object for MeSH DescriptorRecord
  '''
  def __init__(self):
    self.dic = None
    return

  def init_from_XML_DR(self, XML_DR):
    '''
    Parse an XML Element to create members

    Inputs:
      XML_DR: XML Element for MeSH DescriptorRecord 
    '''
    assert self.dic is None

    self.dic = extract_target(XML_DR, MESH_D_TARGETS)
    
    return 

  def init_from_formatted_str(self, st):
    '''
    Parse a formatted str to create members
      based on MESH_D_FORMATTED_TARGETS_LIST

    Inputs:
      st: formatted st
    '''
    assert self.dic is None

    self.dic = {}

    items = [x for x in st.split(':') if x.strip()!='']
    N = len(MESH_D_FORMATTED_TARGETS_LIST)
    assert len(items)==N

    for i in range(N):
      k = MESH_D_FORMATTED_TARGETS_LIST[i]
      v = items[i]
      self.dic[k] = v
    return

  def formatted_str(self):
    '''
    return a formatted string describing self.dic
      e.g. <uniq_id>:<name>
    '''
    st = ''
    for k in MESH_D_FORMATTED_TARGETS_LIST:
      v = self.dic[k]
      if isinstance(v, str):
        v = v.strip()
        if v=='':
          st += 'EMPTY:'
        else:
          st += '%s:'%v
    st = st[:-1]
    return st

  def __str__(self):
    '''
    return a string describing self.dic
    '''
    st = ''
    for k in MESH_D_TARGETS_LIST:
      if k not in self.dic:
        #possible to print a MeSH D
        #init from formatted str
        continue
      v = self.dic[k]
      if isinstance(v, str):
        v = v.strip()
        if v=='':
          st += '%s:EMPTY'%k
        else:
          st += '%s:%s'%(k, v)
        st += '\n'
      elif isinstance(v, list):
        if v==[]:
          st += '%s:EMPTY'%k
        else:
          st += '%s:'%k
          for e in v:
            st += '%s,'%e
        st += '\n'
    return st

  def is_Cell(self):
    '''
    check if MeSH_D belongs to Cell (A11) in MeSH tree
    ''' 
    TNL = self.dic.get('TreeNumberList',[])
    if TNL==[]: return False
    for TN in TNL:
      if TN[0:min(len('A11'), len(TN))]=='A11':
        return True
    return False

  def get_uid(self):
    return self.dic.get(MESH_D_UID, '')

  def get_name(self):
    return self.dic.get(MESH_D_NAME, '')

def load_formatted_MeSH_D(formatted_file):
  '''
  Input:
    formatted_file: a formatted file containing MeSH_D records
      e.g. <uniq_id>:<name>
  Output:
    a dic with key as uniq_id and val as MeSH_D
  '''
  res = {}

  with open(formatted_file, 'r') as fi:
    for line in fi:
      if line!='' and line[0]=='#':
        continue
      dic = MeSH_D()
      dic.init_from_formatted_str(line)
      res[dic.get_uid()] = dic

  return res

def test_load_formatted_MeSH_D():
  #formatted_file = 'export_temp.txt'
  formatted_file = 'MeSH_desc2019_formatted.txt'
  res = load_formatted_MeSH_D(formatted_file)
  pdb.set_trace()
  return

def test_export_cell_MeSH_D():

  '''
  meshFile = 'mesh_desc2019_first20records.xml'
  export_file = 'export_temp.txt'
  debug = True
  '''

  meshFile = '/data1/shunfu/scLit/mesh/desc2019.xml'
  export_file = 'MeSH_desc2019_formatted.txt'
  debug = False

  export_cell_MeSH_D(meshFile, export_file, debug)
  return 

if __name__ == "__main__":
  #test_select_cell_MeSH_D()
  #test_export_cell_MeSH_D()
  test_load_formatted_MeSH_D()
