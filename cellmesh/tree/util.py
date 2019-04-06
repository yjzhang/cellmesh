from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import time
import pdb
import subprocess
import json
import sys
import numpy
import scipy

def get_matrix_data_hist_str(mat, bins=None):
  '''
  Input:
    mat: could be scipy sparse csr_matrix, csc_matrix or np.ndarray
  Output:
    string of hist, bins of mat data
  '''
  if isinstance(mat, numpy.ndarray):
    data = mat
  elif isinstance(mat, scipy.sparse.csc_matrix):
    data = mat.data
  elif isinstance(mat, scipy.sparse.csr_matrix):
    data = mat.data

  if bins is None:
    hist, bins = numpy.histogram(data)
  else:
    hist, bins = numpy.histogram(data, bins)

  st = ''
  for i in range(len(hist)):
    if i==len(hist)-1:
      right_brac = "]"
    else:
      right_brac = ")"
    st += '[%f, %f%s %d\n'%(bins[i], bins[i+1], right_brac, hist[i])
  return st

def run_cmd(cmd, shell=True, show_cmd=True):
  if show_cmd:
    print(cmd)
  subprocess.call(cmd, shell=shell)

class Timer:
  def __init__(self):
    self.prev = time.time()

  def reset(self):
    self.prev = time.time()

  def elapse(self):
    elapse_time = time.time() - self.prev
    self.prev = time.time()
    return elapse_time 

class CoCount:
  '''
  an object to store co-occurrence statistics

  self.dic:
    key: co-occured (k1,k2,...,k_N)
    val: {'cnt': cnt, 'place': [p1,p2,...,p_T]}
      p1 ~ p_T are places (e.g. an article id) where co-occurred terms (k1,...,K_N) appear
        and p is expected to be unique
      due to memory concern, we may not include all places (i.e. T is a threshold, T<=cnt)
  '''
  def __init__(self, T=-1):
    '''
    T is a threshold to store max number of places the co-occurence key appears. -1 no limit.
    '''
    self.dic = {}
    self.timer = Timer()
    self.T = T
    return

  def update(self, k, p):
    '''
    Input:
      k: is a key representing co-occurrence info
      p: is a place info where k appears (e.g. k=(gene_A, cell_B) appears in p=pmid C)
    '''
    if k not in self.dic:
      self.dic[k] = {'cnt':0, 'place':[]}

    self.dic[k]['cnt']+=1

    if self.T==-1 or len(self.dic[k]['place'])<self.T:
      self.dic[k]['place'].append(p)

    return

  def update_cc(self, k, dic):
    '''
    update from a CoCount (cc) item

    Input:
      k: a key representing co-occurrence info
      dic: {'cnt': cnt, 'place': [p1,p2,...,p_T]} wrt k in the cc item
    Output:
      self.dic updated

    Note:
      p is expected to be unique
    '''
    if k not in self.dic:
      self.dic[k] = {'cnt':0, 'place': []}

    self.dic[k]['cnt'] += dic['cnt']
    self.dic[k]['place'] += dic['place']

    if self.T>0 and len(self.dic[k]['place'])>self.T:
      self.dic[k]['place'] = self.dic[k]['place'][0:self.T]

    return

  def merge(self, cc):
    '''
    merge another CoCount into self

    Input:
      cc: a CoCount object to be merged
    Output:
      the self CoCount object
    '''
    for k in cc.dic.keys():
      self.update_cc(k, cc.dic[k])
    return self

  def str_json(self):
    '''
    return a json string

    Note:
      key needs to be a string
    '''
    st = json.dumps(self.dic,indent=2)
    return st


class Count:
  '''
  Basic Usage:
    c = Count()
    ...
    c.update(count_key)
    ...
    print c.get(count_key)
  '''
  def __init__(self):
    self.dic = {}
    self.timer = Timer()
    return

  def update(self, k, v=1):
    if k not in self.dic:
      self.dic[k] = 0
    self.dic[k]+=v
    return

  def get(self, k):
    return self.dic.get(k, 0)

  def merge(self, c):
    '''
    merge Count c into self
    '''
    for k, v in list(c.dic.items()):
      if k not in self.dic:
        self.dic[k] = v
      else:
        self.dic[k] += v
    return

  def __str__(self):
    st = ''
    items = sorted(list(self.dic.items()), key=lambda x:x[0])
    for k, cnt in items:
      st += '%s:%d\n'%(k, cnt)
    return st

  def formatted_str(self, \
    sort_v=True, show_idx=True, restrict_key_set=None, key2description=None):
    '''
    Input:
      sort_v: sort based on v in descreasing order
      show_idx: show idx
      restrict_key_set: if not None, show only k in restrict_key_set
      key2description: {key of self.dic: description string}
        provide additional info. e.g. key of self.dic is gene id,
        we also want to show its name
    '''
    st = ''
    if sort_v is False: #sort by key
      items = sorted(list(self.dic.items()), key=lambda x:x[0])
    else:
      items = sorted(list(self.dic.items()), key=lambda x:-x[1])

    cnt_keys = 0
    for i in range(len(items)):
      k, cnt = items[i]
      if restrict_key_set is not None and k not in restrict_key_set:
        continue

      cnt_keys += 1
      if show_idx is True:
        st1 = '(%d) '%i
      else:
        st1 = ''
      
      if key2description is not None and k in key2description:
        st += st1 + '%s {%s}:%d\n'%(str(k), key2description[k], cnt)
      else:
        st += st1 + '%s:%d\n'%(str(k), cnt)
    st = 'Count.formatted_str: %d keys listed\n'%cnt_keys + st
    return st

  def export_json(self, output_json):
    with open(output_json, 'w') as fo:
      self.timer.reset()
      json.dump(self.dic, fo)

      used_time = self.timer.elapse()
      print('%s written (%d sec)'%(
        output_json, used_time))

  def load_json(self, input_json):
    with open(input_json, 'r') as fi:
      assert self.dic=={}
      self.timer.reset()
      self.dic = json.load(fi)

      used_time = self.timer.elapse()
      print('%s loaded (%d sec)'%(
        input_json, used_time))

def write2json(json_dic, output_json):
  '''
  write json_dic to json_path

  Note:
    indent does not seem to hurt load_json
  '''
  timer = Timer()
  with open(output_json, 'w') as fo:
    json.dump(json_dic, fo, indent=2)

    used_time = timer.elapse()
    print('%s written (%d sec)'%(
      output_json, used_time))
  return

def load_json(input_json):
  timer = Timer()
  with open(input_json, 'r') as fi:
    dic = json.load(fi)

    used_time = timer.elapse()
    print('%s loaded (%d sec)'%(
      input_json, used_time))
    return dic

class Hist:
  def __init__(self, cnt_list):
    '''
    Input:
      cnt_list: list of **int** values
    Output:
      members of Hist
    '''
    self.hist = {} #key: an int number val: freq of k appears in cnt_list
    for v in cnt_list:
      if v not in self.hist:
        self.hist[v] = 0
      self.hist[v] += 1
    return

  def __str__(self):
    '''
    return a string as:
      key:freq_of_key; key sorted from small to large
    '''
    items = sorted(list(self.hist.items()), key=lambda x: x[0])
    st = ''
    for itm in items:
      st += '%d:%d\n'%(itm[0],itm[1])
    return st

  def merge(self, h):
    '''
    merge another Hist object h into self
    '''
    for k, f in list(h.hist.items()):
      if k in self.hist:
        self.hist[k] += f
      else:
        self.hist[k] = f
    return
    
class Log:
  '''
  a log class
    - to be used for debugging. skip log if no log_file specified.

  members:
    log_path: path for logging
    fo: file handler for log_path
    exception_counter: cnt of write exceptions
      if a write exception occurs, an optional tag ##EXCEPTION<idx>## will occur
      before the contents
  '''
  def __init__(self, log_file=None):
    self.log_path = ''
    self.fo = None
    self.exception_counter = 0
    if log_file is not None:
      self.log_path = log_file
      self.fo = open(log_file, 'w')
    return

  def logging(self, msg, verbose=False, exception_tag=False):
    if verbose is True:
      print(msg)
    if self.fo is not None:
      try:
        self.fo.write(msg)
      except:
        self.exception_counter += 1
        if exception_tag is False:
          exception_st = ''
        else:
          exception_st = '##EXCEPTION%d##\n'%self.exception_counter
        self.fo.write(exception_st+str(msg))
    return

  def close(self):
    if self.fo is not None:
      print('%s written'%self.log_path)
      self.fo.close()

def json_print_dic(dic):
  '''
  print a dictionary in json string format
  '''
  #TODO(serialize a customed class)
  st = json.dumps(dic,indent=2)
  print(st)
  return

class iterCounter:
  def __init__(self, N, msg):
    self.N = N
    self.msg = msg
    self.T = N/100
    self.p = 0
    self.q = 0
    return

  def finish(self):
    print('')#print an empty line for display purpose
    return

  def inc(self):
    if self.N<100:
      self.p += 1
      sys.stdout.write('\r');
      sys.stdout.write('[%s] %s: %.2f %%'%(str(time.asctime()), self.msg, float(self.p)*100.0/self.N));
      sys.stdout.flush()
    else:
      self.p += 1
      if self.p >= self.T:
        self.p = 0
        self.q += 1
        sys.stdout.write('\r');
        sys.stdout.write('[%s] %s: %d %%'%(str(time.asctime()), self.msg, self.q));
        sys.stdout.flush()
    return  