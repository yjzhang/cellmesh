from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

import time
import pdb
import json
import sys
import numpy
import scipy

class Timer:
  def __init__(self):
    self.prev = time.time()

  def reset(self):
    self.prev = time.time()

  def elapse(self):
    elapse_time = time.time() - self.prev
    self.prev = time.time()
    return elapse_time 

def load_json(input_json):
  timer = Timer()
  with open(input_json, 'r') as fi:
    dic = json.load(fi)

    used_time = timer.elapse()
    print('%s loaded (%d sec)'%(
      input_json, used_time))
    return dic
  
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