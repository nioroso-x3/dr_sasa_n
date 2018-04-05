#!/usr/bin/python3
"""
argv1: RESI of the target residue
"""
import sys
from glob import glob
import numpy as np

if __name__ == "__main__":
  fnames = glob("*.asa")
  data_atom = {}
  data_res = {}
  resi = sys.argv[1]
  for fname in fnames:
    resdict = {}
    lines = open(fname).readlines()
    for line in lines:
      tokens = line.split("\t")
      if tokens[0] == "ATOM":
        if tokens[5] == resi:
          t = int(tokens[12].strip())
          resid = tuple(tokens[3:6])
          
          if not t in data_atom:
            data_atom[t] = []
          if not resid in resdict:
            resdict[resid] = 0
          data_atom[t].append(float(tokens[10]))
          resdict[resid]+=(float(tokens[10]))
    for key in resdict:
      if not key[0] in data_res:
        data_res[key[0]] = []
      data_res[key[0]].append(resdict[key])
  for t in data_atom:
    avg = np.average(data_atom[t])
    sys.stdout.write("{"+str(t)+","+str(avg)+"},\n")
  print("\n")
  for t in data_res:
    avg = np.average(data_res[t])
    sys.stdout.write("{\""+str(t)+"\","+str(avg)+"},\n")
