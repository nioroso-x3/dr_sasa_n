#!/usr/bin/python
"""
argv1: file to open
argv2: output filename (optional)
"""
import sys
import os
from os.path import basename,splitext
def GetTable(fname):
  lines = open(fname).readlines()
  result = []
  firstline = True
  for line in lines:
    tokens = line.split("\t")
    tokens = [item for item in tokens if not item in ["\n","\r"] ]
    if firstline:
      result.append(tokens)
      firstline = False
    else:
      values = [float(item) for item in tokens[1:]]
      tag = [tokens[0]]
      tag.extend(values)
      result.append(tag)
  return result
if __name__ == "__main__":
  table = GetTable(sys.argv[1])
  outfname = ""
  try:
    outfname = sys.argv[2]
  except:
    outfname = splitext(sys.argv[1])[0]+".cut.tsv"
  cols = len(table[1])
  rows = len(table)
  del_row = []
  del_col = []
  for i in range(1,len(table)):
    linesum = sum(table[i][1:])
    if linesum == 0:
      del_row.append(i)

  for j in range(1,cols):
    curcol = [item[j] for item in table[1:]]
    colsum = sum(curcol)
    if(colsum == 0):
      del_col.append(j)

  #print "Rows deleted:",del_row
  #print "Cols deleted:",del_col

  ntable = []
  for i in range(rows):
    if i in del_row:
      continue
    else:
      ntable.append(table[i])
  ntable2 = []
  for j in range(cols):
    curcol = [item[j] for item in ntable]
    if j in del_col:
      continue
    else:
      ntable2.append(curcol)
  outf = open(outfname,"w")
  for i in range(rows - len(del_row)):
    currow = [item[i] for item in ntable2]
    outf.write("\t".join([str(item) for item in currow])+ "\n" )
  outf.close()