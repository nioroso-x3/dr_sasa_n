#!/usr/bin/python3
"""
argv1: interaction table file name
"""
import sys
import numpy as np
import pandas as pd
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from os import path
def ExtractData(fname,add):
  """
  Receives filename of tab separated table with row names and column names
  Returns Pandas data frame and maxvalue
  """
  maxvalue = 0
  data = []
  colnames = []
  rownames = []
  lines = open(fname).readlines()
  lines = [line.strip() for line in lines]
  tag = lines[0].split("\t")[0]
  buried = ""
  buries = ""

  if "->" in tag:
    tk = tag.split("->")
    buried = tk[1]
    buries = tk[0].split("/")[0]
  if "<-" in tag:
    tk = tag.split("<-")
    buried = tk[0]
    buries = tk[1].split("/")[0]
 
  colnames = [item.strip() for item in lines[0].split("\t")[1:]]
  for line in lines[1:]:
    tokens = line.split("\t")
    rownames.append(tokens[0].strip())
    values = [float(item)+add for item in tokens[1:]]
    m = max(values)
    if(m > maxvalue):
      maxvalue = m
    data.append(values)
  npdata = np.array(data)
  dataframe = pd.DataFrame(npdata,index=rownames, columns=colnames)
  return dataframe,maxvalue,buried,buries


if __name__ == "__main__":
  color="Greys"
  sns.set_style()
  
  
  parser = argparse.ArgumentParser(description='Graphs heatmap from labeled matrix text file')
  parser.add_argument('infile',type=str)
  parser.add_argument('-s','--skip',type=int)
  parser.add_argument('-l','--log',action="store_true")
  parser.add_argument('-S','--scale',type=float)
  args = parser.parse_args()
  if args.scale:
   sns.set_context("poster",font_scale=args.scale)
  else:
   sns.set_context("poster",font_scale=2)
  fig = plt.figure(figsize=(48,36))
  if args.log:
    data,mvalue,buried,buries = ExtractData(args.infile,1)
  else:
    data,mvalue,buried,buries = ExtractData(args.infile,0)

  labeltick = 1
  if args.skip:
    labeltick = int(args.skip)
  
  mvalue = int((mvalue % 10)*0.1  + mvalue/10.0)*10
  #fig.suptitle('ASDF',fontsize= 96, fontweight='bold')
  if labeltick == 1:
    if not args.log:
      ax = sns.heatmap(data,square=True,vmin=0,linecolor="black",linewidth=0,cmap=color)
    else:
      ax = sns.heatmap(data,square=True,linecolor="black",linewidth=0,cmap=color,
                       norm=LogNorm(vmin=0,vmax=data.max()))
  else:
    if not args.log:
      ax = sns.heatmap(data,square=True,vmin=0,linecolor="black",linewidth=0,
      xticklabels=labeltick,yticklabels=labeltick,cmap=color)
    else:
      ax = sns.heatmap(data,square=True,linecolor="black",linewidth=0,
      xticklabels=labeltick,yticklabels=labeltick,cmap=color,
      norm=LogNorm(vmin=0,vmax=data.max()))

  plt.savefig(path.basename(args.infile)+".png",bbox_inches="tight",dpi=150)
