#!/usr/bin/python3
"""
Converts dr_sasa matrixes to relative BSA based on ALA-X-ALA tripeptides and 9 nucleotide fibers
argv1: matrix
argv2: ATOM or RES
"""
import sys

def GetTable(fname):
  lines = open(fname).readlines()
  result = []
  axes = [[],[]]
  firstline = True
  for line in lines:
    tokens = line.strip().split("\t")
    if firstline:
      axes[0] = tokens
      firstline = False
    else:
      axes[1].append(tokens[0])
      values = [float(item) for item in tokens[1:]]
      result.append(values)
  return axes, result

if __name__ == "__main__":
  dna_map = {1:31.6770375,
             2:6.08434,
             3:3.56428,
             4:16.8736,
             5:10.9066141667,
             6:0,
             7:0.02131585,
             8:3.93474,
             9:0.0,
            10:7.153605,
            11:4.12787116667, 
            12:0.0,
            13:0.0,
            14:2.94159,
            15:0.0170527,
            16:10.12169,
            17:0.0,
            18:0.0298422,
            19:3.8547765,
            20:0.0,
            21:15.0149,
            22:3.009805,
            23:0.0,
            24:1.0999,
            25:33.2442}
  prot_map = {1:4.7865965 ,
              2:34.029 ,
              3:1.60982222222 ,
              4:1.73948557895 ,
              5:20.8108789474 , 
              6:55.213675 ,
              7:10.0935666667 ,
              8:26.9460684211 ,
              9:45.4343 ,
              10:0.0 ,
              11:2.39874 ,
              12:28.585325 ,
              13:2.26155 ,
              14:5.44601 ,
              15:44.8111 ,
              16:31.0536 ,
              17:20.0862 ,
              18:58.90245 ,
              19:70.702 ,
              20:48.1078 ,
              21:4.68939 ,
              22:49.59595 ,
              23:2.73548 ,
              24:26.56905 ,
              25:32.9968 ,
              26:43.302 ,
              27:12.081 ,
              28:32.2662 ,
              29:21.69775 ,
              30:75.2253 ,  
              31:8.97137 ,
              32:30.2361 ,
              33:6.634985 ,
              34:19.21515 ,
              35:32.5688 ,
              36:16.723 ,
              37:27.3692 ,
              38:25.1272 ,
              39:29.6739 ,
              40:43.9464 }
  prot_res = {"CYS":128.3296,
             "ASP":137.080194,
             "SER":105.492819,
             "ASN":145.657907,
             "GLN":169.926572,
             "LYS":185.888376,
             "ILE":166.241966,
             "PRO":135.240152,
             "THR":136.265712,
             "ALA":98.977676,
             "GLY":68.34639,
             "HIS":180.073895,
             "LEU":170.801411,
             "ARG":212.080324,
             "TRP":241.290773,
             "VAL":136.25585,
              "GLU":161.727587,
            "TYR":216.70009,
            "MET":189.668219}

  nuc_res = {"A":163.5912717,
             "C":162.86302775,
             "T":173.03036075,
             "G":164.253836,
             "DA":163.5912717,
             "DC":162.86302775,
             "DT":173.03036075,
             "DG":164.253836}
  
  axes, result = GetTable(sys.argv[1])
  a = "n"
  if "<-" in axes[0][0]:
    a = "l"
    t = axes[0][0].split("<-")[0].strip()
    if t  == "DNA":
      d = dna_map
    elif t == "PROTEIN":
      d = prot_map
  elif "->" in axes[0][0]:
    a ="r"
    t = axes[0][0].split("->")[1].strip()
    if t  == "DNA":
      d = dna_map
    elif t == "PROTEIN":
      d = prot_map
  else:
    print("wut")
    
  for j in range(len(axes[1])):
    for i in range(len(axes[0])-1):
      Ytag = axes[1][j].split("/")
      Xtag = axes[0][i+1].split("/")
      if a == "l":
        div = d[(int(Xtag[-1]))]
      elif a == "r":
        div = d[(int(Ytag[-1]))]
      if div == 0:
        if result[j][i] > 0:
          result[j][i] = 1.0
        else:
          result[j][i] = 0.0
      else:
        result[j][i] /= div
  outf = open(sys.argv[1]+".rel.tsv","w")
  outf.write("\t".join(axes[0])+"\n")
  for j in range(len(result)):
    outf.write(axes[1][j]+"\t"+"\t".join([str(item) for item in result[j]])+"\n")
        
