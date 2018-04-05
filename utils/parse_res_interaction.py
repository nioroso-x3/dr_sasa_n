#!/usr/bin/python3
import sys
from os.path import *

#DNA backbone polar atoms
PDI_Viz_dna_polar_bb = {'DA':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DC':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DG':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'DT':['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'U': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'A': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'C': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'G': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\''],
                               'T': ['P','OP1','OP2','OP3','O3\'','O4\'','O5\'']}
#DNA base polar atoms
PDI_Viz_dna_polar_ba =  {'DA':['N1','N3','N6','N7'],
                          'DC':['N3','O2','N4','N1'],
                          'DG':['O6','N1','N2','N3','N7'],
                          'DT':['O2','O4','N3','N1'],
                          'U': ['O2','O4','N3','N1'],
                          'A':['N1','N7','N3','N6'],
                          'C':['N3','O2','N4','N1'],
                          'G':['O6','N3','N7','N2','N1'],
                          'T':['O2','O4','N3','N1']}
#DNA backbone apolar
PDI_Viz_dna_apolar_bb = {'DA':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DC':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DG':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'DT':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'U': ['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'A':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'C':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'G':['C1\'','C2\'','C3\'','C4\'','C5\''],
                               'T':['C1\'','C2\'','C3\'','C4\'','C5\'']}
#DNA base apolar
PDI_Viz_dna_apolar_ba =   {'DA':['C2','C4','C5','C6','C8','N9'],
                                  'DC':['C2','C4','C5','C6'],
                                  'DG':['C2','C4','C5','C6','C8','N9'],
                                  'DT':['C2','C4','C5','C6','C7'],
                                  'U': ['C2','C4','C5','C6',],
                                  'A':['C2','C4','C5','C6','C8','N9'],
                                  'C':['C2','C4','C5','C6'],
                                  'G':['C2','C4','C5','C6','C8','N9'],
                                  'T':['C2','C4','C5','C6','C7'],}

#DNA base H bond acceptors
PDI_Viz_dna_Hbond_A = {'DA':['N1','N7','N3'],
                              'DC':['N3','O2'],
                              'DG':['O6','N3','N7'],
                              'DT':['O2','O4'],
                              'U': ['O2','O4'],
                              'A':['N1','N7','N3'],
                              'C':['N3','O2'],
                              'G':['O6','N3','N7'],
                              'T':['O2','O4']}
#DNA base H bond donors
PDI_Viz_dna_Hbond_D = {'DA':['N6'],
                          'DC':['N4'],
                          'DG':['N2','N1'],
                          'DT':['N3'],
                          'U' :['N3'],
                          'A':['N6'],
                          'C':['N4'],
                          'G':['N2','N1'],
                          'T':['N3']}
#DNA thymine methyl
PDI_Viz_dna_Tmet = {'DA':[],
                       'DC':[],
                       'DG':[],
                       'DT':['C7'],
                        'U':[],
                        'A':[],
                       'C':[],
                       'G':[],
                       'T':['C7']}
#DNA base atoms
PDI_Viz_dna_Hn = {'DA':['C2','C8'],
                     'DC':['C5','C6'],
                     'DG':['C8'],
                     'DT':['C6'],
                     'U': ['C5','C6'],
                     'A':['C2','C8'],
                     'C':['C5','C6'],
                     'G':['C8'],
                     'T':['C6']}
#Protein side chain acceptors
PDI_Viz_prot_Hbond_A = {'GLY':[],
                 'ALA':[],
                 'VAL':[],
                 'LEU':[],
                 'ILE':[],
                 'MET':['SD'],
                 'PRO':[],
                 'PHE':[],
                 'TRP':[],
                 'SER':['OG'],
                 'THR':['OG1'],
                 'GLN':['OE1'],
                 'LYS':[],
                 'TYR':['OH'],
                 'ASN':['OD1'],
                 'CYS':[],
                 'GLU':['OE1','OE2'],
                 'ASP':['OD1','OD2'],
                 'ARG':[],
                 'HIS':['ND1'],
                 'HID':['NE2'],
                 'HIE':['ND1'],
                 'HIP':[]}
#Protein side chain donors
PDI_Viz_prot_Hbond_D = {'GLY':[],
                 'ALA':[],
                 'VAL':[],
                 'LEU':[],
                 'ILE':[],
                 'MET':[],
                 'PRO':[],
                 'PHE':[],
                 'TRP':['NE1'],
                 'SER':['OG'],
                 'THR':['OG1'],
                 'GLN':['NE2'],
                 'LYS':['NZ'],
                 'TYR':['OH'],
                 'ASN':['ND2'],
                 'CYS':['SG'],
                 'GLU':[],
                 'ASP':[],
                 'ARG':['NH1','NH2','NE'],
                 'HIS':['NE2'],
                 'HID':['ND1'],
                 'HIE':['NE2'],
                 'HIP':['NE2','ND1']}

PDI_Viz_prot_polar = {'GLY':['N','O','OXT'],
                             'ALA':['N','O','OXT'],
                             'VAL':['N','O','OXT'],
                             'LEU':['N','O','OXT'],
                             'ILE':['N','O','OXT'],
                             'MET':['SD','N','O','OXT'],
                             'PRO':['N','O','OXT'],
                             'PHE':['N','O','OXT'],
                             'TRP':['NE1','N','O','OXT'],
                             'SER':['OG','N','O','OXT'],
                             'THR':['OG1','N','O','OXT'],
                             'GLN':['OE1','NE2','N','O','OXT'],
                             'LYS':['NZ','N','O','OXT'],
                             'TYR':['OH','N','O','OXT'],
                             'ASN':['OD1','ND2','N','O','OXT'],
                             'CYS':['SG','N','O','OXT'],
                             'GLU':['OE1','OE2','N','O','OXT'],
                             'ASP':['OD1','OD2','N','O','OXT'],
                             'ARG':['NH1','NH2','NE','N','O','OXT'],
                             'HIS':['ND1','NE2','N','O','OXT'],
                             'HID':['ND1','NE2','N','O','OXT'],
                             'HIE':['ND1','NE2','N','O','OXT'],
                             'HIP':['ND1','NE2','N','O','OXT']}

PDI_Viz_prot_apolar = {'GLY':['CA','C',],
                 'ALA':['CA','C','CB'],
                 'VAL':['CA','C','CB','CG1','CG2'],
                 'LEU':['CA','C','CB','CG','CD1','CD2'],
                 'ILE':['CA','C','CB','CG1','CG2','CD1'],
                 'MET':['CA','C','CB','CG','CE'],
                 'PRO':['CA','C','CB','CG','CD'],
                 'PHE':['CA','C','CB','CG','CD1','CD2','CE1','CE2','CZ'],
                 'TRP':['CA','C','CB','CG','CD1','CD2','CE2','CZ2','CH2','CZ3','CE3'],
                 'SER':['CA','C','CB'],
                 'THR':['CA','C','CB','CG2'],
                 'GLN':['CA','C','CB','CG','CD'],
                 'LYS':['CA','C','CB','CG','CD','CE'],
                 'TYR':['CA','C','CB','CG','CD1','CD2','CE1','CE2','CZ'],
                 'ASN':['CA','C','CB','CG'],
                 'CYS':['CA','C','CB'],
                 'GLU':['CA','C','CB','CG','CD'],
                 'ASP':['CA','C','CB','CG'],
                 'ARG':['CA','C','CB','CG','CD','CZ'],
                 'HIS':['CA','C','CB','CG','CD2','CE1'],
                 'HID':['CA','C','CB','CG','CD2','CE1'],
                 'HIE':['CA','C','CB','CG','CD2','CE1'],
                 'HIP':['CA','C','CB','CG','CD2','CE1']}

PDI_Viz_dna_majorgroove = {'DA':['C5', 'C6', 'C8', 'N6', 'N7'],
                                  'DC':['C4', 'C5', 'C6', 'N4'],
                                  'DG':['C5', 'C6', 'C8', 'N7', 'O6'],
                                  'DT':['C4', 'C5', 'C6', 'C7', 'O4'],
                                  'U':['C4', 'C5', 'C6', 'O4'],
                                  'A':['C5', 'C6', 'C8', 'N6', 'N7'],
                                  'C':['C4', 'C5', 'C6', 'N4'],
                                  'G':['C5', 'C6', 'C8', 'N7', 'O6'],
                                  'T':['C4', 'C5', 'C6', 'C7', 'O4']} # RNA
PDI_Viz_dna_minorgroove = {'DA':['C2', 'C4', 'N1', 'N3'],
                                  'DC':['C2', 'O2'],
                                  'DG':['C2', 'C4', 'N2', 'N3'],
                                  'DT':['C2', 'N3', 'O2'],
                                   'U':['C2', 'N3', 'O2'],
                                  'A':['C2', 'C4', 'N1', 'N3'],
                                  'C':['C2', 'O2'],
                                  'G':['C2', 'C4', 'N2', 'N3'],
                                  'T':['C2', 'N3', 'O2']} # RNA
PDI_Viz_dna_ambi = {'DA':['N9'],
                           'DC':['N1', 'N3'],
                           'DG':['N1', 'N9'],
                           'DT':['N1'],
                            'U':['N1'],
                           'A':['N9'],
                           'C':['N1', 'N3'],
                           'G':['N1', 'N9'],
                           'T':['N1']} # RNA

#protein backbone and bases
PDI_Viz_prot_backbone = ['N', 'CA', 'C', 'O', 'OXT']
PDI_Viz_prot_base = ['CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CZ', 'CZ2',
                                          'CZ3', 'CH2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NH1', 'NH2','OD1', 'OD2', 'OG',
                                          'OG1', 'OG2', 'OE1', 'OE2', 'OH', 'SD', 'SG']
#dna backbone and bases
PDI_Viz_dna_backbone = ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'','O2\'', 'O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'OP3', 'P']
PDI_Viz_dna_base = ['C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6']

if __name__ == "__main__":
  lines = [item.strip().split("\t") for item in open(sys.argv[1]).readlines()]
  result = {}
  rowDict = {}
  for row in range(len(lines)):
    for col in range(len(lines[row])):
      if row == 0:
        if col > 0:
          rowDict[col] = lines[row][col].split("/")
      else:
        if col == 0:
          atomid = lines[row][col].split("/")
          resid = tuple(atomid[1:4])
          if not resid in result:
            result[resid] = {"SASA": 0.0,"TOTAL_BSA": 0.0, "BACKBONE_BSA": 0.0, "BASE_BSA":0.0, "MAJORG_BSA":0.0,"MINORG_BSA":0.0,"AMBIGUOUSG_BSA":0.0}
        else:
          Oatomid = rowDict[col]
          bsa = float(lines[row][col])
          if bsa == 0.0:
            continue
          result[resid]["TOTAL_BSA"] += bsa
          #base or backbone
          if Oatomid[0] in PDI_Viz_dna_backbone:
            result[resid]["BACKBONE_BSA"]+=bsa
          elif Oatomid[0] in PDI_Viz_dna_base:
            result[resid]["BASE_BSA"]+=bsa
          if Oatomid[0] in PDI_Viz_dna_majorgroove[Oatomid[1]]:
            result[resid]["MAJORG_BSA"]+=bsa
          elif Oatomid[0] in PDI_Viz_dna_minorgroove[Oatomid[1]]:
            result[resid]["MINORG_BSA"]+=bsa
          elif Oatomid[0] in PDI_Viz_dna_ambi[Oatomid[1]]:
            result[resid]["AMBIGUOUSG_BSA"] += bsa

  lines = [item.strip().split("\t") for item in open(sys.argv[2]).readlines() if item.split("\t")[0] == "ATOM"]
  for line in lines:
    idt = tuple(line[3:6])
    free =float(line[10])
    try:
      result[idt]["SASA"] += free
    except:
      pass
  keys = result.keys()
  keys = sorted(keys,key= lambda x : int(x[2]))
  keys = sorted(keys,key = lambda x : x[1])
  bname = basename(sys.argv[1])
  if len(bname)>4:
    obj = bname[0:8]
  else:
    obj = bname[0:4]
  cols = ["SASA","TOTAL_BSA", "BACKBONE_BSA", "BASE_BSA", "MAJORG_BSA","MINORG_BSA","AMBIGUOUSG_BSA"]
  #print("RESN\tCHAIN\tRESI\t"+"\t".join(cols)+"\tSTRUCT")
  for res in keys:
    if res[0] in ["A","T","G","C","DA","DT","DG","DC"]:
      continue
    if result[res]["TOTAL_BSA"] == 0.0:
      continue
    values = [str(result[res][col]) for col in cols]
    values.append(obj)
    print("\t".join(res)+"\t"+"\t".join(values))


