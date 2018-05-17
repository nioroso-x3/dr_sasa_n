
################################################################################
'''
http://pymolwiki.org/index.php/cgo_arrow

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

from pymol import cmd, cgo, CmdException
import numpy as np

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
       color='blue red', name=''):
  '''
DESCRIPTION

  Create a CGO arrow between two picked atoms.

ARGUMENTS

  atom1 = string: single atom selection or list of 3 floats {default: pk1}

  atom2 = string: single atom selection or list of 3 floats {default: pk2}

  radius = float: arrow radius {default: 0.5}

  gap = float: gap between arrow tips and the two atoms {default: 0.0}

  hlength = float: length of head

  hradius = float: radius of head

  color = string: one or two color names {default: blue red}

  name = string: name of CGO object
  '''
  from chempy import cpv

  radius, gap = float(radius), float(gap)
  hlength, hradius = float(hlength), float(hradius)

  try:
    color1, color2 = color.split()
  except:
    color1 = color2 = color
  color1 = list(cmd.get_color_tuple(color1))
  color2 = list(cmd.get_color_tuple(color2))

  def get_coord(v):
    if not isinstance(v, str):
      return v
    if v.startswith('['):
      return cmd.safe_list_eval(v)
    return cmd.get_atom_coords(v)

  xyz1 = get_coord(atom1)
  xyz2 = get_coord(atom2)
  normal = cpv.normalize(cpv.sub(xyz1, xyz2))

  if hlength < 0:
    hlength = radius * 3.0
  if hradius < 0:
    hradius = hlength * 0.6

  if gap:
    diff = cpv.scale(normal, gap)
    xyz1 = cpv.sub(xyz1, diff)
    xyz2 = cpv.add(xyz2, diff)

  xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

  obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
     [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
     [1.0, 0.0]
  #obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
     #[cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
     #[1.0, 0.0]

  if not name:
    name = cmd.get_unused_name('arrow')

  cmd.load_cgo(obj, name)
stored.cgo_arrow = cgo_arrow
cmd.extend('DRSASA_cgo_arrow', cgo_arrow)
################################################################################
def ReadDRSASAIntTable(fname):
  lines = open(fname).readlines()
  if("<-" in lines[0]):
    direction = "col"
  elif ("->" in lines[0]):
    direction = "row"
  else:
    raise "incompatible file?"
  tag = lines[0].split("\t")[0]
  columns = lines[0].strip().split("\t")[1:]
  sumv = []
  #print columns
  columns = [ item.split("/") for item in columns]
  rows = []
  matrix = []
  for line in lines[1:]:
    line = line.strip()
    tokens = [item.strip() for item in line.split("\t")]
    if len(tokens) < 3:
      continue
    rlbl = tokens[0].split("/")
    #print rlbl
    rows.append(rlbl)
    matrix.append([float(value) for value in tokens[1:]])
  if direction == 0:
    for row in range(len(matrix)):
      sumv.append(sum(matrix[row]))
  elif direction == 1:
    sumv = [0.0 for i in range(len(matrix[0]))]
    for row in range(len(matrix)):
      for col in range(len(matrix[row])):
        sumv[col]+=matrix[row][col]
    
  return columns,rows,matrix,direction,sumv,tag

def GenSeleResn(resn,chain,resi):
  if resn in ["DG","DA","A","G"]:
    atom = "(name N9) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  elif resn in ["DC","DT","C","T","U"]:
    atom = "(name N1) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  else:
    atom = "(name CA) and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  
  resn = "(resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  return atom,resn

def GenSeleAtom(name,resn,chain,resi):
  atom = "(name "+name+") and (resn "+resn+") and (chain "+chain+") and (resi "+str(resi).strip()+")"
  return atom

def DrawDRSASA_ByRes(fname):
  cols,rows,matrix,direc,sumv = ReadDRSASAIntTable(fname)

  mtrx = np.array(matrix).flatten()
  mtrx = mtrx[mtrx != 0]
  mvalue = np.max(mtrx)
  minvalue = np.min(mtrx)
  mv075 = mvalue*0.75
  median = np.median(mtrx)
  avg = np.average(mtrx)
  cutoff = mvalue*0.05
  cmd.color("grey","all")
  print median,avg,cutoff,minvalue
  Rl = []
  Gl = []
  Bl = []
  cmd.set('suspend_updates', 'on')
  view = cmd.get_view()
  for row in range(len(matrix)):
    for col in range(len(matrix[row])):
      v = float(matrix[row][col])
      if v < cutoff:
        continue
      R = 1.8*v/mvalue

      rid = [str(item) for item in rows[row]]
      cid = [str(item) for item in cols[col]]
      sele1,resn1 = GenSeleResn(rid[0],rid[1],rid[2])
      sele2,resn2 = GenSeleResn(cid[0],cid[1],cid[2])
      print(sele1,sele2)
      cmd.select("atom1",sele1)
      cmd.select("atom2",sele2)
      d = cmd.get_distance("atom1","atom2")
      #name = "/".join(rid)+"_"+"/".join(cid)
      name = "a_"+str(row)+"_"+str(col)
      if direc == "row":
        resn_s = resn2
      elif direc == "col":
        resn_s = resn1
      if (v >= mv075):
        color = "red red"
        Rl.append(resn_s)
        cmd.color("red",resn_s)
      elif ( v > median):
        color = "green green"
        Gl.append(resn_s)
        if not resn_s in Rl:
          cmd.color("green",resn_s)
      else:
        color = "blue blue"
        if (not resn_s in Rl) and (not resn_s in Gl):
          cmd.color("blue",resn_s)
      
      if direc == "row":
        cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.1 , R, color, name)
      elif direc == "col":
        cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.1 , R, color, name)
  cmd.set_view(view)
  cmd.set('suspend_updates', 'off')
cmd.extend("DRSASA_DrawInter_ByRes",DrawDRSASA_ByRes)

def DrawDRSASA_ByAtom(fname,mode='0'):
  mode = int(mode)
  cols,rows,matrix,direc,sumv = ReadDRSASAIntTable(fname)
  mtrx = np.array(matrix).flatten()
  mtrx = mtrx[mtrx != 0]
  mvalue = np.max(mtrx)
  minvalue = np.min(mtrx)
  mv075 = mvalue*0.75
  median = np.median(mtrx)
  avg = np.average(mtrx)
  cutoff = mvalue*0.05
  cmd.color("grey","all")
  print median,avg,cutoff,minvalue
  Rl = []
  Gl = []
  Bl = []
  cmd.set('suspend_updates', 'on')
  view = cmd.get_view()
  for row in range(len(matrix)):
    for col in range(len(matrix[row])):
      v = float(matrix[row][col])
      if v < cutoff:
        continue
      R = 1.8*v/mvalue
      nR = v/mvalue
      colorn = "drsasa_"+str(row)+"_"+str(col)
      cmd.set_color(colorn,(int(255*nR),0, 255 - int(255*nR)))
      rid = [str(item) for item in rows[row]]
      cid = [str(item) for item in cols[col]]
      print (rid,cid)
      sele1 = GenSeleAtom(rid[0],rid[1],rid[2],rid[3])
      sele2 = GenSeleAtom(cid[0],cid[1],cid[2],cid[3])
      print(sele1,sele2)
      cmd.select("atom1",sele1)
      cmd.select("atom2",sele2)
      d = cmd.get_distance("atom1","atom2")
      if direc == "row":
        atom_s = sele2
      elif direc == "col":
        atom_s = sele1
      #name = "/".join(rid)+"_"+"/".join(cid)
      name = "a_"+str(row)+"_"+str(col)
      
      if mode == 0:
        if (v >= mv075):
          color = colorn+" "+colorn
          Rl.append(atom_s)
          cmd.color("red",atom_s)
        elif ( v > median):
          color = "green green"
          Gl.append(atom_s)
          if not atom_s in Rl:
            cmd.color("green",atom_s)
        else:
          color = "blue blue"
          if (not atom_s in Rl) and (not atom_s in Gl):
            cmd.color("blue",atom_s)
      elif mode == 1:
        color = colorn+" "+colorn
      if direc == "row":
        if mode == 0:
          cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.1 , R, color, name)
        elif mode == 1:
          cgo_arrow("atom2","atom1",0.1, d*0.25 , d*0.1 , 0.2, color, name)
      elif direc == "col":
        if mode == 0:
          cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.1 , R, color, name)
        elif mode == 1:
          cgo_arrow("atom1","atom2",0.1, d*0.25 , d*0.1 , 0.2, color, name)
  cmd.set_view(view)
  cmd.set('suspend_updates', 'off')
cmd.extend("DRSASA_DrawInter_ByAtom",DrawDRSASA_ByAtom)
