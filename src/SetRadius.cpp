#include "stdafx.h"
#include "atom_struct.h"
#include "SearchFunctions.h"
#include "SetRadius.h"


VDWcontainer::VDWcontainer(string vdwfile){
  vector<vector<string>> fpdb = ReadTabSeparatedString(pdbradii);
  vector<vector<string>> fmol2 = ReadTabSeparatedString(mol2radii);

  load_data(fpdb);
  load_data(fmol2);
  if(vdwfile != ""){
    cerr << "#Loading external VdW radius file.\n";
    try{
      vector<vector<string>> f = ReadTabSeparatedFile(vdwfile);
      load_data(f);
    } 
    catch(exception e){
      cerr << "Invalid VDW file.\n";
    }
  }
}

VDWcontainer::VDWcontainer(){
}

void
VDWcontainer::load_data(vector<vector<string>> f){
  if(f[0][0] == "PDB_RADIUS_FILE"){
    map<string, map<string,float> > newmap;
    map<string, map<string,int> > newmap_pol;
    string cur_res = "";
    for (auto& line : f){

      if (line.size() != 4) continue;
      if (line[0] == "RESIDUE"){
        cur_res = line[2];
        continue; 
      }
      if (line[0] == "ATOM"){
        string name = line[1];
        float vdw = stof(line[2]);
        int   pol = stoi(line[3]);
        if (newmap.count(cur_res) == 0) {
           map<string,float> tmp;
           newmap[cur_res] = tmp;
        }
        if (newmap_pol.count(cur_res) == 0) {
           map<string,int> tmp;
           newmap_pol[cur_res] = tmp;
        }
        newmap[cur_res][name] = vdw;
        newmap_pol[cur_res][name] = pol;
      }
    }
    Map = newmap;
    MapPolarity = newmap_pol;
  }
  else if(f[0][0] == "MOL2_RADIUS_FILE"){
    map<string,float> newmap;
    map<string,int> newmap_pol;
    for (auto& line : f){
      if(line.size() != 4) continue;
      if(line[0] == "SYBYL"){
        newmap[line[1]] = stof(line[2]);
        newmap_pol[line[1]] = stoi(line[3]);
      }
    }        
    u_vdw_radii = newmap;
    u_vdw_radii_polarity = newmap_pol;
  }
  else throw;    
}

string
VDWcontainer::CorrectRESN(string resn){
  extern map<string,string> nucresn;
  try{
    return nucresn.at(resn);
  }
  catch(exception e){
    return resn;
  }
}

void VDWcontainer::SetRadius(vector<atom_struct>& atoms,float probe){

  auto set_unknown = [](atom_struct& a){
      string set_to = "NONE";
    if (a.ELEMENT == "H"){
      a.VDW = 1.0;
      set_to = "H";
      a.POLAR = 0;
      a.ACTIVE=false;
    }
    if (a.ELEMENT == "C"){
      a.VDW = 1.8;
      a.POLAR = 0;
      set_to = "C";
    }
    if (a.ELEMENT == "N"){
      a.VDW = 1.6;
      a.POLAR = 1;
      set_to = "N";
    }
    if (a.ELEMENT == "O"){
      a.VDW = 1.4;
      a.POLAR = 1;
      set_to = "O";
    }
    if (a.ELEMENT == "P"){
      a.VDW = 1.9;
      a.POLAR = 1;
      set_to = "P";
    }
    if (a.ELEMENT == "S"){
      a.VDW = 1.85;
      a.POLAR = 0;
      set_to = "S";
    }
    if (a.ELEMENT == "SE"){
      a.VDW = 1.9;
      a.POLAR = 0;
      set_to = "SE";
    }
    if (a.ELEMENT == "I"){
      a.VDW = 2.094;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
      set_to = "I";
    }
    if (a.ELEMENT == "F"){
      a.VDW = 1.560;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "F";
    }
    if (a.ELEMENT == "Br" || a.ELEMENT == "BR"){
      a.VDW = 1.978;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Br";
    }
    if (a.ELEMENT == "Cl" || a.ELEMENT == "CL"){
      if (a.NAME == a.RESN && (a.RESN == "CL" || a.RESN == "Cl")) a.VDW = 1.81;
      else a.VDW=1.74;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Cl";
    }
    if (a.ELEMENT == "Al" || a.ELEMENT == "AL"){
      a.VDW = 0.54;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Al";
    }
    if (a.ELEMENT == "As" || a.ELEMENT == "AS"){
      a.VDW = 0.58;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "As";
    }
    if (a.ELEMENT == "Au" || a.ELEMENT == "AU"){
      a.VDW = 1.37;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Au";
    }
    if (a.ELEMENT == "Ba" || a.ELEMENT == "BA"){
      a.VDW = 1.35;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ba";
    }
    if (a.ELEMENT == "Be" || a.ELEMENT == "BE"){
      a.VDW = 0.45;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Be";
    }
    if (a.ELEMENT == "Bi" || a.ELEMENT == "BI"){
      a.VDW = 1.03;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Bi";
    }
    if (a.ELEMENT == "Ca" || a.ELEMENT == "CA"){
      a.VDW = 1.00;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ca";
    }
    if (a.ELEMENT == "Cd" || a.ELEMENT == "CD"){
      a.VDW = 0.95;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Cd";
    }
    if (a.ELEMENT == "Cr" || a.ELEMENT == "CR"){
      a.VDW = 0.73;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Cr";
    }
    if (a.ELEMENT == "Cs" || a.ELEMENT == "CS"){
      a.VDW = 1.67;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Cs";
    }
    if (a.ELEMENT == "Cu" || a.ELEMENT == "CU"){
      a.VDW = 0.73;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Cu";
    }
    if (a.ELEMENT == "Fe" || a.ELEMENT == "FE"){
      a.VDW = 0.65;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Fe";
    }
    if (a.ELEMENT == "Ga" || a.ELEMENT == "GA"){
      a.VDW = 0.62;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ga";
    }
    if (a.ELEMENT == "Ge" || a.ELEMENT == "GE"){
      a.VDW = 0.73;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ge";
    }
    if (a.ELEMENT == "Hg" || a.ELEMENT == "HG"){
      a.VDW = 1.02;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Hg";
    }
    if (a.ELEMENT == "K" || a.ELEMENT == "K"){
      a.VDW = 1.38;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "K";
    }
    if (a.ELEMENT == "Li" || a.ELEMENT == "LI"){
      a.VDW = 0.76;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Li";
    }
    if (a.ELEMENT == "Mg" || a.ELEMENT == "MG"){
      a.VDW = 0.72;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Mg";
    }
    if (a.ELEMENT == "Mn" || a.ELEMENT == "MN"){
      a.VDW = 0.83;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Mn";
    }
    if (a.ELEMENT == "Mo" || a.ELEMENT == "MO"){
      a.VDW = 0.69;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Mo";
    }
    if (a.ELEMENT == "Na" || a.ELEMENT == "NA"){
      a.VDW = 1.02;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Na";
    }
    if (a.ELEMENT == "Ni" || a.ELEMENT == "NI"){
      a.VDW = 0.69;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ni";
    }
    if (a.ELEMENT == "Pb" || a.ELEMENT == "PB"){
      a.VDW = 1.19;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Pb";
    }
    if (a.ELEMENT == "Pd" || a.ELEMENT == "PD"){
      a.VDW = 0.86;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Pd";
    }
    if (a.ELEMENT == "Rb" || a.ELEMENT == "RB"){
      a.VDW = 1.52;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Rb";
    }
    if (a.ELEMENT == "Sb" || a.ELEMENT == "SB"){
      a.VDW = 0.76;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Sb";
    }
    if (a.ELEMENT == "Sc" || a.ELEMENT == "SC"){
      a.VDW = 0.75;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Sc";
    }
    if (a.ELEMENT == "Sn" || a.ELEMENT == "SN"){
      a.VDW = 0.69;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Sn";
    }
    if (a.ELEMENT == "Sr" || a.ELEMENT == "SR"){
      a.VDW = 1.18;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Sr";
    }
    if (a.ELEMENT == "Tc" || a.ELEMENT == "TC"){
      a.VDW = 0.65;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Tc";
    }
    if (a.ELEMENT == "Ti" || a.ELEMENT == "TI"){
      a.VDW = 0.86;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Ti";
    }
    if (a.ELEMENT == "V"){
      a.VDW = 0.79;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "V";
    }
    if (a.ELEMENT == "Zn" || a.ELEMENT == "ZN"){
      a.VDW = 0.74;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Zn";
    }
    if (a.ELEMENT == "Zr" || a.ELEMENT == "ZR"){
      a.VDW = 0.72;
      if(a.NAME == a.ELEMENT) a.POLAR=1;
      else a.POLAR =0;
set_to = "Zr";
    }

    if (set_to == "NONE"){
      a.VDW = 0;
      a.ACTIVE = false;
    }else{
      a.ACTIVE = true;
    }
  };

  for (uint32 i = 0;i<atoms.size();++i){
    auto& atom = atoms[i];
    if(atom.DTYPE == "MOL2"){
      string ch_name = "";
      float vdw;
      int pol;
      try{
        vdw = u_vdw_radii.at(atom.CHARGE);
        pol = u_vdw_radii_polarity.at(atom.CHARGE); 

        atom.VDW = vdw;
        atom.POLAR = pol;
        atom.ACTIVE = true;
      }
      catch(exception e){
        cerr << "MOL2_UNKNOWN_TYPE";
        set_unknown(atom);
        if (atom.VDW != 0.0) cerr << "_FIXED ";
        else cerr << " ";
        cerr  <<atom.NAME <<"|" <<atom.RESN << "|" <<atom.RESI<< "|" << atom.CHAIN << "|" << atom.CHARGE << "|" << atom.ELEMENT << "\n";
 
      }
    }
    else{
      if(Map.count(atom.RESN) > 0){
        auto& names = Map.at(atom.RESN);
        auto& namesP = MapPolarity.at(atom.RESN);
        if(names.count(atom.NAME) > 0){
            atom.VDW = names.at(atom.NAME);
            atom.POLAR = namesP.at(atom.NAME);
        }
        else{
          cerr << "UNKNOWN_VDW_NAME";
          set_unknown(atom);
          if(atom.VDW != 0) cerr << "_FIXED ";
          else cerr << " ";
          cerr <<atom.STRUCTURE << "|" << atom.NAME << "|" << atom.RESN <<
           "|" << atom.RESI << "|" <<atom.CHAIN << atom.iCODE<< "|" << atom.ELEMENT<< "\n";
        }
      }
      else {
        string resn = CorrectRESN(atom.RESN);
        if(Map.count(resn) > 0 && Map.at(resn).count(atom.NAME) > 0){
          atom.VDW = Map.at(resn).at(atom.NAME);
        }else {
          cerr << "UNKNOWN_VDW";
          set_unknown(atom);
          if (atom.VDW != 0)cerr << "_FIXED ";
          else cerr << " ";
          cerr <<atom.STRUCTURE << "|" << atom.NAME << "|" << atom.RESN <<
                  "|" << atom.RESI << "|" <<atom.CHAIN << atom.iCODE<< "|" << atom.ELEMENT<< "|\n";
        }
      }
    }
    if(atom.ACTIVE){
      atom.RADIUS = atom.VDW + probe;
      atom.RADIUS2 = atom.RADIUS * atom.RADIUS;
      atom.AREA = atom.RADIUS2*4.0*PI;
      //cout << atom.ID << "\t" << atom.RADIUS << "\t" << atom.RADIUS2 << "\n";
    }
  }
}

void VDWcontainer::SetShell(vector<atom_struct>& pdb){
  auto* p = Points.data();
  auto p_s3 = Points.size();
  auto p_s = p_s3 / 3;
  auto pdb_s = pdb.size();
  for (uint32 i = 0; i < pdb_s; ++i){
    auto& atom_i = pdb[i];
    if (!atom_i.ACTIVE) continue;
    atom_i.SHELL_BURIED.resize(p_s,false);
    atom_i.ACCESSIBLE_P = p_s;
    atom_i.SHELL.resize(p_s3);
    auto R_I = atom_i.RADIUS;
    auto* C = atom_i.COORDS.data();
    auto* s = atom_i.SHELL.data();
    for(uint32 j = 0; j < p_s3; ++j){
      s[j] = p[j] * R_I + C[j%3];
    }
  }
}


void VDWcontainer::GenPoints(string pfile){
/*  //golden spiral
    Result.resize(N * 3);
    uint32 N2 = N * 2;
    float inc = PI*(3.0 - sqrt(5.0));
    float off = 2.0 / (float)N;
  for (uint32 i = 0; i < N; i++){
      float y = i * off - 1.0 + (off / 2.0);
      float r = sqrt(1.0 - y*y);
      float phi = i*inc;
      Result[i] = cos(phi)*r;
      Result[i + N] = y;
      Result[i + N2] = sin(phi)*r;
    }
    return Result;
  }*/

//thomson electrostatic simulation http://thomson.phy.syr.edu/
  ifstream xyz_thomson(pfile);
  uint32 i = 0;
  string line;
  getline(xyz_thomson, line);
  auto N = stoi(line);
  Points.resize(N * 3);
  while (xyz_thomson.good()){
    getline(xyz_thomson, line);
    istringstream buffer(line);
    vector<string> tokens((istream_iterator<string>(buffer)), istream_iterator<string>());
    if (tokens.size() == 3){
      //row order, keeps vector access local [x + y*3].
      Points[i]     = stod(tokens[0]);
      Points[1 + i] = stod(tokens[1]);
      Points[2 + i] = stod(tokens[2]);
      i+=3;
    }
  }
  xyz_thomson.close();
}
void VDWcontainer::GenPoints(){
  Points = p15k;
}
void VDWcontainer::GenPoints(int select){
  if (select == 0)  Points = p15k;
  else Points = p2k;
}

VDWcontainer::~VDWcontainer(void)
{
}
