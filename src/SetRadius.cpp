#include "stdafx.h"
#include "atom_struct.h"
#include "SetRadius.h"


VDWcontainer::VDWcontainer(string vdwfile){
  if(vdwfile != ""){
    ifstream vdw(vdwfile,ifstream::in);
    string current_res;
    map<string, map<string,float> > newmap;  
    while(vdw.good()){
      string line;
      getline(vdw,line);
      if(line == "#vdw_end"){
        break;
      }
      if(line.size() > 0 && line[0] != '#'){
        istringstream buffer(line);
        vector<string> tokens((istream_iterator<string>(buffer)),istream_iterator<string>());
        if(tokens[0] == "RESIDUE"){
           current_res = tokens[2];
           if(newmap.count(current_res)==0){
             newmap[current_res] = map<string,float>();
           }
        }
        if(tokens[0] == "ATOM"){
          string name(line, 5, 5);
          name.erase( std::remove_if( name.begin(), name.end(), ::isspace ), name.end() );
          newmap[current_res][name] = stod(tokens[tokens.size()-2]);
        }
      }
    }
    vdw.close();
    Map = newmap;
  }
}

VDWcontainer::VDWcontainer(){
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
      a.ACTIVE=false;
    }
    if (a.ELEMENT == "C"){
      a.VDW = 1.8;
      set_to = "C";
    }
    if (a.ELEMENT == "N"){
      a.VDW = 1.6;
      set_to = "N";
    }
    if (a.ELEMENT == "O"){
      a.VDW = 1.4;
      set_to = "O";
    }
    if (a.ELEMENT == "P"){
      a.VDW = 1.9;
      set_to = "P";
    }
    if (a.ELEMENT == "S"){
      a.VDW = 1.85;
      set_to = "S";
    }
    if (a.ELEMENT == "SE"){
      a.VDW = 1.9;
      set_to = "SE";
    }
    if (a.ELEMENT == "I"){
      a.VDW = 2.094;
      set_to = "SE";
    }
    if (a.ELEMENT == "F"){
      a.VDW = 1.560;
      set_to = "SE";
    }
    if (a.ELEMENT == "Br"){
      a.VDW = 1.978;
      set_to = "SE";
    }
    if (a.ELEMENT == "Cl"){
      a.VDW = 1.735;
      set_to = "SE";
    }
    if (set_to == "NONE"){
      a.VDW = 0;
      a.ACTIVE = false;
    }};

  for (uint32 i = 0;i<atoms.size();++i){
    auto& atom = atoms[i];
    if(atom.DTYPE == "MOL2"){
      string ch_name = "";
      if(sybyl2chimera.find(atom.CHARGE) != sybyl2chimera.end()) ch_name = sybyl2chimera[atom.CHARGE];
      else{
         
        cerr << "MOL2_UNKNOWN_NAME "  <<atom.NAME <<"|" <<atom.RESN << "|" <<atom.RESI<< "|" << atom.CHAIN << "|" << atom.CHARGE<< "\n";
        continue;
      }
      float vdw = u_vdw_radii[ch_name];
      atom.VDW = vdw;
      atom.ACTIVE = true;
    }
    else{
      if(Map.count(atom.RESN) > 0){
        auto& names = Map.at(atom.RESN);
        if(names.count(atom.NAME) > 0){
            atom.VDW = names.at(atom.NAME);
        }
        else{
          cerr << "UNKNOWN_VDW_NAME";
          set_unknown(atom);
          if(atom.VDW != 0) cerr << "_FIXED";
          cerr << "\t"<<atom.STRUCTURE << "\t" << atom.NAME << "\t" << atom.RESN <<
           "\t" << atom.RESI << "\t" <<atom.CHAIN << atom.iCODE<< "\t" << atom.ELEMENT<< "\n";
        }
      }
      else {
        string resn = CorrectRESN(atom.RESN);
        if(Map.count(resn) > 0 && Map.at(resn).count(atom.NAME) > 0){
          atom.VDW = Map.at(resn).at(atom.NAME);
        }else {
          cerr << "UNKNOWN_VDW";
          set_unknown(atom);
          if (atom.VDW != 0)cerr << "_FIXED";
          cerr <<"\t"<<atom.STRUCTURE << "\t" << atom.NAME << "\t" << atom.RESN <<
                  "\t" << atom.RESI << "\t" <<atom.CHAIN << atom.iCODE<< "\t" << atom.ELEMENT<< "\n";
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
