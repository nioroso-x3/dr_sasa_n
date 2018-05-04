#include "stdafx.h"
#include "atom_struct.h"
#include "SearchFunctions.h"
#include "PDBparser2.h"

string
GetBasename(string path){
  string result;
  stringstream strstream(path);
  string segment;
  vector<string> seglist;
#ifdef _WIN32
  while(std::getline(strstream, segment, '\\')){
#else
  while(std::getline(strstream, segment, '/')){
#endif
    seglist.push_back(segment);
  }
  result = seglist.back();
  return result;
}

map<string, map<string, int> > ANOLEA = {
    { "GLY", { { "N", 3 },
               { "CA", 2 },
               { "C", 4 },
               { "O", 5 },
               { "OXT", 28 } } },

    { "ALA", { { "N",  3 },
               { "CA", 1 },
               { "C", 4  },
               { "O", 5 },
               { "CB", 6 },
               { "OXT", 28 } } },

    { "VAL", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 7 },
               { "CG1", 6 },
               { "CG2", 6 },
               { "OXT", 28 } } },

    { "LEU", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 7 },
               { "CD1", 6 },
               { "CD2", 6 },
               { "OXT", 28 } } },

    { "ILE", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 7 },
               { "CG1", 8 },
               { "CG2", 6 },
               { "CD1", 6 },
               { "OXT", 28 } } },
    { "MET", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 29 },
               { "SD", 9 },
               { "CE", 30 },
               { "OXT", 28 } } },
   { "MSE", { { "N", 3 },
              { "CA", 1 },
              { "C", 4 },
              { "O", 5 },
              { "CB", 8 },
              { "CG", 29 },
              { "SE", 9 },
              { "CE", 30 },
              { "OXT", 28 } } },

    { "PRO", { { "N", 10 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 8 },
               { "CD", 32 },
               { "OXT", 28 } } },

    { "TRP", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 13 },
               { "CD1", 24 },
               { "CD2", 11 },
               { "NE1", 39 },
               { "CE2", 14 },
               { "CE3", 12 },
               { "CZ2", 12 },
               { "CZ3", 12 },
               { "CH2", 12 },
               { "OXT", 28 } } },

    { "SER", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 15 },
               { "OG", 16 },
               { "OXT", 28 } } },

    { "THR", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 17 },
               { "OG1", 16 },
               { "CG2", 6 },
               { "OXT", 28 } } },

    { "PHE", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 11 },
               { "CD1", 12 },
               { "CD2", 12 },
               { "CE1", 12 },
               { "CE2", 12 },
               { "CZ", 12 },
               { "OXT", 28 } } },

    { "ARG", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 8 },
               { "CD", 37 },
               { "NE", 36 },
               { "CZ", 21 },
               { "NH1", 22 },
               { "NH2", 22 },
               { "OXT", 28 } } },

    { "ASN", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 33 },
               { "OD1", 34 },
               { "ND2", 18 },
               { "OXT", 28 } } },

    { "ASP", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 27 },
               { "OD1", 28 },
               { "OD2", 28 },
               { "OXT", 28 } } },

    { "CYS", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 29 },
               { "SG", 19 },
               { "OXT", 28 } } },

    { "GLN", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 8 },
               { "CD", 33 },
               { "OE1", 34 },
               { "NE2", 18 },
               { "OXT", 28 } } },

    { "GLU", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 8 },
               { "CD", 27 },
               { "OE1", 28 },
               { "OE2", 28 },
               { "OXT", 28 } } },

    { "HIS", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 23 },
               { "ND1", 38 },
               { "CD2", 24 },
               { "CE1", 26 },
               { "NE2", 25 },
               { "OXT", 28 } } },

    { "LYS", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 8 },
               { "CD", 8 },
               { "CE", 35 },
               { "NZ", 20 },
               { "OXT", 28 } } },

    { "TYR", { { "N", 3 },
               { "CA", 1 },
               { "C", 4 },
               { "O", 5 },
               { "CB", 8 },
               { "CG", 11 },
               { "CD1", 12 },
               { "CD2", 12 },
               { "CE1", 12 },
               { "CE2", 12 },
               { "CZ", 31 },
               { "OH", 40 },
               { "OXT", 28 } } },


    //Modified Capriotti etal 2010 U -> DT
    { "DA", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 }, 
              { "C5'", 4 }, { "C4'", 5 },
              { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 }, { "C2'", 5 }, { "C1'", 7 }, 
              { "N9", 9 }, { "C8", 10 },
              { "N7", 11 }, { "C5", 12 }, { "C6", 15 }, { "N6", 16 }, { "N1", 11 }, 
              { "C2", 14 }, { "N3", 11 }, { "C4", 13 } } },

    { "DT", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 }, 
              { "C5'", 4 }, { "C4'", 5 },
              { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 }, { "C2'", 5 }, { "C1'", 7 }, 
              { "N1", 9 }, { "C2", 20 },
              { "O2", 19 }, { "N3", 23 }, { "C4", 18 }, { "O4", 19 }, { "C5", 24 }, 
              { "C7", 25 }, { "C6", 22 },{"C5M",25} } }, //types 24 and 25 added

    { "DC", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 }, 
              { "C5'", 4 }, { "C4'", 5 },
              { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 }, { "C2'", 5 }, { "C1'", 7 }, 
              { "N1", 9 }, { "C2", 20 },
              { "O2", 19 }, { "N3", 11 }, { "C4", 15 }, { "N4", 16 }, { "C5", 21 }, 
              { "C6", 22 } } },

    { "DG", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 }, 
              { "C5'", 4 }, { "C4'", 5 },
              { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 }, { "C2'", 5 }, { "C1'", 7 }, 
              { "N9", 9 }, { "C8", 10 },
              { "N7", 11 }, { "C5", 12 }, { "C6", 18 }, { "O6", 19 }, { "N1", 23 }, 
              { "C2", 17 }, { "N2", 16 },
              { "N3", 11 }, { "C4", 13 } } },
    { "A", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 },
             { "C5'", 4 }, { "C4'", 5 }, { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 },
             { "O2'", 6 }, { "C2'", 5 }, { "C1'", 7 }, { "N9", 9 }, { "C8", 10 },
             { "N7", 11 }, { "C5", 12 }, { "C6", 15 }, { "N6", 16 }, { "N1", 11 },
             { "C2", 14 }, { "N3", 11 }, { "C4", 13 } } },

    { "U", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 },
             { "C5'", 4 }, { "C4'", 5 }, { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 },
             { "O2'", 6 }, { "C2'", 5 }, { "C1'", 7 }, { "N1", 9 }, { "C2", 20 },
             { "O2", 19 }, { "N3", 23 }, { "C4", 18 }, { "O4", 19 }, { "C5", 21 },
             { "C6", 22 } } },

    { "C", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 },
             { "C5'", 4 }, { "C4'", 5 }, { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 },
             { "O2'", 6 }, { "C2'", 5 }, { "C1'", 7 }, { "N1", 9 }, { "C2", 20 },
             { "O2", 19 }, { "N3", 11 }, { "C4", 15 }, { "N4", 16 }, { "C5", 21 },
             { "C6", 22 } } },

    { "G", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 },
             { "C5'", 4 }, { "C4'", 5 }, { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 },
             { "O2'", 6 }, { "C2'", 5 }, { "C1'", 7 }, { "N9", 9 }, { "C8", 10 },
             { "N7", 11 }, { "C5", 12 }, { "C6", 18 }, { "O6", 19 }, { "N1", 23 },
             { "C2", 17 }, { "N2", 16 }, { "N3", 11 }, { "C4", 13 } } },
    { "T", { { "P", 2 }, { "OP1", 1 }, { "OP2", 1 }, { "OP3", 1 }, { "O5'", 3 }, 
             { "C5'", 4 }, { "C4'", 5 },
             { "O4'", 8 }, { "C3'", 5 }, { "O3'", 3 },{"O2'",6}, { "C2'", 5 }, 
             { "C1'", 7 }, { "N1", 9 }, { "C2", 20 },
             { "O2", 19 }, { "N3", 23 }, { "C4", 18 }, { "O4", 19 }, { "C5", 24 }, 
             { "C7", 25 }, { "C6", 22 },{"C5M",25} } },
{"PTS", {{"C03",0},  {"C04",0},  {"C05",0}, {"C08",0}, {"C09",0},  {"C10",0}, {"C11",0}, {"C12",0}, {"C13",0}, {"C14",0}, {"C15",0}, {"C17",0}, {"C18",0}, {"C19",0}, {"O02",0}, {"O06",0},{"O16",0} }}
  };


vector<vector<string> > //ordered strings for each category in the atom tags
ReadPDB_ATOM(string fname){ //filename
  vector<vector<string>> result;
  uint32 cols[15][2] = {{0,6},{6,5},{12,4},{16,1},{17,3},{21,1},{22,4},
                        {26,1},{30,8},{38,8},{46,8},{54,6},{60,6},{76,2},{78,2}};

  //   0  1    2      3    4     5    6     7 8 9 10     11    12   13     14
  //ATOM ID NAME ALTLOC RESN CHAIN RESI ICODE X Y  Z  OCCUP TFACT ELEM CHARGE

  ifstream pdb(fname);

  auto Trim = [](string& str){str.erase(std::remove_if(str.begin(),
                                                       str.end(),
                                                       ::isspace),
                                                       str.end());
  };
  while(pdb.good()){
    string line;
    vector<string> cline;
    getline(pdb,line);
    uint32 tag;
    string atomtag;
    bool notATOM = false;
    try{
      for (tag = 0; tag < 15; ++tag){
        string temp(line,cols[tag][0],cols[tag][1]);
        Trim(temp);
        replace(temp.begin(),temp.end(),'*','\'');
        
        if( tag == 0){ 
          if(temp != "ATOM" && temp != "HETATM"){
            notATOM = true;
            break;
          }
        }
        //cerr << "|"<< temp << "\t";
        cline.push_back(temp);
      }
    }
    catch(exception e){
      if (tag < 10){
        cerr << "Line \"" << tag << "\" has insufficient data.\n";
        continue;
      }
    }
    if(notATOM) continue;
    result.push_back(cline);
    //cout << result.back()[1] << "\n";
  }
  if(result.size() == 0) {
    cerr <<fname << " NO_VALID_ATOMS_IN_PDB\n";
    throw -1;
  }
  return result;
}

void
RemoveAltLocs(vector<atom_struct>& pdb){
/*
Steps:
1- Search for resi/chain pairs with altloc atoms
2- Search how many sets are available
3- Select set with highest occupancy/lowest letter
4- If there is pairs, select highest occupancy or if occupancy is equal select lowest position
*/
  vector<atom_struct> new_pdb;
  //chain,resi+icode to altloc set map for normal altloc in same chain
  map<string,map<pair<int,string>, set<string>>> chain_resi_set;
  bool Normal = false;
  vector<int> removed_pos;
  for (auto atom : pdb){
    if(atom.ALTLOC.empty()) continue;
    auto resi = make_pair(atom.RESI,atom.iCODE);
    chain_resi_set[atom.CHAIN][resi].insert(atom.ALTLOC);
  }
  for (auto ch_resi : chain_resi_set){
    string CHAIN = ch_resi.first;
    auto resi_map = ch_resi.second;
    for(auto resi_set : resi_map){
      int RESI = get<0>(resi_set.first);
      string iCODE = get<1>(resi_set.first);
      //variables to find out best altloc
      //altocs
      vector<string> alt_set(resi_set.second.begin(),resi_set.second.end());
      //do not step into per chain altlocs
      if (alt_set.size() < 2){
        Normal = true;
        continue;
      }
      //positions in the pdb of each altloc set
      vector<vector<uint32>> alt_set_atom_pos(alt_set.size());
      //avg occupancy of each altloc set
      vector<float> avg_oc(alt_set.size());

      for (uint32 i = 0; i < alt_set.size(); ++i){
        string altloc = alt_set[i];
        float oc_count = 0.0;
        float oc_sum = 0.0;
        for(uint32 j = 0; j < pdb.size(); ++j){
          if(pdb[j].CHAIN  != CHAIN  ||
             pdb[j].ALTLOC != altloc ||
             pdb[j].RESI   != RESI   ||
             pdb[j].iCODE  != iCODE) continue;
          oc_count += 1.0;
          oc_sum += pdb[j].OCCUPANCY;
          alt_set_atom_pos[i].push_back(j);
        }
        avg_oc[i] = oc_sum / oc_count;
      }
      float cur_avg_oc = avg_oc[0];
      string cur_altloc = alt_set[0];
      int cur_set = 0;
      for (uint32 i = 1; i < alt_set.size(); ++i){
        if((avg_oc[i] > cur_avg_oc)         ||
           (avg_oc[i] == cur_avg_oc &&
            alt_set[i].compare(cur_altloc) < 0 )){
          removed_pos.insert(removed_pos.end(),
                             alt_set_atom_pos[cur_set].begin(),
                             alt_set_atom_pos[cur_set].end());
          cur_set = i;
          cur_avg_oc = avg_oc[i];
          cur_altloc = alt_set[i];
        }
        else{
          removed_pos.insert(removed_pos.end(),
                             alt_set_atom_pos[i].begin(),
                             alt_set_atom_pos[i].end());
        }
      }
    }
  }
  for (uint32 i = 0; i < pdb.size(); ++i){
    auto& atom = pdb[i];
    if(removed_pos.end() != find(removed_pos.begin(),removed_pos.end(),i)) {
//       cerr << "NORMAL_ALT_RM\t" << atom.STRUCTURE << "\t" << atom.ID <<"\t"
//            << atom.NAME << "\t"  << atom.RESN << "\t" << atom.RESI << "\t" 
//            << atom.CHAIN << atom.iCODE << "\n";
      continue;
    }
    new_pdb.push_back(pdb[i]);
  }

  pdb = new_pdb;
  if (Normal){
    for (auto ch_resi : chain_resi_set){
      string CHAIN = ch_resi.first;
      auto resi_map = ch_resi.second;
        for(auto resi_set : resi_map){
          //int RESI = get<0>(resi_set.first);
          string iCODE = get<1>(resi_set.first);
          vector<string> alt_set(resi_set.second.begin(),resi_set.second.end());
//           if(alt_set.size() == 1){
//             cerr << "NON_STD_ALT_CHAIN\t" << pdb[0].STRUCTURE << "\t" <<  alt_set[0] << "\n";
//           }
        }
    }


    //starting chain altlocs

  }
  if(new_pdb.size() == 0) throw -1;
}

bool
IsRESNSolvent(string resn){
  const vector<string> solv = {"SO4", "HOH", "EDO", "SF4", "GOL", "EOH", "DOD",  
                               "CPT", "A71", "EHN", "PTN", "O", "MOE", "ACE", "TEA", 
                               "EDO", "CAC", "PEO", "PG4", "SPK", "DRI", "DAG", "TMR",  
                               "NBU", "TRS","RKT", "NAG"};
  return (find(solv.begin(),solv.end(),resn ) != solv.end());
}

bool
IsRESNIon(string resn){
  const vector<string> ions = {"NA","MG","K", "MN","ZN","BA","CO","3CO","AG","NCO","2HP","ACT",
                               "CS","CU","NI","PT","RB","SR","TL","BR","NRU","NH4"};
  return (find(ions.begin(),ions.end(),resn ) != ions.end());
}

void
GetAtomType(string typefname,  //type dictionary
            vector<atom_struct>&            pdb,            //vector of atoms
            bool                            KeepUnknown){
  map<string, map<string, int> > typemap;
  if (typefname.empty()){
    typemap = ANOLEA;
  }
  else{
    typemap = GetTypeMap(typefname);
  }  
  string c_chain = "NULL";
  vector<atom_struct> new_pdb;

  for (auto& atom : pdb){
    auto name = atom.NAME;
    auto resn = atom.RESN;
    auto chain = atom.CHAIN;
    //cout << atom.NAME << " " << atom.RESN << "\n"; 
    if(IsRESNSolvent(atom.RESN) ){ 
       cerr << "SOLVENT\t" << atom.STRUCTURE << "\t" << atom.ID <<"\t"
            << atom.NAME << "\t"  << atom.RESN << "\t" << atom.RESI << "\t" 
            << atom.CHAIN << atom.iCODE << "\n";
      continue;
    }
    if( IsRESNIon(atom.RESN)){ 
       cerr << "ION\t" << atom.STRUCTURE << "\t" << atom.ID <<"\t"
            << atom.NAME << "\t"  << atom.RESN << "\t" << atom.RESI << "\t" 
            << atom.CHAIN << atom.iCODE << "\n";
      continue;
    }
    if(!typemap.count(resn)) resn = FixRESN(resn);
    if(typemap.count(resn) && typemap.at(resn).count(name)){
      int type = typemap.at(resn).at(name);
      //entering new chain
      if(c_chain != chain && atom.MOL_TYPE == "PROTEIN"){
        c_chain = chain;
        if(type == 3 || type == 10){
          type = 20;
        }
      }
      if(c_chain != chain && atom.MOL_TYPE == "RNA" ){
        c_chain = chain;
        if(type == 3 && atom.NAME == "O3'"){
          type = 6;
        }
      }
      atom.ATOM_TYPE = type;
      new_pdb.push_back(atom);
    }
    else{
      if(KeepUnknown && (atom.ELEMENT != "H")) {
        new_pdb.push_back(atom);
      
        /*cerr << "UNKNOWN_TYPE_POT\t"<< atom.STRUCTURE << "\t" << atom.ID <<"\t"
              << atom.NAME << "\t" << atom.RESN << "\t" << atom.RESI << "\t" 
              << atom.CHAIN << atom.iCODE << "\n";*/
      }
    } 
  }

  auto orgpdb = OrganizePDB(new_pdb);
  for (auto& chain : orgpdb){
    auto& last_res = chain.back();
    for (auto& atom_p : last_res){
      if (atom_p->MOL_TYPE == "RNA" &&
          atom_p->ATOM_TYPE == 3 &&
          atom_p->NAME == "O3'") atom_p->ATOM_TYPE = 6;
    }
  }

  if(new_pdb.size() == 0)  {
    cerr << "NO_VALID_ATOMS_TYPE\n";
    throw -1;
  }
  pdb = new_pdb;
}

string
FixRESN(string resn){
  extern map<string,string> nucresn;
  try{
    return nucresn.at(resn);
  }
  catch(exception e){
    //cerr << "UNKNOWN_NON_STANDARD_RESN " << resn << "\n";
    return resn;
  }
  return resn;
}

string
FixNAME(string name){
  if(name == "O1P")
    return "OP1";
  if(name == "O2P")
    return "OP2";
  if(name == "O3P")
    return "OP3";
  if(name == "O4P")
    return "OP4";
  return name;
}



void
GetMolType(vector<atom_struct>& pdb){
  vector<string> amino_acids = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", 
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"
             };
  vector<string> nucs = {"A","T","G","C","U","DA","DT","DG","DC"};
  auto orgpdb = OrganizePDB(pdb);
  for (auto& chain : orgpdb){
    for (auto& res : chain){
      string type = "";
      bool nuc = false;
      bool o2p = false;
      for (auto& atom_p : res){
        if (atom_p->HETATM) continue;
        //if(atom_p->NAME == "CA" || atom_p->NAME == "C" || atom_p->NAME =="N" || atom_p->NAME == "O" || atom_p->NAME=="OXT"){ 
          //type = "PROTEIN";
        if(find(amino_acids.begin(),amino_acids.end(),atom_p->RESN) != amino_acids.end()){
          type = "PROTEIN";
          //else type ="";
          break;
        }
        if(atom_p->NAME == "P" || atom_p->NAME == "O5'" || atom_p->NAME == "C5'"){ 
          if(find(nucs.begin(),nucs.end(),atom_p->RESN) != nucs.end()) nuc = true;
        }
        if(atom_p->NAME == "O2'" && nuc){
          type = "RNA";
          o2p = true;
          break;
        }
      }
      if(nuc && !o2p) type = "DNA";
      else if (type.empty()) type = "LIGAND";
      //cout << "___" << std::endl;
      for (auto& atom_p : res){
        //if (type == "LIGAND") cout << atom_p->print() << std::endl;
        if (!atom_p->HETATM) atom_p->MOL_TYPE = type;
        else atom_p->MOL_TYPE="LIGAND";
      }
    }
  }
}

map<string, map<string, int> >
GetTypeMap(string filename){

  map<string, map<string, int> > result;
  ifstream potfile(filename, ifstream::in);
  while (potfile.good()){
    string line;
    getline(potfile, line);
    if (line.size() >= 1){
      if (string(line, 0, 1) == string("#"))
        continue;
      istringstream buffer(line);
      vector<string> tokens((istream_iterator<string>(buffer)), istream_iterator<string>());
      if (tokens.size() == 5){

        /*for (auto str : tokens){
        cerr << str << "\t";
        }*/
        auto& atoms = result[tokens[1]];
        atoms[tokens[2]] = stoi(tokens[4]);
        //cerr << "\n";
      }
    }
  }
  return result;
}


vector<atom_struct> //final structure
PDBparser(string pdbfname, //pdb filename or mol2
          string typefname, //type filename
          bool KeepUnknown){
  if(pdbfname.find(".mol2") != std::string::npos){
    return MOL2parser(pdbfname,typefname,KeepUnknown);
  }

  string basename = GetBasename(pdbfname);
  vector<atom_struct> result;


  auto lines = ReadPDB_ATOM(pdbfname);
  //cout << lines.size() << "\n";
  for (auto line : lines){
    int ID = 0;
    string name;
    string altloc;
    string resn;
    string chain;
    int resi = 0;
    string icode = " ";
    float X = 0;
    float Y = 0;
    float Z = 0;
    float oc = 0;
    float tf = 0;
    string ele = "";
    string chg = "";
    bool hetatm = line[0] == "HETATM";
    for (uint32 i = 1; i < line.size(); ++i){
      switch (i){
        case 1: //ID
          try{
            ID = stoi(line[i]);
          }
          catch(exception e){
            cerr << "ID error\n";
          }
          break;
        case 2: //NAME
          name = FixNAME(line[i]);
          break;
        case 3: //ALTLOC
          altloc = line[i];
          break;
        case 4: //RESN
          resn = line[i];
          break;
        case 5: //CHAIN
          if(line[i].empty()) chain = "_";
          else                chain = line[i];
          break;
        case 6: //RESI
          try{
            resi = stoi(line[i]);
          }
          catch(exception e){
            cerr << "RESI error\n";
          }
          break;
        case 7: //ICODE
          icode = line[i];
          break;
        case 8:  //X
          try{
            replace(line[i].begin(),line[i].end(),',','.');
            X = stod(line[i]);
          }
          catch(exception e){
            cerr << "X error\n";
          }
          break;
        case 9:  //Y
          try{
            replace(line[i].begin(),line[i].end(),',','.');
            Y = stod(line[i]);
          }
          catch(exception e){
            cerr << "Y error\n";
          }
          break;
        case 10: //Z
          try{
            replace(line[i].begin(),line[i].end(),',','.');
            Z = stod(line[i]);
          }
          catch(exception e){
            cerr << "Z error\n";
          }
          break;
        case 11: //OCCUP
          try{
            replace(line[i].begin(),line[i].end(),',','.');
            oc = stod(line[i]);
          }
          catch(exception e){
            cerr << "OC error\n";
          }
          break;
        case 12: //TFACT
          try{
            replace(line[i].begin(),line[i].end(),',','.');
            tf = stod(line[i]);
          }
          catch(exception e){
            cerr << "TF error\n";
          }
          break;
        case 13: //ELEM
          ele = line[i];
          break;
        case 14: //CHARGE
          chg = line[i];
          break;
      }
    }
    atom_struct current_atom(ID,
                             resi,
                             icode,
                             name,
                             resn,
                             chain,
                             ele,
                             basename,
                             "",
                             X,
                             Y,
                             Z,
                             altloc,
                             oc,
                             tf,
                             chg);
    current_atom.HETATM = hetatm;
    result.push_back(current_atom);
//     cout << line[1] << "\n";
//     cout << current_atom.print() << "\n";
  }

  RemoveAltLocs(result);
  GetMolType(result);
  GetAtomType(typefname,result,KeepUnknown);
  FixRepeatID(result);
//   for(auto atom : result){
//     cerr << "|"<<atom.print() << "|"<< "\n";
//   }
  if(result.size() == 0)  {
    cerr << "NO_VALID_ATOMS_FINAL\t" << pdbfname;
    throw -1;
  }
  //for (auto& atom : result) atom.DD_DISTANCES.resize(result.size());
  return result;
}

void
ChainSelector(vector<vector<string>>& selection,
              vector<atom_struct>& pdb){
  vector<uint32> removed_pos;
  vector<atom_struct> new_pdb;
  for (uint32 i = 0; i < pdb.size(); ++i){
    bool selected = false;
    for(auto obj : selection){
      for(auto chain : obj){
        if(chain == pdb[i].CHAIN){
          selected = true;
          break;
        }
      }
    }
    if(!selected) removed_pos.push_back(i);
  }
  for (uint32 i = 0; i < pdb.size(); ++i){
    if(removed_pos.end() != find(removed_pos.begin(),removed_pos.end(),i)) {
      //cerr << pdb[i].print() << " *\n";
      continue;
    }
    //cerr << pdb[i].print() << "\n";
    new_pdb.push_back(pdb[i]);
  }
  if(new_pdb.size() == 0) {
    cerr << "NO_VALID_ATOMS_IN_SELECTION\t";
    throw -1;
  }
  pdb = new_pdb;
}

void
FixRepeatID(vector<atom_struct>& pdb){
  map<uint32,uint32> counter;
  for (auto& atom : pdb){
    counter[atom.ID]++;
  }
  bool repeatID = false;
  for (auto it : counter){
    if (it.second != 1){
      repeatID = true;
      break;
    }
  }
  if(repeatID){
    for (uint32 i = 0; i < pdb.size(); ++i){
      pdb[i].ID = i+1; 
    }
  }
}

vector<atom_struct>
MOL2parser(string fname,
          string typefname,
          bool KeepUnknown){
  vector<atom_struct> result;
 

  ifstream mol2(fname);
/*  auto Trim = [](string& str){str.erase(std::remove_if(str.begin(),
                                                       str.end(),
                                                       ::isspace),
                                                       str.end());};*/

 
  string lastsection = "";
  map<string,vector<string>> sections;
  uint CHAINn = 0;
  while(mol2.good()){ 
    string line;
    vector<string> cline;
    getline(mol2,line);
    std::replace(line.begin(),line.end(),'\t',' ');
    if (std::string::npos != line.find("@<TRIPOS>MOLECULE")){
      //new molecule found
      if(sections.size() != 0){ 
        vector<atom_struct> temp;
        MOL2_parse_map(sections,
                     temp,
                     CHAINn);
        result.insert(result.end(),temp.begin(),temp.end());
        sections.clear();
        lastsection = "";
      }
    }

    if (std::string::npos != line.find("@<TRIPOS>")) {
      lastsection = line;
      sections[line] = vector<string>();
      continue;
    }
    
    if (lastsection == "") continue;
    else{
      sections[lastsection].push_back(line);
    }      
  }

  vector<atom_struct> temp;
  
  MOL2_parse_map(sections,
                 temp,
                 CHAINn);
  
  result.insert(result.end(),temp.begin(),temp.end());
  GetMolType(result);
  GetAtomType(typefname,result,KeepUnknown);
  for(auto& atom : result) atom.STRUCTURE = GetBasename(fname);
  return result;
}



void 
MOL2_parse_map(map<string,vector<string>>& sections,
               vector<atom_struct>& result,
               uint& CHAINn){


  auto Split = [](string const &input) { 
    std::istringstream buffer(input);
    vector<string> ret((std::istream_iterator<string>(buffer)),
                        std::istream_iterator<string>());
    return ret;};
 
  map<string,string> sybyl2ele = {{"C.3","C"},{"C.2","C"}, {"C.ar","C"}, {"C.2","C"}, {"C.1","C"}, {"C.cat","C"}, 
                                  {"N.4","N"}, {"N.3","N"}, {"N.2","N"}, {"N.pl3","N"}, {"N.1","N"},{"N.am","N"},{"N.ar","N"},
                                  {"O.3","O"}, {"O.2","O"},{"O.co2","O"},
                                  {"S.3","S"},{"S.2","S"}, {"S.O2","S"},{"S.O","S"}, {"S.o2","S"},{"S.o","S"}, 
	                          {"H","H"},    
                                  {"P.3","P"},
                                  {"CL","CL"},{"Cl","CL"},
                                  {"F","F"}};
  map<int,int> ID2pos;
  set<int> valid_ID;
 
  //get atoms
  if (!sections.count("@<TRIPOS>ATOM")){
    cerr << "NO_VALID_ATOMS_MOL2\n";
    return;
  }


  for (string& line : sections["@<TRIPOS>ATOM"]){
    
    vector<string> tokens = Split(line);
    if(tokens.size() < 6) continue;
    int ID = stoi(tokens[0]);
    string name = FixNAME(tokens[1]);
    float X = stof(tokens[2]);
    float Y = stof(tokens[3]);
    float Z = stof(tokens[4]);
    string type = tokens[5];
    if (type.size() > 5) type = type.substr(0,5);    
    if (type == "H") continue;
    valid_ID.insert(ID);
    string altloc = "";
    string resn = "";
    string chain = "_";
    int resi = 0;
    string icode = "";
    //cout << line << "\n";
    if(tokens.size() >= 6 || line.size() >= 56){ 
      try{
        int resi_ = stoi(line.substr(52,4));
        resi = resi_;
      }catch(exception e){}
      if (resi == 0){
       resi = stoi(tokens[6]);
      }
     
      if(tokens.size() > 7 || line.size() >= 69){
        string resn_;
        if(line.size() >= 69){ 
          try{
            resn_= line.substr(58,9);
          }
          catch(exception e){
            resn_ = tokens[7];
          }
        }
        else resn_ = tokens[7];
        resn_.erase(resn_.find_last_not_of(" \n\r\t")+1);
        resn_ = resn_.substr(0,3);
        for(char c : resn_){
          if (!isdigit(c)){
            if (c == '-') break;
            resn += c;
          }
          else break;
        }
      }
      else resn = "___";
    }   
    float oc = 0;
    float tf = 0;
    string ele = "";
    if(sybyl2ele.find(type) != sybyl2ele.end()) ele = sybyl2ele[type];
    atom_struct current_atom(ID,
                             resi,
                             icode,
                             name,
                             resn,
                             chain,
                             ele,
                             "",
                             "",
                             X,
                             Y,
                             Z,
                             altloc,
                             oc,
                             tf,
                             type);
    current_atom.DTYPE = "MOL2";
    current_atom.HETATM = false;
    ID2pos[ID] = result.size();
     
    result.push_back(current_atom);
  
  }
  if (!sections.count("@<TRIPOS>BOND")){
    cerr << "No bond section, falling back to chain identification by molecule type. (unimplemented)\n";
  }
  else{
    vector<vector<int>> bonded_sets;
    vector<vector<std::pair<int,int>>> plist;
    for(string& line : sections["@<TRIPOS>BOND"]){
      vector<string> tokens = Split(line);
      if (tokens.size() < 3) continue;
      int bID = stoi(tokens[0]);
      int oID = stoi(tokens[1]);
      int tID = stoi(tokens[2]);
      if(!valid_ID.count(oID) || !valid_ID.count(tID)) continue;
      int RoID = ID2pos[oID];
      int RtID = ID2pos[tID];
      vector<std::pair<int,int>> uset;
      uset.push_back(std::make_pair(RoID,RtID));
      plist.push_back(uset);
    }
    bool notdone = true;
    while(notdone){
      notdone = false;
      for(uint i = 0; i < plist.size(); ++i){
        auto& set_list = plist[i];
        for(uint j = 0; j < set_list.size(); ++j){
          int set[2] = {set_list[j].first,set_list[j].second};
          for(int& p : set){
            for(uint k = 0; k < plist.size(); ++k){
              if(k == i) continue;
              auto& o_set_list = plist[k];
              for(uint l = 0; l < o_set_list.size();++l){
                int o_set[2] = {o_set_list[l].first,o_set_list[l].second};
                for(int& o_p : o_set){
                  if(p == o_p){
                    set_list.push_back(o_set_list[l]);
                    o_set_list.pop_back();
                    notdone = true;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }

    
    for (uint i = 0; i < plist.size(); ++i){
      set<int> uID;
      auto& set_list = plist[i];
      for (int j = 0; j < set_list.size(); ++j){
        auto& pair = set_list[j];
        uID.insert(pair.first);
        uID.insert(pair.second);
      }
      if(uID.empty()) continue;
      vector<int> tmp;
      for(const int& id : uID) tmp.push_back(id);
      bonded_sets.push_back(tmp);
    }
    const char* alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_________________________";
    
    for(int i = 0; i < bonded_sets.size() && i < 40 ; ++i) {
      for (int& id : bonded_sets[i]){
        result[id].CHAIN=alphabet[i+CHAINn];
      }
    }
    CHAINn += bonded_sets.size();
  }
}
