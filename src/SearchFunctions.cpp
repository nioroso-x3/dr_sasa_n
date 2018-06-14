#include "stdafx.h"
#include "atom_struct.h"
#include "SearchFunctions.h"
#define MAX_INT pow(2,24)
/*vector<residue_struct>
GetResidues(vector<atom_struct>& pdb){
    map<int,residue_struct> data;
    vector<residue_struct> result;
    for (auto* p = pdb.data(); p != p + pdb.size(); ++p){
        if(data.count(p->RESI) > 0){
            auto& res = data[p->RESI];
            res.RESN = p->RESN;
            res.RESI = p->RESI;
            res.CHAIN = p->CHAIN;
            res.STRUCT_TYPE = p->STRUCT_TYPE
            res.STRUCTURE = p->STRUCTURE;
            res.SASA+=p->SASA;
            res.dSASA+=p->EXT1;
            res.ATOMS.push_back(p);
        }
        else{
            auto& res = data[p->RESI]
            res.SASA+=p->SASA;
            res.dSASA+=p->EXT1;
            res.ATOMS.push_back(p);
        }
    }
    result.reserve(data.size());
    for (auto it = data.begin(); it != data.end(); ++it){
        result.push_back(it->second);
    }
    return result;
}*/

vector<atom_struct*>
GetAtoms(string               resn,
         string               chain,
         int                  resi,
         string               icode,
         vector<atom_struct>& pdb){
  vector<atom_struct*> result;
  auto* p = pdb.data();
  auto* end = p + pdb.size();
  auto searcher = [](string left,string right){
    bool result = right.empty() ? true : left == right;
    return result;
  };
  for (;p != end ; ++p){
    if(searcher(p->RESN,resn) &&
       searcher(p->CHAIN,chain) &&
       ((resi == MAX_INT) ? true : (p->RESI == resi)) && 
       searcher(p->iCODE,icode)){
//      cout << chain << " " << resi << " | " << p->RESN << " " << int(p->CHAIN == chain) << " " << p->RESI << "\n";
      result.push_back(p);
      }
    }
   return result;
}

vector<vector<vector<atom_struct*>>>
OrganizePDB(vector<atom_struct>& pdb){
  vector<vector<vector<atom_struct*>>> result;
  set<string> CHAINs;
  map<string,set<int>> RESIs;
  for(auto& atom : pdb){
    CHAINs.insert(atom.CHAIN);
    RESIs[atom.CHAIN].insert(atom.RESI);
  }
  vector<string> CHAIN(CHAINs.begin(),CHAINs.end());
  sort(CHAIN.begin(),CHAIN.end());
  for (auto& chain : CHAIN){
    vector<vector<atom_struct*>> cur_chain;
    vector<int> RESI(RESIs[chain].begin(),RESIs[chain].end());
    sort(RESI.begin(),RESI.end());
//    cout <<"___\n";
    for (auto& resi : RESI){
       auto atoms = GetAtoms("",chain,resi,"",pdb);
//       for(auto& atom : atoms){
//         cout << atom->print() << "\n";
//       }
       cur_chain.push_back(atoms);
    }
    result.push_back(cur_chain);
  }
  return result;
}


/*vector<atom_struct*>
GetAtomsFromRESI(int,vector<residue_struct>&){
    vector<atom_struct*> r;
    return r;
}*/

/*vector<residue_struct*>
GetRESN(string,vector<residue_struct>&){
    vector<residue_struct*> r;
    return r;
}*/

map<string,int>
GetTypeFromRESI(string resn,
                string chain,
                int resi,
        string icode,
        vector<atom_struct>& pdb){
  map<string,int> result;
  auto* p = pdb.data();
  auto* end = p + pdb.size();
  for (;p != end ; ++p){
    if(p->RESN == resn &&
           p->CHAIN == chain &&
           p->RESI == resi &&
           p->iCODE == icode) {
      result["ATOM_TYPE"] = p->ATOM_TYPE;
      result["ATOM_TYPE_40"] = p->ATOM_TYPE_40;
      break;
    }
  }

  return result;

}

// vector<string> 
// LineSplitter(string line,
//              char   sep){
//   stringstream   linestream(line);
//   string         data;
//   vector<string> thisline;
//   getline(linestream, data, sep);
//   thisline.push_back(data);
//   while(linestream){
//     string cur;
//     linestream >> cur;
//     thisline.push_back(cur);
//     //cout << cur << "\t";
//   }
//     //cout << "\n";
//   return thisline;
// }

vector<string> 
LineSplitter(string line,
             char   sep){
  stringstream   linestream(line);
  string         data;
  vector<string> thisline;
  while(linestream){
    if(getline(linestream, data, sep))
      thisline.push_back(data);
    //cout << cur << "\t";
  }
    //cout << "\n";
  return thisline;
}

vector<vector<string>>
ReadTabSeparatedFile(string fname){
  //cout << fname << "\n";
  ifstream file(fname);
  string   line;
  vector<vector<string>> result;
  while(getline(file, line))
  {
    vector<string> thisline = LineSplitter(line,'\t');
    result.push_back(thisline);
  }
  return result;
}

vector<vector<string>>
ReadTabSeparatedString(string data){
  //cout << fname << "\n";
  istringstream file(data);
  string   line;
  vector<vector<string>> result;
  while(getline(file, line))
  {
    vector<string> thisline = LineSplitter(line,'\t');
    result.push_back(thisline);
  }
  return result;
}
