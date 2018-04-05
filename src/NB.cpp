#include "stdafx.h"
#include "atom_struct.h"
#include "NB.h"

  //ATOM NODE DEFINITIONS#########################################################
  //first value is the center node. residues marked with * use commonCB.
  //common core atoms
  vector<string> commonCA =  {"CA","C","CB","N"};
  vector<string> commonO =   {"O","C"};
  vector<string> commonOXT = {"OXT","C"};
  vector<string> commonCB =  {"CB","CG"};
  //positive charges
  //arginine nodes
  vector<string> argCG = {"CG","CB","CD"};
  vector<string> argNE = {"NE","CD","CZ"};
  vector<string> argNH = {"CZ","NH1","NH2"};
  //histidine nodes*
  vector<string> hisCG =  {"CG","ND1","CD2"};
  vector<string> hisCE1 = {"CE1","ND1","NE2"};
  vector<string> hisCD2 = {"CD2","NE2","CG"};
  //lysine nodes*
  vector<string> lysCG = {"CG","CD","CB"};
  vector<string> lysCE = {"CE","CD","NZ"};
  //negative charges
  //aspartic acid nodes*
  vector<string> aspCG = {"CG","OD1","OD2"};
  //glutamic acid nodes*
  vector<string> gluCG = {"CG","CD"};
  vector<string> gluCD = {"CD","OE1","OE2"};
  //polar uncharged
  //serine
  vector<string> serCB = {"CB","OG"};
  //threonine
  vector<string> thrCB = {"CB","CG2","OG1"};
  //asparagine*
  vector<string> asnCG = {"CG","ND2","OD1"};
  //glutamine*
  vector<string> glnCG = {"CG","CD"};
  vector<string> glnCD = {"CD","OE1","NE2"};
  //specials
  //Cysteine
  vector<string> cysCB = {"CB","SG"};
  //Glycine
  vector<string> glyCA = {"CA","N","C"};
  //Proline*
  vector<string> proCG = {"CG","CB","CD"};
  vector<string> proCD = {"CD","N"};
  //Alanine
  //only common
  //Isoleucine
  vector<string> ileCB =  {"CB","CG2","CG1"};
  vector<string> ileCG1 = {"CG1","CD1"};
  //leucine*
  vector<string> leuCG = {"CG","CD1","CD2"};
  //methionine*
  vector<string> metSD = {"SD","CG","CE"};
  //phenylalanine*
  vector<string> pheCG =  {"CG","CD1","CD2"};
  vector<string> pheCE1 = {"CE1","CD1","CZ"};
  vector<string> pheCE2 = {"CE2","CD2","CZ"};
  //tryptophan*
  vector<string> trpCG =  {"CG","CD1","CD2"};
  vector<string> trpNE1 = {"NE1","CD1","CE2"};
  vector<string> trpCE2 = {"CE2","CD2","CZ2"};
  vector<string> trpCH2 = {"CH2","CZ2","CZ3"};
  vector<string> trpCE3 = {"CE3","CD2","CZ3"};
  //tyrosine*
  vector<string> tyrCG =  {"CG","CD1","CD2"};
  vector<string> tyrCE1 = {"CE1","CD1","CZ"};
  vector<string> tyrCE2 = {"CE2","CD2","CZ"};
  vector<string> tyrCZ =  {"CZ","OH"};
  //valine
  vector<string> valCB = {"CB","CG1","CG2"};

  //DNA/RNA definitions
  //core ribose atoms
  vector<string> commonRC1 = {"C1'","C2'","O4'"};
  vector<string> commonRC2 = {"C2'","C1'","C3'","O2'"};
  vector<string> commonRC3 = {"C3'","C2'","C4'","O3'"};
  vector<string> commonRC4 = {"C4'","C3'","O4'"};
  vector<string> commonRC5 = {"C5'","O5'","C4'"};
  vector<string> commonRP =  {"P","O5'","OP1","OP2","OP3"};
  //Adenine base for RNA and DNA
  vector<string> A_N9 = {"N9","C1'","C4","C8"};
  vector<string> A_N7 = {"N7","C8","C5"};
  vector<string> A_C2 = {"C2","N1","N3"};
  vector<string> A_C4 = {"C4","C5","N3"};
  vector<string> A_C6 = {"C6","C5","N1","N6"};
  //Guanine base for RNA and DNA
  vector<string> G_N9 = {"N9","C1'","C4","C8"};
  vector<string> G_N7 = {"N7","C8","C5"};
  vector<string> G_C2 = {"C2","N1","N2","N3"};
  vector<string> G_C4 = {"C4","N3","N9"};
  vector<string> G_C6 = {"C6","C5","N1","O6"};
  //Thymine base for RNA and DNA
  vector<string> T_N1 = {"N1","C1'","C2","C6"};
  vector<string> T_C5 = {"C5","C4","C6","C7"};
  vector<string> T_C4 = {"C4","N3","C5","O4"};
  vector<string> T_C2 = {"C2","N1","N3","O2"};
  //Cytosine base for RNA and DNA
  vector<string> C_N1 = {"N1","C1'","C2","C6"};
  vector<string> C_C5 = {"C5","C4","C6"};
  vector<string> C_C4 = {"C4","N3","C5","N4"};
  vector<string> C_C2 = {"C2","N1","N3","O2"};
  //END ATOM NODE DEFINITIONS###################################################


vector<vector<string>>
NB_GenPerm(vector<string>& array){
  auto l = array.size();
  vector<vector<string> > result;
  for (uint32 i = 1; i < l; i++){
    vector<string> p1 = {array[0],array[i]};
    vector<string> p2 = {array[i],array[0]};
    result.push_back(p1);
    result.push_back(p2);
  }
  return result;
}

bool
NB_IsBonded(uint32               i,
            uint32               j,
            vector<atom_struct>& pdb){

  auto& atom_i = pdb[i];
  auto& atom_j = pdb[j];
  float S_S_dist = 2.1;
  //##########################################################################
  int res_d = atom_i.RESI - atom_j.RESI;

  //cases between adjacent residues/nucleotides
  if (res_d == 1 && atom_i.CHAIN == atom_j.CHAIN){
    //bond N terminal to C terminal in proteins
    if(atom_i.MOL_TYPE == "PROTEIN" && atom_j.MOL_TYPE == "PROTEIN" &&
       atom_i.NAME == "N" && atom_j.NAME == "C"){
      return true;
    }
    //bond P terminal to O3' terminal in RNA/DNA
    if(((atom_i.MOL_TYPE == "RNA" && atom_j.MOL_TYPE == "RNA") ||
        (atom_i.MOL_TYPE == "DNA" && atom_j.MOL_TYPE == "DNA")) &&
       atom_i.NAME == "P" && atom_j.NAME == "O3'" ){
      return true;
    }
  }

  if (res_d == -1 && atom_i.CHAIN == atom_j.CHAIN){
    //bond N terminal to C terminal in proteins
    if(atom_i.MOL_TYPE == "PROTEIN" && atom_j.MOL_TYPE == "PROTEIN" &&
       atom_i.NAME == "C" && atom_j.NAME == "N"){
      return true;
    }
    //bond P terminal to O3' terminal in RNA/DNA
    if(((atom_i.MOL_TYPE == "RNA" && atom_j.MOL_TYPE == "RNA") ||
       (atom_i.MOL_TYPE == "DNA" && atom_j.MOL_TYPE == "DNA")) &&
       atom_i.NAME == "O3'" && atom_j.NAME == "P"){
      return true;
    }
  }

  //sulfur bonds, will fail if there are many sulfur bonds close together
  if((atom_i.MOL_TYPE == "PROTEIN" && atom_j.MOL_TYPE == "PROTEIN") &&
    atom_i.NAME == "SG" && atom_j.NAME == "SG"){
    if (atom_i.DISTANCES.count(j) && atom_i.DISTANCES[j] <= S_S_dist) return true;
    else return false;
  }
  //cout << "NB0" << atom_i.print() << "\n";
  //cases inside residues/nucleotides
  if (res_d == 0 &&
      atom_i.RESI == atom_j.RESI &&
      atom_i.CHAIN == atom_j.CHAIN){
    if(atom_i.MOL_TYPE == "PROTEIN" && atom_j.MOL_TYPE == "PROTEIN"){
    //common core atoms
      vector<vector<vector<string>>> CommonAtoms = { NB_GenPerm(commonCA),
                                                     NB_GenPerm(commonO),
                                                     NB_GenPerm(commonOXT) };
      for (auto& set : CommonAtoms){
        for(auto& pair : set){
          if(atom_i.NAME == pair[0] && atom_j.NAME == pair[1]) return true;
        }
      }

      if(atom_i.RESN == "ARG" && atom_j.RESN == "ARG"){
        CommonAtoms = { NB_GenPerm(argCG),
                        NB_GenPerm(argNE),
                        NB_GenPerm(argNH) };
      }

      if(atom_i.RESN == "HIS" && atom_j.RESN == "HIS"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(hisCG),
                        NB_GenPerm(hisCE1),
                        NB_GenPerm(hisCD2) };
      }
      if(atom_i.RESN == "LYS" && atom_j.RESN == "LYS"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(lysCG),
                        NB_GenPerm(lysCE)};
      }
      if(atom_i.RESN == "ASP" && atom_j.RESN == "ASP"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(aspCG)};
      }
      if(atom_i.RESN == "GLU" && atom_j.RESN == "GLU"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(gluCG),
                        NB_GenPerm(gluCD)};
      }
      if(atom_i.RESN == "SER" && atom_j.RESN == "SER"){
        CommonAtoms = { NB_GenPerm(serCB)};
      }
      if(atom_i.RESN == "THR" && atom_j.RESN == "THR"){
        CommonAtoms = { NB_GenPerm(thrCB)};
      }
      if(atom_i.RESN == "ASN" && atom_j.RESN == "ASN"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(asnCG)};
      }
      if(atom_i.RESN == "GLN" && atom_j.RESN == "GLN"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(glnCG),
                        NB_GenPerm(glnCD)};
      }
      if(atom_i.RESN == "CYS" && atom_j.RESN == "CYS"){
        CommonAtoms = { NB_GenPerm(cysCB) };
      }
      if(atom_i.RESN == "GLY" && atom_j.RESN == "GLY"){
        CommonAtoms = { NB_GenPerm(glyCA) };
      }
      if(atom_i.RESN == "PRO" && atom_j.RESN == "PRO"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(proCG),
                        NB_GenPerm(proCD)};
      }
      if(atom_i.RESN == "ALA" && atom_j.RESN == "ALA"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(gluCG),
                        NB_GenPerm(gluCD)};
      }
      if(atom_i.RESN == "ILE" && atom_j.RESN == "ILE"){
        CommonAtoms = { NB_GenPerm(ileCB),
                        NB_GenPerm(ileCG1)};
      }

      if(atom_i.RESN == "LEU" && atom_j.RESN == "LEU"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(leuCG)
        };
      }
      if(atom_i.RESN == "MET" && atom_j.RESN == "MET"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(metSD)};
      }
      if(atom_i.RESN == "PHE" && atom_j.RESN == "PHE"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(pheCG),
                        NB_GenPerm(pheCE1),
                        NB_GenPerm(pheCE2)};
      }
      if(atom_i.RESN == "TRP" && atom_j.RESN == "TRP"){
        CommonAtoms = { NB_GenPerm(commonCB),
                        NB_GenPerm(trpCG),
                        NB_GenPerm(trpNE1),
                        NB_GenPerm(trpCE2),
                        NB_GenPerm(trpCH2),
                        NB_GenPerm(trpCE3) };
      }
      if(atom_i.RESN == "TYR" && atom_j.RESN == "TYR"){
        CommonAtoms = {NB_GenPerm(commonCB),
                       NB_GenPerm(tyrCG),
                       NB_GenPerm(tyrCE1),
                       NB_GenPerm(tyrCE2),
                       NB_GenPerm(tyrCZ)};
      }
      if(atom_i.RESN == "VAL" && atom_j.RESN == "VAL"){
        CommonAtoms = { NB_GenPerm(valCB) };
      }

      for (auto& set : CommonAtoms){
        for(auto& pair : set){
          if(atom_i.NAME == pair[0] && atom_j.NAME == pair[1]) return true;
        }
      }
    }
    //cout << "NB1" << atom_i.print() << " " << "|" <<atom_i.MOL_TYPE << "|" << "\n";
    if((atom_i.MOL_TYPE == "RNA" && atom_j.MOL_TYPE == "RNA") ||
       (atom_i.MOL_TYPE == "DNA" && atom_j.MOL_TYPE == "DNA")){
      //common core atoms
      vector<vector<vector<string>>> CommonAtoms = { NB_GenPerm(commonRC1),
                                                     NB_GenPerm(commonRC2),
                                                     NB_GenPerm(commonRC3),
                                                     NB_GenPerm(commonRC4),
                                                     NB_GenPerm(commonRC5),
                                                     NB_GenPerm(commonRP) };
      for (auto& set : CommonAtoms){
        for(auto& pair : set){
          if(atom_i.NAME == pair[0] && atom_j.NAME == pair[1]) return true;
        }
      }
      if((atom_i.RESN == "DA" || atom_i.RESN == "A" )&&
         (atom_j.RESN == "DA" || atom_j.RESN == "A")){
        CommonAtoms = { NB_GenPerm(commonRC1),
                        NB_GenPerm(commonRC2),
                        NB_GenPerm(commonRC3),
                        NB_GenPerm(commonRC4),
                        NB_GenPerm(commonRC5),
                        NB_GenPerm(commonRP),
                        NB_GenPerm(A_N9),
                        NB_GenPerm(A_N7),
                        NB_GenPerm(A_C2),
                        NB_GenPerm(A_C4),
                        NB_GenPerm(A_C6)};
      }
      if((atom_i.RESN == "DG" || atom_i.RESN == "G" )&&
         (atom_j.RESN == "DG" || atom_j.RESN == "G")){
        CommonAtoms = { NB_GenPerm(commonRC1),
                        NB_GenPerm(commonRC2),
                        NB_GenPerm(commonRC3),
                        NB_GenPerm(commonRC4),
                        NB_GenPerm(commonRC5),
                        NB_GenPerm(commonRP),
                        NB_GenPerm(G_N9),
                        NB_GenPerm(G_N7),
                        NB_GenPerm(G_C2),
                        NB_GenPerm(G_C4),
                        NB_GenPerm(G_C6)};
      }
      if((atom_i.RESN == "DT" || atom_i.RESN == "T" )&&
         (atom_j.RESN == "DT" || atom_j.RESN == "T")){
        CommonAtoms = { NB_GenPerm(commonRC1),
                        NB_GenPerm(commonRC2),
                        NB_GenPerm(commonRC3),
                        NB_GenPerm(commonRC4),
                        NB_GenPerm(commonRC5),
                        NB_GenPerm(commonRP),
                        NB_GenPerm(T_N1),
                        NB_GenPerm(T_C5),
                        NB_GenPerm(T_C4),
                        NB_GenPerm(T_C2)};
      }
      if((atom_i.RESN == "DC" || atom_i.RESN == "C" )&&
         (atom_j.RESN == "DC" || atom_j.RESN == "C")){
        CommonAtoms = { NB_GenPerm(commonRC1),
                        NB_GenPerm(commonRC2),
                        NB_GenPerm(commonRC3),
                        NB_GenPerm(commonRC4),
                        NB_GenPerm(commonRC5),
                        NB_GenPerm(commonRP),
                        NB_GenPerm(C_N1),
                        NB_GenPerm(C_C5),
                        NB_GenPerm(C_C4),
                        NB_GenPerm(C_C2)};
      }
      for (auto& set : CommonAtoms){
        for(auto& pair : set){
          if(atom_i.NAME == pair[0] && atom_j.NAME == pair[1]) return true;
        }
      }
    }
  }
  return false;
}

void
NB_LinkPDB(vector<atom_struct>&   pdb){
  auto l = pdb.size();
  for (uint32 i = 0; i < l; i++){
    auto& atom_i = pdb[i];
    for (uint32 j = i + 1; j < l; j++ ){
      auto& atom_j = pdb[j];
      if (NB_IsBonded(i,j,pdb)){
        atom_i.BONDED.push_back(j);
        atom_j.BONDED.push_back(i);
      }
    }
  }
}

vector<vector<uint32>>
NB_GetBonded(uint32               i,
             uint32               D,
             vector<atom_struct>& pdb){
  //cout << "Depth " << D << "\n";
  uint32 depth = 0;
  int depth_inc = 1;
  int n_depth_inc = 0;
  queue<uint32> Nodes;
  vector<vector<uint32>> result;
  vector<char> walked(pdb.size(),false);

  result.resize(D+1);
  Nodes.push(i);
  walked[i] = true;
  while(!Nodes.empty()){
    uint32 cur = Nodes.front();
    Nodes.pop();
    n_depth_inc += pdb[cur].BONDED.size();
    for (uint32 child : pdb[cur].BONDED){
      if(!walked[child]){
        result[depth].push_back(child);
      }
      else n_depth_inc--;
    }

    if ( --depth_inc == 0){
      if(++depth > D) break;
      depth_inc = n_depth_inc;
      n_depth_inc = 0;
    }

    for (uint32 child : pdb[cur].BONDED){
      if(!walked[child]){
        walked[child] = true;
        Nodes.push(child);
      }
    }
  }
  return result;
}

void
NB_RemoveBonded(vector<atom_struct>& pdb,
                int                  nb,
                vector<char>&        Imtrx){
  if (nb == -1) return;
  auto l = pdb.size();
  for (uint32 i = 0; i < l; i++){
    auto list = NB_GetBonded(i,nb,pdb);
    for (auto& depth : list){
      for (uint32 j : depth){
        Imtrx[i + j * l] = false;
      }
    }
  }
}
