#include "stdafx.h"
#include "atom_struct.h"
#include "histogram.h"
#include "SearchFunctions.h"
#include <cstdio>
#include <limits>
#include "SolverDataProcessing.h"


void
RelativeSASA(vector<atom_struct>& pdb){
  map<string,float> reldata = {{"ALA",92.40211452},
                                {"GLN",184.71688},
                                {"TRP",235.3483229},
                                {"SER",122.353095},
                                {"ARG",219.941475},
                                {"ILE",158.07277},
                                {"ASN",146.06073},
                                {"ASP",154.71124},
                                {"HIS",194.617376},
                                {"MET",191.547244},
                                {"LYS",201.689792},
                                {"LEU",169.04496452},
                                {"THR",145.526463}, 
                                {"PHE",201.7065},
                                {"TYR",209.07715},
                                {"GLU",160.482561},
                                {"CYS",130.28853},
                                {"PRO",126.877028},
                                {"GLY",80.23533}, 
                                {"VAL",138.90233},  
                                { "A",163.5912717},
                                { "C",162.86302775},
                                { "T",173.03036075},
                                { "G",164.253836},
                                { "DA",163.5912717},
                                { "DC",162.86302775},
                                { "DT",173.03036075},
                                { "DG",164.253836}};
  for(uint i = 0; i< pdb.size(); ++i){
    auto& atom = pdb[i];
    try{
      float div = reldata.at(atom.RESN);
      atom.SASA /= div;
    }
    catch(exception e){
      continue;
    }
  }
  
}

void
PrintDNA_ProtResultsByAA(vector<atom_struct>& pdb,  // pdb struct
                         string               output){
  vector<pair<uint32,uint32> > pairs;
  ofstream outf(output, ofstream::out);
  map<pair<int,int>,float> resi_data;
  map<int,float> resi_dsasa;

  outf << "<RESI_1>\t<RESI_2>\t<dSASA_by_2>\t<RESI_1_dSASA>\t<RESI_2_dSASA>" << std::endl;


  for (uint32 i = 0; i < pdb.size(); ++i){
    auto& atom_i = pdb[i];
    if(!atom_i.ACTIVE ||
      atom_i.EXT1 == 0) continue;
    if (resi_dsasa.count(atom_i.RESI) == 0) resi_dsasa[atom_i.RESI] = 0.0;
    resi_dsasa[atom_i.RESI] += atom_i.EXT1;
    for (uint32 j = 0; j < atom_i.INTERACTION_SASA_P.size(); ++j){
      auto j_pos = atom_i.INTERACTION_SASA_P[j];
      auto& atom_j = pdb[j_pos];
      pair<uint32,uint32> p1(i,j_pos);
      pair<uint32,uint32> p2(j_pos,i);
      if(!atom_j.ACTIVE ||
        atom_j.EXT1 == 0)
        continue;
      if(pairs.end() != find(pairs.begin(),pairs.end(),p1) ||
        pairs.end() != find(pairs.begin(),pairs.end(),p2) )
        continue;
      pairs.push_back(p1);
      pairs.push_back(p2);
      float ji = atom_i.CONTACT_AREA.at(atom_j.ID);
      float ij = atom_j.CONTACT_AREA.count(atom_i.ID) ? atom_j.CONTACT_AREA.at(atom_i.ID) : 0.0 ;

      auto res_ij = make_pair(atom_i.RESI,atom_j.RESI);
      auto res_ji = make_pair(atom_j.RESI,atom_i.RESI);
      if (resi_data.count(res_ij) == 0) resi_data[res_ij] = 0.0;
      if (resi_data.count(res_ji) == 0) resi_data[res_ji] = 0.0;
      resi_data[res_ij] += ji;
      resi_data[res_ji] += ij;
    }
  }
  for (auto& it : resi_data){
    auto res_pair = it.first;
    auto value = it.second;
    auto resI = get<0>(res_pair);
    auto resJ = get<1>(res_pair);
    auto dsasa1 = resi_dsasa[resI];
    auto dsasa2 = resi_dsasa[resJ];
    outf << resI << "\t" << resJ << "\t" << value << "\t" << dsasa1 << "\t" << dsasa2 << std::endl;
  }
  outf.close();
}

void
GenerateInterBSAMatrix(vector<atom_struct>&                  pdb,
                        map<vector<string>, vector<float>>& matrixIJatom,
                        map<vector<string>, vector<float>>& matrixIJres,
                        map<vector<string>, vector<string>>& COLatom,
                        map<vector<string>, vector<string>>& COLres,
                        map<vector<string>, vector<string>>& ROWatom,
                        map<vector<string>, vector<string>>& ROWres,
                        map<vector<string>, vector<uint32>>& COLatomtype,
                        map<vector<string>, vector<uint32>>& ROWatomtype ){
  CalculateDNA_ProtInteractions(pdb);
  set<string> object_types;
  for (auto& atom : pdb){
    object_types.insert(atom.STRUCT_TYPE);
  }
  vector<string> objs(object_types.begin(),object_types.end());
//pragma omp parallel for schedule(dynamic)
  for(uint32 i = 0; i < objs.size()-1; ++i){
    for (uint32 j = i+1; j < objs.size(); ++j){
      string objA = objs[i];
      string objB = objs[j];
      //cout << "OBJS: "<< objA << " " << objB << "\n";

      vector<string> ColL;
      vector<string> LineL;
      vector<uint32> ColL_type;
      vector<uint32> LineL_type;
      map<string,uint32> ColM;
      map<string,uint32> LineM;
      vector<float> IcJ; //J causes to I
      vector<float> JcI; //I causes to J

      vector<string> ColLres;
      vector<string> LineLres;
      map<string,uint32> ColMres;
      map<string,uint32> LineMres;
      vector<float> IcJres; //J causes to I
      vector<float> JcIres; //I causes to J

      for (auto& atom : pdb){
        string aID = atom.sID();
        string rID = atom.rsID();
        if(atom.STRUCT_TYPE == objA){
          ColL.push_back(aID);
          ColL_type.push_back(atom.ATOM_TYPE);
          ColM[aID] = ColL.size()-1;
          ColLres.push_back(rID);
          //cout << " COL " << aID << " " << ColL.size()-1 << std::endl;
        }
        if(atom.STRUCT_TYPE == objB){
          LineL.push_back(aID);
          LineL_type.push_back(atom.ATOM_TYPE);
          LineM[aID] = LineL.size()-1;
          LineLres.push_back(rID);
          //cout << " LIN " << aID << " " << LineL.size()-1 << std::endl;
        }
      }
      ColLres.erase( unique( ColLres.begin(), ColLres.end() ), ColLres.end() );
      LineLres.erase( unique( LineLres.begin(), LineLres.end() ), LineLres.end() );

      for (uint32 k = 0; k < ColLres.size();++k) ColMres[ColLres[k]] = k;
      for (uint32 k = 0; k < LineLres.size();++k) LineMres[LineLres[k]] = k;

      IcJres.resize(ColLres.size()*LineLres.size(),0.0);
      JcIres.resize(ColLres.size()*LineLres.size(),0.0);

      IcJ.resize(ColL.size()*LineL.size(),0.0);
      JcI.resize(ColL.size()*LineL.size(),0.0);
      
      
      uint32 lm = ColL.size();
      uint32 lm_res = ColLres.size();

      for(auto& atom : pdb){
        string aID = atom.sID();
        string rID = atom.rsID();
        if(atom.STRUCT_TYPE == objA){
          //cout << "ATOM: " << aID << " " << atom.STRUCT_TYPE << "\n";
          for(uint32 pos : atom.INTERACTION_SASA_P){
            if (pdb[pos].STRUCT_TYPE != objB) continue;
            //cout <<"    OBJ_A: " << pdb[pos].sID() << "\n";
            for (uint32 r = 0; r < atom.ov_table.size(); ++r){ //iterate through this atoms overlap table
              auto& ov = atom.ov_table[r];                     //get vector of atom positions that form this overlap
              if(ov.end() != find(ov.begin(),ov.end(),pos)){   //check if this other atom is part of current overlap
                uint32 col = ColM[aID];                        //current atom position in column
                uint32 line = LineM[pdb[pos].sID()];           //other atom position in line
                IcJ[col + line * lm] += atom.ov_norm_area[r];  //sum area to atom matrix
                uint32 col_res = ColMres[rID];                 //column of current residue 
                uint32 line_res = LineMres[pdb[pos].rsID()];   //line of current residue
                IcJres[col_res + line_res * lm_res] += atom.ov_norm_area[r]; //sum area to residue matrix
              }
            }
          }
        }
        //same commentas as object A
        if(atom.STRUCT_TYPE == objB){
          //cout << "ATOM: " << atom.sID() << " " << atom.STRUCT_TYPE << "\n";
          for(uint32 pos : atom.INTERACTION_SASA_P){
            if (pdb[pos].STRUCT_TYPE != objA) continue;
            //cout <<"    OBJ_B: " << pdb[pos].sID() << "\n";
            for (uint32 r = 0; r < atom.ov_table.size(); ++r){
              auto& ov = atom.ov_table[r];
              if(ov.end() != find(ov.begin(),ov.end(),pos)){
                uint32 col = ColM[pdb[pos].sID()];
                uint32 line = LineM[aID];
                JcI[col + line * lm] += atom.ov_norm_area[r];
                uint32 col_res = ColMres[pdb[pos].rsID()];
                uint32 line_res = LineMres[rID];
                JcIres[col_res + line_res * lm_res] += atom.ov_norm_area[r];
              }
            }
          }
        }
      }
      vector<string> keyIJ = { objA,objB };
      vector<string> keyJI = { objB,objA };
      
      matrixIJatom[keyIJ] = IcJ;
      matrixIJatom[keyJI] = JcI;
      
      matrixIJres[keyIJ] = IcJres;
      matrixIJres[keyJI] = JcIres;
      
      COLatom[keyIJ] = ColL;
      COLatom[keyJI] = ColL;
      
      ROWatom[keyIJ] = LineL;
      ROWatom[keyJI] = LineL;
      
      COLres[keyIJ] = ColLres;
      COLres[keyJI] = ColLres;
      
      ROWres[keyIJ] = LineLres;
      ROWres[keyJI] = LineLres;
      
      COLatomtype[keyIJ] = ColL_type;
      COLatomtype[keyJI] = ColL_type;
      
      ROWatomtype[keyIJ] = LineL_type;
      ROWatomtype[keyJI] = LineL_type;
      
    }
  }
}

void
GenerateIntraBSAMatrix(vector<atom_struct>& pdb,
                       vector<float>&      matrixIJatom,
                       vector<float>&      matrixIJres,
                       vector<string>&      COLatom,
                       vector<string>&      COLres,
                       vector<string>&      ROWatom,
                       vector<string>&      ROWres,
                       vector<uint32>&      COLatomtype,
                       vector<uint32>&      ROWatomtype ){
  CalculateDNA_ProtInteractions(pdb);
  set<string> object_types;
  for (auto& atom : pdb){
    object_types.insert(atom.CHAIN);
  }
  vector<string> objs(object_types.begin(),object_types.end());
  
  //variables for atom matrix
  vector<string> ColL;
  vector<string> LineL;
  vector<uint32> ColL_type;
  vector<uint32> LineL_type;
  map<string,uint32> ColM;
  map<string,uint32> LineM;
  vector<float> mtrx; 
  //variables for residue matrix
  vector<string> ColLres;
  vector<string> LineLres;
  map<string,uint32> ColMres;
  map<string,uint32> LineMres;
  vector<float> mtrx_res; 

  for (auto& atom : pdb){
    string aID = atom.sID();
    string rID = atom.rsID();
    ColL_type.push_back(atom.ATOM_TYPE);
    ColL.push_back(aID);
    ColM[aID] = ColL.size()-1;
    ColLres.push_back(rID);
    //cout << " COL " << aID << " " << ColL.size()-1 << std::endl;
    LineL_type.push_back(atom.ATOM_TYPE);
    LineL.push_back(aID);
    LineM[aID] = LineL.size()-1;
    LineLres.push_back(rID);
    //cout << " LIN " << aID << " " << LineL.size()-1 << std::endl;
  }
  ColLres.erase( unique( ColLres.begin(), ColLres.end() ), ColLres.end() );
  LineLres.erase( unique( LineLres.begin(), LineLres.end() ), LineLres.end() );
  for (uint32 k = 0; k < ColLres.size();++k) ColMres[ColLres[k]] = k;
  for (uint32 k = 0; k < LineLres.size();++k) LineMres[LineLres[k]] = k;
  
  mtrx.resize(ColL.size()*LineL.size(),0.0);
  mtrx_res.resize(ColLres.size()*LineLres.size(),0.0);
  uint32 lm = ColL.size();
  uint32 lm_res = ColLres.size();
  
  //pragma omp parallel for schedule(dynamic)
  for(uint32 i = 0; i < pdb.size(); ++i){
    auto& atom = pdb[i];
    string aID = atom.sID();
    string rID = atom.rsID();
   //cout << "ATOM: " << aID << " " << atom.STRUCT_TYPE << "\n";
   for(uint32 pos : atom.INTERACTION_SASA_P){
      //cout <<"    OBJ_A: " << pdb[pos].sID() << "\n";
      for (uint32 r = 0; r < atom.ov_table.size(); ++r){
        auto& ov = atom.ov_table[r];
        if(ov.end() != find(ov.begin(),ov.end(),pos)){
          uint32 col = ColM[aID];
          uint32 line = LineM[pdb[pos].sID()];
          
          uint32 col_res = ColMres[rID];
          uint32 line_res = LineMres[pdb[pos].rsID()];
          //pragma omp critical(intra_bsa_matrix)
          {
            mtrx[col + line * lm] += atom.ov_norm_area[r];
            mtrx_res[col_res + line_res * lm_res] += atom.ov_norm_area[r];
          }
        }
      }
    }
  }
  matrixIJatom= mtrx;
  matrixIJres = mtrx_res;
  COLatom = ColL;
  ROWatom = LineL;
  COLres = ColLres;
  ROWres = LineLres;
  COLatomtype = ColL_type;
  ROWatomtype = LineL_type;
}

void
PrintDDMatrix(vector<atom_struct>& pdb,
              vector<char>&        imtrx,
              string               I_TYPE,
              string               J_TYPE,
              string               output){
  

  ofstream outf(output);
  outf << "x\t";
  if (I_TYPE == J_TYPE){
    for (int i = 0; i < pdb.size(); ++i){
      if(pdb[i].MOL_TYPE != I_TYPE) continue;
      outf << pdb[i].sID();
      if( i != (pdb.size()-1)) outf << "\t";
    }
    outf << "\n";
    for (int i = 0; i < pdb.size(); ++i){
      if (pdb[i].MOL_TYPE != I_TYPE) continue;
      outf << pdb[i].sID() << "\t";
      for (int j = 0; j < pdb.size(); ++j){
        if(pdb[j].MOL_TYPE != I_TYPE) continue;
        if(imtrx[i + j * pdb.size()]) outf << pdb[i].DD_DISTANCES.at(j);
        else outf << "0.0";
        if (j != (pdb.size()-1)) outf << "\t";
      }
      outf << "\n";
    }

  }
  else{
    for (int i = 0; i < pdb.size(); ++i){
      outf << pdb[i].sID();
      if( i != (pdb.size()-1)) outf << "\t";
    }
    outf << "\n";
    for (int i = 0; i < pdb.size(); ++i){
      outf << pdb[i].sID() << "\t";
      for (int j = 0; j < pdb.size(); ++j){
        if(imtrx[i + j * pdb.size()]) outf << pdb[i].DD_DISTANCES.at(j);
        else outf << "0.0";
        if (j != (pdb.size()-1)) outf << "\t";
      }
      outf << "\n";
    }
  }
  outf.close();
}

void 
PrintTMatrix(vector<int>& v,
                      int            x,
                      int            y,
                      int            z,
                      string         output){
  ofstream outf(output);
  outf << "x\t";
  if (x > 0){
    if (y > 0){
      if (z > 0){
        
      }
      else{
        for (int i=0; i < x; ++i){
          outf << i+1;
          if (i == x-1) outf << "\n";
          else outf << "\t";
        }
        for (int j=0; j < y; ++j){
          outf << j+1 << "\t";
          for (int i=0; i < x; ++i){
            outf << v[i + j*x];
            if (i == x-1) outf << "\n";
            else outf << "\t";
          }
        }
      }
    }
  }
  outf.close();
}

void 
PrintTMatrix(vector<float>& v,
                      int            x,
                      int            y,
                      int            z,
                      string         output){
  ofstream outf(output);
  outf << "x\t";
  if (x > 0){
    if (y > 0){
      if (z > 0){
        
      }
      else{
        for (int i=0; i < x; ++i){
          outf << i+1;
          if (i == x-1) outf << "\n";
          else outf << "\t";
        }
        for (int j=0; j < y; ++j){
          outf << j+1 << "\t";
          for (int i=0; i < x; ++i){
            outf << v[i + j*x];
            if (i == x-1) outf << "\n";
            else outf << "\t";
          }
        }
      }
    }
  }
  outf.close();
}


void
PrintDNA_ProtResultsByAtomMatrix(vector<atom_struct>& pdb,     // pdb struct
                                 string               output, // output filename
                                 int                  mode){  
  CalculateDNA_ProtInteractions(pdb);
  set<string> object_types;
  for (auto& atom : pdb){
    object_types.insert(atom.STRUCT_TYPE);
  }
  vector<string> objs(object_types.begin(),object_types.end());

//pragma omp parallel for schedule(dynamic)
  for(uint32 i = 0; i < objs.size()-1; ++i){
    for (uint32 j = i+1; j < objs.size(); ++j){
      float  objA_bsa = 0.0;
      float  objB_bsa = 0.0;
      float  objA_sasa = 0.0;
      float  objB_sasa = 0.0;
      string objA = objs[i];
      string objB = objs[j];
      //cout << "OBJS: "<< objA << " " << objB << "\n";

      vector<string> ColL;
      vector<string> LineL;
      vector<uint32> ColL_type;
      vector<uint32> LineL_type;
      map<string,uint32> ColM;
      map<string,uint32> LineM;
      vector<float> IcJ; //J causes to I
      vector<float> JcI; //I causes to J

      vector<string> ColLres;
      vector<string> LineLres;
      map<string,uint32> ColMres;
      map<string,uint32> LineMres;
      vector<float> IcJres; //J causes to I
      vector<float> JcIres; //I causes to J
      
      map<string,uint32> oriIDX;

      for (uint32 i = 0; i < pdb.size(); ++i){
        auto& atom = pdb[i];
        string aID = atom.sID();
        string rID = atom.rsID();
        oriIDX[aID] = i;
        if(atom.STRUCT_TYPE == objA){
          ColL.push_back(aID);
          ColL_type.push_back(atom.ATOM_TYPE);
          ColM[aID] = ColL.size()-1;
          ColLres.push_back(rID);
          //cout << atom.print() << std::endl;
          objA_sasa += atom.SASA;
          //cout << " COL " << aID << " " << ColL.size()-1 << std::endl;
        }
        if(atom.STRUCT_TYPE == objB){
          LineL.push_back(aID);
          LineL_type.push_back(atom.ATOM_TYPE);
          LineM[aID] = LineL.size()-1;
          LineLres.push_back(rID);
          objB_sasa += atom.SASA ;
          //cout << " LIN " << aID << " " << LineL.size()-1 << std::endl;
        }
      }
      ColLres.erase( unique( ColLres.begin(), ColLres.end() ), ColLres.end() );
      LineLres.erase( unique( LineLres.begin(), LineLres.end() ), LineLres.end() );

      for (uint32 k = 0; k < ColLres.size();++k) ColMres[ColLres[k]] = k;
      for (uint32 k = 0; k < LineLres.size();++k) LineMres[LineLres[k]] = k;

      IcJres.resize(ColLres.size()*LineLres.size(),0.0);
      JcIres.resize(ColLres.size()*LineLres.size(),0.0);

      IcJ.resize(ColL.size()*LineL.size(),0.0);
      JcI.resize(ColL.size()*LineL.size(),0.0);
      
      
      uint32 lm = ColL.size();
      uint32 lm_res = ColLres.size();

      for(auto& atom : pdb){
        string aID = atom.sID();
        string rID = atom.rsID();
        if(atom.STRUCT_TYPE == objA){
          //cout << "ATOM: " << aID << " " << atom.STRUCT_TYPE << "\n";
          for(uint32 pos : atom.INTERACTION_SASA_P){
            if (pdb[pos].STRUCT_TYPE != objB) continue;
            //cout <<"    OBJ_A: " << pdb[pos].sID() << "\n";
            for (uint32 r = 0; r < atom.ov_table.size(); ++r){ //iterate through this atoms overlap table
              auto& ov = atom.ov_table[r];                     //get vector of atom positions that form this overlap
              if(ov.end() != find(ov.begin(),ov.end(),pos)){   //check if this other atom is part of current overlap
                uint32 col = ColM[aID];                        //current atom position in column
                uint32 line = LineM[pdb[pos].sID()];           //other atom position in line
                IcJ[col + line * lm] += atom.ov_norm_area[r];  //sum area to atom matrix
                objA_bsa += atom.ov_norm_area[r];
                uint32 col_res = ColMres[rID];                 //column of current residue 
                uint32 line_res = LineMres[pdb[pos].rsID()];   //line of current residue
                IcJres[col_res + line_res * lm_res] += atom.ov_norm_area[r]; //sum area to residue matrix
              }
            }
          }
        }
        //same comments as object A
        if(atom.STRUCT_TYPE == objB){
          //cout << "ATOM: " << atom.sID() << " " << atom.STRUCT_TYPE << "\n";
          for(uint32 pos : atom.INTERACTION_SASA_P){
            if (pdb[pos].STRUCT_TYPE != objA) continue;
            //cout <<"    OBJ_B: " << pdb[pos].sID() << "\n";
            for (uint32 r = 0; r < atom.ov_table.size(); ++r){
              auto& ov = atom.ov_table[r];
              if(ov.end() != find(ov.begin(),ov.end(),pos)){
                uint32 col = ColM[pdb[pos].sID()];
                uint32 line = LineM[aID];
                JcI[col + line * lm] += atom.ov_norm_area[r];
                objB_bsa += atom.ov_norm_area[r];
                uint32 col_res = ColMres[pdb[pos].rsID()];
                uint32 line_res = LineMres[rID];
                JcIres[col_res + line_res * lm_res] += atom.ov_norm_area[r];
              }
            }
          }
        }
      }
      if((accumulate(IcJ.begin(),IcJ.end(),0.0) == 0.0) &&
         (accumulate(JcI.begin(),JcI.end(),0.0) == 0.0)){
       //pragma omp critical(log)
       // {
          //cerr << "NULL_INTERACTION\tOBJ\t"<< objA << "\t<->\tOBJ\t"<<objB<<"\n";
        //}
        continue;
      }
      //atom files
      stringstream AcB;                //dSASA in A caused by B outputfilename
      AcB << output << "." <<objA << "_vs_" << objB << ".by_atom.tsv";
      stringstream BcA;                //dSASA in B caused by A outputfilename
      BcA << output << "." <<objB << "_vs_" << objA << ".by_atom.tsv";
      ofstream AcBfile(AcB.str());
      ofstream BcAfile(BcA.str());

      AcBfile << objA << " <- " << objB << "\t"; //dSASA B causes to A
      BcAfile << objA << " -> " << objB << "\t"; //dSASA A causes to B
      
      for (uint s = 0; s < ColL.size(); ++s) {

          AcBfile << ColL[s] << "/" << pdb[oriIDX[ColL[s]]].MOL_TYPE << "/" << pdb[oriIDX[ColL[s]]].ATOM_TYPE << "\t";
          BcAfile << ColL[s] << "/" << pdb[oriIDX[ColL[s]]].MOL_TYPE << "/" << pdb[oriIDX[ColL[s]]].ATOM_TYPE << "\t";
      }
      AcBfile << std::endl;
      BcAfile << std::endl;
      for (uint32 l = 0; l < LineL.size(); ++l){
          AcBfile << LineL[l] << "/" << pdb[oriIDX[LineL[l]]].MOL_TYPE << "/" << pdb[oriIDX[LineL[l]]].ATOM_TYPE << "\t";
          BcAfile << LineL[l] << "/" << pdb[oriIDX[LineL[l]]].MOL_TYPE << "/" << pdb[oriIDX[LineL[l]]].ATOM_TYPE << "\t";

        for (uint32 c = 0; c < ColL.size(); ++c){
          AcBfile << IcJ[c + l * ColL.size()] << "\t";
          BcAfile << JcI[c + l * ColL.size()] << "\t";
        }
        AcBfile << std::endl;
        BcAfile << std::endl;
      }
      AcBfile.close();
      BcAfile.close();

      //residuefiles
      stringstream AcBres;                //dSASA in A caused by B outputfilename
      AcBres << output << "." <<objA << "_vs_" << objB << ".by_res.tsv";
      stringstream BcAres;                //dSASA in B caused by A outputfilename
      BcAres << output << "." <<objB << "_vs_" << objA << ".by_res.tsv";
      ofstream AcBfileres(AcBres.str());
      ofstream BcAfileres(BcAres.str());

      AcBfileres << objA << "<-" << objB << "\t"; //dSASA B causes to A
      BcAfileres << objA << "->" << objB << "\t"; //dSASA A causes to B
      for (auto s : ColLres) {
        AcBfileres << s << "\t";
        BcAfileres << s << "\t";
      }
      AcBfileres << std::endl;
      BcAfileres << std::endl;
      for (uint32 l = 0; l < LineLres.size(); ++l){
        AcBfileres << LineLres[l] << "\t";
        BcAfileres << LineLres[l] << "\t";
        for (uint32 c = 0; c < ColLres.size(); ++c){
          AcBfileres << IcJres[c + l * ColLres.size()] << "\t";
          BcAfileres << JcIres[c + l * ColLres.size()] << "\t";
        }
        AcBfileres << std::endl;
        BcAfileres << std::endl;
      }
      AcBfileres.close();
      BcAfileres.close();
      //log to stdout
      cout << "Object " << objA << " complexed surface (A^2):\t" << objA_sasa << std::endl;
      cout << "Object " << objB << " complexed surface (A^2):\t" << objB_sasa << std::endl; 
      cout << "Object " << objA << " uncomplexed surface (A^2):\t" << (objA_sasa+objA_bsa) << std::endl;
      cout << "Object " << objB << " uncomplexed surface (A^2):\t" << (objB_sasa+objB_bsa) << std::endl; 
      cout << objA << " <--- " << objB << " buried surface (A^2):\t" << objA_bsa  << std::endl;
      cout << objA << " ---> " << objB << " buried surface (A^2):\t" << objB_bsa  << std::endl;
      cout << "Interface "<< objA << "/" << objB << " (A^2):\t" << ((objA_bsa+objB_bsa) / 2.0)  << std::endl;
    }
  }
}

void
Print_MatrixInsideAtom(vector<atom_struct>& pdb,
                       string               output){
  CalculateDNA_ProtInteractions(pdb);
  set<string> object_types;
  for (auto& atom : pdb){
    object_types.insert(atom.CHAIN);
  }
  vector<string> objs(object_types.begin(),object_types.end());
  
  
  //variables for atom matrix
  vector<string> ColL;
  vector<string> LineL;
  vector<string> ColL_type;
  vector<string> LineL_type;
  map<string,uint32> ColM;
  map<string,uint32> LineM;
  vector<float> mtrx; 
  //variables for residue matrix
  vector<string> ColLres;
  vector<string> LineLres;
  map<string,uint32> ColMres;
  map<string,uint32> LineMres;
  vector<float> mtrx_res; 

  for (auto& atom : pdb){
    string aID = atom.sID();
    string rID = atom.rsID();
    stringstream atomtype;
    atomtype << aID << "/" << atom.MOL_TYPE << "/" << atom.ATOM_TYPE; 
    ColL_type.push_back(atomtype.str());
    ColL.push_back(aID);
    ColM[aID] = ColL.size()-1;
    ColLres.push_back(rID);
    //cout << " COL " << aID << " " << ColL.size()-1 << std::endl;
    LineL_type.push_back(atomtype.str());
    LineL.push_back(aID);
    LineM[aID] = LineL.size()-1;
    LineLres.push_back(rID);
    //cout << " LIN " << aID << " " << LineL.size()-1 << std::endl;
  }
  ColLres.erase( unique( ColLres.begin(), ColLres.end() ), ColLres.end() );
  LineLres.erase( unique( LineLres.begin(), LineLres.end() ), LineLres.end() );
  for (uint32 k = 0; k < ColLres.size();++k) ColMres[ColLres[k]] = k;
  for (uint32 k = 0; k < LineLres.size();++k) LineMres[LineLres[k]] = k;
  
  mtrx.resize(ColL.size()*LineL.size(),0.0);
  mtrx_res.resize(ColLres.size()*LineLres.size(),0.0);
  uint32 lm = ColL.size();
  uint32 lm_res = ColLres.size();
//pragma omp parallel for schedule(dynamic)
  for(uint32 i = 0; i < pdb.size(); ++i){
    auto& atom = pdb[i];
    string aID = atom.sID();
    string rID = atom.rsID();
   //cout << "ATOM: " << aID << " " << atom.STRUCT_TYPE << "\n";
   for(uint32 pos : atom.INTERACTION_SASA_P){
      //cout <<"    OBJ_A: " << pdb[pos].sID() << "\n";
      for (uint32 r = 0; r < atom.ov_table.size(); ++r){
        auto& ov = atom.ov_table[r];
        if(ov.end() != find(ov.begin(),ov.end(),pos)){
          uint32 col = ColM[aID];
          uint32 line = LineM[pdb[pos].sID()];
          uint32 col_res = ColMres[rID];
          uint32 line_res = LineMres[pdb[pos].rsID()];
          //pragma omp critical(printatommatrix)
          {
          mtrx_res[col_res + line_res * lm_res] += atom.ov_norm_area[r];
          mtrx[col + line * lm] += atom.ov_norm_area[r];
          }
        }
      }
    }
  }
  //atomfiles

/*  stringstream fname;
  fname << output << "." << "internal_chain_";
  for (auto& str : objs) fname << str;
  fname << "_by_atom.tsv";
  ofstream mtrxfile(fname.str());
  mtrxfile << "ATOM-v" << "\t";
  for (uint32 s = 0; s < ColL.size(); ++s) {
    mtrxfile << ColL[s] << "\t";
  }
  mtrxfile << endl;
  for (uint32 l = 0; l < LineL.size(); ++l){
    mtrxfile << LineL[l] << "\t";
    for (uint32 c = 0; c < ColL.size(); ++c){
      mtrxfile << mtrx[c + l * ColL.size()] << "\t";
    }
    mtrxfile << endl;
  }
  mtrxfile.close();*/
  
    //atomfiles by type

  stringstream fname_type;
  fname_type << output << ".matrix.";
  for (auto& str : objs) fname_type << str;
  fname_type << ".by_atom.tsv";
  ofstream mtrxfile_type(fname_type.str());
  mtrxfile_type << "ATOM-v" << "\t";
  for (uint32 s = 0; s < ColL.size(); ++s) {

      mtrxfile_type << ColL_type[s] << "\t";
    
  }
  mtrxfile_type << endl;
  for (uint32 l = 0; l < LineL.size(); ++l){

      mtrxfile_type << LineL_type[l] << "\t";
    
    for (uint32 c = 0; c < ColL.size(); ++c){
      mtrxfile_type << mtrx[c + l * ColL.size()] << "\t";
    }
    mtrxfile_type << endl;
  }
  mtrxfile_type.close();
  
  //residuefiles
  
  stringstream fname_res;
  fname_res << output << ".matrix.";
  for (auto& str : objs) fname_res << str; 
  fname_res << ".by_res.tsv";
  ofstream mtrxfile_res(fname_res.str());
  mtrxfile_res << "RESIDUE-v" << "\t"; 
  for (auto s : ColLres) {
    mtrxfile_res << s << "\t";
  }
  mtrxfile_res << endl;
  for (uint32 l = 0; l < LineLres.size(); ++l){
    mtrxfile_res << LineLres[l] << "\t";
    for (uint32 c = 0; c < ColLres.size(); ++c){
      mtrxfile_res << mtrx_res[c + l * ColLres.size()] << "\t";
    }
    mtrxfile_res << endl;
  }
  mtrxfile_res.close();
}

void
PrintDNA_ProtResults(vector<atom_struct>& pdb,
                     string               output){

  ofstream outf(output, ofstream::out);
  outf << "//ATOM\t<ID>\t<NAME>\t<RESN>\t<CHAIN>\t<RESI>\t<STRUCT_TYPE>\t<X>\t<Y>\t<Z>\t<SASA>\t<dSASA>\t<VDW>" << std::endl;
  //outf << "// ATOM_CONTACT_SURFACE <ID> <NAME> <RESN> <CHAIN> <RESI> <STRUCT TYPE> <AREA BURIED BY OTHER MOL TYPES>" << std::endl;
  outf << "//ATOM_CONTACT_SURFACE_DETAIL\t<ID>\t<NAME>\t<RESN>\t<CHAIN>\t<RESI>\t<STRUCT_TYPE>\t<AREA_BURIED_BY_THIS_ATOM>\t<DISTANCE_FROM_THIS_ATOM>" << std::endl;
  outf << "//ATOM_CONTACT_SURFACE_OVERLAPS\t<COUNT_OF_OVERLAPPING_ATOMS>\t<ID_1>\t<ID_2>\t...\t<OVERLAPPING_AREA>" << std::endl;
  //outf << "// The ATOM_CONTACT_SURFACE* tags are always in order, and their information belongs to the last ATOM tag." << std::endl;

  for (uint32 i = 0; i < pdb.size(); ++i){
    auto& atom_i = pdb[i];
    if(!atom_i.ACTIVE) continue;
    if (atom_i.HETATM) outf << "HETATM\t";
    else outf << "ATOM\t";
   
    outf << atom_i.ID << "\t" << atom_i.NAME << "\t" << atom_i.RESN << "\t" << atom_i.CHAIN << "\t" << atom_i.RESI << "\t" << atom_i.MOL_TYPE << "\t"<<
    atom_i.COORDS[0] <<  "\t" << atom_i.COORDS[1] << "\t" << atom_i.COORDS[2] << "\t" <<
    atom_i.SASA << "\t" << atom_i.EXT1  << "\t" << atom_i.VDW << std::endl;
    if (atom_i.AREA_BURIED_BY_ATOM_area.size() != 0){
      outf << PrintDNA_ProtInteractions(pdb, i);
    }

  }

  outf.close();
}

void
PrintSASAResults(vector<atom_struct>& pdb,
                 string               output){
  float total = 0;
  float mw = 0;
  ofstream outf(output, ofstream::out);
    outf << "REMARK 000 SASA SOLVER" << std::endl;
  map<string,float> aw = {{"C",12.011},{"H",1.0},{"O",15.999},{"N",14.007},{"P",30.973},
                          {"S",32.06},{"F",19.0},{"Br",80.0 },{"BR",80.0 },{"Cl",35.45},{"CL",35.45}};

  for (uint32 i = 0; i<pdb.size(); ++i){
    auto& atom_i = pdb[i];
    if(!atom_i.ACTIVE) continue;
    if(atom_i.HETATM)  outf << "HETATM";
    else               outf << "ATOM  ";
    outf << std::right << std::setw(5) << atom_i.ID;
    if(atom_i.NAME.size() == 4){
      outf << " ";
      outf << std::left << std::setw(4) << atom_i.NAME;
    }else{
      outf << "  ";
      outf << std::left << std::setw(3) << atom_i.NAME;
    }
    
    outf << std::right << std::setw(4) << atom_i.RESN;
    outf << " ";
    outf << atom_i.CHAIN;
    string resi = std::to_string(atom_i.RESI)+ atom_i.iCODE;
    outf << std::setw(4) << resi; 
    outf << "    ";    
    //atom_i.MOL_TYPE  << " ";
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[0];
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[1];
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[2];
    //outf << "  1.00";
    outf << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom_i.VDW;
    outf << std::setw(6) << std::fixed << std::setprecision(2) << atom_i.SASA; 
    //outf << "          ";
    outf <<" ";
    outf << std::left << std::setw(8) << atom_i.MOL_TYPE;
    outf << std::right << std::setw(3) << atom_i.ELEMENT;
    outf << std::setw(2) << atom_i.ATOM_TYPE << std::endl;
    try{
      mw += aw.at(atom_i.ELEMENT);
    }
    catch(exception e){
    }
    total += atom_i.SASA;
  }
  //float Arel = total / (4.44*pow(mw,0.77));

  outf << "REMARK 001 TOTAL SURFACE AREA (A^2)\t" << total << std::endl;
  cout << "Selected complex surface (A^2):\t" << total << std::endl;
  outf << "REMARK 002 MOLECULAR WEIGHT (DALTON)\t" << mw << std::endl;
  //cout << "MW: " << mw << "\n";
//  outf << "REMARK 002  " << Arel << std::endl;
//  cout << "AREL: " << Arel << "\n";
  outf.close();
}
void
PrintDSASAResults(vector<atom_struct>& pdb,
                 string               output){
  float total = 0;
  float mw = 0;
  ofstream outf(output, ofstream::out);
    outf << "REMARK 000 BSA SOLVER" << std::endl;
  map<string,float> aw = {{"C",12.011},{"H",1.0},{"O",15.999},{"N",14.007},{"P",30.973},
                          {"S",32.06},{"F",19.0},{"Br",80.0 },{"BR",80.0 },{"Cl",35.45},{"CL",35.45}};

  for (uint32 i = 0; i<pdb.size(); ++i){
    auto& atom_i = pdb[i];
    if(!atom_i.ACTIVE) continue;
    if(atom_i.HETATM)  outf << "HETATM";
    else               outf << "ATOM  ";
    outf << std::right << std::setw(5) << atom_i.ID;
    if(atom_i.NAME.size() == 4){
      outf << " ";
      outf << std::left << std::setw(4) << atom_i.NAME;
    }else{
      outf << "  ";
      outf << std::left << std::setw(3) << atom_i.NAME;
    }
    
    outf << std::right << std::setw(4) << atom_i.RESN;
    outf << " ";
    outf << atom_i.CHAIN;
    string resi = std::to_string(atom_i.RESI)+ atom_i.iCODE;
    outf << std::setw(4) << resi; 
    outf << "    ";    
    //atom_i.MOL_TYPE  << " ";
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[0];
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[1];
    outf << std::setw(8) << std::fixed << std::setprecision(3) << atom_i.COORDS[2];
    //outf << "  1.00";
    outf << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom_i.VDW;
    outf << std::setw(6) << std::fixed << std::setprecision(2) << atom_i.EXT0; 
    //outf << "          ";
    outf <<" ";
    outf << std::left << std::setw(8) << atom_i.MOL_TYPE;
    outf << std::right << std::setw(3) << atom_i.ELEMENT;
    outf << std::setw(2) << atom_i.ATOM_TYPE << std::endl;
    try{
      mw += aw.at(atom_i.ELEMENT);
    }
    catch(exception e){
    }
    total += atom_i.EXT0;
  }

  //cout << "MW: " << mw << "\n";
//  outf << "REMARK 002  " << Arel << std::endl;
//  cout << "AREL: " << Arel << "\n";
  outf.close();
}

void
PrintSASAResults_per_type(vector<atom_struct>& pdb,         //pdb struct
              string output1,          //protein SASA
              string output2){                  //DNA SASA

  map<uint32,vector<float> > protvalues;
  map<uint32,vector<float> > dnavalues;
  for (auto& atom : pdb){
    if(atom.STRUCT_TYPE == "PROTEIN"){
      if (atom.ATOM_TYPE != 0)
      protvalues[atom.ATOM_TYPE].push_back(atom.SASA);
    }
    else{
      if (atom.ATOM_TYPE != 0)
      dnavalues[atom.ATOM_TYPE].push_back(atom.SASA);
    }
  }

  ofstream protf(output1, ofstream::out);
  for (auto it : protvalues){
    if (it.second.size() == 0) continue;
    protf << it.first << "_" << 0 << "\t";

    for (auto sasa : it.second){
      protf << sasa << "\t";
    }
    protf << std::endl;
  }
  protf.close();

  ofstream dnaf(output2, ofstream::out);
  for (auto it : dnavalues){
    if (it.second.size() == 0) continue;
    dnaf << it.first << "_" << 0 << "\t";
    for (auto sasa : it.second){
      dnaf << sasa << "\t";
    }
    dnaf << std::endl;
  }
  dnaf.close();
}


void
PrintSASAResults_per_type_and_res(vector<atom_struct>& pdb,         //pdb struct
                 string                output1,
                     string                output2){
  ofstream protf(output1, ofstream::out);
  ofstream dnaf(output2, ofstream::out);
  for (auto& ATOM : pdb){
    if(ATOM.STRUCT_TYPE == "PROTEIN"){
      if (ATOM.ATOM_TYPE != 0)
        protf << ATOM.RESI << ATOM.iCODE << "\t" << ATOM.ATOM_TYPE << "\t" << ATOM.SASA << std::endl;
    }
    else{
      if (ATOM.ATOM_TYPE != 0)
        dnaf << ATOM.RESI << ATOM.iCODE << "\t" << ATOM.ATOM_TYPE << "\t" << ATOM.SASA << std::endl;
    }
  }
  protf.close();
  dnaf.close();
}

/*void
PrintSASAone_type_(vector<atom_struct>& pdb,         //pdb struct
           int type,                           //type
           string output){
  ofstream outf(output,ofstream::out);
  for (auto& ATOM : pdb){
    if(ATOM.STRUCT_TYPE == "PROTEIN" && ATOM.ATOM_TYPE == type){
      outf << ATOM.
    }
  }
}*/

void
CalculateDNA_ProtInteractions(vector<atom_struct>& pdb){
//pragma omp parallel for schedule(dynamic)
  for (uint32 pos = 0; pos < pdb.size(); ++pos){
    auto& atom_i = pdb[pos];
    vector<char> valid_overlaps;
    valid_overlaps.resize(atom_i.AREA_BURIED_BY_ATOM_vector.size(), true);
    vector<uint32> o_atoms;
    vector<uint32> interac_pos;
    auto& ov_table = atom_i.ov_table;
    auto& ov_table_area = atom_i.ov_table_area;
    auto& ov_norm_area = atom_i.ov_norm_area;
    ov_table.clear();
    ov_table_area.clear();
    ov_norm_area.clear();
    for (uint32 i = 0; i < atom_i.INTERACTION_P.size(); ++i){
      if (atom_i.STRUCT_TYPE != pdb[atom_i.INTERACTION_P[i]].STRUCT_TYPE){
        interac_pos.push_back(atom_i.INTERACTION_P[i]);
      }
    }
    for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i){
      for (uint32 k = 0; k < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++k){
        uint32 other = atom_i.AREA_BURIED_BY_ATOM_vector[i][k];
        if (interac_pos.end() == find(interac_pos.begin(), interac_pos.end(), other)) valid_overlaps[i] = false;
      }
      if (!valid_overlaps[i]) continue;
      atom_i.AREA_BURIED_BY_ATOM_vector_valid.push_back(i);
      for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++j){
        if (o_atoms.end() == find(o_atoms.begin(), o_atoms.end(), atom_i.AREA_BURIED_BY_ATOM_vector[i][j])){
          o_atoms.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i][j]);
        }
      }
    }
    sort(o_atoms.begin(),o_atoms.end());
    //cout << atom_i.NAME <<  " " << atom_i.RESN << " " << atom_i.RESI << " " << atom_i.CHAIN << " " << atom_i.EXT0 << " " << atom_i.EXT1 << std::endl;
    for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i){
      uint32 overlaps = atom_i.AREA_BURIED_BY_ATOM_vector[i].size();
      if (valid_overlaps[i]){
        vector<uint32> overlap;
        for (uint32 j = 0; j < overlaps; ++j){
          overlap.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i][j]);
          //cout << atom_i.AREA_BURIED_BY_ATOM_vector[i][j] << " ";
        }
        ov_table.push_back(overlap);
        ov_table_area.push_back(atom_i.AREA_BURIED_BY_ATOM_area[i]);
        ov_norm_area.push_back((float)atom_i.AREA_BURIED_BY_ATOM_area[i]/(float)overlap.size());
        //cout << atom_i.AREA_BURIED_BY_ATOM_area[i] << std::endl;
      }
    }
  }
}

string
PrintDNA_ProtInteractions(vector<atom_struct>& pdb,
                          uint32               pos){

  stringstream Result;
  vector<uint32> o_atoms;
  vector<uint32> interac_pos;
  auto& atom_i = pdb[pos];
  vector<char> valid_overlaps(atom_i.AREA_BURIED_BY_ATOM_vector.size(), true);


  for (uint32 i = 0; i < atom_i.INTERACTION_P.size(); ++i){
    if (atom_i.STRUCT_TYPE != pdb[atom_i.INTERACTION_P[i]].STRUCT_TYPE){
      interac_pos.push_back(atom_i.INTERACTION_P[i]);
    }
  }


  for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i){
    for (uint32 k = 0; k < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++k){
      uint32 other = atom_i.AREA_BURIED_BY_ATOM_vector[i][k];
      if (interac_pos.end() == find(interac_pos.begin(), interac_pos.end(), other)) valid_overlaps[i] = false;
    }
    if (!valid_overlaps[i]) continue;
    atom_i.AREA_BURIED_BY_ATOM_vector_valid.push_back(i);
    for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++j){
      if (o_atoms.end() == find(o_atoms.begin(), o_atoms.end(), atom_i.AREA_BURIED_BY_ATOM_vector[i][j])){
        o_atoms.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i][j]);
      }
    }
  }

  sort(o_atoms.begin(),o_atoms.end());

    for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_area.size(); ++i){
    if (valid_overlaps[i])  atom_i.EXT0 += atom_i.AREA_BURIED_BY_ATOM_area[i];
  }
  //Result << "ATOM_CONTACT_SURFACE\t" << atom_i.ID << "\t" << atom_i.NAME << "\t" << atom_i.RESN << "\t" << atom_i.CHAIN << "\t" << atom_i.RESI << "\t"<< atom_i.STRUCT_TYPE << "\t" << atom_i.EXT0 << std::endl;

  for (uint32 i = 0; i < o_atoms.size(); ++i){
    float T_area = 0;
    uint32 ATOM_POS = o_atoms[i];
    auto& atom_j = pdb[ATOM_POS];
    for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++j){
      if (valid_overlaps[j] && atom_i.AREA_BURIED_BY_ATOM_vector[j].end() != find(atom_i.AREA_BURIED_BY_ATOM_vector[j].begin(), atom_i.AREA_BURIED_BY_ATOM_vector[j].end(), ATOM_POS)){
        T_area += atom_i.AREA_BURIED_BY_ATOM_area[j];
      }
    }
    Result << "ATOM_CONTACT_SURFACE_DETAIL\t" << atom_j.ID << "\t" << atom_j.NAME << "\t" << atom_j.RESN << "\t" << atom_j.CHAIN << "\t" << atom_j.RESI << "\t"<< atom_j.STRUCT_TYPE <<"\t" << T_area << "\t" << atom_i.DISTANCES[ATOM_POS] << std::endl;
  }

  for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i){
    uint32 overlaps = atom_i.AREA_BURIED_BY_ATOM_vector[i].size();
    if (valid_overlaps[i]){
      Result << "ATOM_CONTACT_SURFACE_OVERLAP\t" << overlaps << "\t";
      for (uint32 j = 0; j < overlaps; ++j){
        Result << pdb[atom_i.AREA_BURIED_BY_ATOM_vector[i][j]].ID << "\t";
      }
      Result << atom_i.AREA_BURIED_BY_ATOM_area[i] << std::endl;
    }
  }
  return Result.str();
}

void
PrintDNA_ProtResultsTable(vector<atom_struct>& pdb,
                          string               pair_table){

    vector<pair<uint32,uint32> > pairs;
    ofstream outf(pair_table, ofstream::out);

    vector<float> diff_p;
    vector<float> abs_diff;

    outf << "// ATOM INTERACTION PAIRS" << std::endl;
    //outf << "// <2->1> means 'buried surface in atom 1 caused by atom 2'" << std::endl;
    //outf << "// <1->2> means 'buried surface in atom 2 caused by atom 1'" << std::endl;
  outf << "// <ID_1>\t<NAME_1>\t<RESN_1>\t<CHAIN_1>\t<RESI_1>\t<ID_2>\t<NAME_2>\t<RESN_2>\t<CHAIN_2>\t<RESI_2>" << std::endl;


  for (uint32 i = 0; i < pdb.size(); ++i){
      auto& atom_i = pdb[i];
      if(!atom_i.ACTIVE ||
          atom_i.EXT1 == 0) continue;
      for (uint32 j = 0; j < atom_i.INTERACTION_SASA_P.size(); ++j){
          auto j_pos = atom_i.INTERACTION_SASA_P[j];
          auto& atom_j = pdb[j_pos];
          pair<uint32,uint32> p1(i,j_pos);
          pair<uint32,uint32> p2(j_pos,i);
          if(!atom_j.ACTIVE ||
             atom_j.EXT1 == 0)
             continue;
          if(pairs.end() != find(pairs.begin(),pairs.end(),p1) ||
             pairs.end() != find(pairs.begin(),pairs.end(),p2) )
             continue;
          pairs.push_back(p1);
          pairs.push_back(p2);
          float j_i = atom_i.CONTACT_AREA.at(atom_j.ID);
          float i_j = atom_j.CONTACT_AREA.count(atom_i.ID) ? atom_j.CONTACT_AREA.at(atom_i.ID) : 0.0 ;
          float diff = abs(j_i - i_j);
          float percent = 100 * diff / (i_j < j_i ? i_j : j_i);

          diff_p.push_back(percent);

          abs_diff.push_back(diff);

          outf << atom_i.ID << "\t" << atom_i.NAME << "\t" << atom_i.RESN << "\t" << atom_i.CHAIN << "\t" << atom_i.RESI << atom_i.iCODE << "\t" <<
                  atom_j.ID << "\t" << atom_j.NAME << "\t" << atom_j.RESN << "\t" << atom_j.CHAIN << "\t" << atom_j.RESI << atom_j.iCODE << "\t" << std::endl;

      }
  }

  /*ofstream stat_f(stats);
  //stat_f << "<abs_diff> <percentage_diff>" << std::endl;
  stat_f << "<abs_diff>" << std::endl;
  for(uint32 i = 0; i < diff_p.size(); ++i){
     // stat_f << abs_diff[i] << "\t" << diff_p[i] << std::endl;
      stat_f << abs_diff[i] << std::endl;
  }
    stat_f.close();*/
  outf.close();

}

void
GeneratePDistribution(map<pair<int,int>,vector<float> >&     data,
                        map<pair<int,int>,histogram>&         hst,
                        uint32                                bins,
                        float                                xmax,
                        string                                output){
    vector<float> gdsasa;
    for(auto it = data.begin(); it != data.end(); ++it){
        auto p = it->first;
        auto& v = it->second;
        gdsasa.insert(gdsasa.end(),v.begin(),v.end());
    histogram c_hst(bins,0,xmax);
        for (auto i : v){
            c_hst.inc(i);
        }
        c_hst.Pdist();
    hst[p] = c_hst;
    }
    histogram Pr(bins,0,xmax);
    for (auto i : gdsasa){
         Pr.inc(i);
    }
    Pr.Pdist();
    hst[make_pair(0,0)] = Pr;

    ofstream outf(output,ofstream::out);
    for (auto& it : hst){
        auto& p = it.first;
        auto& hst = it.second;
        outf << p.first << "_" << p.second << "_range: ";
        for(auto rng : hst.range){
           outf << rng << " ";
        }
        outf << std::endl << p.first << "_" << p.second << "_bins: ";
        for(auto p : hst.p){
           outf << p << " ";
        }
        outf << std::endl;
    }
    outf.close();
}

// void
// AddRawAAData(map<string,vector<float> >& table,
//              vector<atom_struct>&         pdb ){
// 
//     map<uint32,residue_struct> res_dsasa;
// 
//     uint32 pdb_size = pdb.size();
//     for (uint32 i = 0; i < pdb_size; ++i){
//         auto& atom_i = pdb[i];
//         if(!atom_i.ACTIVE) continue;
//         if(res_dsasa.find(atom_i.RESI) == res_dsasa.end()){
//             res_dsasa[atom_i.RESI] = residue_struct(atom_i.RESN,atom_i.CHAIN,atom_i.STRUCTURE,atom_i.STRUCT_TYPE,atom_i.RESI);
//         }
//         res_dsasa.at(atom_i.RESI).SASA+=atom_i.SASA;
//         res_dsasa.at(atom_i.RESI).dSASA+= atom_i.EXT1;
//     }
// 
//   for (auto it = res_dsasa.begin(); it != res_dsasa.end(); ++it) {
//     auto res = it->second;
//     if(res.dSASA > 0) table[res.RESN].push_back(res.dSASA);
//   }
// }

void
AddRawAAxNUCData(map<pair<int,int>,vector<float> >& data,
         map<pair<int,int>,vector<float> >& data_w,
                 vector<atom_struct>&                pdb,
         vector<int>&                     interac){

  map<string, int > RES_TYPE = {
    {"ALA",  1 },
    {"ARG",  2 },
    {"ASN",  3 },
    {"ASP",  4 },
    {"CYS",  5 },
    {"GLN",  6 },
    {"GLU",  7 },
    {"GLY",  8 },
    {"HIS",  9 },
    {"ILE", 10 },
    {"LEU", 11 },
    {"LYS", 12 },
    {"MET", 13 },
    {"PHE", 14 },
    {"PRO", 15 },
    {"SER", 16 },
    {"THR", 17 },
    {"TRP", 18 },
    {"TYR", 19 },
    {"VAL", 20 },
    {"DA", 1 },
    {"DT", 2 },
    {"DC", 3 },
    {"DG", 4 },
    {"A" , 5},
    {"U" , 6},
    {"C" , 7},
    {"G" , 8}};

  auto GetNucFromPosVec = [&](vector<uint32> atoms){
                            set<tuple<string,string,int,string> > nucs;
                for(auto pos : atoms){
                  auto& atom = pdb[pos];
                  auto nuc_tuple = make_tuple(atom.RESN,atom.CHAIN,atom.RESI,atom.iCODE);
                  nucs.insert(nuc_tuple);
                    }
                return nucs;
               };
  //residues and nucleotides are identified by a (resn,chain,resi,icode) tuple

    map<tuple<string,string,int,string>, vector<uint32> > res_vec;                           //residue id and vector for interacting atom positions
  map<tuple<string,string,int,string>, vector<tuple<string,string,int,string> > > res_nuc; //residue id and vector of interacting nucleotides
  //map<pair<int,string>, vector<pair<string,float> > > res_sasa;                         //residue id and vector interacting nucleotides and dsasa caused in the protein

  for (auto& it : pdb){
    if(!it.ACTIVE ||
      it.EXT1 == 0 ||
      it.STRUCT_TYPE != "PROTEIN") continue;

    auto id = make_tuple(it.RESN,it.CHAIN,it.RESI,it.iCODE);
    auto inter = res_vec[id];                          //atom positions in the pdb that cause dSASA

    sort(inter.begin(),inter.end());
    auto atom_interac = it.INTERACTION_SASA_P;
    sort(atom_interac.begin(),atom_interac.end());

    //vector of all unique atoms that interact with this residue
    vector<uint32> union_result(inter.size() + it.INTERACTION_SASA_P.size());
    auto union_it = set_union(inter.begin(),inter.end(),
                        atom_interac.begin(),atom_interac.end(),
                  union_result.begin());
    union_result.resize(union_it - union_result.begin());
    res_vec[id] = union_result;
  }

    for (auto it : res_vec){
    auto id = it.first;            //identifier tuple
    auto& interac_p = it.second;    //vector of interacting atoms
    vector<tuple<string,string,int,string> > nucs; //vector of nucleotides
    for (auto pos : interac_p){
        auto& atom = pdb[pos]; //atom in position in pdb vector
      auto atom_nuc = make_tuple(atom.RESN, atom.CHAIN, atom.RESI, atom.iCODE);
      auto p = find(nucs.begin(),nucs.end(),atom_nuc); //search for unique nucleotides
      if (p == nucs.end()) nucs.push_back(atom_nuc); //adds all unique nucleotides
    }
    res_nuc[id] = nucs;
   }

  for (auto it : res_nuc){
    auto resid = it.first;
    auto nuc_vec = it.second;
    auto res_atoms = GetAtoms(get<0>(resid), //pointers to pdb
                                      get<1>(resid),
                                      get<2>(resid),
                                      get<3>(resid),
                                      pdb);
    auto res_type = RES_TYPE[get<0>(resid)];
    map<tuple<string,string,int,string>,float> dsasa;
    map<tuple<string,string,int,string>, int> nuc_types;
    interac.push_back(nuc_vec.size());
    for (auto nuc : nuc_vec){

      dsasa[nuc] = 0;
      nuc_types[nuc] = RES_TYPE[get<0>(nuc)];
    }
    for(auto res_p : res_atoms){
      for (auto i : res_p->AREA_BURIED_BY_ATOM_vector_valid){
        auto nucs = GetNucFromPosVec(res_p->AREA_BURIED_BY_ATOM_vector[i]);
        for (auto nuc : nucs){
          dsasa[nuc] += res_p->AREA_BURIED_BY_ATOM_area[i];
        }
      }
    }
    for(auto nuc : nuc_vec){
      auto data_id = make_pair(res_type,nuc_types[nuc]);
      auto A = dsasa[nuc];

      if (A > 0) data[data_id].push_back(dsasa[nuc]);

    }
  }
}

void
AddRawAAxNUCDataWithExcluding(map<pair<int,int>,vector<float> >& data,
                              map<pair<int,int>,vector<float> >& data_w,
                vector<atom_struct>&                pdb,
                              vector<int>&                        interac){
    vector<string> aa_bb = {"N", "CA", "C", "O", "OXT"};
    vector<string> dna_bb = {"C1\"", "C2\"", "C3\"", "C4\"", "C5\"", "O3\"",
               "O4\"", "O5\"", "OP1", "OP2", "OP3", "P", "O2\""};

    map<string, int > RES_TYPE = {
        {"ALA",  1 },
        {"ARG",  2 },
        {"ASN",  3 },
        {"ASP",  4 },
        {"CYS",  5 },
        {"GLN",  6 },
        {"GLU",  7 },
        {"GLY",  8 },
        {"HIS",  9 },
        {"ILE", 10 },
        {"LEU", 11 },
        {"LYS", 12 },
        {"MET", 13 },
        {"PHE", 14 },
        {"PRO", 15 },
        {"SER", 16 },
        {"THR", 17 },
        {"TRP", 18 },
        {"TYR", 19 },
        {"VAL", 20 },
        {"DA", 1 },
        {"DT", 2 },
        {"DC", 3 },
        {"DG", 4 },
        {"A" , 5},
        {"U" , 6},
        {"C" , 7},
        {"G" , 8}};

        auto GetNucFromPosVec = [&](vector<uint32> atoms){
            set<tuple<string,string,int,string> > nucs;
            for(auto pos : atoms){
                auto& atom = pdb[pos];
                auto nuc_tuple = make_tuple(atom.RESN,atom.CHAIN,atom.RESI,atom.iCODE);
                nucs.insert(nuc_tuple);
            }
            return nucs;
        };
        //residues and nucleotides are identified by a (resn,chain,resi,icode) tuple

        map<tuple<string,string,int,string>, vector<uint32> > res_vec;                           //residue id and vector for interacting atom positions
        map<tuple<string,string,int,string>, vector<tuple<string,string,int,string> > > res_nuc; //residue id and vector of interacting nucleotides
        //map<pair<int,string>, vector<pair<string,float> > > res_sasa;                         //residue id and vector interacting nucleotides and dsasa caused in the protein

        for (auto& it : pdb){
            if(!it.ACTIVE ||
                it.EXT1 == 0 ||
                it.STRUCT_TYPE != "PROTEIN") continue;

            auto id = make_tuple(it.RESN,it.CHAIN,it.RESI,it.iCODE);
            auto inter = res_vec[id];                          //atom positions in the pdb that cause dSASA

            sort(inter.begin(),inter.end());
            auto atom_interac = it.INTERACTION_SASA_P;
            sort(atom_interac.begin(),atom_interac.end());

            //vector of all unique atoms that interact with this residue
            vector<uint32> union_result(inter.size() + it.INTERACTION_SASA_P.size());
            auto union_it = set_union(inter.begin(),inter.end(),
                                      atom_interac.begin(),atom_interac.end(),
                                      union_result.begin());
            union_result.resize(union_it - union_result.begin());
            res_vec[id] = union_result;
        }

        for (auto it : res_vec){
            auto id = it.first;            //identifier tuple
            auto& interac_p = it.second;    //vector of interacting atoms
            vector<tuple<string,string,int,string> > nucs; //vector of nucleotides
            for (auto pos : interac_p){
                auto& atom = pdb[pos]; //atom in position in pdb vector
                auto atom_nuc = make_tuple(atom.RESN, atom.CHAIN, atom.RESI, atom.iCODE);
                auto p = find(nucs.begin(),nucs.end(),atom_nuc); //search for unique nucleotides
                if (p == nucs.end()) nucs.push_back(atom_nuc); //adds all unique nucleotides
            }
            res_nuc[id] = nucs;
        }

        for (auto it : res_nuc){
            auto resid = it.first;
            auto nuc_vec = it.second;
            auto res_atoms = GetAtoms(get<0>(resid), //pointers to pdb
                                              get<1>(resid),
                                              get<2>(resid),
                                              get<3>(resid),
                                              pdb);
            auto res_type = RES_TYPE[get<0>(resid)];
            map<tuple<string,string,int,string>,float> dsasa;

            map<tuple<string,string,int,string>, int> nuc_types;
            interac.push_back(nuc_vec.size());
            for (auto nuc : nuc_vec){

                dsasa[nuc] = 0;

                nuc_types[nuc] = RES_TYPE[get<0>(nuc)];
            }
            for(auto res_p : res_atoms){
                if (find(aa_bb.begin(),aa_bb.end(),res_p->NAME) != aa_bb.end()) continue;

                for (auto i : res_p->AREA_BURIED_BY_ATOM_vector_valid){
                    auto nucs = GetNucFromPosVec(res_p->AREA_BURIED_BY_ATOM_vector[i]);
                    for (auto nuc : nucs){
                        if (find(dna_bb.begin(),dna_bb.end(),get<0>(nuc)) != dna_bb.end()) continue;
                        dsasa[nuc] += res_p->AREA_BURIED_BY_ATOM_area[i];

                    }
                }
            }
            for(auto nuc : nuc_vec){
                auto data_id = make_pair(res_type,nuc_types[nuc]);
                auto A = dsasa[nuc];

                if (A > 0) data[data_id].push_back(dsasa[nuc]);

            }
        }
}


void
AddRawAAxNUCOverlapData(map<pair<int,int>,vector<float> >& data,
            map<pair<int,int>,vector<float> >& data_w,
            vector<atom_struct>&                pdb ){

  map<string, int > RES_TYPE = {
    {"ALA",  1 },
    {"ARG",  2 },
    {"ASN",  3 },
    {"ASP",  4 },
    {"CYS",  5 },
    {"GLN",  6 },
    {"GLU",  7 },
    {"GLY",  8 },
    {"HIS",  9 },
    {"ILE", 10 },
    {"LEU", 11 },
    {"LYS", 12 },
    {"MET", 13 },
    {"PHE", 14 },
    {"PRO", 15 },
    {"SER", 16 },
    {"THR", 17 },
    {"TRP", 18 },
    {"TYR", 19 },
    {"VAL", 20 },
    {"DA", 1 },
    {"DT", 2 },
    {"DC", 3 },
    {"DG", 4 },
    {"DA DA", 1},
    {"DA DT", 2},
    {"DA DC", 3},
    {"DA DG", 4},
    {"DT DG", 5},
    {"DT DC", 6},
    {"DT DT", 7},
    {"DG DG", 8},
    {"DG DC", 9},
    {"DC DC", 10},
    {"DA DA DA", 11},
    {"DA DA DT", 12},
    {"DA DA DG", 13},
    {"DA DA DC", 14},
    {"DA DT DA", 13},
    {"DA DT DG", 13},
    {"DA DT DC", 13},
    {"DA DT DT", 13},
    /*{"A" , 5},
    {"U" , 6},
    {"C" , 7},
    {"G" , 8}*/};
  map<tuple<int,int>, int> OVERLAP_TYPE;

  auto GetNucFromPosVec = [&](vector<uint32> atoms){
                            set<tuple<string,string,int,string> > nucs;
                for(auto pos : atoms){
                  auto& atom = pdb[pos];
                  auto nuc_tuple = make_tuple(atom.RESN,atom.CHAIN,atom.RESI,atom.iCODE);
                  nucs.insert(nuc_tuple);
                    }
                return nucs;
               };



  //residues and nucleotides are identified by a (resn,chain,resi,icode) tuple

    map<tuple<string,string,int,string>, vector<uint32> > res_vec;                           //residue id and vector for interacting atom positions
  map<tuple<string,string,int,string>, vector<tuple<string,string,int,string> > > res_nuc; //residue id and vector of interacting nucleotides
  //map<pair<int,string>, vector<pair<string,float> > > res_sasa;                         //residue id and vector interacting nucleotides and dsasa caused in the protein

  for(auto& it : pdb){
    if(!it.ACTIVE ||
      it.EXT1 == 0 ||
      it.STRUCT_TYPE != "PROTEIN") continue;

    auto id = make_tuple(it.RESN,it.CHAIN,it.RESI,it.iCODE);
    auto inter = res_vec[id];                          //atom positions in the pdb that cause dSASA

    sort(inter.begin(),inter.end());
    auto atom_interac = it.INTERACTION_SASA_P;
    sort(atom_interac.begin(),atom_interac.end());

    //vector of all unique atoms that interact with this residue
    vector<uint32> union_result(inter.size() + it.INTERACTION_SASA_P.size());
    auto union_it = set_union(inter.begin(),inter.end(),
                        atom_interac.begin(),atom_interac.end(),
                  union_result.begin());
    union_result.resize(union_it - union_result.begin());
    res_vec[id] = union_result;
  }

    for (auto it : res_vec){
    auto id = it.first;            //identifier tuple
    auto& interac_p = it.second;    //vector of interacting atoms
    vector<tuple<string,string,int,string> > nucs; //vector of nucleotides
    for (auto pos : interac_p){
        auto& atom = pdb[pos]; //atom in position in pdb vector
      auto atom_nuc = make_tuple(atom.RESN, atom.CHAIN, atom.RESI, atom.iCODE);
      auto p = find(nucs.begin(),nucs.end(),atom_nuc); //search for unique nucleotides
      if (p == nucs.end()) nucs.push_back(atom_nuc); //adds all unique nucleotides
    }
    res_nuc[id] = nucs;
   }

  for (auto it : res_nuc){
      auto resid = it.first;
    auto nuc_vec = it.second;
    auto res_atoms = GetAtoms(get<0>(resid), //pointers to pdb
                                      get<1>(resid),
                                      get<2>(resid),
                                      get<3>(resid),
                                      pdb);
    auto res_type = RES_TYPE[get<0>(resid)];
    map<tuple<string,string,int,string>,float> dsasa;

    map<tuple<string,string,int,string>, int> nuc_types;
    for (auto nuc : nuc_vec){
      dsasa[nuc] = 0;

      nuc_types[nuc] = RES_TYPE[get<0>(nuc)];
    }
    for(auto res_p : res_atoms){
      for (auto i : res_p->AREA_BURIED_BY_ATOM_vector_valid){
        auto nucs = GetNucFromPosVec(res_p->AREA_BURIED_BY_ATOM_vector[i]);
        for (auto nuc : nucs){
          dsasa[nuc] += res_p->AREA_BURIED_BY_ATOM_area[i];

        }
      }
    }
    for(auto nuc : nuc_vec){
      auto data_id = make_pair(res_type,nuc_types[nuc]);
      auto A = dsasa[nuc];

      if (A > 0) data[data_id].push_back(dsasa[nuc]);

    }
  }
}


void
AddRawAtomxNucData(map<pair<int,int>,vector<float> >& data_p,
                   map<pair<int,int>,vector<float> >& data_p_w,
                   vector<atom_struct>&                 pdb){

  for (auto& atom_i : pdb){
      if(!atom_i.ACTIVE ||
          atom_i.EXT1 == 0 ||
          atom_i.ATOM_TYPE == 0) continue;
      for (auto pos : atom_i.INTERACTION_SASA_P){
          auto& atom_j = pdb[pos];
/*          if(!atom_j.ACTIVE ||
             atom_j.EXT1 == 0 ||
             atom_j.ATOM_TYPE == 0)
             continue;*/
          auto dsasa   = atom_i.CONTACT_AREA.at(atom_j.ID);

          if(dsasa > 0){
        if (atom_i.STRUCT_TYPE == "PROTEIN" || atom_i.STRUCT_TYPE == "A"){
                  data_p[make_pair(atom_i.ATOM_TYPE,atom_j.ATOM_TYPE)].push_back(dsasa);
            }
      }
    }
  }
}

void
GeneratePairInteractionData(vector<atom_struct>& pdb){

//pragma omp parallel for
  for (uint32 i = 0; i<pdb.size(); ++i){
    auto& atom_i = pdb[i];
    if(!atom_i.ACTIVE) continue;

    vector<uint32> o_atoms;
      vector<uint32> interac_pos;
      vector<char> valid_overlaps(atom_i.AREA_BURIED_BY_ATOM_vector.size(), true);

      for (uint32 i = 0; i < atom_i.INTERACTION_P.size(); ++i){
        if (atom_i.STRUCT_TYPE != pdb[atom_i.INTERACTION_P[i]].STRUCT_TYPE){
          interac_pos.push_back(atom_i.INTERACTION_P[i]);
        }
      }


      for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i){

        for (uint32 k = 0; k < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++k){
          uint32 other = atom_i.AREA_BURIED_BY_ATOM_vector[i][k];
          if (interac_pos.end() == find(interac_pos.begin(), interac_pos.end(), other)) valid_overlaps[i] = false;
        }
        if (!valid_overlaps[i]) continue;
      atom_i.AREA_BURIED_BY_ATOM_vector_valid.push_back(i);
        for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++j){
          if (o_atoms.end() == find(o_atoms.begin(), o_atoms.end(), atom_i.AREA_BURIED_BY_ATOM_vector[i][j])){
            o_atoms.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i][j]);
          }
        }
      }
      atom_i.INTERACTION_SASA_P = o_atoms;

      for (uint32 i = 0; i < o_atoms.size(); ++i){
        float T_area = 0;

        auto ATOM_POS = o_atoms[i];
        auto& atom_j = pdb[ATOM_POS];

        for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++j){

          if (valid_overlaps[j] &&
              atom_i.AREA_BURIED_BY_ATOM_vector[j].end() != find(atom_i.AREA_BURIED_BY_ATOM_vector[j].begin(), atom_i.AREA_BURIED_BY_ATOM_vector[j].end(), ATOM_POS)){
            T_area += atom_i.AREA_BURIED_BY_ATOM_area[j];

          }
        }
        atom_i.CONTACT_AREA[atom_j.ID] = T_area;

      }
  }
}

void
PrintAARawData(map<string,vector<float> >& data,
               string                       output){
    ofstream outf(output, ofstream::out);
    map<string,float> max_dsasa;
    for (auto it = data.begin(); it != data.end(); ++it){
        auto resn = it->first;
        auto vector = it->second;
        auto size = vector.size();
        max_dsasa[resn] = *max_element(vector.begin(),vector.end());
        outf << resn << "\t" << size <<"\t";
        for (uint32 i = 0; i < size; ++i){
            outf << vector[i] << "\t";
        }
        outf << max_dsasa[resn] << std::endl;
    }
    outf.close();
}

void
PrintAAxNUCRawData(map<pair<int,int>,vector<float> >& data,
                   string                              output){

  ofstream outf(output, ofstream::out);
    map<string,float> max_dsasa;
    for (auto it = data.begin(); it != data.end(); ++it){
        auto resn = it->first.first;
        auto nuc = it->first.second;
        auto vector = it->second;
        outf << resn<<"_"<< nuc << "\t" ;
        for (uint32 i = 0; i < vector.size(); ++i){
            outf << vector[i] << "\t";
        }
        outf << std::endl;
    }
    outf.close();
}

// void
// PrintAAdSASASequence(vector<atom_struct>& pdb,
//                      string               output){
// 
//     ofstream outf(output,ofstream::out);
//     map<int,residue_struct> res;
//     vector<residue_struct> v;
//     for (auto it = pdb.begin(); it != pdb.end(); ++it){
//         if ( it->ACTIVE) {
//             if (res.count(it->RESI) == 0) {
//                 res[it->RESI] = residue_struct(it->RESN,it->CHAIN,it->STRUCTURE,it->STRUCT_TYPE,it->RESI);
//                 res[it->RESI].SASA += it->SASA;
//                 res[it->RESI].dSASA += it->EXT1;
//                 }
//             else {
//                 res[it->RESI].SASA += it->SASA;
//                 res[it->RESI].dSASA += it->EXT1;
//             }
//         }
//     }
//     for (auto it = res.begin(); it != res.end(); ++it){
//         v.push_back(it->second);
//     }
//     auto sorter = [](const residue_struct& i, const residue_struct& j){return i.RESI < j.RESI;};
//     sort(v.begin(),v.end(),sorter);
//     int i = 0;
//     for (auto it = v.begin(); it != v.end(); ++it){
//         outf << it->RESN << "\t" << i <<"\t" << it->RESI <<"\t" <<  it->dSASA << "\t" << it->SASA << std::endl;
//         i++;
//     }
//     outf.close();
// }

void
PrintAtomxNucData(map<pair<int,int>, vector<float>>& data,
                  string                              output){
  ofstream outf(output, ofstream::out);
  for (auto it = data.begin(); it != data.end(); ++it){
    auto t1 = it->first.first;
    auto t2 = it->first.second;
    auto* vector = it->second.data();
    outf << t1 <<"_" << t2 << "\t";
    for (uint32 i = 0; i < it->second.size(); ++i){
      outf << vector[i] << " ";
    }
    outf << std::endl;
  }
  outf.close();
}

void CountInteractionsP(float&               prot_atoms,
                        float&               nuc_atoms,
                        vector<atom_struct>&  pdb){
    for (auto it = pdb.begin();it != pdb.end(); ++it){
        if (it->STRUCT_TYPE != "PROTEIN") continue;
        if (it->CONTACT_AREA.size() > 0){
            prot_atoms++;
            nuc_atoms += it->CONTACT_AREA.size();
        }
    }
}

void 
PrintSplitAsaAtom(vector<atom_struct>& pdb,
                  string outfolder){
  vector<string> amino_acids = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", 
                "PRO", "SER", "THR", "TRP", "TYR", "VAL"
             };
             
  vector<string> bases = { "DA", "DC", "DG", "DT", "A", "C", "G", "T", "U" }; // DNA, RNA

  map<string,string> deoxybase_lookup = { {"A", "DA"}, {"C", "DC"}, {"G", "DG"}, {"T", "DT"} };

  vector<string> aa_backbone = {"N", "CA", "C", "O", "OXT"};
  vector<string> aa_polar_backbone = {"N", "O", "OXT"};
  vector<string> aa_hyd_backbone = {"CA", "C"};
  vector<string> aa_sidechain = {"CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3",
                "CZ", "CZ2", "CZ3", "CH2", "ND1", "ND2", "NE", "NE1", "NE2", "NZ", "NH1", "NH2",
                "OD1", "OD2", "OG", "OG1", "OG2", "OE1", "OE2", "OH", "SD", "SG"};
  vector<string> aa_polar_sidechain = {"ND1", "ND2", "NE", "NE1", "NE2", "NZ", "NH1", "NH2",
                      "OD1", "OD2", "OG", "OG1", "OG2", "OE1", "OE2", "OH"};
  vector<string> aa_hyd_sidechain = {"CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3","CZ", "CZ2", "CZ3", "CH2"};

  vector<string> dna_backbone = {"C1\'", "C2\'", "C3\'", "C4\'", "C5\'", "O3\'", "O4\'", "O5\'", "OP1", "OP2", "OP3", "P", "O2\'"}; // O2" of RNA
  vector<string> dna_polar_backbone = {"O3\'", "O4\'", "O5\'", "OP1", "OP2", "OP3", "O2\'"}; // O2" of RNA
  vector<string> dna_hyd_backbone = {"C1\'", "C2\'", "C3\'", "C4\'", "C5\'","P"};
  vector<string> dna_base = {"C2", "C4", "C5", "C6", "C7", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
  map<string,vector<string>> dna_majorgroove = {{"DA",{"C5", "C6", "C8", "N6", "N7"}},
                   {"DC",{"C4", "C5", "C6", "N4"}},
                   {"DG",{"C5", "C6", "C8", "N7", "O6"}},
                   {"DT",{"C4", "C5", "C6", "C7", "O4"}},
                   {"U",{"C4", "C5", "C6", "O4"}}}; // RNA
  map<string,vector<string>> dna_minorgroove = {{"DA",{"C2", "C4", "N1", "N3"}},
                   {"DC",{"C2", "O2"}},
                   {"DG",{"C2", "C4", "N2", "N3"}},
                   {"DT",{"C2", "N3", "O2"}},
                   {"U",{"C2", "N3", "O2"}}}; // RNA
  map<string,vector<string>> dna_ambiguous = {{"DA",{"N9"}},
                {"DC",{"N1", "N3"}},
                {"DG",{"N1", "N9"}},
                {"DT",{"N1"}},
                {"U",{"N1"}}}; // RNA

  map<string,vector<string>> dna_polar_majorgroove = {{"DA",{"N6", "N7"}},
                         {"DC",{"N4"}},
                         {"DG",{"N7", "O6"}},
                         {"DT",{"O4"}},
                         {"U",{"O4"}}}; // RNA
  map<string,vector<string>> dna_polar_minorgroove = {{"DA",{"N1", "N3"}},
                         {"DC",{"O2"}},
                         {"DG",{"N2", "N3"}},
                         {"DT",{"N3", "O2"}},
                         {"U",{"N3", "O2"}}};
  map<string,vector<string>> dna_polar_ambiguous = {{"DA",{"N9"}},
                      {"DC",{"N1", "N3"}},
                      {"DG",{"N1", "N9"}},
                      {"DT",{"N1"}},
                      {"U",{"N1"}}}; // RNA

  map<string,vector<string>>  dna_hyd_majorgroove = {{"DA",{"C5", "C6", "C8"}},
                       {"DC",{"C4", "C5", "C6"}},
                       {"DG",{"C5", "C6", "C8"}},
                       {"DT",{"C4", "C5", "C6", "C7"}},
                       {"U",{"C4", "C5", "C6"}}}; // RNA
  map<string,vector<string>> dna_hyd_minorgroove = {{"DA",{"C2", "C4"}},
                       {"DC",{"C2"}},
                       {"DG",{"C2", "C4"}},
                       {"DT",{"C2"}},
                       {"U",{"C2"}}}; // RNA
  map<string,vector<string>>dna_hyd_ambiguous = {{"DA",{}},
                    {"DC",{}},
                    {"DG",{}},
                    {"DT",{}},
                    {"U",{}}}; // RNA

  const char* formatS = "%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f";
  string header = "atmname\tresname\tchain\tresnum\ttotal_ASA\tbb_ASA\tsc_ASA\tmajorgroove_ASA\tminorgroove_ASA\tnogroove_ASA\tpolar_ASA\tpolar_bb_ASA\tpolar_sc_ASA\tpolar_majorgroove_ASA\tpolar_minorgroove_ASA\tpolar_nogroove_ASA\thyd_ASA\thyd_bb_ASA\thyd_sc_ASA\thyd_majorgroove_ASA\thyd_minorgroove_ASA\thyd_nogroove_ASA\tlig_ASA\tlig_polar_ASA\tlig_hyd_ASA\tX\tY\tZ";
  //ligands total, polar(O,N) and hydrophobic(everything else).
  string obj_name = pdb[0].STRUCTURE;
#ifdef WIN32__
  string fS = "\\";
#else
  string fS = "/";
#endif
  vector<string> fnames = {outfolder,
                           /*outfolder+fS+obj_name+"_protein.atmasa",
                           outfolder+fS+obj_name+"_dna.atmasa",
                           outfolder+fS+obj_name+"_rna.atmasa",
                           outfolder+fS+obj_name+"_bsa.atmasa",
                           outfolder+fS+obj_name+"_ligand.atmasa"*/};
  ofstream outf(fnames[0]);
  outf << header << std::endl;
/*  vector<atom_struct*> prot;
  vector<atom_struct*> dna;
  vector<atom_struct*> rna;
  vector<atom_struct*> ligd;
  vector<atom_struct*> bsa;
  
  for (uint32 i; i < pdb.size(); ++i){
    atom_struct* cur_atom;
    atom_struct* data = pdb.data();
    if(data[i].MOL_TYPE == "PROTEIN") prot.push_back(data+i);
    else if (data[i].MOL_TYPE == "DNA") dna.push_back(data+i);
    else if (data[i].MOL_TYPE == "RNA") rna.push_back(data+i);
    else if (data[i].MOL_TYPE == "LIGAND") ligd.push_back(data+i);
    else if (data[i].EXT1 > 0.0) bsa.push_back(data+i);
  }*/
    auto opdb = OrganizePDB(pdb);
  for(auto& chain : opdb){
    for(auto& res : chain){
      for(auto& atom : res){
        float total_ASA,bb_ASA,sc_ASA,majorgroove_ASA,minorgroove_ASA,nogroove_ASA,polar_ASA,polar_bb_ASA,polar_sc_ASA,polar_majorgroove_ASA,polar_minorgroove_ASA,polar_nogroove_ASA,hyd_ASA,hyd_bb_ASA,hyd_sc_ASA,hyd_majorgroove_ASA,hyd_minorgroove_ASA,hyd_nogroove_ASA,lig_ASA,lig_polar_ASA,lig_hyd_ASA;

        total_ASA = NAN;
        bb_ASA = NAN;
        sc_ASA = NAN;
        majorgroove_ASA = NAN;
        minorgroove_ASA = NAN;
        nogroove_ASA = NAN;
        polar_ASA = NAN;
        polar_bb_ASA = NAN;
        polar_sc_ASA = NAN;
        polar_majorgroove_ASA = NAN;
        polar_minorgroove_ASA = NAN;
        polar_nogroove_ASA = NAN;
        hyd_ASA = NAN;
        hyd_bb_ASA = NAN;
        hyd_sc_ASA = NAN;
        hyd_majorgroove_ASA = NAN;
        hyd_minorgroove_ASA = NAN;
        hyd_nogroove_ASA = NAN;
        lig_ASA = NAN;
        lig_polar_ASA = NAN;
        lig_hyd_ASA = NAN;
        if(find(amino_acids.begin(),amino_acids.end(),atom->RESN) != amino_acids.end()){
          if(find(aa_backbone.begin(),aa_backbone.end(),atom->NAME) != aa_backbone.end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(bb_ASA)) bb_ASA = 0.0;
            
            total_ASA += atom->SASA;
            bb_ASA += atom->SASA;
          }
          else if(find(aa_sidechain.begin(),aa_sidechain.end(),atom->NAME) != aa_sidechain.end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(sc_ASA)) sc_ASA = 0.0;
            
            total_ASA += atom->SASA;
            sc_ASA += atom->SASA;
          }
          if(find(aa_polar_backbone.begin(),aa_polar_backbone.end(),atom->NAME) != aa_polar_backbone.end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_bb_ASA)) polar_bb_ASA = 0.0;
 
            polar_ASA += atom->SASA;
            polar_bb_ASA += atom->SASA;
          }
          else if(find(aa_polar_sidechain.begin(),aa_polar_sidechain.end(),atom->NAME) != aa_polar_sidechain.end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_sc_ASA)) polar_sc_ASA = 0.0;

            polar_ASA += atom->SASA;
            polar_sc_ASA += atom->SASA;
          }
          if(find(aa_hyd_backbone.begin(),aa_hyd_backbone.end(),atom->NAME) != aa_hyd_backbone.end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_bb_ASA)) hyd_bb_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_bb_ASA += atom->SASA;
          }
          else if(find(aa_hyd_sidechain.begin(),aa_hyd_sidechain.end(),atom->NAME) != aa_hyd_sidechain.end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_sc_ASA)) hyd_sc_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_sc_ASA += atom->SASA;
          }
        }
        else if(find(bases.begin(),bases.end(),atom->RESN) != bases.end()){
          string resname2;
          if (deoxybase_lookup.count(atom->RESN) != 0) resname2 = deoxybase_lookup.at(atom->RESN);
          else resname2 = atom->RESN;

          if(find(dna_backbone.begin(),dna_backbone.end(),atom->NAME) != dna_backbone.end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(bb_ASA)) bb_ASA = 0.0;
            total_ASA += atom->SASA;
            bb_ASA+= atom->SASA;
          }
          else if(find(dna_majorgroove[resname2].begin(),dna_majorgroove[resname2].end(),atom->NAME) != dna_majorgroove[resname2].end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(sc_ASA)) sc_ASA = 0.0;
            if(std::isnan(majorgroove_ASA)) majorgroove_ASA = 0.0;

            total_ASA += atom->SASA;
            sc_ASA += atom->SASA;
            majorgroove_ASA += atom->SASA;
          }
          else if(find(dna_minorgroove[resname2].begin(),dna_minorgroove[resname2].end(),atom->NAME) != dna_minorgroove[resname2].end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(sc_ASA)) sc_ASA = 0.0;
            if(std::isnan(minorgroove_ASA)) minorgroove_ASA = 0.0;
            total_ASA += atom->SASA;
            sc_ASA += atom->SASA;
            minorgroove_ASA += atom->SASA;
          }
          else if(find(dna_ambiguous[resname2].begin(),dna_ambiguous[resname2].end(),atom->NAME) != dna_ambiguous[resname2].end()){
            if(std::isnan(total_ASA)) total_ASA = 0.0;
            if(std::isnan(sc_ASA)) sc_ASA = 0.0;
            if(std::isnan(nogroove_ASA)) nogroove_ASA = 0.0;
            total_ASA += atom->SASA;
            sc_ASA += atom->SASA;
            nogroove_ASA += atom->SASA;
          }
          if(find(dna_polar_backbone.begin(),dna_polar_backbone.end(),atom->NAME) != dna_polar_backbone.end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_bb_ASA)) polar_bb_ASA = 0.0;
            polar_ASA += atom->SASA;
            polar_bb_ASA += atom->SASA;
          }
          else if(find(dna_polar_majorgroove[resname2].begin(),dna_polar_majorgroove[resname2].end(),atom->NAME) != dna_polar_majorgroove[resname2].end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_sc_ASA)) polar_sc_ASA = 0.0;
            if(std::isnan(polar_majorgroove_ASA)) polar_majorgroove_ASA = 0.0;
            polar_ASA += atom->SASA;
            polar_sc_ASA += atom->SASA;
            polar_majorgroove_ASA += atom->SASA;
          }
          else if(find(dna_polar_minorgroove[resname2].begin(),dna_polar_minorgroove[resname2].end(),atom->NAME) != dna_polar_minorgroove[resname2].end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_sc_ASA)) polar_sc_ASA = 0.0;
            if(std::isnan(polar_minorgroove_ASA)) polar_minorgroove_ASA = 0.0;
            polar_ASA += atom->SASA;
            polar_sc_ASA += atom->SASA;
            polar_minorgroove_ASA += atom->SASA;
          }
          else if(find(dna_polar_ambiguous[resname2].begin(),dna_polar_ambiguous[resname2].end(),atom->NAME) != dna_polar_ambiguous[resname2].end()){
            if(std::isnan(polar_ASA)) polar_ASA = 0.0;
            if(std::isnan(polar_sc_ASA)) polar_sc_ASA = 0.0;
            if(std::isnan(polar_nogroove_ASA)) polar_nogroove_ASA = 0.0;
            polar_ASA += atom->SASA;
            polar_sc_ASA += atom->SASA;
            polar_nogroove_ASA += atom->SASA;
          }
          
          if(find(dna_hyd_backbone.begin(),dna_hyd_backbone.end(),atom->NAME) != dna_hyd_backbone.end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_bb_ASA)) hyd_bb_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_bb_ASA += atom->SASA;
          }
          else if(find(dna_hyd_majorgroove[resname2].begin(),dna_hyd_majorgroove[resname2].end(),atom->NAME) != dna_hyd_majorgroove[resname2].end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_sc_ASA)) hyd_sc_ASA = 0.0;
            if(std::isnan(hyd_majorgroove_ASA)) hyd_majorgroove_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_sc_ASA += atom->SASA;
            hyd_majorgroove_ASA += atom->SASA;
          }
          else if(find(dna_hyd_minorgroove[resname2].begin(),dna_hyd_minorgroove[resname2].end(),atom->NAME) != dna_hyd_minorgroove[resname2].end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_sc_ASA)) hyd_sc_ASA = 0.0;
            if(std::isnan(hyd_minorgroove_ASA)) hyd_minorgroove_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_sc_ASA += atom->SASA;
            hyd_minorgroove_ASA += atom->SASA;
          }
          else if(find(dna_hyd_ambiguous[resname2].begin(),dna_hyd_ambiguous[resname2].end(),atom->NAME) != dna_hyd_ambiguous[resname2].end()){
            if(std::isnan(hyd_ASA)) hyd_ASA = 0.0;
            if(std::isnan(hyd_sc_ASA)) hyd_sc_ASA = 0.0;
            if(std::isnan(hyd_nogroove_ASA)) hyd_nogroove_ASA = 0.0;
            hyd_ASA += atom->SASA;
            hyd_sc_ASA += atom->SASA;
            hyd_nogroove_ASA += atom->SASA;
          }
        }
        else if (atom->MOL_TYPE == "LIGAND"){
          if(std::isnan(total_ASA)) total_ASA = 0.0;
          total_ASA += atom->SASA;
          if(atom->ELEMENT == "C")
          {
            if(std::isnan(lig_ASA)) lig_ASA = 0.0;
            if(std::isnan(lig_hyd_ASA)) lig_hyd_ASA = 0.0;
            lig_ASA += atom->SASA;
            lig_hyd_ASA += atom->SASA;
          }
          else{
            if(std::isnan(lig_ASA)) lig_ASA = 0.0;
            if(std::isnan(lig_polar_ASA)) lig_polar_ASA = 0.0;
            lig_ASA += atom->SASA;
            lig_polar_ASA += atom->SASA;
          }
        }
        outf << atom->NAME << "\t";
        char buffer[8192];
        snprintf(buffer,8192,formatS,atom->RESN.c_str(),atom->CHAIN.c_str(),atom->RESI,total_ASA,bb_ASA,sc_ASA,majorgroove_ASA,minorgroove_ASA,nogroove_ASA,polar_ASA,polar_bb_ASA,polar_sc_ASA,polar_majorgroove_ASA,polar_minorgroove_ASA,polar_nogroove_ASA,hyd_ASA,hyd_bb_ASA,hyd_sc_ASA,hyd_majorgroove_ASA,hyd_minorgroove_ASA,hyd_nogroove_ASA,lig_ASA,lig_polar_ASA,lig_hyd_ASA,atom->COORDS[0],atom->COORDS[1],atom->COORDS[2]);
        string data = string(buffer);
        outf << data << std::endl;
      }
    }
  }
  outf.close();
}
