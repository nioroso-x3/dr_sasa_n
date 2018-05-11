#include "stdafx.h"
#include "atom_struct.h"
#include "SurfaceSolverCL.h"
#include "SurfaceSolverOnTheFly.h"

void
SolveInteractions(vector<atom_struct>& pdb,
                  uint32               mode){


  auto pdb_size = pdb.size();

  function<bool(const atom_struct&,const atom_struct&)> selector[] = {
      [](const atom_struct& i,const atom_struct& j){return true;},
      [](const atom_struct& i,const atom_struct& j){return i.MOL_TYPE == j.MOL_TYPE;},
      [](const atom_struct& i,const atom_struct& j){return i.MOL_TYPE != j.MOL_TYPE;},
      [](const atom_struct& i,const atom_struct& j){return i.CHAIN    != j.CHAIN;},
  };
#pragma omp parallel for schedule(dynamic)
  for (uint32 i = 0; i < pdb_size; ++i){
    for (uint32 j = i + 1; j < pdb_size; ++j){
      auto& atom_i = pdb[i];
      auto* C_I = atom_i.COORDS.data();
      auto& atom_j = pdb[j];
      auto* C_J = atom_j.COORDS.data();
      float dist = ((C_I[0] - C_J[0])*(C_I[0] - C_J[0]) +
                     (C_I[1] - C_J[1])*(C_I[1] - C_J[1]) +
                     (C_I[2] - C_J[2])*(C_I[2] - C_J[2]));
      float radsum = atom_i.RADIUS + atom_j.RADIUS;
      radsum*=radsum;
      if (dist <= radsum  && selector[mode](atom_i,atom_j) ){
        float rdist = sqrt(dist);
        #pragma omp critical(sasa_simple)
        {
          atom_i.INTERACTION_P.push_back(j);
          atom_j.INTERACTION_P.push_back(i);
          atom_i.DISTANCES[j] = rdist;
          atom_j.DISTANCES[i] = rdist;
        }
      }
    }
  }
}

void
SolveInteractionsALL(vector<atom_struct>&            pdb,
                     int                             mode,
                     vector<vector<string>>          objects){


  auto pdb_size = pdb.size();

  auto selector = [&](atom_struct& atom_i,atom_struct& atom_j){

    vector<char> isIobj(objects.size(),false);
    vector<char> isJobj(objects.size(),false);
    if(objects.size() > 1 && (mode == 1 || mode == 0)){
    for (uint32 i = 0; i < objects.size(); i++){
        stringstream ss;
        auto chainset = objects[i];
        for (auto c : chainset){
          ss << c;
        }
        string obj = ss.str();
        for (auto C1 : chainset){
          if (atom_i.CHAIN == C1){
            isIobj[i] = true;
            atom_i.STRUCT_TYPE = obj;
          }
          if (atom_j.CHAIN == C1){
            isJobj[i] = true;
            atom_j.STRUCT_TYPE = obj;
          }
        }
      }
    } else if (mode == 2){
      stringstream objI;
      stringstream objJ;
      objI << atom_i.RESN << "/" << atom_i.RESI << "/" << atom_i.CHAIN;
      objJ << atom_j.RESN << "/" << atom_j.RESI << "/" << atom_j.CHAIN;
      atom_i.STRUCT_TYPE = objI.str();
      atom_j.STRUCT_TYPE = objJ.str();
    } else if (mode == 3){
      stringstream objI;
      stringstream objJ;
      objI << atom_i.NAME << "/" << atom_i.RESN << "/" << atom_i.RESI << "/" << atom_i.CHAIN;
      objJ << atom_j.NAME << "/" << atom_j.RESN << "/" << atom_j.RESI << "/" << atom_j.CHAIN;
      atom_i.STRUCT_TYPE = objI.str();
      atom_j.STRUCT_TYPE = objJ.str();
      //cout << objI.str() << " " << objJ.str() << "\n";
    }else if (mode == 4){
      atom_i.STRUCT_TYPE = atom_i.MOL_TYPE;
      atom_j.STRUCT_TYPE = atom_j.MOL_TYPE;
    }
    
    if(mode == 2 || mode == 3 || mode == 4){
      bool samestruct = atom_i.STRUCT_TYPE == atom_j.STRUCT_TYPE;
      
      return samestruct;
    }

    //cout << "m: "<<mode << "\t" << atom_i.STRUCT_TYPE << "\t" << atom_j.STRUCT_TYPE << "\n";
    if (mode == 0) return true;
    if (mode == 1){
      for (uint32 i = 0; i < objects.size(); i++){
        if (isIobj[i] && isJobj[i]){
          //cout<< atom_i.ID<< "\t" << atom_i.CHAIN << "\t" << atom_i.STRUCT_TYPE << "\t|\t"<<atom_j.ID<< "\t" <<atom_j.CHAIN<<"\t"<< atom_j.STRUCT_TYPE << "\tT" <<"\n";
          return true;
        }
      }
    }
    return false;
  };

#pragma omp parallel for schedule(dynamic)
  for (uint32 i = 0; i < pdb_size; ++i){
    for (uint32 j = i + 1; j < pdb_size; ++j){
      auto& atom_i = pdb[i];
      auto* C_I = atom_i.COORDS.data();
      auto& atom_j = pdb[j];
      auto* C_J = atom_j.COORDS.data();
      float dist = ((C_I[0] - C_J[0])*(C_I[0] - C_J[0]) +
                     (C_I[1] - C_J[1])*(C_I[1] - C_J[1]) +
                     (C_I[2] - C_J[2])*(C_I[2] - C_J[2]));
      float radsum = atom_i.RADIUS + atom_j.RADIUS;
      radsum*=radsum;

      if (dist <= radsum){
        #pragma omp critical(sasa_complex)
        {
          if(selector(atom_i,atom_j)){
            float rdist = sqrt(dist);
            //if (mode == 0)cout << atom_i.print() << " / " << atom_j.print() << "\n";
            atom_i.INTERACTION_P.push_back(j);
            atom_j.INTERACTION_P.push_back(i);
            atom_i.DISTANCES[j] = rdist;
            atom_j.DISTANCES[i] = rdist;
          }
        }
      }
    }
  }
}

vector<char> 
SURF_GenerateImtrx(vector<atom_struct>& pdb){
  vector<char> result;
  result.resize(pdb.size()*pdb.size(),false);
  for (uint i = 0; i < pdb.size(); i++){
    auto& atom = pdb[i];
    for (uint& j : atom.INTERACTION_P ){
      result[i + j*pdb.size()] = true;
    }
  }
  return result;
}


void
SimpleSolver_2(vector<atom_struct>& pdb,
               vector<float>&      points){
  auto shell_s3 = points.size();
  auto shell_s = shell_s3 / 3;
  auto* p = points.data();
  uint32 pdbs = pdb.size();
#pragma omp parallel for
  for(uint32 i = 0; i < pdbs; ++i){
    auto& atom_i = pdb[i];

    if (!atom_i.ACTIVE) continue;

    atom_i.SHELL_BURIED.resize(shell_s,false);
    atom_i.ACCESSIBLE_P = shell_s;
    auto* b = atom_i.SHELL_BURIED.data();
    auto* C_I = atom_i.COORDS.data();
    auto R_I = atom_i.RADIUS;

    for (uint32 j = 0; j < atom_i.INTERACTION_P.size();++j){
      auto& atom_j = pdb[atom_i.INTERACTION_P[j]];
      if (!atom_j.ACTIVE) continue;
      auto* C_J = atom_j.COORDS.data();
      auto R_J2 = atom_j.RADIUS2;
      uint counter = 0;
      for (uint32 k = 0; k < shell_s; ++k){
        if(!b[k]){
          counter++;
          auto k3 = k * 3;
          auto Xj = p[k3]     * R_I + C_I[0] - C_J[0];
          auto Yj = p[k3 + 1] * R_I + C_I[1] - C_J[1];
          auto Zj = p[k3 + 2] * R_I + C_I[2] - C_J[2];
          float dist_j = Xj*Xj + Yj*Yj + Zj*Zj;
          b[k] = (dist_j <= R_J2);
        }
      }
    }
    atom_i.ACCESSIBLE_P -= count(atom_i.SHELL_BURIED.begin(), atom_i.SHELL_BURIED.end(), true);
    atom_i.SASA = atom_i.AREA*((float)atom_i.ACCESSIBLE_P / (float)shell_s);
  }
}

void
Generic_Solver(vector<atom_struct>&       pdb,
               vector<float>&         points,
               vector<vector<string>>  obj1,
               int                     mode,
               int                     cl_mode){

  vector<atom_struct> pdb_n = pdb;
  uint32 shell_s3 = points.size();
  uint32 shell_s = shell_s3 / 3;
  auto* p = points.data();
  SolveInteractionsALL(pdb,0,obj1);
  SolveInteractionsALL(pdb_n,mode,obj1);
  SimpleSolverCL(pdb,points,cl_mode);
  SimpleSolverCL(pdb_n,points,cl_mode);


  vector<uint32> dsasa_pos;
  for (uint32 i = 0; i < pdb.size();++i){
    pdb[i].EXT1 = pdb_n[i].SASA - pdb[i].SASA;
    pdb[i].STRUCT_TYPE = pdb_n[i].STRUCT_TYPE;
    if(pdb[i].EXT1 > 0){
      dsasa_pos.push_back(i);
    }
  }
  auto* primes = PRIMES.data();
  auto perm_comp= [&](vector<uint32>& a,vector<uint32>& b){
    if(a.size() != b.size())
      return false;
    else{
      long long a1 = 1;
      long long b1 = 1;
      for(uint32 i = 0; i < a.size();++i){
        a1*=primes[a[i]];
        b1*=primes[b[i]];
      }
      return (a1 == b1);
    }
  };

#pragma omp parallel for
  for (uint32 idx = 0; idx < dsasa_pos.size();++idx){
    uint32 i = dsasa_pos[idx];
    map<uint32, vector<uint32> > SHELL_BURIED_MAP;       //stores which atoms bury a certain point
    vector<vector<uint32> > SHELL_BURIED_COUNT;          //stores vectors of atom positions
    vector<vector<uint32> > AREA_BURIED_BY_ATOM_vector;  //stores all unique sets of atoms that bury points
    vector<float> AREA_BURIED_BY_ATOM_area;             //stores the area that each unique set of atoms buries

    auto& atom_i = pdb[i];
    auto* C_I = atom_i.COORDS.data();
    auto R_I = atom_i.RADIUS;
// #pragma omp critical(debug)
//     {
//     cout << atom_i.ID  <<"\t" << atom_i.NAME <<  "\t" << atom_i.EXT1 << "\t" <<
//             pdb_n[i].SASA <<"\t" << atom_i.SASA << "\t" << atom_i.STRUCT_TYPE <<"\n";
//   }
    if(!atom_i.ACTIVE) continue;
    for (uint32 j = 0; j < atom_i.INTERACTION_P.size(); ++j){
      auto& atom_j = pdb[atom_i.INTERACTION_P[j]];

      if(!atom_j.ACTIVE) continue;
      auto* C_J = atom_j.COORDS.data();
      //if((atom_i.EXT1 == 0.0) && (atom_j.EXT1 == 0.0)) continue;


      auto R_J2 = atom_j.RADIUS2;

      for (uint32 k = 0; k < shell_s; ++k){
        auto k3= k*3;
        auto Xj = p[k3]     * R_I + C_I[0] - C_J[0];
        auto Yj = p[k3 + 1] * R_I + C_I[1] - C_J[1];
        auto Zj = p[k3 + 2] * R_I + C_I[2] - C_J[2];

        float dist_j = Xj*Xj + Yj*Yj + Zj*Zj;
        if (dist_j <= R_J2){
          SHELL_BURIED_MAP[k].push_back(atom_i.INTERACTION_P[j]);
        }
      }
    }
    for (auto it = SHELL_BURIED_MAP.begin(); it != SHELL_BURIED_MAP.end(); ++it) {
      //sort(it->second.begin(), it->second.end());
      SHELL_BURIED_COUNT.push_back(it->second);
    }
    if(SHELL_BURIED_COUNT.size() == 0) continue;

    //count which overlaps are the same
    auto s_size = SHELL_BURIED_COUNT.size();
    vector<char> COUNT_CHECK(s_size, false);
    auto* c = COUNT_CHECK.data();
    for (uint32 k = 0; k < s_size - 1; ++k) {
      if (c[k]) continue;
      c[k] = true;
      uint32 count = 1;
      for (uint32 l = k + 1; l < s_size; ++l){
        if (c[l]) continue;
        if (perm_comp(SHELL_BURIED_COUNT[k],SHELL_BURIED_COUNT[l])){
        //if (SHELL_BURIED_COUNT[k] == SHELL_BURIED_COUNT[l]){
          c[l] = true;
          count++;
        }
      }
      sort(SHELL_BURIED_COUNT[k].begin(),SHELL_BURIED_COUNT[k].end());
      AREA_BURIED_BY_ATOM_vector.push_back(SHELL_BURIED_COUNT[k]);
      AREA_BURIED_BY_ATOM_area.push_back(atom_i.AREA*(float)count / (float)shell_s);
    }
    atom_i.AREA_BURIED_BY_ATOM_vector = AREA_BURIED_BY_ATOM_vector;
    atom_i.AREA_BURIED_BY_ATOM_area = AREA_BURIED_BY_ATOM_area;
  }
}

void 
DecoupledSolver(vector<atom_struct>&       pdb,
                vector<float>&             points){
  uint32 shell_s = points.size()/3.0;
  auto* primes = PRIMES.data();
  auto* p = points.data();
  auto perm_comp= [&](vector<uint32>& a,vector<uint32>& b){
    if(a.size() != b.size())
      return false;
    else{
      long long a1 = 1;
      long long b1 = 1;
      for(uint32 i = 0; i < a.size();++i){
        a1*=primes[a[i]];
        b1*=primes[b[i]];
      }
      return (a1 == b1);
    }
  };


  #pragma omp parallel for
  for (uint32 i = 0; i < pdb.size();++i){
    map<uint32, vector<uint32> > SHELL_BURIED_MAP;       //stores which atoms bury a certain point
    vector<vector<uint32> > SHELL_BURIED_COUNT;          //stores vectors of atom positions
    vector<vector<uint32> > AREA_BURIED_BY_ATOM_vector;  //stores all unique sets of atoms that bury points
    vector<float> AREA_BURIED_BY_ATOM_area;             //stores the area that each unique set of atoms buries

    auto& atom_i = pdb[i];
    auto* C_I = atom_i.COORDS.data();
    auto R_I = atom_i.RADIUS;
// #pragma omp critical(debug)
//     {
//     cout << atom_i.ID  <<"\t" << atom_i.NAME <<  "\t" << atom_i.EXT1 << "\t" <<
//             pdb_n[i].SASA <<"\t" << atom_i.SASA << "\t" << atom_i.STRUCT_TYPE <<"\n";
//   }
    if(!atom_i.ACTIVE) continue;
    for (uint32 j = 0; j < atom_i.INTERACTION_P.size(); ++j){
      auto& atom_j = pdb[atom_i.INTERACTION_P[j]];

      if(!atom_j.ACTIVE) continue;
      auto* C_J = atom_j.COORDS.data();
      //if((atom_i.EXT1 == 0.0) && (atom_j.EXT1 == 0.0)) continue;


      auto R_J2 = atom_j.RADIUS2;

      for (uint32 k = 0; k < shell_s; ++k){
        auto k3= k*3;
        auto Xj = p[k3]     * R_I + C_I[0] - C_J[0];
        auto Yj = p[k3 + 1] * R_I + C_I[1] - C_J[1];
        auto Zj = p[k3 + 2] * R_I + C_I[2] - C_J[2];

        float dist_j = Xj*Xj + Yj*Yj + Zj*Zj;
        if (dist_j <= R_J2){
          SHELL_BURIED_MAP[k].push_back(atom_i.INTERACTION_P[j]);
        }
      }
    }
    for (auto it = SHELL_BURIED_MAP.begin(); it != SHELL_BURIED_MAP.end(); ++it) {
      //sort(it->second.begin(), it->second.end());
      SHELL_BURIED_COUNT.push_back(it->second);
    }
    if(SHELL_BURIED_COUNT.size() == 0) continue;

    //count which overlaps are the same
    auto s_size = SHELL_BURIED_COUNT.size();
    vector<char> COUNT_CHECK(s_size, false);
    auto* c = COUNT_CHECK.data();
    for (uint32 k = 0; k < s_size - 1; ++k) {
      if (c[k]) continue;
      c[k] = true;
      uint32 count = 1;
      for (uint32 l = k + 1; l < s_size; ++l){
        if (c[l]) continue;
        if (perm_comp(SHELL_BURIED_COUNT[k],SHELL_BURIED_COUNT[l])){
        //if (SHELL_BURIED_COUNT[k] == SHELL_BURIED_COUNT[l]){
          c[l] = true;
          count++;
        }
      }
      sort(SHELL_BURIED_COUNT[k].begin(),SHELL_BURIED_COUNT[k].end());
      AREA_BURIED_BY_ATOM_vector.push_back(SHELL_BURIED_COUNT[k]);
      AREA_BURIED_BY_ATOM_area.push_back(atom_i.AREA*(float)count / (float)shell_s);
    }
    atom_i.AREA_BURIED_BY_ATOM_vector = AREA_BURIED_BY_ATOM_vector;
    atom_i.AREA_BURIED_BY_ATOM_area = AREA_BURIED_BY_ATOM_area;
    int l = AREA_BURIED_BY_ATOM_vector.size();
    for(uint i = 0; i < l; ++i){
      auto v = AREA_BURIED_BY_ATOM_vector[i];
      float a = AREA_BURIED_BY_ATOM_area[i];
      atom_i.EXT0 += (a / v.size()); 
    }

  }
}

