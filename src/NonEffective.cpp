#include "stdafx.h"
#include "atom_struct.h"
#include "NonEffective.h"

void
NE_RemoveNonEffective(vector<atom_struct>& pdb,   //structure
                      float               o,     //angle
                      vector<char>&        Imtrx){//interaction matrix
  if (o == 180.0) return;
  o = PI * (o / 180.0);
  uint64 l = pdb.size();
  auto* pdbv = pdb.data();
  #pragma omp parallel for schedule(dynamic)
  for (uint64 i = 0; i < l; i++){
    auto& atom_i = pdbv[i];
    auto* C_I = atom_i.COORDS.data();
    for (uint64 j : atom_i.DD_INTERACTING){
      auto* C_J = pdbv[j].COORDS.data();
      for (uint64 k : atom_i.DD_INTERACTING){
        if (k == j)continue;
        auto* C_K = pdbv[k].COORDS.data();
        float Vik[] = {C_I[0] - C_K[0],
                        C_I[1] - C_K[1],
                        C_I[2] - C_K[2]};
        float Vjk[] = {C_J[0] - C_K[0],
                        C_J[1] - C_K[1],
                        C_J[2] - C_K[2]};
        float mVik =   sqrt(Vik[0]*Vik[0] + Vik[1]*Vik[1] + Vik[2]*Vik[2]);
        float mVjk =   sqrt(Vjk[0]*Vjk[0] + Vjk[1]*Vjk[1] + Vjk[2]*Vjk[2]);
        float nVik[] = { Vik[0]/mVik, Vik[1]/mVik, Vik[2]/mVik };
        float nVjk[] = { Vjk[0]/mVjk, Vjk[1]/mVjk, Vjk[2]/mVjk };
        float dot =    nVik[0]*nVjk[0] + nVik[1]*nVjk[1] + nVik[2]*nVjk[2];
        float angle =  acos(dot);
        if (angle > o){
          Imtrx[i + j * l] = false;
          break;
        }
      }
    }
  }
}