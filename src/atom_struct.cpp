#include "stdafx.h"
#include "atom_struct.h"


atom_struct::atom_struct(uint32          id,
                         int             resi,
                         string          icode,
                         string          name,
                         string          resn,
                         string          chain,
                         string          elem,
                         string          obj,
                         string          mol,
                         float           x,
                         float           y,
                         float           z,
                         string          altloc,
                         float           occu,
                         float           tfactor,
                         string          charge){
  ID = id;
  RESI = resi;
  iCODE = icode;
  NAME = name;
  RESN = resn;
  CHAIN = chain;
  ELEMENT = elem;
  STRUCTURE = obj;
  vector<float> coord = {x,y,z};
  COORDS = coord;
  VDW = 0;
  RADIUS = 0;
  RADIUS2 = 0;
  AREA = 0;
  SASA = 0;
  ACCESSIBLE_P = 0;
  EXT0 = 0;
  EXT1 = 0;
  ACTIVE = true;
  STRUCT_TYPE = mol;
  ALTLOC = altloc;
  OCCUPANCY = occu;
  TFACTOR = tfactor;
  CHARGE = charge;
  ENERGY = 0;
  ATOM_TYPE = 0;
  ATOM_TYPE_40 = 0;
  INT_NUM = 0;
  DD_ENERGY = 0;
  AS_ENERGY = 0;
  CS_ENERGY = 0;
  sDD_ENERGY = 0;
  BSA_ENERGY = 0;
  EF_INT_NUM = 0;
  TERMINAL = false;
}

atom_struct&
atom_struct::operator=(const atom_struct& other){
  STRUCTURE = other.STRUCTURE;
  ID = other.ID;
  NAME = other.NAME;
  RESN = other.RESN;
  CHAIN = other.CHAIN;
  ELEMENT = other.ELEMENT;
  STRUCT_TYPE = other.STRUCT_TYPE;
  RESI = other.RESI;
  iCODE = other.iCODE;
  VDW = other.VDW;
  RADIUS = other.RADIUS;
  RADIUS2 = other.RADIUS2;
  AREA = other.AREA;                                        //Full sphere area, vdw+probe radius
  SASA = other.SASA;                                        //exposed area
  EXT0 = other.EXT0;                                        //currently used to store the dSASA calculated from overlap discovering
  EXT1 = other.EXT1;                                        //stores fast dSASA
  COORDS = other.COORDS;
  SHELL_BURIED = other.SHELL_BURIED;
  AREA_BURIED_BY_ATOM_vector = other.AREA_BURIED_BY_ATOM_vector; //holds set of atoms that cause a certain buried area, uses position in pdb vector
  AREA_BURIED_BY_ATOM_vector_valid = other.AREA_BURIED_BY_ATOM_vector_valid;    //holds the valid sets
  AREA_BURIED_BY_ATOM_area = other.AREA_BURIED_BY_ATOM_area;                    //area corresponding to each set

  INTERACTION_P = other.INTERACTION_P;                                          //holds interacting atoms as positions in the pdb vector
  INTERACTION_SASA_P = other.INTERACTION_SASA_P;                                //holds interacting atoms that cause sasa loss as positions in the pdb vector
  CONTACT_AREA = other.CONTACT_AREA;                                            //holds the contact area caused by atom ID

  DISTANCES = other.DISTANCES;                                                  //holds atom distances, uses position in pdb vector
  ACCESSIBLE_P = other.ACCESSIBLE_P;
  ACTIVE = other.ACTIVE;
  SHELL = other.SHELL;
  ov_table = other.ov_table;                                                    //contains all valid overlaps and the atom positions that make them
  ov_table_area = other.ov_table_area;                                          //holds areas of overlaps in ov_table
  ov_norm_area = other.ov_norm_area;                                            //holds areas of normalized overlaps (divided by number of participating atom
  ATOM_TYPE = other.ATOM_TYPE;
  ATOM_TYPE_40 = other.ATOM_TYPE_40;
  OCCUPANCY = other.OCCUPANCY;
  TFACTOR = other.TFACTOR;
  ALTLOC = other.ALTLOC;
  CHARGE = other.CHARGE;
  BONDED = other.BONDED;
  ENERGY = other.ENERGY;
  MOL_TYPE = other.MOL_TYPE;
  INT_NUM = other.INT_NUM;
  DD_ENERGY = other.DD_ENERGY;
  AS_ENERGY = other.AS_ENERGY;
  CS_ENERGY = other.CS_ENERGY;
  sDD_ENERGY = other.sDD_ENERGY;
  BSA_ENERGY = other.BSA_ENERGY;
  DD_DISTANCES = other.DD_DISTANCES;
  TERMINAL = other.TERMINAL;
  DD_INTERACTING = other.DD_INTERACTING;
  EF_INT_NUM = other.EF_INT_NUM;
  HETATM = other.HETATM;
  DTYPE = other.DTYPE;
  return *this;
}

atom_struct::atom_struct(const atom_struct& other){
  STRUCTURE = other.STRUCTURE;
  ID = other.ID;
  NAME = other.NAME;
  RESN = other.RESN;
  CHAIN = other.CHAIN;
  ELEMENT = other.ELEMENT;
  STRUCT_TYPE = other.STRUCT_TYPE;
  RESI = other.RESI;
  iCODE = other.iCODE;
  VDW = other.VDW;
  RADIUS = other.RADIUS;
  RADIUS2 = other.RADIUS2;
  AREA = other.AREA;                                                            //Full sphere area, vdw+probe radius
  SASA = other.SASA;                                                            //exposed area
  EXT0 = other.EXT0;                                                            //currently used to store the dSASA calculated from overlap discovering
  EXT1 = other.EXT1;                                                            //stores fast dSASA
  COORDS = other.COORDS;
  SHELL_BURIED = other.SHELL_BURIED;
  AREA_BURIED_BY_ATOM_vector = other.AREA_BURIED_BY_ATOM_vector;                //holds set of atoms that cause a certain buried area, uses position in pdb vector
  AREA_BURIED_BY_ATOM_vector_valid = other.AREA_BURIED_BY_ATOM_vector_valid;    //holds the valid sets
  AREA_BURIED_BY_ATOM_area = other.AREA_BURIED_BY_ATOM_area;                    //area corresponding to each set
  INTERACTION_P = other.INTERACTION_P;                                          //holds interacting atoms as positions in the pdb vector
  INTERACTION_SASA_P = other.INTERACTION_SASA_P;                                //holds interacting atoms that cause sasa loss as positions in the pdb vector
  CONTACT_AREA = other.CONTACT_AREA;                                            //holds the contact area caused by atom ID

  DISTANCES = other.DISTANCES;                                                  //holds atom distances, uses position in pdb vector
  ACCESSIBLE_P = other.ACCESSIBLE_P;
  ACTIVE = other.ACTIVE;
  SHELL = other.SHELL;
  ov_table = other.ov_table;                                                    //contains all valid overlaps and the atom positions that make them
  ov_table_area = other.ov_table_area;                                          //holds areas of overlaps in ov_table
  ov_norm_area = other.ov_norm_area;
  ATOM_TYPE = other.ATOM_TYPE;
  ATOM_TYPE_40 = other.ATOM_TYPE_40;
  OCCUPANCY = other.OCCUPANCY;
  TFACTOR = other.TFACTOR;
  ALTLOC = other.ALTLOC;
  CHARGE = other.CHARGE;
  BONDED = other.BONDED;
  ENERGY = other.ENERGY;
  MOL_TYPE = other.MOL_TYPE;
  INT_NUM = other.INT_NUM;
  DD_ENERGY = other.DD_ENERGY;
  AS_ENERGY = other.AS_ENERGY;
  CS_ENERGY = other.CS_ENERGY;
  sDD_ENERGY = other.sDD_ENERGY;
  BSA_ENERGY = other.BSA_ENERGY;
  DD_DISTANCES = other.DD_DISTANCES;
  TERMINAL = other.TERMINAL;
  DD_INTERACTING = other.DD_INTERACTING;
  EF_INT_NUM = other.EF_INT_NUM;
  HETATM = other.HETATM;
  DTYPE= other.DTYPE;
}

string
atom_struct::sID(){
  stringstream sID;
  sID << NAME << ALTLOC << "/" << RESN << "/" << CHAIN << "/" << RESI;
  if (!iCODE.empty()){ sID << iCODE;}
  sID << "/" << ID;
  string sid = sID.str();
  return sid;
}

string
atom_struct::rsID(){
  stringstream rsID;
  rsID << RESN << "/" << CHAIN << "/" << RESI;
  if (!iCODE.empty())  rsID << iCODE;
  string rsid = rsID.str();
  return rsid;
}

string atom_struct::print(){
  stringstream tmp;
  string chain = CHAIN;
  string result;
  if( CHAIN.empty()) chain = "_";

  string t = "\t";
  tmp    << STRUCTURE << t<< ID << t << NAME << t << ALTLOC << t << RESN << iCODE << t
         << chain << t << RESI << t << SASA << t <<MOL_TYPE << t <<EXT0 <<t  << EXT1 << t << ATOM_TYPE<< t  << INT_NUM << t << ENERGY;
  result = tmp.str();
  return result;
}
string atom_struct::print2(){
  stringstream tmp;
  tmp << print();
  for( auto i : BONDED){
    tmp << "\t" << i+1;
  }
  string result = tmp.str();
  return result;
}


atom_struct::~atom_struct()
{
}
