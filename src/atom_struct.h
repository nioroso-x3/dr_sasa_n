class atom_struct
{
public:
  atom_struct(uint32,        //id
              int,            //resi
              string,         //icode
              string,         //name
              string,         //resn
              string,          //chain
              string,          //elem
              string,           //obj name
              string,          //mol type
              float x,
              float y,
              float z,
              string,            //altloc
              float,           //occu
              float,             //tf
              string);              //charge

  string STRUCTURE;
  uint32 ID;
  string NAME;
  string ALTLOC;
  string RESN;
  string CHAIN;
  string ELEMENT;
  string STRUCT_TYPE;
  string MOL_TYPE;
  int RESI;
  string iCODE;
  float VDW;
  float RADIUS;
  float RADIUS2;
  float AREA;                                        //Full sphere area, pi*(vdw+probe)**2
  float SASA;                                        //exposed area
  float EXT0;                                        //currently used to store the dSASA calculated from overlap discovering
  float EXT1;                                        //stores fast dSASA
  vector<float> COORDS;
  vector<char> SHELL_BURIED;
  vector<vector<uint32> > AREA_BURIED_BY_ATOM_vector; //holds set of atoms that cause a certain buried area, uses position in pdb vector
  vector<uint32> AREA_BURIED_BY_ATOM_vector_valid;    //holds the valid sets
  vector<float> AREA_BURIED_BY_ATOM_area;            //area corresponding to each set
  vector<uint32> INTERACTION_P;                       //holds interacting atoms as positions in the pdb vector
  vector<uint32> INTERACTION_SASA_P;                  //holds interacting atoms that cause sasa loss as positions in the pdb vector
  map<uint32,float> CONTACT_AREA;                    //holds the contact area caused by atom ID
  map<uint32,float> DISTANCES;                       //holds atom distances, uses position in pdb vector
  map<uint32,float> DD_DISTANCES;                    //holds atom distances, uses position in pdb vector, for distance dependent potential
  uint32 ACCESSIBLE_P;
  bool ACTIVE;
  vector<float> SHELL;
  vector<vector<uint32>> ov_table;                    //contains all valid overlaps and the atom positions that make them
  vector<float> ov_table_area;                       //holds areas of overlaps in ov_table
  vector<float> ov_norm_area;
  int ATOM_TYPE;
  int ATOM_TYPE_40;
  float OCCUPANCY;
  float TFACTOR;
  string CHARGE;
  int INT_NUM;
  int EF_INT_NUM;
  float ENERGY;
  float DD_ENERGY;
  float AS_ENERGY;
  float CS_ENERGY;
  float sDD_ENERGY;
  float BSA_ENERGY;
  bool TERMINAL;
  bool HETATM;
  int POLAR;
  string DTYPE;
  vector<uint32> BONDED;
  vector<uint32> DD_INTERACTING;
  ~atom_struct();
  void clear();
  string sID();
  string rsID();
  string print();
  string print2();
  atom_struct(const atom_struct& other);
  atom_struct& operator=(const atom_struct& other);
};
