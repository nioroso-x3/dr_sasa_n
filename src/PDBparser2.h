string
GetBasename(string);

vector<vector<string> > //ordered strings for each category in the atom tags
ReadPDB_ATOM(string); //filename

void
RemoveAltLocs(vector<atom_struct>&); //removes altlocs from intermediary list

void
GetMolType(vector<atom_struct>&);

void
GetAtomType(string,  //type dictionary fname
            vector<atom_struct>&,            //vector of atoms
            bool);                          //keep unknown type atoms

vector<atom_struct> //final structure
PDBparser(string, //pdb filename
          string, //type filename
          bool);  //keep unknown type atoms
//read type map from text file
map<string, map<string, int> >
GetTypeMap(string);

void
ChainSelector(vector<vector<string>>&,
              vector<atom_struct>&);
string
FixRESN(string);

string
FixNAME(string);

void 
FixRepeatID(vector<atom_struct>&);

bool
IsRESNSolvent(string);

bool
IsRESNIon(string);

int
InputType(string);

vector<atom_struct>
MOL2parser(string,string,bool); //mol2 filename

void 
MOL2_parse_map(map<string,vector<string>>& sections,
               vector<atom_struct>& result,
               unsigned int& CHAINn);
