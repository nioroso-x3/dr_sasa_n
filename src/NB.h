vector<vector<string>>
NB_GenPerm(vector<string>&); //Generates permutations of a N length array, with the
                             //first element being the center

//checks if atoms are bonded
bool
NB_IsBonded(uint32,                // atom i
            uint32,                // atom j
            vector<atom_struct>&); // atom vector

//Links bonded atoms in a atom vector
void
NB_LinkPDB(vector<atom_struct>&);

//Gets all bonded and non-bonded atoms up to a certain depth.
vector<vector<uint32>>
NB_GetBonded(uint32,                // root atom
             uint32,                // max depth 0: directly bonded atoms
             vector<atom_struct>&); // structure


void
NB_RemoveBonded(vector<atom_struct>&,   //pdb
                int,                    //select bonding depth
                vector<char>&);         //interaction matrix