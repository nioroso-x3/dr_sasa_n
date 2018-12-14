vector<atom_struct*>
GetAtoms(string,
         string,
         int,
				 string,
				 vector<atom_struct>&);

vector<vector<vector<atom_struct*>>>
OrganizePDB(vector<atom_struct>&);

vector<atom_struct>
ReorderPDB(vector<atom_struct>&);

/*vector<atom_struct*>
GetAtomsFromRESI(int,
				 vector<residue_struct>&);*/

/*vector<residue_struct>
GetRESN(string,
		vector<atom_struct>&);*/

/*vector<residue_struct>
GetRESN(string,
		vector<residue_struct>&);*/

map<string,int>
GetTypeFromRESI(string,
                string,
                int,
								string,
								vector<atom_struct>);

vector<vector<string>>
ReadTabSeparatedFile(string);

vector<vector<string>>
ReadTabSeparatedString(string);


vector<string> 
LineSplitter(string,
             char);
