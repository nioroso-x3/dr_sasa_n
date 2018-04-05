class residue_struct
{
public:
	residue_struct();
	residue_struct(string, string, string, string,int);
	string STRUCTURE;
	string RESN;
	string CHAIN;
	string STRUCT_TYPE;
	string sID;
	int RESI;
	float SASA;
    float dSASA;
    vector<uint32> ATOMS;
	vector<vector<uint32>> ov_table;       //overlaps by atom
	vector<float> ov_table_area;         
	vector<float> ov_table_area_norm;
	vector<vector<tuple<string,string,int>>> res_ov_table;
	vector<tuple<string,string,int>> INTERACTING;
	vector<float> INTERACTING_RES_dSASA;
	~residue_struct();
};
