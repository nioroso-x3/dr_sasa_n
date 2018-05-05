void
RelativeSASA(vector<atom_struct>&);

void
PrintDNA_ProtResults(vector<atom_struct>&,       //pdb struct
           string);             //filename

void
PrintDNA_ProtResultsTable(vector<atom_struct>&,  // pdb struct
              string);               //filename

void
PrintDNA_ProtResultsByAA(vector<atom_struct>&,  // pdb struct
             string);

void
PrintDDMatrix(vector<atom_struct>&,
              vector<char>&,
              string, 
              string,
              string);


void 
PrintTMatrix(vector<int>&,int,int,int,string);
void 
PrintTMatrix(vector<float>&,int,int,int,string);



// void
// PrintDNA_ProtResultsMatrix(vector<atom_struct>&,  // pdb struct
//                           string);               // output filename

void
GenerateInterBSAMatrix(vector<atom_struct>&                 ,
                        map<vector<string>, vector<float>>& ,
                        map<vector<string>, vector<float>>& ,
                        map<vector<string>, vector<string>>& ,
                        map<vector<string>, vector<string>>& ,
                        map<vector<string>, vector<string>>& ,
                        map<vector<string>, vector<string>>& ,
                        map<vector<string>, vector<uint32>>&    ,
                        map<vector<string>, vector<uint32>>&   );

void
GenerateIntraBSAMatrix(vector<atom_struct>& ,
                        vector<float>&      ,
                        vector<float>&      ,
                        vector<string>&      ,
                        vector<string>&      ,
                        vector<string>&      ,
                        vector<string>&      ,
                        vector<uint32>&      ,
                        vector<uint32>&     );

void
PrintDNA_ProtResultsByAtomMatrix(vector<atom_struct>&,  // pdb struct
                                 string,                // output filename
                                 int);
                                 

void
PrintSASAResults(vector<atom_struct>&,         //pdb struct
         string);             //filename
void
PrintDSASAResults(vector<atom_struct>&,
                 string               );



void
PrintSASAResults_per_type(vector<atom_struct>&,         //pdb struct
              string,            //protein sasa
              string);                     //dna sasa

void
PrintSASAResults_per_type_and_res(vector<atom_struct>&,         //pdb struct
              string,
              string);                     //filename

void
Print_MatrixInsideAtom(vector<atom_struct>& pdb,
                       string);

/*void
PrintSASAone_type_(vector<atom_struct>&,         //pdb struct
           int,                           //type
           string);                       //output file */


void
CalculateDNA_ProtInteractions(vector<atom_struct>&);  //pdb struct

string                                         //returns string with data
PrintDNA_ProtInteractions(vector<atom_struct>&,  //pdb struct
                          uint32);          //position of atom to query

void
PrintAARawData (map<string,vector<float> >&,  //data to print
                string);             //filename

void
PrintAAxNUCRawData(map<pair<int,int>, vector<float> >&, //data to print
                   string);                 //filename

// void
// PrintAAdSASASequence(vector<atom_struct>&,    //data to pring
//            string);          //filename

void
PrintAtomxNucData(map<pair<int,int>, vector<float> >&, //data to print
                  string);                //filename

void
GeneratePDistribution(map<pair<int,int>,vector<float> >&,  //dsasa type1 <- type2
                      map<pair<int,int>,histogram>&,   //histogram
                      uint32,                           //Bins
                      float,                           //Xmax
                      string);                          //Output file

// void
// AddRawAAData(map<string,vector<float> >&,
//        vector<atom_struct>& );

void
AddRawAAxNUCData(map<pair<int,int>,vector<float> >&, //standard dsasa
                 map<pair<int,int>,vector<float> >&, //weighted dsasa
                 vector<atom_struct>&,                 //pdb structure
                 vector<int>&);

void
AddRawAAxNUCDataWithExcluding(map<pair<int,int>,vector<float> >&, //standard dsasa
                              map<pair<int,int>,vector<float> >&, //weighted dsasa
                              vector<atom_struct>&,                 //pdb structure
                              vector<int>&);
void
AddRawAAxNUCOverlapData(map<pair<int,int>,vector<float> >&, //standard dsasa
            map<pair<int,int>,vector<float> >&, //weighted dsasa
            vector<atom_struct>& );         //pdb struct

//deprecated
void
AddRawAtomxNucData(map<pair<int,int>,vector<float> >&, //protein data
                   map<pair<int,int>,vector<float> >&, //dna data
                   vector <atom_struct>&);        //pdb struct

void
GeneratePairInteractionData(vector<atom_struct>&);     //Generates pair interactions in pdb struct

void
CountInteractionsP(float&,                  //reference to value holding amount of interacting protein atoms
                   float&,           //reference to value holding amount of interacting xna atoms
                   vector<atom_struct>&);    //pdb struct
void
PrintSplitAsaAtom(vector<atom_struct>&,string);

void
PrintSplitAsaRes(vector<atom_struct>&,string);

