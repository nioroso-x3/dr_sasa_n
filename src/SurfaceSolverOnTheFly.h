extern vector<uint32> PRIMES;
void
SimpleSolver_2(vector<atom_struct>&,
               vector<float>&);
void
DNA_Prot_Solver_2(vector<atom_struct>&,
                  vector<float>&);
void
SolveInteractions(vector<atom_struct>&,
                  uint32);
void
Generic_Solver(vector<atom_struct>&,
               vector<float>&,
               vector<vector<string>>,
               int,
               int);
void
SolveInteractionsALL(vector<atom_struct>&,
                     int,
                     vector<vector<string>>);

vector<char> 
SURF_GenerateImtrx(vector<atom_struct>&);