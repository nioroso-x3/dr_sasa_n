// sr_sasa.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include "atom_struct.h"
#include "histogram.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverOnTheFly.h"
#include "SurfaceSolverCL.h"
#include "SolverDataProcessing.h"
#include "NB.h"
#include "NonEffective.h"
#include "SearchFunctions.h"
//#include "ShapeComplementarity.h"

string help = "Usage:\n"
" *Simple SASA solver: (mode 0, default)\n"
"Calculates SASA and outputs a NACCESS PDB like file and a second file with categorized atoms.\n"
"Outputs just add a .asa and .atmasa to the input filename if not defined by the user.\n"
"The .asa file has the SASA values in the B-factor column, and a internal atom typing in the charge column reserved for future use.\n"
"EXAMPLE:\n"
"  ./dr_sasa -m 0 -i 4ins.pdb -o 4ins.asa\n\n "
" *Delta SASA by chain ID or automatic: (mode 1)\n"
"Calculates the delta SASA in various objects contained in a single pdb or mol2 file.\n"
"Objects are currently defined only by their chains.\n"
"Outputs an interaction table, all surface overlaps, and a dSASA matrix for each\n"
"permutation of defined objects.\n"
"If no chains are selected the interactions will be defined by molecular type\n"
"EXAMPLE:\n"
"  ./dr_sasa -m 1 -i 4ins.pdb -chain AB -chain CD\n\n"
"  ./dr_sasa -m 1 -i 1bl0.pdb\n\n"
" *Aminoacid dSASA mode: (mode 2)\n"
"Calculates the delta SASA of all aminoacids inside a single object. (intramolecular contacts)\n"
"Outpus an interaction table and overlap table.\n"
"EXAMPLE:\n"
"  ./dr_sasa -m 2 -i 4ins.pdb -chain ABCD\n\n"
" *Atom dSASA mode: (mode 3)\n"
"Calculates the delta SASA of all atoms inside a single object. (intramolecular contacts)\n"
"Outpus an interaction table and overlap table.\n"
"EXAMPLE:\n"
"  ./dr_sasa -m 3 -i 4ins.pdb -chain ABCD\n\n"
" *Switches\n\n"
"-nomatrix\tswitch will disable matrix output.\n\n"
"-r float\tswitch will set the water probe radius in Angstroms. Default value is 1.4. Setting to 0 is equal to using the molecular surface.\n\n"
"-v\tAllows the user to define their own VdW radii for PDBs or MOL2 files. Examples are distributed under the utils folder.\n\n"
;

vector<char>
float2bool(vector<float> mtrx){
  vector<char> result(mtrx.size());
  for (uint32 i = 0; i < mtrx.size(); ++i){
    result[i] = mtrx[i] > 0.0 ? true : false; 
  }
  return result;
}

void
Fixrna2dna(vector<atom_struct>& pdb){
  for(auto& atom: pdb){
    if(atom.MOL_TYPE == "RNA") {
      atom.MOL_TYPE = "DNA";
    }
  }
}

string trim(string value){
    value.erase(value.find_last_not_of(" \n\r\t")+1);
    return value;
}

vector<vector<float> > GetHSTparam(string file){
  vector<vector<float> > result;
  ifstream input(file);
  string line;
  while(getline(input,line)){
    stringstream iss(line);
    float bins,xmax,xmin;
    if(!(iss >> bins >> xmax >> xmin)) break;
    vector<float> p = {bins,xmax,xmin};
    result.push_back(p);
  }
  return result;
}

vector<string> GetLines(string file){
  ifstream input(file);
  vector<string> result;
  string line;
  while(getline(input,line)){
    stringstream iss(line);
    string pdb;
    if(!(iss >> pdb)) break;
    pdb = trim(pdb);
    if (pdb.size()) result.push_back(pdb);
  }
  return result;
}

vector<char>
MatrixAND(vector<char>& Imtrx1,
      vector<char>& Imtrx2){
  vector<char> result;
  if (Imtrx1.size() == Imtrx2.size()){
    auto l = Imtrx1.size();
    result.resize(l);
    for (uint32 i = 0; i < l; i++){
      result[i] = (Imtrx1[i] && Imtrx2[i]);
    }
  }
  return result;
}

int main(int argc, char* argv[])
{
  //control variables
  uint32 mode = 0;
  uint32 maxin = 8;
  uint32 maxout = 8;
  uint32 outs = 0;
  float probe = 1.4;
  vector<string> inputs;
  vector<string> outputs;
  string vdwfile = "";
  string p_file;
  string types;
  vector<vector<string>> chain_sep;
  bool mtrx = true;
  int cl_mode = 0;
  bool keepunknown = true;
  int mtrxtype = 0;
  
  for (uint32 i = 0; i < maxout; ++i){
  #ifdef __WIN32
  outputs.push_back("NUL");
  #else
  outputs.push_back("/dev/null");
  #endif
    //cout << outputs.back() << "\n";
  }

  for (int i = 0; i < argc; ++i){
    string c = string(argv[i]);
    //cout << c << "\n";
      //mode
    try{
      if(c == "-m"){
        if(i+1 < argc) mode = stoi(string(argv[i+1]));
      }
      //input filenames
      if(c == "-i"){
        if(i+1 < argc && inputs.size() < maxin){
          try{
            inputs.push_back(string(argv[i+1]));
          }
          catch(exception e){
            cout << "input error";
          }
          //cout << inputs.back() << "\n";
        }
      }
      //output filenames
      if (c == "-o"){
        if(i+1 < argc && outs < maxout){
          outputs[outs++] = (string(argv[i+1]));
          //cout << outputs.back()<< "\n";;
        }
      }
      //vdw radii filename
      if (c == "-v"){
        if(i+1 < argc && vdwfile.empty()){
          vdwfile = string(argv[i+1]);
          //cout << vdwfile << "\n";
        }
      }
      //unit sphere points file
      if(c == "-s"){
        if(i+1 < argc && p_file.empty()){
          p_file = string(argv[i+1]);
        }
      }
      //probe size radius
      if(c == "-r"){
        if(i+1 < argc) probe = stod(string(argv[i+1]));
      }
      //atom type definitions
      if (c == "-t"){
        if(i+1 < argc && types.empty()) types = argv[i+1];
      }
      //object separator by chain
      if(c == "-chain"){
        if(i+1 < argc){
          string str = string(argv[i+1]);
          vector<string> result;
          for (uint32 c = 0; c < str.length();c++){
            //cout <<i << " "  <<str.substr(c,1) << "\n";
            result.push_back(str.substr(c,1));
          }
          chain_sep.push_back(result);
        }
      }
      if(c == "-h" || c == "--help" || argc == 1){
        cout << help;
        return 0;
      }
      //disable matrix output
      if (c == "-nomatrix") mtrx = false;
      if (c == "-usegpu"){
        if (i + 1 < argc) cl_mode = stoi(string(argv[i + 1])) + 1;
      }
      if (c == "-nohetatm"){
        keepunknown = false;
      }
      if (c == "-outtype"){
        mtrxtype = 1;
      }
    }
    catch(exception e){
      cerr <<"Exception in parsing argument \""<< c<< "\"\n";
      return -1;
    }
  }
  if (inputs.size() == 0) {
    cerr << "No input.\n";
    return -1;
  }

  //basic sasa solver

  if (mode != 100){
    cout << inputs[0] << "\n";
    if(mode == 0){
      if (chain_sep.size()>1){
        cerr << "Please define a single object.\n";
        return 0;
      }
    }
    string input = inputs[0];
    string output = inputs[0]+".asa";
    string splitasa = inputs[0]+".atmasa";
    if (mode == 0){
      if(outs >= 1) output = outputs[0];
      if(outs >= 2) splitasa = outputs[1];
    }
    VDWcontainer rad(vdwfile);
    rad.GenPoints();
    auto pdb = PDBparser(input,"",keepunknown);
    if(mode == 0){
      if(chain_sep.size() == 1) ChainSelector(chain_sep,pdb);
    }
    else{
      if(chain_sep.size() == 1) ChainSelector(chain_sep,pdb);
      else if(chain_sep.size() > 1){
        vector<vector<string>> sele;
        vector<string> tmp;
        for(auto& item : chain_sep){
          for (auto& c : item) tmp.push_back(c);
        } 
        sele.push_back(tmp);
        ChainSelector(sele,pdb);
      }
    }
    rad.SetRadius(pdb, probe);
    SolveInteractions(pdb,0);
    SimpleSolverCL(pdb,rad.Points,cl_mode);
    PrintSASAResults(pdb,output);
    PrintSplitAsaAtom(pdb,splitasa);
    if (mode == 0) return 0;
  }
  //relative ASA
  if(mode==100){
    if (chain_sep.size()>1){
      cerr << "Please define a single object.\n";
      return 0;
    }
    string input = inputs[0];
    string output = outputs[0];
    string splitasa = outputs[1];
    VDWcontainer rad;
    rad.GenPoints();
    auto pdb = PDBparser(input,"",keepunknown);
    if(chain_sep.size() == 1)
      ChainSelector(chain_sep,pdb);
    rad.SetRadius(pdb, probe);
    SolveInteractions(pdb,0);
    SimpleSolverCL(pdb,rad.Points,cl_mode);
    RelativeSASA(pdb);
    cout << input << " done:\n";
    PrintSASAResults(pdb,output);
    PrintSplitAsaAtom(pdb,splitasa); 
    return 0;
  }

  //Simple ASA solver, output is by atom type definition
  if(mode==106){
    if (chain_sep.size()!=1){
      cerr << "Please define a single object.\n";
      return 0;
    }
    string input = inputs[0];
    string output1 = outputs[0];
    string output2 = outputs[1];

    VDWcontainer rad(vdwfile);
    rad.GenPoints();
    auto pdb = PDBparser(input,types,keepunknown);
    ChainSelector(chain_sep,pdb);
    rad.SetRadius(pdb, probe);
    SolveInteractions(pdb,0);
    SimpleSolverCL(pdb,rad.Points,0);
    PrintSASAResults_per_type(pdb,
                              output1,
                              output2);
    return 0;
  }

 
  //Generic structure dSASA
  if(mode==1){
    int Imode = 0;
    if (chain_sep.size() <= 1){
      Imode = 4;
      cout << "#Automatic interaction solver selected.\n";
    }
    else{
      Imode = 1;
      cout << "#Manual chain interaction solver selected.\n";
    } 
    stringstream stdinput;
    stdinput << inputs[0];
    if(chain_sep.size() == 1){
      stdinput << ".";
      for (uint32 i = 0; i < chain_sep.size(); ++i){
        for (auto c : chain_sep[i]){
          stdinput << c;
        }
        if (i < chain_sep.size() - 1) stdinput << "_vs_";
      }
    }
    else if (chain_sep.size() >= 2){
      stdinput << ".";
      for (uint32 i = 0; i < chain_sep.size(); ++i){
        for (auto c : chain_sep[i]){
          stdinput << c;
        }
        if (i < chain_sep.size() - 1) stdinput << "_vs_";
      }
    }
    string vsinput = stdinput.str();
    string input = inputs[0];
    string output1 = vsinput+".int_table";
    string output2 = vsinput+".overlaps";
    string output3 = input;
    string output4 = input+".dsasa";
    VDWcontainer rad(vdwfile);
    rad.GenPoints();
    auto pdb = PDBparser(input,types,keepunknown);
    //
    //for(auto item : pdb) cout << item.print() << "\n";
    //
    if(chain_sep.size() >= 1) ChainSelector(chain_sep,pdb);

    rad.SetRadius(pdb, probe);

    Generic_Solver(pdb,rad.Points,chain_sep,Imode,cl_mode);

    GeneratePairInteractionData(pdb);

    //PrintDNA_ProtResultsTable(pdb, output1);
    PrintDNA_ProtResults(pdb, output2);
    PrintDSASAResults(pdb,output4);
    if(mtrx){ 
      PrintDNA_ProtResultsByAtomMatrix(pdb, output3,mtrxtype);
      Print_MatrixInsideAtom(pdb,vsinput);
    }

    return 0;
  }
 

//dSASA inside aminoacids mode
  if(mode == 2){
    if (chain_sep.size() > 1){
      cerr << "Please define a single object.\n";
      return 0;
    }
    stringstream stdinput;
    stdinput << inputs[0] << ".internal.";
    if (chain_sep.size() == 1){
      for (auto c : chain_sep[0]){
        stdinput << c;
      }
    }
    else stdinput << "all";
    stdinput << ".by_res";
    string vsinput = stdinput.str();
    string input = inputs[0];
    string output1 = vsinput+".int_table";
    string output2 = vsinput+".overlaps";

    VDWcontainer rad(vdwfile);
    rad.GenPoints();
    auto pdb = PDBparser(input,types,keepunknown);
    if (chain_sep.size() == 1) ChainSelector(chain_sep,pdb);

    rad.SetRadius(pdb, probe);

    Generic_Solver(pdb,rad.Points,chain_sep,2,cl_mode);

    GeneratePairInteractionData(pdb);

    //PrintDNA_ProtResultsTable(pdb, output1);
   
    if(mtrx)Print_MatrixInsideAtom(pdb,vsinput);
    else PrintDNA_ProtResults(pdb, output2);
    return 0;
  }
  //dSASA inside atoms mode
  if(mode == 3){
    if (chain_sep.size() > 1){
      cerr << "Please define a single object.\n";
      return 0;
    }
    stringstream stdinput;
    stdinput << inputs[0] << ".internal.";
    if (chain_sep.size() == 1){
      for (auto c : chain_sep[0]){
        stdinput << c;
      }
    }
    else stdinput << "all";
    stdinput << ".by_atom";
    string vsinput = stdinput.str();
    string input = inputs[0];
    string output1 = vsinput+".int_table";
    string output2 = vsinput+".overlaps";

    VDWcontainer rad(vdwfile);
    rad.GenPoints();
    auto pdb = PDBparser(input,types,keepunknown);
    if (chain_sep.size() == 1) ChainSelector(chain_sep,pdb);

    rad.SetRadius(pdb, probe);

    Generic_Solver(pdb,rad.Points,chain_sep,3,cl_mode);

    GeneratePairInteractionData(pdb);

    //PrintDNA_ProtResultsTable(pdb, output1);
    if(mtrx)Print_MatrixInsideAtom(pdb,vsinput);
    else PrintDNA_ProtResults(pdb, output2);
    return 0;
  }
  cout << "Invalid mode\n"; 
  return 0;
}
