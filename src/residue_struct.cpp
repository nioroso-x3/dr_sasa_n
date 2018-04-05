#include "stdafx.h"
#include "atom_struct.h"

#include "residue_struct.h"


residue_struct::residue_struct(string resn, string chain, string obj,string mol, int resi){

	RESI = resi;
	RESN = resn;
	CHAIN = chain;
	STRUCT_TYPE = mol;
	STRUCTURE = obj;
	SASA = 0;
	dSASA = 0;
	stringstream temp;
	temp << RESN << "/" << CHAIN << "/" << RESI;
	sID = temp.str();
}
residue_struct::residue_struct()
{
	SASA = 0;
	dSASA = 0;
}
residue_struct::~residue_struct()
{
}
