# dr-sasa - Solvent Accessible Surface Area calculation software for biological molecules

dr-sasa is a solvent accessible surface area calculation software for biological molecules, that supports proteins, DNA, RNA and ligand inputs. The input files can be in the PDB or MOL2 format. PDB format files will use a NACCESS compatible VdW radii lookup table, while MOL2 formats will use the same VdW radii used in Chimera.

If you use this software in your research, please acknowledge it by citing the following:

*Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) “Calculation of accurate contact surface areas
between atoms for the quantitative analysis of non-bonded molecular interactions”. Bioinformatics (submitted).*

## Operation Modes

There are five operation modes. The first does common ASA (Accessible Surface Area) calculations, the second calculates the BSA (Buried Surface Area) between various chains or molecular objects, and third and fourth calculate CSA (Contact Surface Area) between residues and atoms as if they were independent objects.  
The fifth mode does not use the SASA or BSA for its calculations. It simply calculates the raw contact surfaces between chains or molecular objects.

## Usage

##### -Simple SASA solver: (mode 0, default)  
Calculates SASA (Solvent Accessible Surface Area) and outputs a NACCESS PDB like file and a second file with categorized atoms.
Output files just add a .asa and .atmasa between the input filename extesion if not defined by the user.
The .asa file has the SASA values in the B-factor column, and a internal atom typing in the charge column reserved for future use.  
EXAMPLE:
###
```
  ./dr_sasa -m 0 -i 4ins.pdb -o 4ins.asa
```
###### -Delta SASA by chain ID or automatic: (mode 1)  
Calculates the delta SASA (change in Solvent Accessible Surface Area) or BSA in various objects contained in a single pdb or mol2 file.
Objects are currently defined only by their chain IDs.
Outputs an interaction table, all surface overlaps, and a delta SASA matrix for each
permutation of defined objects.
If no chains are selected the interactions will be defined by molecular type (ie. protein,DNA,RNA,ligand).
EXAMPLE:
###
```
  ./dr_sasa -m 1 -i 4ins.pdb -chain AB -chain CD

  ./dr_sasa -m 1 -i 1bl0.pdb
```
###### -Residue dSASA mode: (mode 2)  
Calculates the delta SASA of all residues (ie. aminoacis, nucleotides or ligands) inside a single object. (intramolecular contacts)
Outputs an interaction table and an overlap table.  
EXAMPLE:
###
```
  ./dr_sasa -m 2 -i 4ins.pdb -chain ABCD
```
###### -Atom dSASA mode: (mode 3)  
Calculates the delta SASA of all atoms inside a single object (ie. intramolecular contacts).
Outputs an interaction table and an overlap table.  
EXAMPLE:
###
```
  ./dr_sasa -m 3 -i 4ins.pdb -chain ABCD
```
###### -Atom contact surface mode: (mode 4)  
Calculates the contact surface between all the chains of a defined object or, if no chains are defined, between different molecular types.
Outputs interaction tables.  
EXAMPLE:
###
```
  ./dr_sasa -m 4 -i 1bl0.pdb
```
 -Switches  

-nomatrix switch will disable matrix output.

-r float  switch will set the water probe radius in Angstroms. Default value is 1.4. Setting to 0 is equal to using the molecular surface.

-v  Allows the user to define their own VdW radii for PDBs or MOL2 files. Examples are distributed under the utils folder, as vdw.radii and vdiw.radii.mol2.

-no_atmasa_autosort Special setting for atmasa output files. Disables the autosorter, useful for malformed mol2 or pdbs with atoms with missing chain identifers.
## Compiling
After cloning the source, enter the directory and create a build directory:

###
```
$ cd $(drsasa_dir)
$ mkdir build
$ cd build
$ cmake ../
$ make
```

Then copy the binary to your local or system binaries folder.
