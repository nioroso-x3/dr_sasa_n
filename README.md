# dr-sasa Solvent Accessible Surface Area calculation software for biological molecules

dr-sasa is a solvent accessible surface area calculation software for biological molecules, that supports proteins, dna and rna inputs. The input files can be in the PDB or MOL2 format. PDB format files will use a NACCESS compatible VdW radii lookup table, while MOL2 formats will use the same VdW radii used in Chimera.

## Operation Modes

There are three operation modes. The first does common ASA calculations, the second calculates 

## Usage

placeholder

## Compiling
After cloning the source, enter the directory and create a build directory:

###
```
cd $(drsasa_dir)
mkdir build
cd build
cmake ../
make
```

Then copy the binary to your local or system binaries folder.
