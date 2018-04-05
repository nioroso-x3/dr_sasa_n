#!/bin/bash
for item in ./structures_110/*.pdb
do
    sr_sasa.bin 2 $item vdw.radii 1.4 thomson15092.xyz /dev/null /dev/null $item.dat $item.wdat
done
