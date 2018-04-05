#!/usr/bin/python
import sys
#import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
if __name__ == "__main__":
	max_col = 6
	inputf = open(sys.argv[1])
	lines = [line.split() for line in inputf.readlines()]
	inputf.close()
	#overlap dictionary, key is ids that are part of overlap
	#structure stored is : (dSASA,
	#                       ((RESN1,RESI1,CHAIN1),(RESN2,RESI2,CHAIN2),...))
	overlaps = {}
	#atom dict, key is atom ID, structure stored is (RESN,RESI,CHAIN)
	atoms = {}
	for line in lines:
		if line[0] == "ATOM":
			atoms[int(line[1])] = (line[3],line[4],int(line[5]))
	
	for line in lines:
		if line[0] == "ATOM_CONTACT_SURFACE_OVERLAP":
			ids = (int(item) for item in line[2:-1])
			dSASA = float(line[-1])
			residues = set()
			for i in ids:
				residues.add(atoms[i])
			overlaps[ids] = (dSASA,residues)
	columns = [[] for i in range(max_col)]
	for item in overlaps.values():
		count = len(item[1])
		columns[count].append(item[0])
	for i in range(1,max_col):
		print str(i)+"\t"+ "\t".join([str(value) for value in columns[i]])
	
	#tries to draw a ugly 2d histogram

#	x = [len(item[1]) for item in overlaps.values()]
#	y = [item[0] for item in overlaps.values()]	
#	mybins = [[1,2,3,4,5,6],[i*0.2 for i in range(36)]]
#	H, xedges, yedges = np.histogram2d(x,y,bins=mybins)
#	plt.figure()
#	myextent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#	plt.imshow(H.T,origin='low',extent=myextent,interpolation='nearest',aspect='auto')
#	plt.colorbar()
#	plt.show()