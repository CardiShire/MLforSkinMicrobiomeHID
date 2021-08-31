#!/usr/bin/env python3
import sys
from MyPopGen_04292021 import FST

args = sys.argv[1:]

print("My arguments are: ", args)

for line in sys.stdin:
	line_split = line.rstrip().split('\t')
	marker = line_split[0]
	base_coordinates = line_split[1]
	reads1 = line_split[3]
	reads2 = line_split[6]
	p1 = line_split[4]
	p2 = line_split[7]
	Fst = FST(p1, p2)
	if Fst != None:
		print(marker, base_coordinates, reads1, reads2, Fst)
