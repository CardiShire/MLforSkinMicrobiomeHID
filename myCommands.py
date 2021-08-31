"""This file can take all files in a directory.
Input list of filenames.
Output is a txt file of all possible combinations of unordered pairs of file input
in the form of a unix command line to analyze with samtools mpileup
and output an Fst file in .gz."""


import sys
import itertools
import os
import datetime

def simplifyFilename(filename):
	'''
	This function takes in a filename
	with a path (or not)
	and it returns the filename
	without this path AND
	without the file externsion(s)
	'''
	whatevs = os.path.basename(filename) #remove the path.
	i = whatevs.find(".") #find finds the "." from the left
	#rfind finds the "." from the right
	if i < 0: #if there is no "." then it returns the file name as is
		return whatevs
	return whatevs [:i] #returns file name until the location of "."
	
#ignoring args 0 which is arguments
arguments = sys.argv[1:]

d = datetime.datetime.today()


#for all combinations of arguments of length 2
for combo in itertools.combinations(arguments, 2):
	f1 = simplifyFilename(combo[0])
	f2 = simplifyFilename(combo[1])
	
	print("samtools mpileup ", combo[0], " ", combo[1], 
		" | ./fstiminator.py | gzip -9 > " , f1 , "_" , f2 , 
		"_", str(d).split()[0], ".fst.gz" , sep="")



