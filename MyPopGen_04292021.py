"""
This function (FST) takes in two mpileup strings and computes Hudson's 
estimator of Fst. 
The two strings must be from two (microbiome) populations.
It returns a tuple of:
The heterozygosity within population 1 (Hw_pop1)
The heterozygosity within population 2 (Hw_pop2)
The heterozygosity between the two populations (Hb)
The Fst estimate (FST_est)
"""
import sys
# p1 = "AAAAA+13XXXXXXXXXXXXX-12RRRRRRRRRRRR^A"print
# p2 = "TTTTTTTT"



#Calculating the number of pairwise differences within a population. 
	#Computes heterozygosity within each population. 
#Defining a function that determines the pi for pop1 and pop2.
def FST(p1, p2):
	pop1 = p1.lower()
	pop2 = p2.lower()
	nucCounter1 = {}
	i1 = 0
  
	#while i (idex) is less than the length of the 
	#string pop1 the loop will continue
	while i1 < len(pop1):
		#c (character) is equal to the index of the string (pop1) 
		c1 = pop1[i1]
		if c1 == '+' or c1 == '-':
			#at this point, pop1[i] is guaranteed to be numeric
			i1 += 1
			#make a list to add each character between 0-9 
			list_int1 = []
			#if c is greater than '0' (character in ord() that is equal to 0)
			#or equal to or less than '9' (equal to 9) then...
			while pop1[i1] >= '0' and pop1[i1] <='9': #this statement should be true at least once
				list_int1.append(pop1[i1]) #keep track of the integers we have seen
				i1 += 1 #and now go one more index in, is the next character an integer too?
			c1_int = int(''.join(list_int1))
				#add the integer (c_int) to the index number plus 1 to skip 
				#the indel and not add it to the dictionary
			i1 += c1_int
		#if c for a particular index equals ^, $, +, or - then add one to 
		#the index and do not add the character to the dictionary nucCounter1
		elif c1 == "^":
			i1 += 2
		elif c1 == "$":
			i1 += 1
		#if c is not equal to any of the above and is already in the 
		#dictionary then add 1 to the value associated with the c's key
		elif c1 in nucCounter1:
			nucCounter1[c1] += 1
			i1 += 1
		#if c is not equal to any of the above and it is NOT in the 
		#dictionary then add a new key (c) to the dictionary with a value
		#of one   
		else:
			nucCounter1[c1] = 1
			i1 += 1
	#Turn the nucCounter1 dictionary into a list.		
	Listpop1 = list(nucCounter1.values())
	
	#Now do the same for pop2
	nucCounter2 = {}
	i2 = 0

	#while i (idex) is less than the length of the 
	#string pop1 the loop will continue
	while i2 < len(pop2):
		#c (character) is equal to the index of the string (pop1) 
		c2 = pop2[i2]
		if c2 == '+' or c2 == '-':
			#at this point, pop1[i] is guaranteed to be numeric
			i2 += 1
			#make a list to add each character between 0-9 
			list_int2 = []
			#if c is greater than '0' (character in ord() that is equal to 0)
			#or equal to or less than '9' (equal to 9) then...
			while pop2[i2] >= '0' and pop2[i2] <='9': #this statement should be true at least once
				list_int2.append(pop2[i2]) #keep track of the integers we have seen
				i2 += 1 #and now go one more index in, is the next character an integer too?
			c2_int = int(''.join(list_int2))
				#add the integer (c_int) to the index number plus 1 to skip 
				#the indel and not add it to the dictionary
			i2 += c2_int
		#if c for a particular index equals ^, $, +, or - then add one to 
		#the index and do not add the character to the dictionary nucCounter1
		elif c2 == "^":
			i2 += 2
		elif c2 == "$":
			i2 += 1
		#if c is not equal to any of the above and is already in the 
		#dictionary then add 1 to the value associated with the c's key
		elif c2 in nucCounter2:
			nucCounter2[c2] += 1
			i2 += 1
		#if c is not equal to any of the above and it is NOT in the 
		#dictionary then add a new key (c) to the dictionary with a value
		#of one   
		else:
			nucCounter2[c2] = 1
			i2 += 1
	#Turn the nucCounter2 dictionary into a list. 
	Listpop2 = list(nucCounter2.values())

	#Check to make sure the count for each populations is equal to or greater than 1.
	if len(nucCounter1) == 1 and len(nucCounter2) == 1:
	  #print("Hi")
	  if nucCounter1.keys() == nucCounter2.keys():
	    return None
		
	if len(nucCounter2) == 0:
		#print("The dictionary for pop2 is empty or only has one nucleotide.", file=sys.stderr)
		return None
	if len(nucCounter1) == 0:
	  return None
		
	Spop1 = 0
	Spop2 = 0
	
	#Total n for each population.
	Npop1 = sum(nucCounter1.values())
	Npop2 = sum(nucCounter2.values())
	
	if sum(nucCounter1.values()) <= 1 or sum(nucCounter2.values()) <= 1:
	  return None
	
	#Calculate site homozygosity
	for j in Listpop1:
		Spop1 += j*(j-1)/(Npop1*(Npop1-1))
	for i in Listpop2:
		Spop2 += i*(i-1)/(Npop2*(Npop2-1))
	#pi - the mean number of differences within a population
	#Calculate heterozygosity, 1 - homozygosity
	Hw_pop1 = 1 - Spop1
	Hw_pop2 = 1 - Spop2
	
	mean_Hw = (Hw_pop1 + Hw_pop2)/2
	
	#Calculate Hb - the mean number of differences between a population
	similarity = 0 
	for k, v in nucCounter1.items():
		if k in nucCounter2:
			similarity += v * nucCounter2[k]

	Hb = 1 - (similarity/(Npop1*Npop2))
	
	FST_est = 1 - (mean_Hw/Hb)
	#print(FST_est)
	return FST_est
  #return Hw_pop1, Hw_pop2, Hb, FST_est
	#print(FST_est)	

#Hudson FST
#FST = 1 - (Hw/Hb)
#Hw is mean number of differences between different sequences sampled 
	#from the same subpopulation
#Hb is the mean number of differences between sequences sampled from the
	#two different subpopulations sampled
#In other words:
	#Hw is the average pairwise comparison of sequences within subpopulations
	#Hb is the average pairwise comparison of sequences from two different
		#subpopulations sampled	
