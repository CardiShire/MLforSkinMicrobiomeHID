# command line
# nohup python3  SNPtotals_nesteddict_up.py 'S*gz' > fstsAll.tsv 2> fstErrs &

import gzip
import sys
import glob

def print_inner_dictionary(dict_obj, marker):
    """
    This takes in a dictionary and prints it
    ...
    the inner dictionary of dOfD in particular
    it will print it as a TSV
    """
    for snp in dict_obj:
      print(marker, snp, dict_obj[snp][0], dict_obj[snp][1], sep="\t")

def append_value(dict_obj, snp, fst):
    """
    This takes in a dictionary, it associates a SNP (position, integer)
    to a tuple (sum, count) which gives the sum of all FSTs at this site
    and the count (for easy averaging and coverage estimation)
    """
  
    # Check if marker exist in dict or not
    # breakpoint()
    
    if not snp in dict_obj:
          dict_obj[snp] = (0,0)
          
    tup = dict_obj[snp]
    s = tup[0]
    c = tup[1]
    
    s += max(0, fst)
    c += 1
    
    dict_obj[snp] = (s,c)

        
Path = "~/src/Fstiminator/FstCalculations/AllBamFiles/"

if len(sys.argv) < 2:
    print("I need a pattern!", file=sys.stderr)
    exit(1)
    
pattern = sys.argv[1]
filelist = glob.glob(pattern)

dOfD = {} # this associates a MARKER to a dictionary (key to value)

counter=0

for i in filelist:
  if i.endswith(".fst.gz"):
    with gzip.open( i, 'rt') as f:
        first_line = f.readline()

        for line in f:
          line = line.strip().split(' ')
          # telling where to pull information for marker, snp, and fst from the bam files
          marker, snp, fst = line[0],line[1], line[4]

          snp = int(snp)
          fst = float(fst)
          # if the marker is not in the dictionary add the marker
          if not marker in dOfD:
            dOfD[ marker ] = {}
          # add countst to the times snp was seen and add fst estimate
          append_value(dOfD[ marker ], snp, fst)
          counter += 1
          if counter%10000==0:
            print(".", file=sys.stderr, end="")
            
# add column headers
print("Marker", "SNP", "SumFST", "Count", sep="\t")
# print the dictionary with marker, snp, and Fst mean estimate
for marker in dOfD:
  print_inner_dictionary(dOfD[marker] , marker)

