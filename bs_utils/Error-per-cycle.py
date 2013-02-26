import string, bz2, fileinput, shelve, time, subprocess, random, gzip
from fileinput import FileInput
from subprocess import Popen
import operator
from optparse import OptionParser
parser = OptionParser()

parser.set_defaults(filename="WA019L1")
parser.add_option("-l", "--lane", type="str",dest="filename",help="input lane name as WB002L1, WB002L2, etc", metavar="FILE")

parser.set_defaults(cyclenumber=100)
parser.add_option("-n", "--number", dest="cyclenumber",help="input read length", metavar="NUMB")


(options, args) = parser.parse_args()
                  
lane_name=options.filename
lane_name=str(lane_name)

length=options.cyclenumber
length=int(length)

#----------------------------------------------------------------
path="/NFS/Users/chen/BulkLabMy/Mouse/Mila-RRBS/mapping/"

file=lane_name+".mapping.gz"
inf=gzip.open(path+file, 'rb')

#file=lane_name+".mapping"
#inf=open(path+file, 'r')

cycles=[[0,0] for x in range(length)] #<-- [error, all]

mat=[[0,0,0,0,0] for x in range(length)] #<---- ATCGN for each cycle
atcg_map={"A":0,"T":1,"C":2,"G":3,"N":4}
#----------------------------------------------------------------

for line in inf:
	l=line.split()
	gseq=str(l[4][3:-3])
	read=str(l[5])
	read_length=min(len(read),length)
	for i in range(read_length):
		r=read[i]
		g=gseq[i]
		cycles[i][1]+=1
		if r!=g:
			if not (r=="T" and g=="C"):
				cycles[i][0]+=1
				mat[i][atcg_map[r]]+=1
inf.close()		

#----------------------------------------------------------------
	
outfile="error-per-cycle-"+lane_name+".txt"
outf=open(outfile,"w")
'''	
for k in range(2):
	for j in range(length):
		outf.write("%d\t"%(cycles[j][k]))
	outf.write("\n")
'''	

outf.write("%s\t"%(lane_name))

for j in range(length):
	error_rate=float(cycles[j][0])/sum(cycles[j])
	outf.write("%1.5f\t"%(error_rate))
outf.write("\n")
	
'''	
	
for k in range(5):
	for j in range(length):
		outf.write("%d\t"%(mat[j][k]))
	outf.write("\n")
	
for j in range(length):
	outf.write("%d\t"%(sum(mat[j])))
outf.write("\n")

for k in range(5):
	for j in range(length):
		if sum(mat[j]) !=0:
			outf.write("%1.6f\t"%(float(mat[j][k])/(sum(mat[j]))))
		else:
			outf.write("NA\t")
	outf.write("\n")
'''	
	
outf.close()

