import string, bz2, fileinput, shelve, time, subprocess, random
from fileinput import FileInput
from subprocess import Popen
import operator
from optparse import OptionParser
parser = OptionParser()

parser.set_defaults(lane_name="WA016L8")
parser.add_option("-l", "--lane", type="string", dest="lane_name",help="Input your lane (U006L1...)", metavar="FILE")
                  
(options, args) = parser.parse_args()
lane_name=options.lane_name
myfile_name=lane_name+".mapping"

#----------------------------------------------------------------
def uniq_pos(f1,f2):
	nn=0
	nn1=0
	nn2=0
	outff=open(f2,"w")
	pos10=""
	pos20=""
	for line in fileinput.input(f1):
		nn+=1
		l=line.split()
		pos=l[0]
		sign=pos[4]
		if sign=="+":
			if pos!=pos10:
				nn1+=1
				outff.write("%s"%(line))
				pos10=pos
		elif sign=="-":
			if pos!=pos20:
				nn2+=1
				outff.write("%s"%(line))
				pos20=pos
	fileinput.close()
	outff.close()
	print "-- sorted : %d"%(nn);
	print "-- uniq+sorted : %d"%(nn1);
	print "-- uniq-sorted : %d"%(nn2);
	print "-- uniq_sorted : %d"%(nn1+nn2);
		
#----------------------------------------------------------------
path="./"
#----------------------------------------------------------------
outfile1="reduced_"+myfile_name
outf1=open(path+outfile1,"w")

n=0
for line in fileinput.input(path+myfile_name):
	n+=1
	l=line.split()
	coordinate=str(l[3])
	outf1.write("%s\t%s\t%s"%(coordinate, l[5], str(l[-1]))+"\n")

fileinput.close()
outf1.close()
print "%s : %d"%(myfile_name,n);
#----------------------------------------------------------------
rid=str(random.randint(1000000,9999999))
folder1="./TMP_"+lane_name+"_"+rid

Popen('nohup mkdir %s '%(folder1),shell=True)

sorting1=Popen('nohup python sorting.py -b 5000000 -t %s -k "str(line[0:4]),int(line[5:16])" %s sorted_%s'%(folder1, outfile1, myfile_name),shell=True)

sorting1.wait()

uniq_pos("sorted_%s"%(myfile_name),"Uniq_sorted_%s"%(myfile_name))

Popen('nohup rm -r %s '%(folder1),shell=True)
Popen('nohup rm reduced_%s sorted_%s'%(myfile_name,myfile_name),shell=True)
#----------------------------------------------------------------
