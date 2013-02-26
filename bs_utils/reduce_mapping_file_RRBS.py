import string, bz2, fileinput, shelve, time, subprocess, random
from fileinput import FileInput
from subprocess import Popen
import operator
from optparse import OptionParser
parser = OptionParser()

parser.set_defaults(filename="WB002L1.mapping")
parser.add_option("-f", "--file", type="string", dest="filename",help="inputfile as FILE", metavar="FILE")

parser.set_defaults(pathname="./")
parser.add_option("-p", "--path", type="string", dest="pathname",help="input path", metavar="PATH")

(options, args) = parser.parse_args()
                  
myfile_name=options.filename
path=options.pathname

#----------------------------------------------------------------
outfile="reduced-"+myfile_name
outf=open(path+outfile,"w")
for line in fileinput.input(path+myfile_name):
	l=line.split()
	outf.write("%s\t%s\t%s"%(str(l[3]),l[5],str(l[-1]))+"\n")
fileinput.close()
outf.close()

#----------------------------------------------------------------
rid=str(random.randint(1000000,9999999))
folder1="TMP_"+myfile_name+"_"+rid

Popen('nohup mkdir %s '%(path+folder1),shell=True)

sorting1=Popen('nohup python sorting.py -b 500000 -t ./%s -k "str(line[0:4]),int(line[5:16])" %s sorted_%s'%(path+folder1, path+outfile, myfile_name),shell=True)

sorting1.wait()

Popen('nohup rm -r %s '%(path+folder1),shell=True)
Popen('nohup rm %s &'%(path+outfile),shell=True)
Popen('nohup gzip %s &'%(path+myfile_name),shell=True)
			