import string, bz2, fileinput, shelve, time, subprocess 
from fileinput import FileInput
from subprocess import Popen
import operator
from optparse import OptionParser
parser = OptionParser()

#----------------------------------------------------------------
def rev_string(x):
	rev=x[::-1]
	return rev;


def context_calling(seq_lst,position):
	word=seq_lst[position]
	word=word.upper()
	context="--"
	context_CH="--"
	if word=="C":
		word2=seq_lst[position+1]
		context_CH=word+word2
		if word2=="G":
			context="CG"
		elif word2 in ['A','C','T']:
			word3=seq_lst[position+2]
			if word3=="G":
				context="CHG"
			elif word3 in ['A','C','T']:
				context="CHH"
	elif word=="G":
		word2=seq_lst[position-1]
		context_CH=word+word2
		context_CH=context_CH.translate(string.maketrans("ATCG", "TAGC"))
		if word2=="C":
			context="CG"
		elif word2 in ['A','G','T']:
			word3=seq_lst[position-2]
			if word3=="C":
				context="CHG"
			elif word3 in ['A','G','T']:
				context="CHH"
	return word,context,context_CH;

def myout(atcg,olist):
	mm="na"
	strmc="na"
	d=0
	if atcg=="C":
		mc=olist[2]
		uc=olist[1]
		d=mc+uc
		if d>0:
			mm="%1.3f"%(float(mc)/d)
			mm=str(mm)
			strmc=str(mc)
	elif atcg=="G":
		mc=olist[7]
		uc=olist[6]
		d=mc+uc
		if d>0:
			mm="%1.3f"%(float(mc)/d)
			mm=str(mm)
			strmc=str(mc)
	return d,strmc,mm;		

def add_up(s,dmer,lmer,xyz):
	if s=="+":
		for i in range(lmer):
			if xyz[i]=="A":
				dmer[lmer+i][0]+=1
			elif xyz[i]=="T":
				dmer[lmer+i][1]+=1
			elif xyz[i]=="C":
				dmer[lmer+i][2]+=1
			elif xyz[i]=="G":
				dmer[lmer+i][3]+=1
			elif xyz[i]=="N":
				dmer[lmer+i][4]+=1
	elif s=="-":
		for i in range(lmer):
			if xyz[i]=="a":
				dmer[lmer+i][5]+=1
			elif xyz[i]=="t":
				dmer[lmer+i][6]+=1
			elif xyz[i]=="c":
				dmer[lmer+i][7]+=1
			elif xyz[i]=="g":
				dmer[lmer+i][8]+=1
			elif xyz[i]=="n":
				dmer[lmer+i][9]+=1
	return dmer;
	
#----------------------------------------------------------------

start = time.time()
#genome_path='/Scratch/chen/Mouse/reference_genome/'
#genome_path='/Scratch/chen/Mouse/reads/RRBS_reference_genome_160_260/'
#genome_path="/u/home/mcdb/paoyang/genome/mouse/RRBS_reference_genome_100_300/"

parser.set_defaults(genome_path_info="/u/home/mcdb/paoyang/genome/mouse/RRBS_reference_genome_100_300/")
parser.add_option("-p", "--path", type="string", dest="genome_path_info",help="Path to RRBS refgeneme folder [/u/home/mcdb/paoyang/genome/mouse/RRBS_reference_genome_100_300/]", metavar="PATH")

genome_file='ref.shelve'

parser.add_option("-l", "--lane", type="string",dest="sampleid",help="input sample id (WB002L3)", metavar="FILE")

parser.add_option("-e", "--cycle", type="int",dest="read_length",help="input read length", metavar="LENG")

(options, args) = parser.parse_args()

genome_path=options.genome_path_info

sid=options.sampleid

lmer=options.read_length

#----------------------------------------------------------------
genome_seq=""
d = shelve.open(genome_path+genome_file,'r')
all_chrs=d.keys()
all_chrs.sort()
print all_chrs;
#----------------------------------------------------------------
filepath="./"
file="sorted_"+sid+".mapping"

f=open(filepath+file,"r")
line1=f.readline()
ll=line1.split()
lmer=lmer+10
dmer=[[0,0,0,0,0,0,0,0,0,0] for x in range(lmer*2)]  #<===== A,T,C,G,N,a,t,c,g,n
p0=str(ll[0][5:])

outfile = sid+".CGmap"

opath ="../CGmap/"
outf=open(opath+outfile,'w')		

chr0=""
for line in fileinput.input(filepath+file):
	l=line.split()
	SteveFilter=str(l[-1])
	coordinate=l[0]
	chr=str(coordinate[:4])
	#-----------------------------------------------
	if chr!=chr0:
		genome_seq=d[chr]
		genome_seq=genome_seq.upper()
		chr0=chr
	#-----------------------------------------------
	if SteveFilter=='0':
		xyz=l[1] # dna mode
		strand=coordinate[4]
		int_pos=int(coordinate[5:])-1
		if strand=="+":
			xyz=xyz.ljust(lmer,"@")
		elif strand=="-":
			xyz=rev_string(xyz)
			xyz=xyz.lower()
			xyz=xyz.ljust(lmer,"@")
		c=chr
		c0=str(c)
		p=str(int_pos)
		if p==p0:
			dmer=add_up(strand,dmer,lmer,xyz)
			#------- update p0 -------
			p0=int(p0)
			p0=str(p0)
		elif p!=p0:
			lag=int(p)-int(p0)
			if lag <= len(dmer):
				for i in range(lag):
					outlst=dmer[i]
					if sum(outlst)>0:
						ATCG,ctx,ctx2=context_calling(genome_seq,int(p0)-lmer+i)
						#pa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate=myout(ATCG,outlst)
						#outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,pa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate)+'\n')
						depth,mc,mrate=myout(ATCG,outlst)
						if ATCG in ["C","G"] and mrate !="na":
							outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%s\t%s\t%d'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,mrate,mc,depth)+'\n')
						
				for i in range(lag):
					del dmer[0]
					dmer.append([0,0,0,0,0,0,0,0,0,0])
				dmer=add_up(strand,dmer,lmer,xyz)
				#------- update p0 -------
				p0=int(p)
				p0=str(p0)
			elif lag > len(dmer):
				for i in range(len(dmer)):
					outlst=dmer[i]
					if sum(outlst)>0:
						ATCG,ctx,ctx2=context_calling(genome_seq,int(p0)-lmer+i)
						#pa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate=myout(ATCG,outlst)
						#outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,pa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate)+'\n')
						depth,mc,mrate=myout(ATCG,outlst)
						if ATCG in ["C","G"] and mrate !="na":
							outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%s\t%s\t%d'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,mrate,mc,depth)+'\n')
				dmer=[[0,0,0,0,0,0,0,0,0,0] for x in range(lmer*2)]
				dmer=add_up(strand,dmer,lmer,xyz)
				p0=int(p)
				p0=str(p0)
#-------- last mer --------
for i in range(len(dmer)):
	outlst=dmer[i]
	if sum(outlst)>0:
		ATCG,ctx,ctx2=context_calling(genome_seq,int(p0)-lmer+i)
		#qpa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate=myout(ATCG,outlst)
		#outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,pa,pt,pc,pg,pn,na,nt,nc,ng,nn,mrate)+'\n')
		depth,mc,mrate=myout(ATCG,outlst)
		if ATCG in ["C","G"] and mrate !="na":
			outf.write('%2s\t%1s\t%s\t%3s\t%2s\t%s\t%s\t%d'%(c0,ATCG,str(int(p0)-lmer+i),ctx,ctx2,mrate,mc,depth)+'\n')

#Popen('nohup gzip %s &'%(opath+outfile),shell=True)			
fileinput.close()
print "END CGMAP -- %s"%(sid)
elapsed = (time.time() - start)
print "# %s Time: %s"%(outfile,elapsed);
outf.close()			
d.close()
