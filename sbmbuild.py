#########################################
## Last revision on September 12, 2017 ##
#########################################

import numpy as np
import math
import linecache
import sys
import os


#######################
## GLOBAL VARIABLES ###
#######################

try:
	name = sys.argv[1]
	dca = sys.argv[2]
	xl = sys.argv[3]
except: 
	print("\n\tIncomplete arguments given!\n\tUsage: sbmbuild input_sequences DCA_contacts XL_pairs\n\tFiles for restriction pairs should be given in a two-column format\n\tSupport can be found at: rikchicfb@gmail.com\n")
	quit(0)


dcal = open(dca,"r")
dcas = len(dcal.readlines())
xll = open(xl,"r")
xls = len(xll.readlines())
seq = linecache.getline(name,1).split()[0]
ss = linecache.getline(name,2).split()[0]


#######################
## GLOBAL CONSTANTS ###
#######################

place=os.path.dirname(os.path.realpath(__file__))
codes = {'G':'GLY','P':'PRO','A': 'ALA','V':'VAL','L':'LEU','I':'ILE','M':'MET','C':'CYS','F':'PHE','Y':'TYR','W':'TRP','H':'HIS','K':'LYS','R':'ARG','Q':'GLN','N':'ASN','E':'GLU','D':'ASP','S':'SER','T':'THR'}
xlcodes= {'KK':2.21,'KS':1.82,'SK':1.82,'SS':1.43,'EE':1.69,'ED':1.56,'DE':1.56,'DD':1.43,'DK':1.14,'KD':1.14,'EK':1.27,'KE':1.27,'DS':0.75,'SD':0.75,'ES':0.88,'SE':0.88}
angle= {'H':'0.908000000E+02','E':'0.130500000E+03','C':'0.745280000E+02','-':'0.745280000E+02'}
dihedral1= {'H':'0.230900000E+03','E':'0.358700000E+03','C':'0.360000000E+03','-':'0.360000000E+03'}
dihedral2= {'H':'0.692750000E+03','E':'0.107600000E+04','C':'0.108000000E+04','-':'0.108000000E+04'}
factor = {'H':1,'E':1,'C':0.01,'-':0.01}
mg = 10
e = math.pow(0.4,12)
f=1.0 #this force factor can be modified as desired to balance the number of restrictions. Default value works fine when the number of couplings is comparable to sequence length

################################################################################
## HEADERS								                                      ##
################################################################################

def headers(head):
	if head=="begin":
		header=" ; Structure-based  topology file for Gromacs with Gaussian (SBM) support: rikchicfb@gmail.com - v-0.1\n [ defaults ]\n ;nbfunc comb-rule gen-pairs\n\t1\t1 no\n\n [ atomtypes ]\n ;name\tmass\tcharge \tptype\tc10\tc12\n CA\t1.000\t0.000\tA\t0.000\t0.167772160E-04\n\n [ moleculetype ]\n ;name   nrexcl\n Macromolecule \t3\n\n [ atoms ]\n ;nr  type  resnr residue atom  cgnr charge  mass\n"
	if head=="pair":
		header="\n [ pairs ]\n ; i j type and weights\n"
	if head=="excl":
		header="\n [ exclusions ]\n ;ai     aj\n"    
	if head=="bond":
		header="\n [ bonds ]\n ;ai     aj      func    r0(nm)  Kb\n"
	if head=="ang":
		header="\n [ angles ]\n ;ai  aj   ak  func  th0(deg)   Ka\n"
	if head=="dih":
		header="\n [ dihedrals ]\n ;ai  aj  ak  al  func  phi0(deg)  Kd  mult\n"
	if head=="end":
		header="\n [ system ]\n ;name\n Macromolecule\n\n [ molecules ]\n ;name   #molec\n Macromolecule   1\n"
	return header

################################################################################
## GENERATE GRO FILE FOR SBM SIMULATIONS				                      ##
################################################################################

def printgro(r,seq):
	outputgro = open(name.split('.')[0]+"_calpha.gro","w")
	outputgro.write("Structure-Based gro file for SBM simulation - rikchicfb@gmail.com\n"+str(len(seq)).rjust(5)+"\n")
	for i in range(len(r[:,0])):
		outputgro.write(str(i+1).rjust(5)+codes[seq[i]].ljust(6)+"CA".rjust(4)+str(i+1).rjust(5))
		outputgro.write(str("{0:.3f}".format(r[i,0])).rjust(8)+str("{0:.3f}".format(round(r[i,1],3))).rjust(8)+str("{0:.3f}".format(round(r[i,2],3))).rjust(8)+"\n")			
	outputgro.write(str("{0:.5f}".format(max(r[:,2])+mg)).rjust(10)+str("{0:.5f}".format(max(r[:,2])+mg)).rjust(10)+str("{0:.5f}".format(max(r[:,2])+mg)).rjust(10)+"\n")
	outputgro.close()
	print "\n\tA protein coordinates file was saved in this folder as "+name.split('.')[0]+"_calpha.gro\n"
	return 0

################################################################################
## GENERATE DCA GAUSSIAN-LIKE POTENTIALS				                      ##
################################################################################

def distances(seq,dca,f):
	table = open(place+"/distavg","r")
	l = len(table.readlines())
	## get distance for each pair combination
	dic = {}
	for i in range(1,l+1):
		line = linecache.getline(place+"/distavg",i).split()
		dic[line[1]+line[2]]=line[0]
		dic[line[2]+line[1]]=line[0]
	## get distances to DCA pairs
	if (len(dca)!=2):
		print "\n\tWrong format for DCA list! Topology file could not be generated!\n\tA two-column format is mandatory! Check your DCA file!\n\n"
		quit(0)
	try:
		d = float(dic[seq[int(dca[0])-1]+seq[int(dca[1])-1]])*0.1
	except: 
			print("\n\tInconsistent DCA pair given at residue pair: "+dca[0]+" and "+dca[1]+"! Topology file could not be generated!\n\tCheck your list of coevolution pairs!\n\n")
			quit(0)
	sigma = 0.5
	return dca[0].rjust(6)+dca[1].rjust(6)+" 6 ".rjust(3)+str("{0:.3f}".format(f)).ljust(7)+str("{0:.8f}".format(d)).rjust(10)+str("{0:.8f}".format(sigma)).rjust(13)+str(e).rjust(15)+"\n"

#########################################################################################################################################
## GENERATE TOPOLOGY FILE WITH DIHEDRALS ASSIGNED FROM THE SS PREDICTION AND COEVOLUTION DCA CONTACTS AS GAUSSIAN ATRACTION POTENTIALS ##
#########################################################################################################################################

def topology(r,seq,ss):
	outputtop = open(name.split('.')[0]+"_calpha.top","w")
	outputtop.write(headers("begin"))
	for i in range(len(r[:,0])):
		outputtop.write(str(i+1).rjust(6)+"CA".rjust(4)+str(i+1).rjust(8)+codes[seq[i]].rjust(5)+"CA".rjust(4)+str(i+1).rjust(8)+"{0:.3f}".format(0).rjust(8)+"{0:.3f}".format(1).rjust(8)+"\n")
	outputtop.write(headers("pair"))
	for i in range(1,dcas+1):
		outputtop.write(distances(seq,linecache.getline(dca,i).split(),f))
	for i in range(1,xls+1):
		px=linecache.getline(xl,i).split()
		if (len(px)!=2):
			print "\n\tWrong format for cross-linking list! Topology file could not be generated!\n\tA two-column format is mandatory! Check your cross-linking  file!\n\n"
			quit(0)
		px1=px[0]
		px2=px[1]
		try:
			d=xlcodes[seq[int(px1)-1]+seq[int(px2)-1]]
		except: 
			print("\n\tInconsistent XL pair given at residue pair: "+px1+" and "+px2+"! Topology file could not be generated!\n\tCheck your list of chemical cross-links!\n\tList of accepted residue pairs:\n\t\t for DSS linker: KK, KS, SS\n\t\t for Hexadiamine linker: EE, ED, DD\n\t\t for zero-length: DK, EK, DS, ES\n")
			quit(0)
		sigmaxl=round((d*0.05),8)
		outputtop.write(px1.rjust(6)+px2.rjust(6)+" 6 ".rjust(3)+str("{0:.3f}".format(f)).ljust(7)+str("{0:.8f}".format(round((d/11)*1,8))).rjust(10)+str("{0:.8f}".format(sigmaxl)).rjust(13)+str(e).rjust(15)+"\n")
		for j in range(2,11):
			outputtop.write(px1.rjust(6)+px2.rjust(6)+" 5 ".rjust(3)+str("{0:.3f}".format(f)).ljust(7)+str("{0:.8f}".format(round((d/11)*j,8))).rjust(10)+str("{0:.8f}".format(sigmaxl)).rjust(13)+"\n")
	outputtop.write(headers("bond"))
	for i in range(len(r[:,0])-1):
		outputtop.write(str(i+1).rjust(6)+str(i+2).rjust(7)+"1".rjust(2)+"0.380200000E+00".rjust(18)+"0.200000000E+05".rjust(18)+"\n")
	outputtop.write(headers("excl"))
	for i in range(1,dcas+1):
		outputtop.write(linecache.getline(dca,i).split()[0].rjust(6)+linecache.getline(dca,i).split()[1].rjust(6)+"\n")
	for i in range(1,xls+1):
		outputtop.write(linecache.getline(xl,i).split()[0].rjust(6)+linecache.getline(xl,i).split()[1].rjust(6)+"\n")
	outputtop.write(headers("ang"))
	for i in range(len(r[:,0])-2):
		try:
			outputtop.write(str(i+1).rjust(6)+str(i+2).rjust(7)+str(i+3).rjust(7)+"1".rjust(2)+angle[ss[i+1]].rjust(18)+str("{0:.3f}".format(0.4*factor[ss[i+1]])).rjust(8)+"000000E+02\n")
		except:
			print "\n\tSomething is wrong with your secondary structure prediction! Please check your file!\n\tList of accepted tags:\n\t\t alpha-helix: H\n\t\t beta-strand: E\n\t\t coiled-coil: C or -\n\n"
			quit(0)
	outputtop.write(headers("dih"))
	for i in range(len(r[:,0])-3):
		outputtop.write(str(i+1).rjust(6)+str(i+2).rjust(7)+str(i+3).rjust(7)+str(i+4).rjust(7)+"1".rjust(2))
		outputtop.write("{:.9E}".format((float(dihedral1[ss[i+1]])+float(dihedral1[ss[i+2]]))/2).rjust(18)+str("{0:.3f}".format(0.1*factor[ss[i+2]])).rjust(8)+"000000E+01 1\n")
		outputtop.write(str(i+1).rjust(6)+str(i+2).rjust(7)+str(i+3).rjust(7)+str(i+4).rjust(7)+"1".rjust(2))
		outputtop.write("{:.9E}".format((float(dihedral2[ss[i+1]])+float(dihedral2[ss[i+2]]))/2).rjust(18)+str("{0:.3f}".format(0.5*factor[ss[i+2]])).rjust(8)+"000000E+00 3\n")
	outputtop.write(headers("end"))
	outputtop.close()
	print "\tA protein topology file was saved in this folder as "+name.split('.')[0]+"_calpha.top\n"
	return 0

################################################################################
## MAIN CODE								      							  ##
################################################################################

distca = 0.3802
ang = math.asin(0.60547)

print "\n\n\tProtein sequence in selected file:"
for j in range(0,len(seq),80):
	print "\t"+seq[j:j+80]

r = np.zeros((len(seq),3))

print "\n\tTotal number of residues: "+str(len(r[:,0]))

print "\n\n\tnumber of DCA restrictions: "+str(dcas)
print "\n\tnumber of XL restrictions: "+str(xls)+"\n"

d = len(seq)*0.115

for i in range(0,len(r[:,0])):
	if(i==0):
		r[i,:] = [d,d,mg]
	else:	
		x=r[i-1,0]
		y=r[i-1,1]+math.cos(ang)*distca*math.pow(-1,i)
		z=r[i-1,2]+math.sin(ang)*distca
		r[i,:] = [x,y,z]

printgro(r,seq)
topology(r,seq,ss)
print "\n\tYou are ready to run a fold simulation.\n\n"


