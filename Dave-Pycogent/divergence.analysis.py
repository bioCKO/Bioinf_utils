#!/usr/bin/python
#David Clarke - 160843
#
#This script takes a list of "Bombyx peptides","Helicoverpa contigs","Heliothis scaffolds" seperated by tabs
#exonerates the files against each other and generates a genbank file and divergence estimates.
#
#


import math
import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

from warnings import filterwarnings

filterwarnings("ignore", "Not using MPI as mpi4py not found")

from cogent import LoadSeqs, LoadTree
from cogent.evolve.models import CNFGTR
from cogent.maths.stats import chisqprob

GBKOUT = 1 #Output Genbank Files


NOFRAMESHIFT = 1 #Removes exons with HAHV frameshifts from paml analysis
NOBMHAFS = 1 #Removes exons with BMHA frameshifts from analysis
NOOVERLAP = 1 #Selects only one overlapping exon

DNDSANALYSIS = 1 #Enables DN/DS analysis

PAMLOUT = 1 #Output PAML sequence Files
PYCOGENTOUT = 1 #Output PYCOGENT analysis
MYDNDSOUT = 1 #Output MYDNDS analysis

USECTLFILE = 1#Use control file rather than defults

#location of sequence database files to be used
HeViESTfile = "/home/dclarke/HeVi19.05.11nc_assembly/HeVi19.05.11nc_d_results/HeVi19.05.11nc_out.unpadded.fasta"  
HeAmGenomicfile = "/var/www/html/blast/db/Helicoverpa.contigs.fa-04-29-2010"
BoMoPepfile = "/var/www/html/blast/db/silkpep"
projectprefex = ""

if PAMLOUT == 1:
	os.system("mkdir paml")
os.system("mkdir data")	
os.system("mkdir DNDSres")
os.system("mkdir temp")
os.system("mkdir seq")
os.system("mkdir genbank")


#This idf statment imports the settings from the control file
if USECTLFILE == 1:
	try:
		ctlfile = open('../DNDSan.ctl','r')
		ctllines = ctlfile.readlines()
		ctlfile.close()
		if ctllines[0][:-1]=="#control file for DNDSanalysis":
	
			try:
				if ctllines[3].split("'")[0] == "EST location = ":
					HeViESTfile = ctllines[3].split("'")[1]
					print "EST location read"		
				else:
					print "EST location not read"
			except:	
				print "ERROR reading EST location"

			try:
				if ctllines[4].split("'")[0] == "Genomic location = ":
					HeAmGenomicfile = ctllines[4].split("'")[1]
					print "Genomic location read"
				else:
					print "Genomic location not read"
			except:	
				print "ERROR reading Genomic location"
		
			try:
				if ctllines[5].split("'")[0] == "pep location = ":
					BoMoPepfile = ctllines[5].split("'")[1]
					print "pep location read"
				else:
					print "pep location not read"
			except:	
				print "ERROR reading pep location"
		
			try:
				if ctllines[6].split("'")[0] == "projectprefex = ":
					projectprefex = ctllines[6].split("'")[1]
					print "project prefex read"
				else:
					print "project prefex not read"
			except:	
				print "ERROR reading pep location"
			
			if PAMLOUT == 1:
				try:
					if ctllines[9].split("/")[0] == "-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----" and ctllines[67].split("/")[0] == "-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----":
						pamlctl = open('paml/codemlHAseq.ctl','w')
						for x in ctllines[10:66]:
							pamlctl.write(x)
						pamlctl.close()
						print "paml control file written"
					else:
						print "paml control file not read"
						print ctllines[9]
						print ctllines[67]
				except:
					print "ERROR reading PAML crl file"
				
		else:
			print "control file line 1 does not match"
	except:
		print "ERROR reading control file"
		
oldcontig = 'x'
oldest = 'x'
oldpep = 'x'
oldChr = 'x' 
HAGenes = ""
HVGenes = ""
HAHVGene = ""
numgenes = 0
totHAGenes = ""
totHVGenes = ""
totnumgenes = 0
totHAautoGenes = ""
totHVautoGenes = ""
totnumautogenes = 0
totHAno0Genes = ""
totHVno0Genes = ""
totnumno0genes = 0

tab1 = 0
tab2 = 0
list = ""
list2 = ""
mode = 0
num = 0
HAcodons = ""
HVcodons = ""
HAnofs = ""
HVnofs = ""
HAnofseq = ""
HVnofseq = ""
seqnum = 0
codonlens = []
listcodonlens = ""
TotalPep = []

codonstop = 0
codonstart = 0
ecodonstop = 0
ecodonstart = 0
startaln = 0
stopaln = 0
HAstart = 0
HAstop = 0
y = 0
D = 1
seqlen = 0

f = open('temp/makefile.gbk','w')
He = open('data/make.data.txt','w')
ex = open('data/exon.data.txt','w')
Pa = open('data/paml.data.txt','w')
HAS = open('data/HAseq.fasta','w')
if PAMLOUT == 1:
	PAMLres = open('DNDSres/PAMLres.txt','w')
PYCOGENTres = open('DNDSres/PYCOGENTres.txt','w')
MYDNDSres = open('DNDSres/MYDNDSres.txt','w')
MYDNDSlongres = open('DNDSres/MYDNDSlongres.txt','w')
exonerateres = open('data/exoner.data.txt','w')
sequences = open('data/sequences.data.txt','w')
nofseq = open('data/nofseq.data.txt','w')
exonmix = open('data/exonmix.txt','w')
exonorder = open('data/exonorder.txt','w')
framecheck = open('data/framecheck.txt','w')

#Method for aligning genomic and EST data
def mRNA(newest,newcontig,newpep,pepstart,pepstop,TotalPep):
	aligntype = 0
	fastacmd = ("fastacmd -d " + HeViESTfile + " -s " + newest + " -o temp/Heliothis_contig.out")
	os.system(fastacmd)

	
	os.system("exonerate --ryo '>%ti(%qi)\n%tcs' -n 1 --showalignment FALSE --showvulgar TRUE --model coding2genome temp/Heliothis_contig.out temp/Helicoverpa_contig.out >temp/exonerate.txt")
	os.system("exonerate --ryo '>%ti(%qi)\n%tcs' -n 1 --showalignment TRUE --showvulgar TRUE --model coding2genome temp/Heliothis_contig.out temp/Helicoverpa_contig.out >temp/exonerateextra.txt")	
	e = open('temp/exonerate.txt','r')
	elist = e.readlines()
	e.close()
	
	if elist[3][:8]=="draw_hsp":
		He.write(newcontig + "(" + newest + ") HA.HV exonerate of bad quality\n")

				
	elen = len(elist)
	if elen <= 4:
		He.write(newcontig + "(" + newest + ") HA.HV exonerate empty\n")
	else:
		addexonerate = open('temp/exonerateextra.txt','r')
		readexonerate = addexonerate.readlines()
		addexonerate.close()
		for x in readexonerate:
			exonerateres.write(x)
		
	endseq = elen-1
	HAseq = elist[3:endseq]
	vulgar = HAseq[0]
		
	V = vulgar.split(' ')
	
	HValnst = V[2]
	HValnsp = V[3]
	startaln = V[6]
	stopaln = V[7]
	HVcon = V[1][4:]
	HVdir = 1
	if int(HValnst) > int(HValnsp):
		HVdir = -1	
	DifSt = (int(HValnst) - int(startaln))
	
	
	Pa.write(newest + " - " + newcontig + ";DifSt: " + str(DifSt) + " HValignstart: " + HValnst + ", HValignstop: " + HValnsp + "\n")

		
	try:
		if DNDSANALYSIS == 1:
			readHA = open("temp/Helicoverpa_contig.out", "rU")
			HAlines = readHA.readlines()
			readHA.close()
			HAstring = ""
			HAtitle = HAlines[0].rstrip()
			for x in HAlines[1:]:
				HAstring = HAstring + x.rstrip()
			HAseq = Seq(HAstring)
		
	
		
			readHV = open("temp/Heliothis_contig.out", "rU")
			HVlines = readHV.readlines()
			readHV.close()
			HVstring = ""
			HVtitle = HVlines[0].rstrip()
			for x in HVlines[1:]:
				HVstring = HVstring + x.rstrip()
			HVseq = Seq(HVstring)
	except:
		print "read sequences error " + HVcon
			
	if GBKOUT == 1:
		f.write('     mRNA            ' + startaln + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
	
	HAstart = 0
	HAstop = 0
	HAcode = ""
	HVcode = ""
	HAcodons = ""
	HVcodons = ""
	HAnofs = ""
	HVnofs = ""
	HAnofseq = ""
	HVnofseq = ""
	
	D = 1
	Ex = 1
	if int(startaln) > int(stopaln):
		D = -1
		Ex = 0
	revrev = 1
	if D & HVdir == -1:
		revrev = -1
		
	exons = ""
	cds = ""
	codonstop = int(startaln)
	codonstart = int(startaln) + 1*Ex
	ecodonstop = int(startaln)
	ecodonstart = int(startaln) + 1*Ex
	pcodonstop = int(startaln)
	pcodonstart = int(startaln) + 1*Ex
	G = 0
	F = 0
	E = 0
	H = 0
	y = 0
	frameshift = 0
	PST = int(pepstart)
	PSP = int(pepstop)
	startfix = 1
	exonY = {}
	HAstart = 0
	HAstop = 0
	frame = 0
	for x in V:			
		if x in ("C","M","F","G"):
			codonstop = codonstop + (int(V[y+2]))*D
			ecodonstop = ecodonstop + (int(V[y+2]))*D
			pcodonstop = pcodonstop + (int(V[y+2]))*D
			if x == "F":
				F = F + (int(V[y+2]))
				E = E + (int(V[y+1]))
				frameshift = 1
			if x == "G":
				G = G + (int(V[y+2]))
				if (int(V[y+2]))%3 != 0:
					frameshift = 1
					print"Gap frameshift"
				H = H + (int(V[y+1]))
				if (int(V[y+1]))%3 != 0:
					frameshift = 1
					print"Gap frameshift"

			HAHVcode = SEQMAKE(HAcode,HAseq,HVcode,HVseq,PST,PSP,pcodonstop,pcodonstart,x,(int(V[y+2])),(int(V[y+1])),DifSt,HVdir,startfix,TotalPep,aligntype)
				
			HAcode = HAHVcode[0]
			HVcode = HAHVcode[1]
			try:
				newHAnofs = HAHVcode[2].tostring()
			except:
				newHAnofs = HAHVcode[2]
				
			try:
				newHVnofs = HAHVcode[3].tostring()
			except:
				newHVnofs = HAHVcode[3]

			if aligntype == 1:
				print newest + " - " + newcontig + "\nHelicoverpa  \n" + newHAnofs + "\nHeliothis  \n" + newHVnofs
			
			
			if (HAHVcode[4] and HAHVcode[5]) > 0:
				if HAstart == 0:
					HAstart = HAHVcode[4] + 1*Ex
				HAstop = HAHVcode[5] + (1-Ex)
			
				
			HAnofs = HAnofs + newHAnofs
			HVnofs = HVnofs + newHVnofs
					
			DifSt = DifSt + ((int(V[y+1]))-(int(V[y+2])*HVdir*D))*HVdir
			pcodonstart = pcodonstop
			startfix = startfix - startfix		
		if x == 'S':
			codonstop = codonstop + (int(V[y+2]))*D
			ecodonstop = ecodonstop + (int(V[y+2]))*D
			pcodonstop = pcodonstop + (int(V[y+2]))*D
			
			DifSt = DifSt + ((int(V[y+1]))-(int(V[y+2])*HVdir*D))*HVdir
			pcodonstart = pcodonstop
		if x == 'I':
			ex.write(newcontig + '\t' + HVcon + '\t' + str(codonstart) + '\t' + str(codonstop) + '\tHV-'+ newpep + '\t' + V[y+2] + '\t' + str(F) + '\t' + str(G) + '\t' + str(E) + '\t' + str(H) + '\n')

			exons = exons + str(ecodonstart) + '..' + str(ecodonstop) + ','
			cds = cds + str(ecodonstart) + '..' + str(ecodonstop) + ','
			
			DifSt = DifSt + ((int(V[y+1]))-(int(V[y+2]))-4)*D	
			codonstart = codonstop + (int(V[y+2]) + 5)*D
			codonstop = codonstop + (int(V[y+2]) + 4)*D
			ecodonstart = ecodonstop + (int(V[y+2]) + 5)*D
			ecodonstop = ecodonstop + (int(V[y+2]) + 4)*D
			pcodonstart = pcodonstop + (int(V[y+2]) + 4)*D
			pcodonstop = pcodonstop + (int(V[y+2]) + 4)*D
			G = 0
			F = 0
			E = 0
			H = 0
			
			if (frameshift == 1) and (len(HAnofs) != 0):
				if "-" in HAnofs or HVnofs:
					frameshift = 1
				else:
					frameshift = 0
					print "exon re-entered"
					print newest + " - " + newcontig + "\nHelicoverpa  \n" + HAnofs + "\nHeliothis  \n" + HVnofs

			
			if frameshift == 0:
				if len(HAnofs)%3==0:
					HAnofseq = HAnofseq + HAnofs
					HVnofseq = HVnofseq + HVnofs

					nexonY = len(exonY)+1
					
					doesstop = CODONSTOP(HAnofs,HVnofs)
					stopcodon = int(doesstop[2])

					if HAstart < HAstop:
						frame = HAstart%3 + 1
					elif HAstop <HAstart:
						frame = (HAstart%3 + 1)*-1
				
					newexon = (HAstart,HAstop,HAnofs,HVnofs,frame)
					if ((HAstart + HAstop) != 0) and (stopcodon == 1):
						exonY[nexonY] = newexon

			HAstart = 0
			HAstop = 0		
			frame = 0
			HAnofs = ""
			HVnofs = ""
			frameshift = 0
			
		if y == (len(V)-1):
			if (frameshift == 1) and (len(HAnofs) != 0):
				if "-" in HAnofs or HVnofs:
					frameshift = 1
				else:
					frameshift = 0
					print "exon re-entered"
					print newest + " - " + newcontig + "\nHelicoverpa  \n" + HAnofs + "\nHeliothis  \n" + HVnofs
			if frameshift == 0:
				if len(HAnofs)%3==0:
					HAnofseq = HAnofseq + HAnofs
					HVnofseq = HVnofseq + HVnofs

					nexonY = len(exonY)+1
					
					doesstop = CODONSTOP(HAnofs,HVnofs)
					stopcodon = int(doesstop[2])

					if HAstart < HAstop:
						frame = HAstart%3 + 1
					elif HAstop <HAstart:
						frame = (HAstart%3 + 1)*-1
				
					newexon = (HAstart,HAstop,HAnofs,HVnofs,frame)
					if ((HAstart + HAstop) != 0) and (stopcodon == 1):
						exonY[nexonY] = newexon
						
			HAnofs = ""
			HVnofs = ""
			frameshift = 0
			try:
				nofseq.write(newest + " - " + newcontig + "\nHAnofseq\n" + HAnofseq + "\nHVnofseq\n" + HVnofseq + "\n")
			except:
				nofseq.write("write HAnofseq error" + newest + " - " + newcontig + "\n")
				
			ex.write(newcontig + '\t' + HVcon + '\t' + str(codonstart) + '\t' + str(stopaln) + '\tHV-' + newpep + '\t' + "0" + '\t' + str(F) + '\t' + str(G) + '\t' + str(E) + '\t' + str(H) + '\n')
			if GBKOUT == 1:
				f.write('     exon             ' + exons + str(ecodonstart) + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
				if int(startaln) > int(stopaln):
					f.write('     CDS             complement(' + cds + str(codonstart) + '..' + stopaln + ')\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
				else:
					f.write('     CDS             ' + cds + str(codonstart) + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
			
			if NOFRAMESHIFT == 1:
				HAcode = HAnofseq
				HVcode = HVnofseq
			
			
			try:				
				sequences.write(newest + " - " + newcontig + "\nHelicoverpa  \n" + HAcode + "\nHeliothis  \n" + HVcode + '\n')

				if (len(HAcode) == len(HVcode)):
					codelen = len(HAcode)
					if codelen == 0:
						if NOOVERLAP == 0:
							if PYCOGENTOUT == 1:
								PYCOGENTres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "FALSE" + '\t' + 'noseq'  + '\t' + "FALSE\n")
							if PAMLOUT == 1:
								PAMLres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "FALSE" + '\t' + 'noseq'  + '\t' + "FALSE\n")
					
					try:
						doesstop = CODONSTOP(HAcode,HVcode)
						HAcodons = str(doesstop[0])
						HVcodons = str(doesstop[1])
						stopcodon = int(doesstop[2])
						z = int(doesstop[3])	
					except:
						print "error looking for stop codons  " + HVcon

					
					if PAMLOUT == 1 and NOOVERLAP == 0:
						RUNPAML(newest,newcontig,newpep,HAcode,HVcode,NOOVERLAP,oldChr)
					
					if PYCOGENTOUT == 1 and NOOVERLAP == 0:
						RUNPYCOGENT(newest,newcontig,newpep,HAcode,HVcode,NOOVERLAP,oldChr)
					
					if MYDNDSOUT == 1 and NOOVERLAP == 0:
						MYDNDS(newest,newcontig,newpep,HAcode,HVcode,NOOVERLAP,oldChr)
					
					if stopcodon == 0:

						if NOOVERLAP == 0:
							PYCOGENTres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "FALSE" + '\t' + 'stop/s' + '\t' + str(codelen) + '\n') 
						try:
							Pa.write("sequences contain stop codons - " + newest + " - " + newcontig + "\n")						
						except:
							Pa.write("sequences contain stop codons - " + newest + " - " + newcontig + "\n")
							
				else:
					try:
						PAMLres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "FALSE" + '\t' + 'len!=' + '\t' + str(len(HAcode)) + '\n') 
						if NOOVERLAP == 0:
							PYCOGENTres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "FALSE" + '\t' + 'len!=' + '\t' + str(len(HAcode)) + '\n')
						Pa.write("sequences are not of equal length - " + newest + " - " + newcontig + "\n")
						return 0									
					except:
						print "sequence lengeth error error " + HVcon
			except:
				print "write sequence error " + HVcon
			G = 0
			F = 0
			E = 0
			H = 0
		y += 1
	return z,HAcodons,HVcodons,stopcodon,exonY

#Method for alignment of genomic sequence and peptides	
def protein(newpep,newcontig,start="0",stop="0"):
	fastacmd = ("fastacmd -d " + BoMoPepfile + " -s " + newpep + " -o temp/Bombyx_pep.out -L " + start + "," + stop)
	os.system(fastacmd)
	fastacmd = ("fastacmd -d " + HeAmGenomicfile + " -s " + newcontig + " -o temp/Helicoverpa_contig.out")
	os.system(fastacmd)
	
	os.system("exonerate --score 500 --ryo '>%ti(%qi)\n%tcs' -n 1 --showalignment FALSE --showvulgar TRUE --model protein2genome temp/Bombyx_pep.out temp/Helicoverpa_contig.out >temp/exonerate.txt")
	os.system("exonerate --score 500 --ryo '>%ti(%qi)\n%tcs' -n 1 --showalignment TRUE --showvulgar TRUE --model protein2genome temp/Bombyx_pep.out temp/Helicoverpa_contig.out >temp/exonerateextra.txt")
	
	e = open('temp/exonerate.txt','r')
	elist = e.readlines()
	elen = len(elist)
	if elen <= 4:
		He.write(newcontig + "(" + newpep + ")" + start + ',' + stop + ' ' + "exonerate empty\n")
	else:
		addexonerate = open('temp/exonerateextra.txt','r')
		readexonerate = addexonerate.readlines()
		addexonerate.close()
		for x in readexonerate:
			exonerateres.write(x)
	
	endseq = elen-1
	HAseq = elist[3:endseq]
	vulgar = HAseq[0]
	e.close()
			
	V = vulgar.split(' ')
	
	fullalignstart = ""
	fullalignend = ""
	startaln = V[6]
	stopaln = V[7]
	HVcon = V[1][4:]
	ALNproS = V[2]
	ALNproE = V[3]
	exonP = {}
	frameshift = 0
	
	try:
		if ALNproS > 1:
			if ALNproS != "0":
				fullalign = protein(newpep,newcontig,"0",ALNproS)				
				fullalignstart = fullalign[0]
				
				newexonp = fullalign[4]
				for x in newexonp:
					exonP[x] = newexonp[x]
	except:
		He.write(newcontig + "(" + newpep + ") " + "start" + ALNproS + "\n")
		
	try:
		if stop == start == "0":
			fullalign = protein(newpep,newcontig,ALNproE,"0")
			fullalignend = fullalign[1]

			newexonp = fullalign[4]
			for x in newexonp:
				exonP[x] = newexonp[x]
			
	except:
		He.write(newcontig + "(" + newpep + ") " + "end" + ALNproE + "\n")
	
	if fullalignstart == "":
		fullalignstart = startaln
	if fullalignend == "":
		fullalignend = stopaln	
			
	He.write(newcontig + '\t' + HVcon + '\t' + startaln + '\t' + stopaln + '\n')
			
	if GBKOUT == 1:
		f.write('     CDS_AFTER            ' + startaln + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
	D = 1
	Ex = 1
	if int(startaln) > int(stopaln):
		D = -1
		Ex = 0
	exons = ""
	cds = ""
	codonstop = int(startaln)
	codonstart = int(startaln) + 1*Ex
	ecodonstop = int(startaln)
	ecodonstart = int(startaln) + 1*Ex
	y = 0
	G = 0
	F = 0
	E = 0
	H = 0

	BMstart = 0
	BMstop = 0
	frame = 0
	for x in V:	
		if x in ("C","M","F","G"):
			if BMstart == 0:
				BMstart = ecodonstop + 1*Ex
			codonstop = codonstop + (int(V[y+2]))*D
			ecodonstop = ecodonstop + (int(V[y+2]))*D
			if x == 'F':
				F = F + (int(V[y+2]))
				E = E + (int(V[y+1]))
				if NOBMHAFS == 1:
					frameshift = 1
			if x == 'G':
				G = G + (int(V[y+2]))
				H = H + (int(V[y+1]))

			BMstop = ecodonstop	+ (1-Ex)
				

		if x == 'S':
			codonstop = codonstop + (int(V[y+2]))*D
			ecodonstop = ecodonstop + (int(V[y+2]))*D
		if x == 'I':
			ex.write(newcontig + '\t' + HVcon + '\t' + str(codonstart) + '\t' + str(codonstop) + '\tBM\t' + V[y+2] +  '\t' + str(F) + '\t' + str(G) + '\t' + str(E) + '\t' + str(H) + '\n')
			exons = exons + str(ecodonstart) + '..' + str(ecodonstop) + ','

			cds = cds + str(ecodonstart) + '..' + str(ecodonstop) + ','

			codonstart = codonstop + (int(V[y+2]) + 5)*D
			codonstop = codonstop + (int(V[y+2]) + 4)*D
			ecodonstart = ecodonstop + (int(V[y+2]) + 5)*D
			ecodonstop = ecodonstop + (int(V[y+2]) + 4)*D
			G = 0
			F = 0	
			E = 0
			H = 0
			
			if frameshift == 0:
				if BMstart < BMstop:
					frame = BMstart%3 + 1
				elif BMstop <BMstart:
					frame = (BMstart%3 + 1)*-1
				nexonP = len(exonP)+1
				newpexon = (BMstart,BMstop,frame)
				exonP[HVcon + "$" + str(nexonP)] = newpexon
			
			BMstart = 0
			BMstop = 0
			frame = 0
			frameshift = 0
			
		if y == (len(V)-1):
			ex.write(newcontig + '\t' + HVcon + '\t' + str(codonstart) + '\t' + str(stopaln) + '\tBM\t' + "0" +  '\t' + str(F) + '\t' + str(G) + '\t' + str(E) + '\t' + str(H) + '\n')
			if GBKOUT == 1:
				f.write('     exon             ' + exons + str(ecodonstart) + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
				if int(startaln) > int(stopaln):
					f.write('     BLASTCDS             complement(' + cds + str(codonstart) + '..' + stopaln + ')\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
				else:
					f.write('     BLASTCDS             ' + cds + str(codonstart) + '..' + stopaln + '\n                     /shared_id="' + HVcon + '"\n                     /note="' + HVcon + '"\n')
			
			
			G = 0
			F = 0
			E = 0
			H = 0
			if frameshift == 0:
				if BMstart < BMstop:
					frame = BMstart%3 + 1
				elif BMstop <BMstart:
					frame = (BMstart%3 + 1)*-1
				nexonP = len(exonP)+1
				newpexon = (BMstart,BMstop,frame)
			
				exonP[HVcon + "$" + str(nexonP)] = newpexon
		y += 1
	
	return startaln,stopaln,fullalignstart,fullalignend,exonP


#method for converting exonerate alignmnet into aligned exons	
def SEQMAKE(HAcode,HAseq,HVcode,HVseq,PST,PSP,pcodonstop,pcodonstart,x,A,B,DifSt,HVdir,startfix,TotalPep,aligntype):
	try:
		CSP = int(pcodonstop)
		CST = int(pcodonstart)
		HAstart = 0
		HAstop = 0
		HAnofs = ""
		HVnofs = ""
		M = 1
		ESTframe = 0
		
		pepdir = 0
		if TotalPep[0][3] > 0:
			pepdir = 1
		elif TotalPep[0][3] < 0:
			pepdir = -1
		
		
		ESTdir = 0
		if CSP > CST:
			ESTdir = 1
			ESTframe = (CST + 1)%3 + 1
		elif CSP < CST:
			ESTdir = -1
			ESTframe = ((CST+1)%3 + 1)*-1
			
		
		conflct = ESTdir * pepdir
		if conflct < 0:
			print "ESTdir and pepdir conflict"
			return
		
		if PST>PSP:		#BM Protein is aligned to the HA sequence in reverse
			if HVdir == 1:	#HV sequence is aligned to the HA sequence in same direction
				Pa.write("making BM rev, HV same " + str(HVdir) + "\n")
				M = -1
				Diff = 0
				if (CSP<PST and CST>PSP):
					if PST<CST:
						Diff = CST - PST
						HAstart = PST
						
					else: HAstart = CST
					if PSP<CSP:
						HAstop = CSP
					else: HAstop = PSP
					if x in ("F","G"):
						ReA = 0
						if ((A + B)%3)>0:
							ReA = 3-((A + B)%3) 
						if A == 0:
							HAcode = HAcode + (B+ReA) * "-"
							HVcode = HVcode + HVseq[(HAstart+DifSt+Diff*2):((HAstart+(HAstart-HAstop))+DifSt)+B+Diff*2].tostring() + (ReA * "-")				
							HAnofs = (B+ReA) * "-"
							HVnofs = HVseq[(HAstart+DifSt+Diff*2):((HAstart+(HAstart-HAstop))+DifSt)+B+Diff*2].tostring() + (ReA * "-")
						else:
							HVcode = HVcode + (A+ReA) * "-"
							HAcode = HAcode + HAseq[HAstop:HAstart].reverse_complement() + ReA * "-"
							HVnofs = (A+ReA) * "-"
							HAnofs = HAseq[HAstop:HAstart].reverse_complement() + ReA * "-"
					elif x in ("M","C"):
						HAcode = HAcode + HAseq[HAstop:HAstart].reverse_complement()
						HVcode = HVcode + HVseq[(HAstart+DifSt+Diff*2):((HAstart+(HAstart-HAstop))+DifSt+Diff*2)].tostring()
						HAnofs = HAseq[HAstop:HAstart].reverse_complement()
						HVnofs = HVseq[(HAstart+DifSt+Diff*2):((HAstart+(HAstart-HAstop))+DifSt+Diff*2)].tostring()
					Pa.write("DifSt: " + str(DifSt) + "\n")
					Pa.write("HVstart " + str((HAstart+DifSt)+Diff*2) + " HVstop " + str(((HAstart+(HAstart-HAstop))+DifSt)+Diff*2) + "\n")
					Pa.write("HAstart " + str(HAstart) + " HAstop " + str(HAstop) + "\n")

			elif HVdir == -1:	#HV sequence is aligned to the HA sequence in reverse
				Pa.write("making BM rev, HV rev " + str(HVdir) + "\n")
				M = 1
				if (CSP<PST and CST>PSP):
					if PST<CST:
						HAstart = PST
					else: HAstart = CST
					if PSP<CSP:
						HAstop = CSP
					else: HAstop = PSP
					if x in ("F","G"):
						ReA = 0
						if ((A + B)%3)>0:
							ReA = 3-((A + B)%3) 
						if A == 0:
							HAcode = HAcode + (B+ReA) * "-"
							HVcode = HVcode + HVseq[(((HAstop)+DifSt)*M):(HAstart+DifSt)*M+B].reverse_complement() + (ReA * "-")				
							HAnofs = (B+ReA) * "-"
							HVnofs = HVseq[(((HAstop)+DifSt)*M):(HAstart+DifSt)*M+B].reverse_complement() + (ReA * "-")
						else:
							HVcode = HVcode + (A+ReA) * "-"
							HAcode = HAcode + HAseq[HAstop:HAstart].reverse_complement() + ReA * "-"
							HVnofs = (A+ReA) * "-"
							HAnofs = HAseq[HAstop:HAstart].reverse_complement() + ReA * "-"
					elif x in ("M","C"):
						HAcode = HAcode + HAseq[HAstop:HAstart].reverse_complement()
						HVcode = HVcode + HVseq[(((HAstop)+DifSt)*M):(HAstart+DifSt)*M].reverse_complement()
						HAnofs = HAseq[HAstop:HAstart].reverse_complement()
						HVnofs = HVseq[(((HAstop)+DifSt)*M):(HAstart+DifSt)*M].reverse_complement()
					Pa.write("DifSt: " + str(DifSt) + "\n")
					Pa.write("HVstart " + str(((HAstop)+DifSt)*M) + " HVstop " + str((HAstart+DifSt)*M) + "\n")					
					Pa.write("HAstart " + str(HAstart) + " HAstop " + str(HAstop) + "\n")

		
		else:
			if HVdir == 1:	#HV sequence is aligned to the HA sequence in same direction
				M = 1
				if (CSP>PST and CST<PSP):	#BM Protein is aligned to the HA sequence in same direction
					Pa.write("making BM same, HV same " + str(HVdir) + "\n")
					if PST>CST:
						HAstart = PST
					else: HAstart = CST-startfix
					if PSP>CSP:
						HAstop = CSP
					else: HAstop = PSP
					if x in ("F","G"):
						ReA = 0
						if ((A + B)%3)>0:
							ReA = 3-((A + B)%3) 
						if A == 0:
							HAcode = HAcode + (B+ReA) * "-"
							HVcode = HVcode + HVseq[(HAstart+DifSt)*M:(HAstop+DifSt)*M+B].tostring() + (ReA * "-")	
							HAnofs = (B+ReA) * "-"
							HVnofs = HVseq[(HAstart+DifSt)*M:(HAstop+DifSt)*M+B].tostring() + (ReA * "-")	
						else:
							HVcode = HVcode + (A+ReA) * "-"
							HAcode = HAcode + HAseq[HAstart:HAstop].tostring() + ReA * "-"
							HVnofs = (A+ReA) * "-"
							HAnofs = HAseq[HAstart:HAstop].tostring() + ReA * "-"
					elif x in ("M","C"):
						HAcode = HAcode + HAseq[HAstart:HAstop].tostring()
						HVcode = HVcode + HVseq[(HAstart+DifSt)*M:(HAstop+DifSt)*M].tostring()
						HAnofs = HAseq[HAstart:HAstop].tostring()
						HVnofs = HVseq[(HAstart+DifSt)*M:(HAstop+DifSt)*M].tostring()
					Pa.write("DifSt: " + str(DifSt) + "\n")
					Pa.write("HVstart " + str((HAstart+DifSt)*M) + " HVstop " + str((HAstop+DifSt)*M) + "\n")
					Pa.write("HAstart " + str(HAstart) + " HAstop " + str(HAstop) + "\n")


			elif HVdir == -1:	#HV sequence is aligned to the HA sequence in reverse
#				M = -1
				Diff = 0
				if (CSP>PST and CST<PSP):	#BM Protein is aligned to the HA sequence in same direction
					Pa.write("making BM same, HV rev " + str(HVdir) + "\n")
					if PST>CST:
						Diff = CST - PST - startfix
						HAstart = PST
					else: HAstart = CST-startfix
					if PSP>CSP:
						HAstop = CSP
					else: HAstop = PSP
					if x in ("F","G"):
						ReA = 0
						if ((A + B)%3)>0:
							ReA = 3-((A + B)%3) 
						if A == 0:
							HAcode = HAcode + (B+ReA) * "-"
							HVcode = HVcode + HVseq[((HAstart-(HAstop-HAstart))+DifSt)*M+2*Diff:(HAstart+DifSt)*M+B+2*Diff].reverse_complement() + (ReA * "-")					
							HAnofs = (B+ReA) * "-"
							HVnofs = HVseq[((HAstart-(HAstop-HAstart))+DifSt)*M+2*Diff:(HAstart+DifSt)*M+B+2*Diff].reverse_complement() + (ReA * "-")
						else:
							HVcode = HVcode + (A+ReA) * "-"
							HAcode = HAcode + HAseq[HAstart:HAstop].tostring() + ReA * "-"
							HVnofs = (A+ReA) * "-"
							HAnofs = HAseq[HAstart:HAstop].tostring() + ReA * "-"

					elif x in ("M","C"):
						HAcode = HAcode + HAseq[HAstart:HAstop].tostring()
						HVcode = HVcode + HVseq[((HAstart-(HAstop-HAstart))+DifSt)*M+2*Diff:(HAstart+DifSt)*M+2*Diff].reverse_complement()
						HAnofs = HAseq[HAstart:HAstop].tostring()
						HVnofs = HVseq[((HAstart-(HAstop-HAstart))+DifSt)*M+2*Diff:(HAstart+DifSt)*M+2*Diff].reverse_complement()
					Pa.write("DifSt: " + str(DifSt) + "\n")
					Pa.write("HVstart " + str(((HAstart-(HAstop-HAstart))+DifSt)*M+2*Diff) + " HVstop " + str((HAstart+DifSt)*M+2*Diff) + "\n")
					Pa.write("HAstart " + str(HAstart) + " HAstop " + str(HAstop) + "\n")

	except:
		Pa.write("make sequence error ")
	return HAcode,HVcode,HAnofs,HVnofs,HAstart,HAstop

#method for adding details to genbank file
def addHAend(writecontig,newChr):
	gbk2 = open("temp/Helicoverpa_contig.out",'r')
	gbk2list = gbk2.readlines()
	gbk3list = gbk2list[1:]

	wholeseq = ""
	for x in gbk3list:
		wholeseq = wholeseq + x.rstrip()
	sA = 0
	sT = 0
	sG = 0
	sC = 0
	for x in wholeseq:
		if x == 'A':
			sA += 1
		if x == 'T':
			sT += 1
		if x == 'G':
			sG += 1
		if x == 'C':
			sC += 1
				
	f.write('BASE COUNT    ' + str(sA) + ' a   ' + str(sC) + ' c   ' + str(sG) + ' g  ' + str(sT) + ' t\nORIGIN\n')
	wholeseq = wholeseq.lower()
	lenseq = len(wholeseq)
	B = ' '
	place = 1
	while place <= len(wholeseq):
		strplace = str((place - 10))
		LP = len(strplace)
		if ((place - 10) % 60) == 0:				
			f.write('\n'+ (B*(9-LP))+str(place-9)+B+wholeseq[(place-10):(place)])
		elif ((place - 10) % 10) == 0:
			f.write(B+wholeseq[(place-10):(place)])
		elif (place) == (len(wholeseq)):
			if (place % 60) <= 10:
				lastseq = place % 10
				f.write('\n'+ (B*(9-LP))+str(place-9)+B+wholeseq[-lastseq:])
			else:
				lastseq = place % 10
				f.write(B+wholeseq[-lastseq:])
		place += 1
	f.write('\n//')
			
	f.close()
	
	
	newfile = open("genbank/" + "Chr" + str(newChr) + "/" + writecontig + ".gbk" , 'w')
	readfile = open('temp/makefile.gbk','r')
	filelines = readfile.readlines()
	for x in filelines:
		newfile.write(x)
	newfile.close()
	readfile.close()
	return

#method for running pycogent divergence analysis	
def RUNPYCOGENT(newest,newcontig,newpep,HAcodons,HVcodons,NOOVERLAP,chr):
	try:
		codelen = len(HAcodons)
		if len(HAcodons) == len(HVcodons):
		
		
			dna = {'Helicoverpa':str(HAcodons),'Heliothis':str(HVcodons)}
			aln = LoadSeqs(data=dna)
	
			assert len(aln.Names) <= 3, "Need to modify the script to handle case of more than 3 taxa"

			tree = LoadTree(tip_names=aln.Names[:])

			sm = CNFGTR()
			lf = sm.makeLikelihoodFunction(tree)
			lf.setAlignment(aln)
	
			# following required for pairs
			lf.setParamRule('length', is_independent=True)
	
			# limit_action='raise' will cause an error if the functiuon cannot be optimise
			# you can set that to 'warn', but I suggest raise
			opt_args = dict(local=True, max_restarts=5, max_evaluations=100000, limit_action='raise') 


			lf.setParamRule('omega', is_constant=False)
			lf.optimise(**opt_args)
			alt_lnL = lf.getLogLikelihood()
			alt_nfp = lf.getNumFreeParams()
	
			# get the results tables
			#stats = lf.getStatistics(with_titles=True)
			#for stats_table in stats:
			#    print stats_table

			# demo just getting out omega
			omega = lf.getParamValue('omega')

			omega_lo, omega_mle, omega_hi = lf.getParamInterval('omega')
			
			if omega_lo == None:
				omega_lo = 0
			if omega_mle == None:
				omega_mle = 0
			if omega_hi == None:
				omega_hi = 0		
			
			if NOOVERLAP == 0:
				PYCOGENTres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "%.4f\t%.4f\t%.4f" % (omega_lo, omega_mle, omega_hi) + '\t' + str(codelen) + '\t' + str(chr) + '\n')
			elif NOOVERLAP == 1:
				PYCOGENTres.write(newcontig + '\t' + newpep + '\t' + "%.4f\t%.4f\t%.4f" % (omega_lo, omega_mle, omega_hi) + '\t' + str(codelen) + '\t' + str(chr) + '\n')
		else:
			if NOOVERLAP == 0:
				PYCOGENTres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "FALSE" + '\t' + "len!=" + '\t' + str(len(HAcodons) - len(HVcodons)) + '\t' + str(codelen) + '\t' + str(chr) + '\n')
			elif NOOVERLAP == 1:
				PYCOGENTres.write(newcontig + '\t' + newpep + '\t' + "FALSE" + '\t' + "len!=" + '\t' + str(len(HAcodons) - len(HVcodons)) + '\t' + str(codelen) + '\t' + str(chr) + '\n')
	except:
		Pa.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "Error in RUNPYCOGENT" + '\n\n')
		print newest + '\t' + newcontig + '\t' + newpep + '\t' + "Error in RUNPYCOGENT"
		PYCOGENTres.write(newcontig + '\t' + newpep + '\t' + "FALSE" + '\t' + "PYerr" + '\t'  + "FALSE" + '\t' + str(len(HAcodons)) + '\t' + str(chr) + '\n')

#method for running PAML divergence analysis
def RUNPAML(newest,newcontig,newpep,HAcode,HVcode,NOOVERLAP,chr):
		try:
			doesstop = CODONSTOP(HAcode,HVcode)
			HAcodons = str(doesstop[0])
			HVcodons = str(doesstop[1])
			stopcodon = int(doesstop[2])
			z = int(doesstop[3])	
		except:
			print "error looking for stop codons  " + newest
						
		try:
			PAMLseq = open('paml/PAMLseq.fasta','w')
			PAMLseq.write("\t2 " + str(z) + "\nHelicoverpa  \n" + HAcodons + "\nHeliothis  \n" + HVcodons)
			PAMLseq.close()	
		except:
			print "error writing PAML sequence " + newest
		
		try:
			PAMLcmd = "cd paml ; codeml codemlHAseq.ctl"
			os.system(PAMLcmd)
			
			rstfile = open('paml/rst','r')
			rstlines = rstfile.readlines()
			rstdata = rstlines[5].split()
			rstdS = rstdata[5]
			rstdN = rstdata[4]
			rstdNdS = rstdata[6]
			print "dS: " + rstdS
			print "dN: " + rstdN	
			print "dN/dS: " + rstdNdS			
			
			if NOOVERLAP == 0:
				PAMLres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + rstdS  + '\t' + rstdN + '\t' + rstdNdS + '\t' + str(z) + '\t' + str(chr) + '\n') 
			elif NOOVERLAP == 1:
				PAMLres.write(newcontig + '\t' + newpep + '\t' + rstdS  + '\t' + rstdN + '\t' + rstdNdS + '\t' + str(z) + '\t' + str(chr) + '\n') 

			
		except:
			Pa.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "Error in RUNPAML" + '\n\n')
			print newest + '\t' + newcontig + '\t' + newpep + '\t' + "Error in RUNPAML"

#method for consolidating multipe exonetate genomic-EST alignments			
def CONPEPSEQ(pepexons,newcontig,newpep,dir,TotalPep):
	try:
		NOOVERLAP = 1
		oldseqstart = -1
		oldseqstop = -1
		oldHAexon = 'x'
		oldHVexon = 'x'
		HAgene = ""
		HVgene = ""
		Exonorder = []
		exonnums = ""
		exonmix.write(str(pepexons))
		
		INFRAME = 1
		
		if dir == 1:
			for x in sorted(pepexons):
				seqstart = pepexons[x][0]
				seqstop = pepexons[x][1]
				HAexon = pepexons[x][2]
				HVexon = pepexons[x][3]
				frame = pepexons[x][4]

				estinfo = x.split('-')
				est = estinfo[1]
				
				try:				
					for PE in TotalPep:
						if (PE[1] >= seqstart) and (PE[2] <= seqstop):
							framecheck.write("prot exon is within EST exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) )
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1
							
						elif (PE[1] <= seqstop) and (PE[2] >= seqstop):
							framecheck.write("prot exon start is within EST exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) + "\n")
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1
							
						elif (PE[1] <= seqstart) and (PE[2] >= seqstart):
							framecheck.write("EST exon start is within prot exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) + "\n")
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1						
				except:
					He.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + " error in forward frame finder")
				
				if INFRAME == 1:
					if (oldseqstart == -1) or ((seqstart > oldseqstart) and (seqstart > oldseqstop)):
						HAgene = HAgene + HAexon
						HVgene = HVgene + HVexon
						Exonorder.append((est,seqstart,seqstop,frame))
						exonnums = exonnums + str(seqstart) + '..' + str(seqstop) + ','
					
					elif (seqstart < oldseqstop) and (seqstop > oldseqstop):
						exonorder.write(newcontig + '\t' + newpep + '\t' + 'exons overlap' + '\n')
						framecheck.write(est + ' exons overlap' + "\n")
					
					oldseqstart = seqstart
					oldseqstop = seqstop
					oldHAexon = HAexon
					oldHVexon = HVexon
								
				
		if dir == -1:
			for x in sorted(pepexons, reverse=True):
				seqstart = pepexons[x][0]
				seqstop = pepexons[x][1]
				HAexon = pepexons[x][2]
				HVexon = pepexons[x][3]
				frame = pepexons[x][4]

				estinfo = x.split('-')
				est = estinfo[1]
				
				try:	
					for PE in TotalPep:
						if (PE[1] <= seqstart) and (PE[2] >= seqstop):
							framecheck.write("prot exon is within EST exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) + "\n")
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1
						
						elif (PE[2] <= seqstop) and (PE[1] >= seqstop):
							framecheck.write("prot exon start is within EST exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) + "\n")
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1
						
						elif (PE[2] <= seqstart) and (PE[1] >= seqstart):
							framecheck.write("EST exon start is within prot exon" + "\n")
							framecheck.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + str(PE) + "\n")
							if frame != PE[3]:
								framecheck.write( "Out of FRAME" + "\n" + "\n")
								INFRAME = 0
							else:
								framecheck.write( "IN FRAME" + "\n" + "\n")
								INFRAME = 1
					
				except:
					He.write(est + ", " + str(seqstart) + ", " + str(seqstop) + ", " + str(frame) + " error in reverse frame finder")
				
				if INFRAME == 1:
					if (oldseqstart == -1) or ((seqstart < oldseqstart) and (seqstart < oldseqstop)):
						HAgene = HAgene + HAexon
						HVgene = HVgene + HVexon
						Exonorder.append((est,seqstart,seqstop,frame))
						exonnums = exonnums + str(seqstart) + '..' + str(seqstop) + ','
					
					elif (seqstart > oldseqstop) and (seqstop < oldseqstop):
						exonorder.write(newcontig + '\t' + newpep + '\t' + 'exons overlap' + '\n')
						framecheck.write(est + ' exons overlap' + "\n")
					
					oldseqstart = seqstart
					oldseqstop = seqstop
					oldHAexon = HAexon
					oldHVexon = HVexon
				
		exonmix.write(newcontig + '-' + newpep + '\n>Helicoverpa\n' + HAgene + '\n>Heliothis\n' + HVgene + '\n')# + str(Exonorder) + '\n'	
		exonorder.write(newcontig + '-' + newpep + 	'\n')
		n = 0
		for y in Exonorder:
			n += 1
			exonorder.write(str(n) + '\t' +y[0] + '\t' + str(y[1]) + '\t' + str(y[2]) + '\t' + str(y[3]) + '\n')
			
		if GBKOUT == 1:
			if dir == 1:
				f.write('     CDS              ' + exonnums[:-1] + '\n                     /shared_id="' + newcontig + '-' + newpep + '"\n                     /note="' + newcontig + '-' + newpep + '"\n')
			if dir == -1:
				f.write('     CDS             complement(' + exonnums[:-1] + '\n                     /shared_id="' + newcontig + '-' + newpep + '"\n                     /note="' + newcontig + '-' + newpep + '"\n')
		
		if PYCOGENTOUT == 1:
			if len(HAgene) > 0:
				RUNPYCOGENT(str(newcontig + '-' + newpep),newcontig,newpep,HAgene,HVgene,NOOVERLAP,oldChr)
			else:
				PYCOGENTres.write(newcontig + '\t' + newpep + '\t' + "FALSE" + '\t' + "noseq" + '\t'  + "FALSE" + '\t' + str(0) + '\t' + str(oldChr) + '\n')

		if PAMLOUT == 1:
			if len(HAgene) > 0:
				RUNPAML(str(newcontig + '-' + newpep),newcontig,newpep,HAgene,HVgene,NOOVERLAP,oldChr)
			else:	
				PAMLres.write(newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "noseq" + '\t' + "FALSE" + '\t' + str(0) + '\t' + str(oldChr) + '\n') 

		if MYDNDSOUT == 1:
			if len(HAgene) > 0:
				MYDNDS(str(newcontig + '-' + newpep),newcontig,newpep,HAgene,HVgene,NOOVERLAP,oldChr)
			else:	
				MYDNDSres.write(newcontig + '\t' + newpep + '\t' + "FALSE"  + '\t' + "noseq" + '\t' + "FALSE" + '\t' + "FALSE" + '\t' + str(0) + '\t' + str(oldChr) + '\n') 

	except:
		He.write(newcontig + '\t' + newpep + '\t' + "conpepseq creation error" + '\n')
	
	return HAgene,HVgene

#Method for checking for stop codon in seqence alignments
def CODONSTOP(HAcode,HVcode):
	try:	
		z = 0
		stopcodon = 1
		current = ""	
		HAcodons = ""
		HVcodons = ""
		e = len(HAcode) - (len(HAcode)%3)
		for x in HAcode:
			if z < e:
				HAcodons = HAcodons + x
				z += 1
				current = current + x								
				if z%3==0:
					HAcodons = HAcodons + " "
					if current in ("TAA","TAG","TGA"):
						stopcodon = 0
					current = ""									
								
		z = 0
		current = ""						
		for x in HVcode:
			if z < e:
				HVcodons = HVcodons + x
				z += 1
				current = current + x
				if z%3==0:
					HVcodons = HVcodons + " "
					if current in ("TAA","TAG","TGA"):
						stopcodon = 0
					current = ""	
	
	except:
		print'error in codonstop'	
		
	return HAcodons,HVcodons,stopcodon,z

#Method for running basic divergence analysis
def MYDNDS(newest,newcontig,newpep,HAcode,HVcode,NOOVERLAP,chr):
	try:	
		z = 0
		stopcodon = 1
		HAcurrent = ""	
		HVcurrent = ""	
		HAS =0.000
		HVS =0.000
		HAN	=0.000
		HVN =0.000
		
		HACNF = 0
		HVCNF = 0
		NonCon = 0
		CNF = 0	
	
		SC = 0.000
		NC = 0.000
		Changes = 0.000
		
		gaps = 0
		e = len(HAcode)
		for x in range(1,(len(HAcode)+1)):
				HAcurrent = HAcurrent + HAcode[x-1]	
				HVcurrent = HVcurrent + HVcode[x-1]
				if x%3==0:
					if ("-") in (HAcurrent[0],HAcurrent[1],HAcurrent[2],HVcurrent[0],HVcurrent[1],HVcurrent[2]):
						gaps += 3
					elif ("N") in (HAcurrent[0],HAcurrent[1],HAcurrent[2],HVcurrent[0],HVcurrent[1],HVcurrent[2]):
						NonCon += 3
					elif (HAcurrent[0]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					elif (HAcurrent[1]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					elif (HAcurrent[2]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					elif (HVcurrent[0]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					elif (HVcurrent[1]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					elif (HVcurrent[2]) in ("B","D","E","F","H","I","J","K","L","M","O","P","Q","R","S","U","V","W","X","Y","Z"):
						CNF += 3
					
					else:
						HACDNDS = CODONDNDS(HAcurrent)
						if HACDNDS[0] == 'NULL':
							stopcodon = 0
						else:
							HAS += HACDNDS[0]
							HAN += HACDNDS[1]
							HACNF += HACDNDS[2]
					
						HVCDNDS = CODONDNDS(HVcurrent)
						if HVCDNDS[0] == 'NULL':
							stopcodon = 0
						else:
							HVS += HVCDNDS[0]
							HVN += HVCDNDS[1]
							HVCNF += HVCDNDS[2]
							
						NCSC = COMPCODNS(HAcurrent,HVcurrent)
						NC += NCSC[0]
						SC += NCSC[1]
						Changes += NCSC[2]

					HAcurrent = ""		
					HVcurrent = ""	
		
		Ts = (HAS+HVS)/2.000
		Tn = (HAN+HVN)/2.000
		Ps = SC/Ts
		Pn = NC/Tn
		
		DS = -3/4*math.log((1-(4*Ps/3)))
		DN = -3/4*math.log((1-(4*Pn/3)))

		
		if NOOVERLAP == 0:
			if stopcodon == 0:
				MYDNDSlongres.write(newest + '\t' + newcontig + '\t' + newpep + '\tFALSE\tstop\tFALSE\tFALSE\t' + "%.0f\t%.0f\t%.0f\t%.0f" % (gaps,NonCon,CNF,HVCNF+HACNF) + '\t' + "%.2f\t%.2f\t%.2f" % (NC,SC,Changes) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
				MYDNDSres.write(newest + '\t' + newcontig + '\t' + newpep + '\tFALSE\tstop\tFALSE\tFALSE\t' +  str(len(HAcode)) + '\t' + str(chr) + '\n') 

			else:
				MYDNDSlongres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "%.2f\t%.2f\t%.2f\t%.2f" % (HAS,HVS,HAN,HVN)  + '\t' + "%.0f\t%.0f\t%.0f\t%.0f" % (gaps,NonCon,CNF,HVCNF+HACNF) + '\t' + "%.2f\t%.2f\t%.2f" % (NC,SC,Changes) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
				MYDNDSres.write(newest + '\t' + newcontig + '\t' + newpep + '\t' + "%.4f\t%.4f\t%.4f\t%.4f" % (Ts,Tn,Ps,Pn) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
		elif NOOVERLAP == 1:
			if stopcodon == 0:
				MYDNDSlongres.write(newcontig + '\t' + newpep + '\tFALSE\tstop\tFALSE\t' + "%.0f\t%.0f\t%.0f\t%.0f" % (gaps,NonCon,CNF,HVCNF+HACNF) + '\t' + "%.2f\t%.2f\t%.2f" % (NC,SC,Changes) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
				MYDNDSlongres.write(newcontig + '\t' + newpep + '\tFALSE\tstop\tFALSE\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
			else:
				MYDNDSlongres.write(newcontig + '\t' + newpep + '\t' + "%.2f\t%.2f\t%.2f\t%.2f" % (HAS,HVS,HAN,HVN) + '\t' + "%.0f\t%.0f\t%.0f\t%.0f" % (gaps,NonCon,CNF,HVCNF+HACNF) + '\t' + "%.2f\t%.2f\t%.2f" % (NC,SC,Changes) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 
				MYDNDSres.write(newcontig + '\t' + newpep + '\t' + "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (Ts,Tn,Ps,Pn,DS,DN) + '\t' + str(len(HAcode)) + '\t' + str(chr) + '\n') 


	except:
		print'error in MYDNDS'	

#Method	for comparing codons for NS and S changes	
def COMPCODNS(HAcurrent,HVcurrent):
	try:
		C = [0,0,0]
		N = 0.000
		S = 0.000
		changes = 0.000
	
		for x in (0,1,2):
				if HAcurrent[x] != HVcurrent[x]:
					C[x] += 1
					changes += 1.000
		
		if sum(C) == 0:
			return N,S,changes
		
		else:
			HAAA = CODONAA(HAcurrent)
			HVAA = CODONAA(HVcurrent)
			
			if sum(C) == 1:
				if HAAA == HVAA:
					S = 1.000
				else:	
					N = 1.000
				return N,S,changes
		
			elif sum(C) in (2,3):
				try:
					N2 = 0.000
					S2 = 0.000
					for x in (0,1,2):
						if C[x] == 1:
							if x == 0:
								AN2S2 = COMPCODNS(HAcurrent,HAcurrent[0]+HVcurrent[1]+HVcurrent[2])
								N2 += AN2S2[0]/changes
								S2 += AN2S2[1]/changes
								BN2S2 = COMPCODNS(HAcurrent[0]+HVcurrent[1]+HVcurrent[2],HVcurrent)
								N2 += BN2S2[0]/changes
								S2 += BN2S2[1]/changes
							if x == 1:
								AN2S2 = COMPCODNS(HAcurrent,HVcurrent[0]+HAcurrent[1]+HVcurrent[2])
								N2 += AN2S2[0]/changes
								S2 += AN2S2[1]/changes
								BN2S2 = COMPCODNS(HVcurrent[0]+HAcurrent[1]+HVcurrent[2],HVcurrent)
								N2 += BN2S2[0]/changes
								S2 += BN2S2[1]/changes
							if x == 2:
								AN2S2 = COMPCODNS(HAcurrent,HVcurrent[0]+HVcurrent[1]+HAcurrent[2])
								N2 += AN2S2[0]/changes
								S2 += AN2S2[1]/changes
								BN2S2 = COMPCODNS(HVcurrent[0]+HVcurrent[1]+HAcurrent[2],HVcurrent)
								N2 += BN2S2[0]/changes
								S2 += BN2S2[1]/changes
				except:
					print "Error in double change"
				return N2,S2,changes

	
	except:
		print'error in COMPCODNS'	

#Method for determining coding amino acid of codon
def CODONAA(codon):
	try:
		amino = 'NULL'
		if codon in ("TTT","TTC"):
			amino = 'Phe'
		elif codon in ("TTA","TTG","CTT","CTC","CTA","CTG"):
			amino = 'Leu'
		elif codon in ("ATT","ATC","ATA"):	
			amino = 'Ile'
		elif codon in ("ATG"):
			amino = 'Met'
		elif codon in ("GTT","GTC","GTA","GTG"):
			amino = 'Val'		
		elif codon in ("TCT","TCC","TCA","TCG","AGT","AGC"):
			amino = 'Ser'		
		elif codon in ("CCT","CCC","CCA","CCG"):
			amino = 'Pro'			
		elif codon in ("ACT","ACC","ACA","ACG"):
			amino = 'Thr'	
		elif codon in ("GCT","GCC","GCA","GCG"):
			amino = 'Ala'	
		elif codon in ("TAT","TAC"):
			amino = 'Tyr'				
		elif codon in ("CAT","CAC"):
			amino = 'His'	
		elif codon in ("CAA","CAG"):
			amino = 'Gln'			
		elif codon in ("AAT","AAC"):
			amino = 'Asn'			
		elif codon in ("AAA","AAG"):
			amino = 'Lys'			
		elif codon in ("GAT","GAC"):
			amino = 'Asp'	
		elif codon in ("GAA","GAG"):
			amino = 'Glu'			
		elif codon in ("TGT","TGC"):
			amino = 'Cys'		
		elif codon in ("TGG"):
			amino = 'Trp'	
		elif codon in ("CGT","CGC","CGA","CGG","AGA","AGG"):
			amino = 'Arg'		
		elif codon in ("GGT","GGC","GGA","GGG"):
			amino = 'Gly'	
		elif codon in ("TAA","TAG","TGA"):
			amino = 'Stop'	
			print "Stop codon found"
		else:
			amino = 'NULL'
			print "NULL codon found"
	except:
		print'error in CODONAA'	

	return amino

#Method for determing NS of S sites for codon	
def CODONDNDS(codon):
	S = 0.000
	N = 0.000

	CNF = 0
	if (codon) in ("TTA","TTG","AGG"):
		S = (2.000/3.000)
		N = (7.000/3.000)

	elif (codon) in ("TTT","TTC","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","AGT","AGC"):#"TAT","TAC","TGT","TGC",
		S = (1.000/3.000)
		N = (8.000/3.000)	

	elif (codon) in ("CTT","CTC","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","CGT","CGC","GGT","GGC","GGA","GGG"):
		S = (1.000)
		N = (2.000)	

	elif (codon) in ("CTA","CTG","CGG"):
		S = (4.000/3.000)
		N = (5.000/3.000)

	elif (codon) in ("ATT","ATC","ATA"):
		S = (2.000/3.000)
		N = (7.000/3.000)

	elif (codon) in ("ATG","TGG"):
		S = 0.000
		N = (3.000)

	elif (codon) in ("TAA","TAG","TGA"):
		S = 'NULL'
		N = 'NULL'

	elif (codon) in ("TAT","TAC"):
		S = 1.000
		N = 2.000
	elif (codon) in ("TGT","TGC"):
		S = 1.000/2.000
		N = 5.000/2.000
	elif (codon) in ("CGA"):
		S = 3.000/2.000
		N = 3.000/2.000
	elif (codon) in ("AGA"):
		S = 5.000/6.000
		N = 13.000/6.000
	else:
		print "Codon not found " + codon
		CNF = 1

	return S,N,CNF

#Method for consolidating multiple exonerate peptide alignmnets	
def CONSOLPRO(exonP):
	Exonorder = []
	oldseqstart = -1
	oldseqstop = -1
	oldframe = 0
	exonnums = ""

	for x in exonP:
		
		num0s = 6-len(str(exonP[x][0]))
		BMexons[num0s*str(0) + str(str(exonP[x][0]) + '$' + str(x))] = exonP[x]
		if (exonP[x][0] < exonP[x][1]):
			BMdir = 1
		elif (exonP[x][0] > exonP[x][1]):
			BMdir = -1
		else:
			exonorder.write(newcontig + '\t' + newpep + '\t' + 'BMexons in different directions' + '\n')
	
	
	if BMdir == 1:
		for x in sorted(BMexons):
	
			pepinfo = x.split('$')
			pep = pepinfo[1]
			
			if (oldseqstart == -1) or ((BMexons[x][0] > oldseqstart) and (BMexons[x][0] > oldseqstop)):

				Exonorder.append((pep,BMexons[x][0],BMexons[x][1],BMexons[x][2]))
				exonnums = exonnums + str(BMexons[x][0]) + '..' + str(BMexons[x][1]) + ','
					
			elif (BMexons[x][0] < oldseqstop) and (BMexons[x][1] > oldseqstop):
				exonorder.write(pep + '\t' + 'BM exons overlap' + '\n')
					
			oldseqstart = BMexons[x][0]
			oldseqstop = BMexons[x][1]
			oldframe = BMexons[x][2]

				
	if BMdir == -1:
		for x in sorted(BMexons, reverse = True):

			pepinfo = x.split('$')
			pep = pepinfo[1]
			
			if (oldseqstart == -1) or ((BMexons[x][0] < oldseqstart) and (BMexons[x][0] < oldseqstop)):

				Exonorder.append((pep,BMexons[x][0],BMexons[x][1],BMexons[x][2]))
				exonnums = exonnums + str(BMexons[x][0]) + '..' + str(BMexons[x][1]) + ','
					
			elif (BMexons[x][0] > oldseqstop) and (BMexons[x][2] < oldseqstop):
				exonorder.write(pep + '\t' + 'BM exons overlap' + '\n')
					
			oldseqstart = BMexons[x][0]
			oldseqstop = BMexons[x][1]
			oldframe = BMexons[x][2]
			
	n = 0
	for y in Exonorder:
		n += 1
		exonorder.write(str(n) + '\t' +y[0] + '\t' + str(y[1]) + '\t' + str(y[2]) + '\t' + str(y[3]) + '\n')

	if GBKOUT == 1:
		if BMdir == 1:
			f.write('     BLASTCDS             ' + exonnums[:-1] + '\n                     /shared_id="' + 'Consol-' + pep + '"\n                     /note="' + 'Consol-' + pep + '"\n')
		if BMdir == -1:
			f.write('     BLASTCDS             complement(' + exonnums[:-1] + '\n                     /shared_id="' + 'Consol-' + pep + '"\n                     /note="' + 'Consol-' + pep + '"\n')

	return Exonorder

#Main part of divergence analysis analysis program
	
while True:

	try:
	
		# raise an exception if there is an error
		line = raw_input()
			
	except:
		mode = 10
#		break	#Exits if no output must be printed
		
	data = line.split('\t')
	newest = data[2][1:-1]
	newcontig = data[1][1:-1]
	newpep = data[0][1:-1]
	newChr = str(data[4])

	#Create new chromosomal folder for genbank files 
	
	if newChr != oldChr:
		os.system("mkdir genbank/Chr" + str(newChr))
	
	#finnish off old and create new scaffold genbank file
	
	try:
		if GBKOUT == 1:
			if newcontig != oldcontig != 'x':
				if NOOVERLAP == 1 and oldpep != 'x':
					HAHVGene = CONPEPSEQ(pepexons,oldcontig,oldpep,dir,TotalPep)
					HAGenes = HAGenes + HAHVGene[0]
					HVGenes = HVGenes + HAHVGene[1]
					if len(HAHVGene[0]) > 0:
						numgenes += 1
					pepexons = {}
					dir = 0
					oldpep = 'x'
				addHAend(oldcontig,oldChr)
				f = open('temp/makefile.gbk','w')
	except:
			He.write(newcontig + "\n")
	
#Create new peptide
	
	try:
		if newpep != oldpep:
			if NOOVERLAP == 1 and oldpep != 'x':
				HAHVGene = CONPEPSEQ(pepexons,oldcontig,oldpep,dir,TotalPep)
				HAGenes = HAGenes + HAHVGene[0]
				HVGenes = HVGenes + HAHVGene[1]				
				if len(HAHVGene[0]) > 0:
					numgenes += 1
			pepexons = {}
			BMexons = {}
			BMdir = 0
			dir = 0

			
			print(newcontig + "(" + newpep + ")")
			He.write(newcontig + "(" + newpep + ") start\n")
			peppos = []
			peppos = protein(newpep,newcontig)
			exonP = []
			exonP = peppos[4]
			TotalPep = []
			TotalPep = CONSOLPRO(exonP)
	
	except:
			He.write(newcontig + "(" + newpep + ") outbreak\n")

#Finnish old and create new Chromosomal file.
			
	try:
		if HAHVGene != "":
			if (newChr != oldChr) and (oldChr != 'x'):
				if PYCOGENTOUT == 1:
					RUNPYCOGENT(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
				if PAMLOUT == 1:
					RUNPAML(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
				if MYDNDSOUT == 1:
					MYDNDS(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
				WriteChrSeq = open(str("seq/" + "Chr" + str(oldChr) + "Seq.fasta"),'w')
				WriteChrSeq.write("Chr" + str(oldChr) + "Seq" + "\nHelicoverpa  \n" + HAGenes + "\nHeliothis  \n" + HVGenes)
				WriteChrSeq.close()
				totHAGenes = totHAGenes + HAGenes
				totHVGenes = totHVGenes + HVGenes
				totnumgenes = totnumgenes + numgenes
				if oldChr in ("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"):
					totHAautoGenes = totHAautoGenes + HAGenes
					totHVautoGenes = totHVautoGenes + HVGenes
					totnumautogenes = totnumautogenes + numgenes
				if oldChr in ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"):
					totHAno0Genes = totHAno0Genes + HAGenes
					totHVno0Genes = totHVno0Genes + HVGenes
					totnumno0genes = totnumno0genes + numgenes
				HAGenes = ""
				HVGenes = ""
				HAHVGene = ""
				numgenes = 0
	except:
		print "Error in chromosomal PYCOGENT"

#Process new EST sequence
		
	try:
		if newest != oldest:
			newseq = mRNA(newest,newcontig,newpep,peppos[2],peppos[3],TotalPep)
			exonY = newseq[4]
			for x in exonY:
				num0s = 6-len(str(exonY[x][0]))
				pepexons[num0s*str(0) + str(str(exonY[x][0]) + '-' + newest + '-' + str(x))] = exonY[x]
				if (exonY[x][0] < exonY[x][1]) and (dir != -1):
					dir = 1
				elif (exonY[x][0] > exonY[x][1]) and (dir != 1):
					dir = -1
				else:
					exonorder.write(newcontig + '\t' + newpep + '\t' + 'HVexons in different directions' + '\n')
					He.write(newcontig + "(" + newest + ")- error sorting exons\n")
			
			if newseq[0]*newseq[3] != 0:
				seqlen = seqlen + int(newseq[0])
				HAcodons = HAcodons + newseq[1]
				HVcodons = HVcodons + newseq[2]
				seqnum += 1
				codonlens.append(newseq[0]/3)
				
	except:
			He.write(newcontig + "(" + newest + ")\n")

#Update old parameters
			
	oldcontig=newcontig	
	oldest = newest
	oldpep = newpep
	oldChr = newChr

#Enter final false
	
	if mode == 10:


	
		if NOOVERLAP == 1:
			HAHVGene = CONPEPSEQ(pepexons,oldcontig,oldpep,dir,TotalPep)
			HAGenes = HAGenes + HAHVGene[0]
			HVGenes = HVGenes + HAHVGene[1]
			if len(HAHVGene[0]) > 0:
				numgenes += 1

		print "Finishing last file"
		if PAMLOUT == 1:
			for x in codonlens:
				listcodonlens = listcodonlens + str(x) + " "
			HAS.write("\t2 " + str(seqlen) + "\nHelicoverpa  \n" + HAcodons + "\nHeliothis  \n" + HVcodons)

#write last genbank file	
	
		if GBKOUT == 1:
			addHAend(oldcontig,newChr)

#Perform divergece analysis on concatinated sequence
			
		try:
			if PYCOGENTOUT == 1:
				RUNPYCOGENT(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
			if PAMLOUT == 1:
				RUNPAML(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
			if MYDNDSOUT == 1:
				MYDNDS(str("Chromosome " + str(oldChr)),str("Chr-" + str(oldChr)),str("ngenes = " + str(numgenes)),HAGenes,HVGenes,1,oldChr)
			WriteChrSeq = open(str("seq/" + "Chr" + str(oldChr) + "Seq.fasta"),'w')
			WriteChrSeq.write("Chr" + str(oldChr) + "Seq" + "\nHelicoverpa  \n" + HAGenes + "\nHeliothis  \n" + HVGenes)
			WriteChrSeq.close()
				
			totHAGenes = totHAGenes + HAGenes
			totHVGenes = totHVGenes + HVGenes
			totnumgenes = totnumgenes + numgenes
			if oldChr in ("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"):
				totHAautoGenes = totHAautoGenes + HAGenes
				totHVautoGenes = totHVautoGenes + HVGenes
				totnumautogenes = totnumautogenes + numgenes
			if oldChr in ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"):
				totHAno0Genes = totHAno0Genes + HAGenes
				totHVno0Genes = totHVno0Genes + HVGenes
				totnumno0genes = totnumno0genes + numgenes
			
			if PYCOGENTOUT == 1:			
				RUNPYCOGENT(str("TotalAutosomes"),str("Chr-auto"),str("ngenes = " + str(totnumautogenes)),totHAautoGenes,totHVautoGenes,1,29)
			if PAMLOUT == 1:
				RUNPAML(str("TotalAutosomes"),str("Chr-auto"),str("ngenes = " + str(totnumautogenes)),totHAautoGenes,totHVautoGenes,1,29)
			if MYDNDSOUT == 1:
				MYDNDS(str("TotalAutosomes"),str("Chr-auto"),str("ngenes = " + str(totnumautogenes)),totHAautoGenes,totHVautoGenes,1,29)
			WriteChrSeq = open(str("seq/" + "Chr-autoSeq.fasta"),'w')
			WriteChrSeq.write("Chr-autoSeq" + "\nHelicoverpa  \n" + totHAautoGenes + "\nHeliothis  \n" + totHVautoGenes)
			WriteChrSeq.close()
			
			if PYCOGENTOUT == 1:
				RUNPYCOGENT(str("TotalChromosomes-no0"),str("Chr-no0"),str("ngenes = " + str(totnumno0genes)),totHAno0Genes,totHVno0Genes,1,30)
			if PAMLOUT == 1:
				RUNPAML(str("TotalChromosomes-no0"),str("Chr-no0"),str("ngenes = " + str(totnumno0genes)),totHAno0Genes,totHVno0Genes,1,30)
			if MYDNDSOUT == 1:
				MYDNDS(str("TotalChromosomes-no0"),str("Chr-no0"),str("ngenes = " + str(totnumno0genes)),totHAno0Genes,totHVno0Genes,1,30)
			WriteChrSeq = open(str("seq/" + "Chr-no0Seq.fasta"),'w')
			WriteChrSeq.write("Chr-no0Seq" + "\nHelicoverpa  \n" + totHAno0Genes + "\nHeliothis  \n" + totHVno0Genes)
			WriteChrSeq.close()
			
			if PYCOGENTOUT == 1:			
				RUNPYCOGENT(str("TotalChromosomes"),str("Chr-all"),str("ngenes = " + str(totnumgenes)),totHAGenes,totHVGenes,1,31)
			if PAMLOUT == 1:
				RUNPAML(str("TotalChromosomes"),str("Chr-all"),str("ngenes = " + str(totnumgenes)),totHAGenes,totHVGenes,1,31)
			if MYDNDSOUT == 1:
				MYDNDS(str("TotalChromosomes"),str("Chr-all"),str("ngenes = " + str(totnumgenes)),totHAGenes,totHVGenes,1,31)
			WriteChrSeq = open(str("seq/" + "Chr-allSeq.fasta"),'w')
			WriteChrSeq.write("Chr-allSeq" + "\nHelicoverpa  \n" + totHAGenes + "\nHeliothis  \n" + totHVGenes)
			WriteChrSeq.close()
				
		except:
			print "Error in final chromosomal PYCOGENT"
			
		HAS.close()
		Pa.close()
		He.close()
		ex.close()
		PAMLres.close()
		PYCOGENTres.close()
		MYDNDSres.close()
		MYDNDSlongres.close()
		exonerateres.close()
		sequences.close()
		nofseq.close()
		exonmix.close()
		exonorder.close()
		framecheck.close()
		
		print "Genbank files created"
		break



