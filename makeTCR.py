#!/usr/bin/python

'''

Authors: Allyson Byrd and Oliver Harrison

##############################################
#03/16/2021
##############################################
Updates: 
	- Modified to work with python3
	- Removed print colored and use of modules
	- Added usage of argparse and checking of inputs

Usage: python makeTCR.py -q input.fasta --alpha -o output_folder_name

##############################################
#06/27/2016
##############################################

Purpose: to generate the sequence for TCR transgenics
make sure the "sequences" folder is somewhere and that location is indicated on line 47
this does require biopython and blast, so before you can run jobs, you have to
#module load BLAST+/2.2.31-goolf-1.7.20-Python-2.7.9
#module load Biopython/1.65-goolf-1.7.20-Python-2.7.9

Usage: python makeTCR.py input.fasta alpha output_folder_name

'''


#current implementation is for a single alpha or beta sequence 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Blast.Applications import NcbitblastnCommandline  
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_protein
#from Bio.Alphabet import generic_dna
from Bio import Entrez
#from termcolor import colored

import argparse, os, glob, sys 
#https://docs.python.org/2/howto/argparse.html

parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description="Write a description here.")
#Input parameters
parser.add_argument("-q", "--queryFile", help="Fasta file of the input sequence", action="store", default="", required=True)
parser.add_argument("--alpha", help="Sequence is alpha", action="store_true")
parser.add_argument("--beta", help="Sequence is beta", action="store_true")
parser.add_argument("-o", "--output_dir", help="Name of the output folder", action="store", default="./", required=False)
parser.add_argument("-db", "--db_dir", help="Directory where the databases are located", action="store", default="/Users/oharrison/Desktop/Terminal_TCR_program/sequences/", required=False)

#db_dir = "/hpcdata/lpd_mis/Ollie/TCR_program/sequences/" #directory where the databases are located 	
#db_dir = "/Users/oharrison/Desktop/Terminal_TCR_program/sequences/" #directory where the databases are located

#Parsing the input parameters 
args = parser.parse_args()

queryFile = args.queryFile #fasta file of the input sequence
output_dir = args.output_dir #name of the output folder
db_dir = args.db_dir #location of databases

#queryFile = sys.argv[1] #fasta file of the input sequence
#type = sys.argv[2] #alpha or beta 
#output_dir = sys.argv[3] #name of the output folder

print ("\n##### Step 0: Checking inputs are correct.")

#checking to make sure the input query file exists
if not os.path.isfile(queryFile):
	print("Error: %s does not exist. Please input a new fasta file. Exitting!" % queryFile)
	sys.exit()
	
#checking to make sure the database directory exists
if not os.path.isdir(db_dir):
	print("Error: %s does not exist. Please input a new database directory. Exitting!" % db_dir)
	sys.exit()

if args.alpha and args.beta:
	print("Must pick alpha OR beta, not both. Exitting!")
	sys.exit()
elif args.alpha:
	alpha = True
	print("\nYour query is an ALPHA sequence!\n")
elif args.beta:
	alpha = False
	print ("\nYour query is a BETA sequence!\n")
else:
	print("Must include --alpha or --beta.")
	#print ("Second input argument needs to be 'alpha' or 'beta' not '%s" % type, 'red')
	print("Exiting!")
	sys.exit()
	
#Verifying the primary output directory exists
if not os.path.isdir(output_dir):	
	print ("Creating the directory %s that will store the output" % output_dir)
	cline = "mkdir %s" % output_dir
	os.system(cline)
elif output_dir == "./":
	print("No output directory was specified so outputs will be written to cwd.")

#output files

BLAST_out = "%s/BLAST_out.txt" % output_dir
BLAST2_out = "%s/BLAST_out_alignments.txt" % output_dir
seqs_out = "%s/sequences.fa" % output_dir

FILE_BLAST = open(BLAST_out, "w")
FILE_BLAST2 = open(BLAST2_out, "w")
FILE_seqs = open(seqs_out, "w")

##################### Setting Up #####################

#making a reverse complement version of the input file and saving the sequence to work with later 
print ("\n##### Step 1: Creating the reverse complement of your input sequence")

rc_query = "%s/%s_reverseComplement.fa" % (output_dir, queryFile)

FILE_in = open(queryFile, 'r')
FILE_out = open(rc_query, 'w')
seq = ''
for line in FILE_in:

	if line.startswith(">"):
		if seq != '':
			seq = Seq(line)
			reverse_c = seq.reverse_complement()
			FILE_out.write("%s\n" % reverse_c)
		FILE_out.write(line)
		seq_name = line.strip()[1:]
		seq = ''
	else:
		seq += line.strip()

seq = Seq(line)
reverse_c = seq.reverse_complement()
FILE_out.write("%s\n" % reverse_c)

FILE_in.close()
FILE_out.close()

FILE_seqs.write(">Sanger_reverseComplement\n")
FILE_seqs.write("%s\n" % reverse_c)

#acquiring stats about your input sequence 
print ("\nStats for input sequence '%s':" % seq_name)
print ("\tNucleotide length:\t%s" % len(seq))
print("\t# N bases:\t%s" % seq.count('N'))
if seq.count('N') > 0:
	print ("\nWARNING there are ambigous bases in your sequence, you should check them more closely!")

##################### Reporting the V region #####################

#blasting against database of V regions 
print ("\n##### Step 2: Searching for the V sequence in the IMGT database")

if alpha:
	V_db = "%s/TRAV.fa" % db_dir
else:
	V_db = "%s/TRBV.fa" % db_dir 

#https://www.biostars.org/p/88944/
blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=V_db, evalue=1, outfmt="7 sseqid length sstart send qstart qend qlen slen evalue pident", max_target_seqs="10", num_threads = 10)
#print blastcline
FILE_BLAST.write("%s\n\n" % blastcline)
stdout, stderr = blastcline()
#print stdout
FILE_BLAST.write("%s\n" % stdout)
stdout = stdout.split("\n")
count = 0
for line in stdout:
	if "# 0 hits" in line:
		print ("There were no matching V sequences... Trouble shooting necessary :(")
		print ("Exiting!")
		sys.exit()
	elif not line.startswith("#") and not line == "":
		line = line.split("\t")
		subject = line[0]
		align_length = int(line[1])
		s_end = int(line[3])
		s_start = int(line[2])
		id = float(line[9])
		if count == 0:
			print ("\nTop V hit(s) in the IMGT database:")
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
			
			if id < 100:
				print ("\t\tWarning: Your top hit has <100.0 identify to a sequence in the database")
			
		elif id > 99.9:
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
		else:
			break
		
		count += 1
		

print ("\n##### Step 3: Searching for the V sequence in the NCBI coding sequences")

if alpha:
	V2_db = "%s/alpha_coding_sequences.fa" % db_dir
else:
	V2_db = "%s/beta_coding_sequences.fa" % db_dir 

#blasting against the alpha coding sequences and for now taking the first hit 
blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=V2_db, evalue=1, outfmt="7 sseqid length sstart send qstart qend qlen slen evalue pident", max_target_seqs="10", num_threads = 10)
#print blastcline
FILE_BLAST.write("\n%s\n\n" % blastcline)
stdout, stderr = blastcline()
#print stdout
FILE_BLAST.write("%s\n" % stdout)
stdout = stdout.split("\n")
count = 0

V_seq_ids = {} #keys: first BLAST hit to a V seq + all hits with perfect sequence id  value: list of pertinent information for this hit [subject, s_start, 'full seq', 'desired proportion']
V_seq_ids_order = [] #stores the order seq_ids occurred in the BLAST, important so write them out in the same order

for line in stdout:
	if "# 0 hits" in line:
		print ("There were no matching V sequences... Trouble shooting necessary :(")
		print ("Exiting!")
		sys.exit()
	elif not line.startswith("#") and not line == "":
		line = line.split("\t")
		subject = line[0]
		align_length = int(line[1])
		subject_protein_id = subject.split("_")[2]
		s_end = int(line[3])
		s_start = int(line[2])
		id = float(line[9])
		if count == 0:
			print ("\nTop V hit(s) in the NCBI database:")
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
			
			V_seq_ids[subject_protein_id] = [subject, s_start, '', '']
			V_seq_ids_order.append(subject_protein_id)
			
			if id < 100:
				print ("\t\tWarning: Your top hit has <100.0 identify to a sequence in the database")
				print ("\t\tPlease see detailed BLAST output in %s" % BLAST2_out)
				perfectFound = False
			else:
				perfectFound = True 
			
		elif id > 99.9:
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
			V_seq_ids[subject_protein_id] = [subject, s_start, '', '']
			V_seq_ids_order.append(subject_protein_id)
		else:
			break
		
		count += 1

#when perfect blast not found
#if not perfectFound:
if True:
	blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=V2_db, evalue=1, outfmt="0", max_target_seqs="5", num_threads = 10)
	FILE_BLAST2.write("\n%s\n\n" % blastcline)
	stdout, stderr = blastcline()
	FILE_BLAST2.write("%s\n" % stdout)
		
		
print ("\n##### Step 3B: Pulling coding sequences for matching V sequences from NCBI")
#pulling the coding sequence from V2b that match seqids in V_seq_ids

if len(V_seq_ids) > 1:
	print ("\nMultiple V sequences had 100% identity to the input query")
	print ("\tSequences for all options will be reported in %s/sequences.fa" % output_dir)


coding_File = open(V2_db, "r")
for line in coding_File:
	if line.startswith('>'):
		header = line.strip().replace(">","")
		protein_id = header.split("_")[2]
	else:
		if protein_id in V_seq_ids:
		#if V_seq_ids.has_key(protein_id):
			V_seq_ids[protein_id][2] += line.strip()
			
for seq_id in V_seq_ids_order:
	
	alpha_seq = V_seq_ids[seq_id][2]
	s_start = V_seq_ids[seq_id][1]
	
	desired_alpha = alpha_seq[0:s_start-1]
	
	V_seq_ids[seq_id][3] = desired_alpha
	if len(desired_alpha) > 10:
		FILE_seqs.write(">%s | 0:%s\n" %  (V_seq_ids[seq_id][0], s_start) )
		FILE_seqs.write("%s\n" % desired_alpha)
	else:
		print ("\n\tThe sequence pulled for %s is <10 bp, likely not a real hit therefore excluding" % seq_id)

#####################Reporting the D region #####################
# if not alpha:
# 	print "\nStep 3.5: Searching for the perfect D"
# 	
# 	TRBD1 = "gggacagggggc"
# 	TRBD2 = "gggactgggggggc"
# 	
# 	lower_rc = str.lower(str(reverse_c))
# 	
# 	print lower_rc
# 	
# 	if TRBD1 in lower_rc:
# 		print "\tTRBD1 was found in the input sequence"
# 		if TRBD2 in lower_rc:
# 			print "\tTRBD2 was also found in the input sequence...Worry"
# 	elif TRBD2 in reverse_c: 
# 		print "\tTRBD2 was found in the input sequence"
# 	else:
# 		print "\tThere was NOT a perfect match to TRBD1 or 2 in your input sequence"
#####################Reporting the J region #####################
 
print ("\n##### Step 4: Searching for the J sequence in the IMGT database")

if alpha:
	J_db = "%s/TRAJ.fa" % db_dir
else:
	J_db = "%s/TRBJ.fa" % db_dir 

blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=J_db, evalue=1, outfmt="7 sseqid length sstart send qstart qend qlen slen evalue pident", max_target_seqs="10", num_threads = 10)
#print blastcline
FILE_BLAST.write("\n%s\n\n" % blastcline)
stdout, stderr = blastcline()
#print stdout
FILE_BLAST.write("%s\n" % stdout)
stdout = stdout.split("\n")
count = 0

for line in stdout:
	if "# 0 hits" in line:
		
		if alpha:
			print ("\tThere were no matching J sequences... Maybe concerned?")
		else:
			print ("\tThere were no matching J sequences, so cannot correctly select a Constant sequence")
			print ("\tExiting!")
			sys.exit()
		
	elif not line.startswith("#") and not line == "":
		line = line.split("\t")
		subject = line[0]
		align_length = int(line[1])
		s_end = int(line[3])
		s_start = int(line[2])
		id = float(line[9])
		if count == 0:
			print ("\nTop J hit(s) in the IMGT database:")
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
			
			j_type = subject.split("|")[1].split("-")[0]
			
			if id < 100:
				print ("\t\tWarning: Your top hit has <100.0 identify to a sequence in the database")
				print ("\t\tPlease see detailed BLAST output in %s" % BLAST2_out)
				perfectFound = False
			else:
				perfectFound = True 
			
		elif id > 99.9:
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
		else:
			break
		
		count += 1

#when perfect blast not found
if not perfectFound:
	blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=J_db, evalue=1, outfmt="0", max_target_seqs="5", num_threads = 10)
	FILE_BLAST2.write("\n%s\n\n" % blastcline)
	stdout, stderr = blastcline()
	FILE_BLAST2.write("%s\n" % stdout)


#####################Finding the necessary portion of the constant region#####################

print ("\n##### Step 5: Searching for the constant region")

if alpha:
	C_db = "%s/a_constant.fa" % db_dir
else:
	if j_type == "TRBJ1":
		C_db = "%s/CB1.fa" % db_dir 
	else:
		assert j_type == "TRBJ2"
		C_db = "%s/CB2.fa" % db_dir 
	print ("\nThe J region is type %s therefore database %s will be used" % (j_type, C_db))

#had to decrease the word size to be able to detect smaller hits 
blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=C_db , evalue=100, outfmt="7 sseqid length sstart send qstart qend qlen slen evalue pident", max_target_seqs="10", num_threads = 10, word_size=10)
#print blastcline
FILE_BLAST.write("\n%s\n\n" % blastcline)
stdout, stderr = blastcline()
#print stdout
FILE_BLAST.write("%s\n" % stdout)
stdout = stdout.split("\n")
count = 0

#this shouldn't really be necessary since there should only really be one C sequence hit 
C_seq_ids = {} #keys: first BLAST hit to a C seq + all hits with perfect sequence id  value: list of pertinent information for this hit [subject, s_start, 'full seq', 'desired proportion']
C_seq_ids_order = [] #stores the order seq_ids occurred in the BLAST, important so write them out in the same order

for line in stdout:
	if "# 0 hits" in line:
		print ("\tThere were no matching C sequences... Trouble shooting necessary :(")
		print ("\tExiting!")
		sys.exit()
	elif not line.startswith("#") and not line == "":
		line = line.split("\t")
		subject = line[0]
		subject_protein_id = subject
		align_length = int(line[1])
		s_end = int(line[3])
		s_start = int(line[2])
		id = float(line[9])
		if count == 0:
			print ("\nTop C hit(s) in the IMGT database:")
			print ("\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length))
			
			C_seq_ids[subject_protein_id] = [subject, s_end, '', '']
			C_seq_ids_order.append(subject_protein_id)
			
			if s_start != 1:
				print ("\tWarning: The constant sequence did not start aligning at basepair 1, but basepair %s" % s_start)
			
			if id < 100:
				print ("\tWarning: Your top hit has <100.0 identify to a sequence in the database")
				print ("\t\tPlease see detailed BLAST output in %s" % BLAST2_out)
				perfectFound = False
			else:
				perfectFound = True 
			
		elif id > 99.9:
			#only want the first sequence
			pass
			#print "\t%s\t%s identity\t%s nt overlap" % (subject, id, align_length)
			#C_seq_ids[subject_protein_id] = [subject, s_end, '', '']
			#C_seq_ids_order.append(subject_protein_id)
		else:
			break
		
		count += 1

#when perfect blast not found
if not perfectFound:
	blastcline = NcbitblastnCommandline(cmd="blastn", query=rc_query, db=C_db, evalue=1, outfmt="'0'", max_target_seqs="5", num_threads = 10, word_size=10)
	FILE_BLAST2.write("\n%s\n\n" % blastcline)
	stdout, stderr = blastcline()
	FILE_BLAST2.write("%s\n" % stdout)


print ("\n##### Step 5B: Pulling coding sequences for matching C sequences")
#pulling the coding sequence from V2b that match seqids in V_seq_ids

if len(C_seq_ids) > 1:
	print ("\nMultiple C sequences had 100% identity to the input query")
	print ("\tSequences for all options will be reported in %s/sequences.fa" % output_dir)


constant_FILE = open(C_db, "r")
for line in constant_FILE:
	if line.startswith('>'):
		header = line.strip().replace(">","")
		protein_id = header.split(" ")[0]
	else:
		if protein_id in C_seq_ids:
		#if C_seq_ids.has_key(protein_id):
			C_seq_ids[protein_id][2] += line.strip()

			
for seq_id in C_seq_ids_order:
	
	alpha_seq = C_seq_ids[seq_id][2]
	s_end = C_seq_ids[seq_id][1]
	
	desired_alpha = alpha_seq[s_end:]
	
	C_seq_ids[seq_id][3] = desired_alpha
	if len(desired_alpha) > 10:
		FILE_seqs.write(">%s | %s:%s\n" %  (C_seq_ids[seq_id][0], s_end, len(desired_alpha)) )
		FILE_seqs.write("%s\n" % desired_alpha)
	else:
		print ("\n\tThe sequence pulled for %s is <10 bp, likely not a real hit therefore excluding" % seq_id)

FILE_BLAST.close()
FILE_BLAST2.close()

#combing the sequences that were the first in the blast output file 

print ("\n##### Step 6: Combing sequences that were the first hit in the BLAST output files")


V_seq = V_seq_ids[V_seq_ids_order[0]][3]
C_seq = C_seq_ids[C_seq_ids_order[0]][3]

combined = "%s%s%s" % (V_seq, reverse_c, C_seq)

print ("\n>combined:%s|ReverseComplementSanger|%s" % (V_seq_ids_order[0], C_seq_ids_order[0]))
print (combined)

FILE_seqs.write(">combined:%s|ReverseComplementSanger|%s\n" % (V_seq_ids_order[0], C_seq_ids_order[0]) )
FILE_seqs.write("%s\n" % combined)


coding_dna = Seq(combined) #, generic_dna
protein = coding_dna.translate()

print ("\n>combined_protein")
print (protein)

FILE_seqs.write(">combined_protein:%s|ReverseComplementSanger|%s\n" % (V_seq_ids_order[0], C_seq_ids_order[0]) )
FILE_seqs.write("%s\n" % protein)

FILE_seqs.close()


print ("\nFinished!")
