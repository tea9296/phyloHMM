#!/opt/local/bin/python2.7
import sys, os, operator
import numpy as np
from HMM import *
from Bio import SearchIO, AlignIO 

###################################################################################
#################################   PARAMETERS   ##################################
###################################################################################

#Input arguments
pfam_output_file = sys.argv[1]  #Pfam output file
query_id = sys.argv[2]  #Query seed id
alignment_file = sys.argv[3]  #Trimmed alignment file with full-length seed seq
aa_freqs_output = sys.argv[3]+".aa_freq"  #aa_freq output file name

#File directories
hmm_file_dir = '/Users/ding/Desktop/for_ding/'
#hmm_file_dir = '/Users/dinghe/Dropbox/Projects/mito_origin/Simon/pfam/hmm/'
hmm_file_basenames = os.listdir(hmm_file_dir)
alignment_file_dir = '/Users/ding/Desktop/for_ding/'
#alignment_file_dir = '/Users/dinghe/Dropbox/Projects/mito_origin/Data/SG_alignment/20160519/'

#Empirical aa_freqs from LG model ordered as ARNDCQEGHILKMFPSTWYV
LG_aa_freqs = '2.34595 2.91338 3.49262 3.01205 4.45606 3.52632 2.83988 2.41146 3.78460 2.67538 2.34684 3.06526 3.47661 3.10027 3.17541 2.95790 2.93381 4.53976 3.46891 2.49357'

###################################################################################
##################################   FUNCTIONS   ##################################
###################################################################################

#Parse the hmmscan output file(s) and return the information needed for alignment-hmm sites mapping
def pfam_parser(out_file, qurey_id):
	
	domain_info = {}  # store domain(s) info
	
	for QueryResult in SearchIO.parse(out_file,'hmmer3-text'):
		if qurey_id in str(QueryResult.id):
			#print QueryResult.description
			for hit in QueryResult.hits:
				for HSPFragments in hit.hsps:
					query_hit_start = HSPFragments.query_start+1 # parsed start position appears to be 1 site short
					query_hit_seq = HSPFragments.query.seq
					query_hit_end = HSPFragments.query_end
					domain_hit_start = HSPFragments.hit_start+1
					domain_hit_end = HSPFragments.hit_end
					
				domain_info[hit.id] = [query_hit_start, query_hit_end, query_hit_seq, domain_hit_start, domain_hit_end]
	
	return domain_info

#Parse the the original alignments (with seed seq)
def alignment_parser(alignment_file):
	alignment = AlignIO.read(alignment_file, "fasta")
	align_array = np.array([list(rec) for rec in alignment], np.character, order="F")  #Convert alignment object to np_array object
	align_columns = np.transpose(align_array)  #align_columns.shape == (#col, #raw), i.e. (#site, #seq)
	#Loop the columns to define the site mapping boundaries (using the full seed seq to map onto the trimmed alignment)
	align_mapping = []
	for col in align_columns:
		#Check if only the seed seq in the column
		col_check = col.tolist()
		col_check = filter(lambda a: a != '-', col_check)
		if len(col_check) == 1:  #Only seed seq present
			align_mapping.append(0)
		elif len(col_check) > 1:
			align_mapping.append(1)

	return align_mapping  #The structure looks like [0,0,0,1,1,1,0,0,1,...0] where 1 indicates mapped onto trimmed alignment

###################################################################################
###############################   MAIN PROCESSES   ################################
###################################################################################

seed_to_alignment_mapping = alignment_parser(alignment_file_dir+alignment_file)
domain = pfam_parser(pfam_output_file, query_id)
seed_seq_master = []  #The master list storing aa_freqs for the entire seed_seq
#Getting the aa_freqs from respective hmms, storing structure see HMM.py
for hmm_name, hit_info in sorted(domain.items(), key=operator.itemgetter([1][0])):  #The aa_seqs should be recorded in from the N-terminal to C-terminal of the seed seq
	if hmm_name+'.hmm' in hmm_file_basenames: hmm = HMM(hmm_file_dir+hmm_name+'.hmm')  #Instantiate a HMM class from a hmm file
	seq_start = int(hit_info[0])  #Seed seq match start
	seq_end = int(hit_info[1])  #Seed seq match end
	seed_seq = str(hit_info[2])  #Domain matching seed seq
	dm_start = int(hit_info[3])  #Doamin match start
	dm_end = int(hit_info[4])  #Doamin match end
	
	#print hmm_name  #Check
	#print len(hmm.hmm_aa_freq['A'][0])  #Check
	
	#Mapping the aa_freqs to seed_seq sites according to the hmmscan output
	hmmscan_seed_seq = list(seed_seq)  #A list storing hmmscan_seed_seq site-by-site
	j = 0  #Check for insert state so the pointer maps correctly to the domain position
	for i, seed_seq_site in enumerate(hmmscan_seed_seq):
		if seed_seq_site.isalpha() and seed_seq_site.isupper():
			#Replace the site with a list of matching-state aa_freqs in the order of 'ARNDCQEGHILKMFPSTWYV'
			#print i, seed_seq_site, hmm.hmm_aa_freq['A'][0][i+dm_start-1-j]
			hmmscan_seed_seq[i] = [hmm.hmm_aa_freq['A'][0][i+dm_start-1-j]]
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['R'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['N'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['D'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['C'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['Q'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['E'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['G'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['H'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['I'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['L'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['K'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['M'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['F'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['P'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['S'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['T'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['W'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['Y'][0][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['V'][0][i+dm_start-1-j])
		elif seed_seq_site.isalpha() and seed_seq_site.islower():
			#Replace the site with a list of insert-state aa_freqs in the order of 'ARNDCQEGHILKMFPSTWYV'
			j += 1
			#print i, seed_seq_site  #Check
			hmmscan_seed_seq[i] = [hmm.hmm_aa_freq['A'][1][i+dm_start-1-j]]
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['R'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['N'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['D'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['C'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['Q'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['E'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['G'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['H'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['I'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['L'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['K'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['M'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['F'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['P'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['S'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['T'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['W'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['Y'][1][i+dm_start-1-j])
			hmmscan_seed_seq[i].append(hmm.hmm_aa_freq['V'][1][i+dm_start-1-j])
		elif seed_seq_site == '-':
			continue  #Skipping the deletion state
		else:
			print "Something is not right!!!"
			break
	
	seed_seq_master += hmmscan_seed_seq

for i, aa_freqs in enumerate(seed_seq_master):
	if isinstance(aa_freqs, list):
		seed_seq_master[i] = " ".join(aa_freqs)  #Construct the aa_freq string

seed_seq_master = filter(lambda a: a != "-", seed_seq_master)  #Discard deletion states "-"

#################################################################################################
######   IMPORTANT PROCESSES TO CORRECTLY MAP THE PFAM DOMAIN ONTO THE ORIGINAL ALIGNMENT  ######
#################################################################################################
seed_to_alignment_mapping_check = seed_to_alignment_mapping[:]  #Duplicate the seed_to_alignment_mapping for checking the sites within the found domains but not in the trimmed alignment
#Loop the domain hits in the order of from the N-terminal to C-terminal of the seed seq
for id, info in sorted(domain.items(), key=operator.itemgetter([1][0])):
	seq_start = int(info[0])  #Seed seq match start
	seq_end = int(info[1])  #Seed seq match end
	#This loop maps aa_freqs onto the positions where the domains were found based on the seed seq
	#So the seed_to_alignment_mapping should have aa_freqs regardless of if sites present in the trimmed alignment
	for i in range(seq_start, seq_end+1):
		seed_to_alignment_mapping[i-1] = seed_seq_master.pop(0)  #First position is 0

#This loop dicards the domain aa_freqs where not present in the trimmed alignment
for i, aa_freq in enumerate(seed_to_alignment_mapping):
	if seed_to_alignment_mapping_check[i] == 0:
		seed_to_alignment_mapping[i] = 0

#################################################################################################
#################################################################################################

#Output the aa_freq
aa_freqs_output_file = open(aa_freqs_output, 'w')

seed_to_alignment_mapping = filter(lambda a: a != 0, seed_to_alignment_mapping)  #Discard sites not in trimmed alignment: 0 

for aa_freq in seed_to_alignment_mapping:
	if aa_freq == 1:
		aa_freqs_output_file.write(LG_aa_freqs+"\n")
	else:
		aa_freqs_output_file.write(str(aa_freq)+"\n")

aa_freqs_output_file.close()

