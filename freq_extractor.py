#!/anaconda2/bin/python
import operator, argparse, subprocess, re, mmap
import numpy as np
from Bio import SearchIO, AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

###################################################################################
#################################   PARAMETERS   ##################################
###################################################################################

# Input arguments
parser = argparse.ArgumentParser(description='Commandline options for freq_extractor')
parser.add_argument('--ali', help='input alignment')
parser.add_argument('--database', help='hmm database')
parser.add_argument('--hmmer_path', help='HMMER executable file path')
args = parser.parse_args()  # e.g. input alignment file: args.ali; output freq_file: args.out

#pfam_output_file = sys.argv[1]  #Pfam output file
#query_id = sys.argv[2]  #Query seed id
#alignment_file = sys.argv[3]  #Trimmed alignment file with full-length seed seq
#aa_freqs_output = sys.argv[3]+".aa_freq"  #aa_freq output file name

# Directories
hmmscan_exe = str(args.hmmer_path)+'hmmscan'
hmmfetch_exe = str(args.hmmer_path)+'hmmfetch'
#hmm_file_dir = '/Users/ding/Desktop/for_ding/'
#hmm_file_dir = '/Users/dinghe/Dropbox/Projects/mito_origin/Simon/pfam/hmm/'
#hmm_file_basenames = os.listdir(hmm_file_dir)
#alignment_file_dir = '/Users/ding/Desktop/for_ding/'
#alignment_file_dir = '/Users/dinghe/Dropbox/Projects/mito_origin/Data/SG_alignment/20160519/'

# Empirical aa_freqs from LG model ordered as ARNDCQEGHILKMFPSTWYV
LG_aa_freqs = '2.34595 2.91338 3.49262 3.01205 4.45606 3.52632 2.83988 2.41146 3.78460 2.67538 2.34684 3.06526 3.47661 3.10027 3.17541 2.95790 2.93381 4.53976 3.46891 2.49357'

# Global variables
global seed_seq_m_start
global seed_seq_m_end
###################################################################################
####################################   CLASS   ####################################
###################################################################################
class HMM:
    
    def __init__(self, hmm_file, seed_seq_id):
        # creating a HMM structure for the target hmm profile(s) -> {profile_name_1:[{hmm_aa_freq_1},profile_length_1], profile_name_2:[{hmm_aa_freq_2},profile_length_2],...}
        self.hmm_database = {}
        match_freq_lines_pattern = re.compile(r"\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\d+\s+[A-Za-z].+")
        insert_freq_lines_pattern = re.compile(r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")

        # separate hmm_profile into blocks in the hmm_file, and parse profile information: NAME, LENG and aa_freq
        with open(hmm_file, 'r+') as hmmF:
            mmap_hmmF = mmap.mmap(hmmF.fileno(), 0)
            target_hmm_profile_match = re.search('NAME\s+(.*)', mmap_hmmF)
            if target_hmm_profile_match:
                profile_name = target_hmm_profile_match.group(1).strip('\r')
                print 'Found target hmm profile: %s' % (profile_name)
                # running hmmfetch to retrieve the target hmm profile
                target_hmm_lines = subprocess.check_output([hmmfetch_exe, args.database, profile_name])
                target_hmm_lines = target_hmm_lines.split("//")[0]
                profile_len = re.split("LENG\s+", target_hmm_lines)[1].split("\n")[0]
                match_freq_lines = re.split("HMM\s+", target_hmm_lines)[1].split("\n")[5:-1][::3]  # first match starts from the 5th index line and then from it every 3rd line, the last line is empty
                insert_freq_lines = re.split("HMM\s+", target_hmm_lines)[1].split("\n")[6:-1][::3]  # first match starts from the 6th index line and then from it every 3rd line, the last line is empty
                                
                # creating aa_freq structure -> (match_freq, insert_freq)
                hmm_aa_freq = {'A': [[],[]], 'R': [[],[]], 'N': [[],[]], 'D': [[],[]], 'C': [[],[]], 'Q': [[],[]], \
					'E': [[],[]], 'G': [[],[]], 'H': [[],[]], 'I': [[],[]], 'L': [[],[]], 'K': [[],[]], \
					'M': [[],[]], 'F': [[],[]], 'P': [[],[]], 'S': [[],[]], 'T': [[],[]], 'W': [[],[]], \
					'Y': [[],[]], 'V': [[],[]]}

                # start processing profile_feq
                for line in match_freq_lines:
                    m_freq = re.findall(match_freq_lines_pattern, line)
                    if m_freq is not None:
                        hmm_aa_freq['A'][0].append(m_freq[0][0])
                        hmm_aa_freq['C'][0].append(m_freq[0][1])
                        hmm_aa_freq['D'][0].append(m_freq[0][2])
                        hmm_aa_freq['E'][0].append(m_freq[0][3])
                        hmm_aa_freq['F'][0].append(m_freq[0][4])
                        hmm_aa_freq['G'][0].append(m_freq[0][5])
                        hmm_aa_freq['H'][0].append(m_freq[0][6])
                        hmm_aa_freq['I'][0].append(m_freq[0][7])
                        hmm_aa_freq['K'][0].append(m_freq[0][8])
                        hmm_aa_freq['L'][0].append(m_freq[0][9])
                        hmm_aa_freq['M'][0].append(m_freq[0][10])
                        hmm_aa_freq['N'][0].append(m_freq[0][11])
                        hmm_aa_freq['P'][0].append(m_freq[0][12])
                        hmm_aa_freq['Q'][0].append(m_freq[0][13])
                        hmm_aa_freq['R'][0].append(m_freq[0][14])
                        hmm_aa_freq['S'][0].append(m_freq[0][15])
                        hmm_aa_freq['T'][0].append(m_freq[0][16])
                        hmm_aa_freq['V'][0].append(m_freq[0][17])
                        hmm_aa_freq['W'][0].append(m_freq[0][18])
                        hmm_aa_freq['Y'][0].append(m_freq[0][19])

                for line in insert_freq_lines:
                    i_freq = re.findall(insert_freq_lines_pattern, line)
                    if i_freq is not None:
                        hmm_aa_freq['A'][1].append(i_freq[0][0])
                        hmm_aa_freq['C'][1].append(i_freq[0][1])
                        hmm_aa_freq['D'][1].append(i_freq[0][2])
                        hmm_aa_freq['E'][1].append(i_freq[0][3])
                        hmm_aa_freq['F'][1].append(i_freq[0][4])
                        hmm_aa_freq['G'][1].append(i_freq[0][5])
                        hmm_aa_freq['H'][1].append(i_freq[0][6])
                        hmm_aa_freq['I'][1].append(i_freq[0][7])
                        hmm_aa_freq['K'][1].append(i_freq[0][8])
                        hmm_aa_freq['L'][1].append(i_freq[0][9])
                        hmm_aa_freq['M'][1].append(i_freq[0][10])
                        hmm_aa_freq['N'][1].append(i_freq[0][11])
                        hmm_aa_freq['P'][1].append(i_freq[0][12])
                        hmm_aa_freq['Q'][1].append(i_freq[0][13])
                        hmm_aa_freq['R'][1].append(i_freq[0][14])
                        hmm_aa_freq['S'][1].append(i_freq[0][15])
                        hmm_aa_freq['T'][1].append(i_freq[0][16])
                        hmm_aa_freq['V'][1].append(i_freq[0][17])
                        hmm_aa_freq['W'][1].append(i_freq[0][18])
                        hmm_aa_freq['Y'][1].append(i_freq[0][19])
                        
                self.hmm_database[profile_name] = [hmm_aa_freq, profile_len]

            else:
                print "No hmm_profile found in the database!"
            
        # checking...
#        for k, v in self.hmm_database.items():
#            print k
#            for aa, f in v[0].items():
#                print '%s : match_poi_1 %s / insert_poi_1 %s' % (aa, f[0][0], f[1][0])
#                print '%s : match_poi_-1 %s / insert_poi_-1 %s' % (aa, f[0][-1], f[1][-1])


###################################################################################
##################################   FUNCTIONS   ##################################
###################################################################################

#Parse the hmmscan output file(s) and return the information needed for alignment-hmm sites mapping
def ali_parser(alignment_file):
    ali_obj = AlignIO.read(alignment_file, "fasta")
    for record in ali_obj:
        record.seq = str(record.seq).replace('X','-')  # in case of 'X', replacing it with '-'
    
    raw_array = np.array([list(rec) for rec in ali_obj], np.character, order="F")  # Convert alignment object to np_array object
    col_array = np.transpose(raw_array)  # align_columns.shape == (#col, #raw), i.e. (#site, #seq)

    # compute the columne scores (as the proxy for seq coverage percentage)
    col_scores = []
    for col in col_array:
        unique, counts = np.unique(col, return_counts=True)
        if '-' in unique:
            col_scores.append(1.0 - float(dict(zip(unique, counts))['-'])/float(len(col)))

    # obtain the seq with the highest coverage in the alignment
    array_score = np.copy(raw_array).astype(object)
    j = 0  # seq count
    while j < len(array_score):
        i = 0  # col count
        while i < len(col_scores):
            if array_score[j][i] == '-': array_score[j][i] = 0
            else:
                array_score[j][i] = col_scores[i]
            i += 1
        j += 1
    seq_score = np.sum(array_score, axis=1)
    
    #obtaining seed_seq related info
    seed_seq_array = raw_array[np.argmax(seq_score)]
    seed_seq = "".join(seed_seq_array)
    for record in ali_obj:
        if str(record.seq) == seed_seq:
            seed_id = str(record.id)

    seed_to_ali_mapping = [0 if x == '-' else 1 for x in seed_seq_array]  # constructing a mapping list of '-'== 0 & amino acids == 1,2,3,... seed_seq indicating position)
    i = 0  # track the alignment position
    j = 1  # track the seed_seq position
    #print len(seed_seq)
    while i < len(seed_seq):
        if seed_to_ali_mapping[i] == 1: 
            seed_to_ali_mapping[i] = j
            i += 1
            j += 1
        else:
            i += 1
    
    #print '%d \n%s %s\n%s' % (np.amax(seq_score), seed_id, seed_seq, seed_to_ali_mapping)  # check
    
    seed_seq_record = SeqRecord(Seq(seed_seq.replace("-","")), id=seed_id)  # remove the 'gap' so it won't affect hmm_aa_freq mapping in the Main Processes
    SeqIO.write(seed_seq_record, "temp_seed_seq.fasta", "fasta")
    
    # run hmmscan
    subprocess.check_output([hmmscan_exe, '-o', 'temp_hmmscan.out', args.database, 'temp_seed_seq.fasta'])

    return (seed_seq_record, seed_to_ali_mapping)
   
# parse the hmmscan output (TO-DO: need some sort of loop to handle multiple profile hits)
def hmmscan_parser(out_file, qurey_id):
    domain_info = {}  # store domain(s) info
    for QueryResult in SearchIO.parse(out_file,'hmmer3-text'):
        if qurey_id in str(QueryResult.id):
            #print QueryResult.description
            for hit in QueryResult.hits:  # TO-DO: need to control for multiple hits, comparing e-values?
                domain_info[hit.id] = []
                for HSPFragments in hit.hsps:  # multiple domains
                    query_hit_start = HSPFragments.query_start+1 # parsed start position appears to be 1 site short
                    query_hit_seq = HSPFragments.query.seq
                    query_hit_end = HSPFragments.query_end
                    domain_hit_start = HSPFragments.hit_start+1
                    domain_hit_end = HSPFragments.hit_end
                
                    domain_info[hit.id].append([query_hit_start, query_hit_end, query_hit_seq, domain_hit_start, domain_hit_end])
	
	return domain_info

#def pfam_parser(out_file, qurey_id):
#	
#	domain_info = {}  # store domain(s) info
#	
#	for QueryResult in SearchIO.parse(out_file,'hmmer3-text'):
#		if qurey_id in str(QueryResult.id):
#			#print QueryResult.description
#			for hit in QueryResult.hits:
#				for HSPFragments in hit.hsps:
#					query_hit_start = HSPFragments.query_start+1 # parsed start position appears to be 1 site short
#					query_hit_seq = HSPFragments.query.seq
#					query_hit_end = HSPFragments.query_end
#					domain_hit_start = HSPFragments.hit_start+1
#					domain_hit_end = HSPFragments.hit_end
#					
#				domain_info[hit.id] = [query_hit_start, query_hit_end, query_hit_seq, domain_hit_start, domain_hit_end]
#	
#	return domain_info

#Parse the the original alignments (with seed seq)
#def alignment_parser(alignment_file):
#	alignment = AlignIO.read(alignment_file, "fasta")
#	align_array = np.array([list(rec) for rec in alignment], np.character, order="F")  #Convert alignment object to np_array object
#	align_columns = np.transpose(align_array)  #align_columns.shape == (#col, #raw), i.e. (#site, #seq)
#	#Loop the columns to define the site mapping boundaries (using the full seed seq to map onto the trimmed alignment)
#	align_mapping = []
#	for col in align_columns:
#		#Check if only the seed seq in the column
#		col_check = col.tolist()
#		col_check = filter(lambda a: a != '-', col_check)
#		if len(col_check) == 1:  #Only seed seq present
#			align_mapping.append(0)
#		elif len(col_check) > 1:
#			align_mapping.append(1)
#
#	return align_mapping  #The structure looks like [0,0,0,1,1,1,0,0,1,...0] where 1 indicates mapped onto trimmed alignment

###################################################################################
###############################   MAIN PROCESSES   ################################
###################################################################################

seed_seq, seed_to_ali_mapping = ali_parser(args.ali)[0], ali_parser(args.ali)[1]
domains = hmmscan_parser('temp_hmmscan.out', seed_seq.id)
hmm_database = HMM(args.database, seed_seq.id).hmm_database  # Storing target hmm profile into a HMM Class which is structured as a dictionary
seed_seq_master = []  #The master list storing aa_freqs for the entire seed_seq

#Getting the aa_freqs from respective hmms, storing structure see HMM.py
for hmm_name, hits in sorted(domains.items(), key=operator.itemgetter([1][0])):  #The aa_seqs should be recorded in from the N-terminal to C-terminal of the seed seq
    if hmm_database[hmm_name] is not None:
        for hit_info in hits:    
            seed_seq_m_start = int(hit_info[0])  #Seed seq match start
            seed_seq_m_end = int(hit_info[1])  #Seed seq match end
            dm_seed_seq = str(hit_info[2])  #Domain matching seed seq
            dm_start = int(hit_info[3])  #Doamin match start
            dm_end = int(hit_info[4])  #Doamin match end
        
#           print hmm_name  #Check
#           print '%d %d' % (len(hmm_database[hmm_name][0]['A'][0]), len(hmm_database[hmm_name][0]['A'][1]))   #Check
	
            # Mapping the aa_freqs to seed_seq sites according to the hmmscan output
            hmmscan_seed_seq = list(dm_seed_seq)  #A list storing dm_seed_seq site-by-site

            j = 0  #Check for insert state so the pointer maps correctly to the domain position
            for i, seed_seq_site in enumerate(hmmscan_seed_seq):
                if seed_seq_site.isalpha() and seed_seq_site.isupper():
                    #Replace the site with a list of matching-state aa_freqs in the order of 'ARNDCQEGHILKMFPSTWYV'
                    #print i, seed_seq_site, i+dm_start-1-j, hmm_database[hmm_name][0]['A'][0][i+dm_start-1-j]  # Check
                    hmmscan_seed_seq[i] = [hmm_database[hmm_name][0]['A'][0][i+dm_start-1-j]]
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['R'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['N'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['D'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['C'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['Q'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['E'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['G'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['H'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['I'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['L'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['K'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['M'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['F'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['P'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['S'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['T'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['W'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['Y'][0][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['V'][0][i+dm_start-1-j])
                elif seed_seq_site.isalpha() and seed_seq_site.islower():
                    #Replace the site with a list of insert-state aa_freqs in the order of 'ARNDCQEGHILKMFPSTWYV'
                    j += 1
                    #print i, seed_seq_site  #Check
                    hmmscan_seed_seq[i] = [hmm_database[hmm_name][0]['A'][1][i+dm_start-1-j]]
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['R'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['N'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['D'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['C'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['Q'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['E'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['G'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['H'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['I'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['L'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['K'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['M'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['F'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['P'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['S'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['T'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['W'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['Y'][1][i+dm_start-1-j])
                    hmmscan_seed_seq[i].append(hmm_database[hmm_name][0]['V'][1][i+dm_start-1-j])
                elif seed_seq_site == '-':
                    continue  #Skipping the deletion state
                else:
                    print "Something is not right!!!"
                    break

            seed_seq_master += hmmscan_seed_seq

    else:
        print "something is not right..." # Check

for i, aa_freqs in enumerate(seed_seq_master):
    if isinstance(aa_freqs, list):
        seed_seq_master[i] = " ".join(aa_freqs)  #Construct the aa_freq string

seed_seq_master = filter(lambda x: x != "-", seed_seq_master)  #Discard deletion states "-"

# Mapping the aa_freqs to original alignment sites
while seed_seq_m_start <= seed_seq_m_end:
    for site_aa_freq in seed_seq_master:
        seed_to_ali_mapping = [site_aa_freq if x == seed_seq_m_start else x for x in seed_to_ali_mapping]
        seed_seq_m_start += 1

##################################################################################################
#########   IMPORTANT PROCESSES TO CORRECTLY MAP THE DOMAIN ONTO THE ORIGINAL ALIGNMENT  #########
##################################################################################################
#seed_to_alignment_mapping_check = seed_to_alignment_mapping[:]  #Duplicate the seed_to_alignment_mapping for checking the sites within the found domains but not in the trimmed alignment
##Loop the domain hits in the order of from the N-terminal to C-terminal of the seed seq
#for id, info in sorted(domains.items(), key=operator.itemgetter([1][0])):
#	seq_start = int(info[0])  #Seed seq match start
#	seq_end = int(info[1])  #Seed seq match end
#	#This loop maps aa_freqs onto the positions where the domains were found based on the seed seq
#	#So the seed_to_alignment_mapping should have aa_freqs regardless of if sites present in the trimmed alignment
#	for i in range(seq_start, seq_end+1):
#		seed_to_alignment_mapping[i-1] = seed_seq_master.pop(0)  #First position is 0
#
##This loop dicards the domain aa_freqs where not present in the trimmed alignment
#for i, aa_freq in enumerate(seed_to_alignment_mapping):
#	if seed_to_alignment_mapping_check[i] == 0:
#		seed_to_alignment_mapping[i] = 0

##################################################################################################
##################################################################################################

#Output the aa_freq
aa_freqs_output_file = open(args.ali+".freq", 'w')

#seed_to_alignment_mapping = filter(lambda a: a != 0, seed_to_alignment_mapping)  #Discard sites not in trimmed alignment: 0 

for aa_freq in seed_to_ali_mapping:
	if type(aa_freq) is str:
		aa_freqs_output_file.write(aa_freq+"\n")
	else:
		aa_freqs_output_file.write("NA\n")

aa_freqs_output_file.close()

