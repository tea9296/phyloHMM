#!/anaconda2/bin/python
import warnings, sys
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
from Bio import SearchIO
import operator, argparse, subprocess, re, mmap, math
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__copyright__ = """
  ################################################################################################
  #                                                                                              #
  # ExCAPE: ExtraCting site-wise Amino acid frequency Profiles for a multiple sEquence alignment #
  #                                                                                              #
  # Ding He, Copyright 2017                                                                      #
  # derrick.he@gmail.com                                                                         #
  #                                                                                              #
  ################################################################################################
"""

print __copyright__
###################################################################################
#################################   PARAMETERS   ##################################
###################################################################################

# Commandline arguments
parser = argparse.ArgumentParser(description='Command-line options for freq_extractor')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('--ali', help='input alignment')
required.add_argument('--database', help='hmm database')
required.add_argument('--hmmer_path', help='HMMER executable directory')
required.add_argument('--gb_path', help='Gblocks executable directory')
optional.add_argument('--iqtree', action='store_true', help='flag to output frequency file in IQ-TREE supported format')
args = parser.parse_args()  # e.g. input alignment file: args.ali; output freq_file: args.out
# Demanding required input
if not (args.ali) or not (args.database) or not (args.hmmer_path) or not (args.gb_path): parser.error(' Please provide necessary input, "python '+sys.argv[0]+' -h" to see more info')

# 3rd-party program directories
hmmscan_exe = str(args.hmmer_path)+'hmmscan'
hmmfetch_exe = str(args.hmmer_path)+'hmmfetch'
gblocks_exe = str(args.gb_path)+'Gblocks'

# Empirical aa_freqs from LG model ordered as ARNDCQEGHILKMFPSTWYV
# Probably do not need it now...
#LG_aa_freqs = '2.34595 2.91338 3.49262 3.01205 4.45606 3.52632 2.83988 2.41146 3.78460 2.67538 2.34684 3.06526 3.47661 3.10027 3.17541 2.95790 2.93381 4.53976 3.46891 2.49357'

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
    full_ali_obj = AlignIO.read(alignment_file, "fasta")
    for record in full_ali_obj:
        record.seq = str(record.seq).replace('X','-')  # in case of 'X', replacing it with '-'
    
    raw_array = np.array([list(rec) for rec in full_ali_obj], np.character, order="F")  # Convert alignment object to np_array object
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
    for record in full_ali_obj:
        if str(record.seq) == seed_seq:
            seed_id = str(record.id)

    seed_to_ali_mapping = [0 if x == '-' else 1 for x in seed_seq_array]  # constructing a mapping list of '-'== 0 & amino acids == 1,2,3,... seed_seq indicating position in the original alignment)
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
        aa_freqs = map(lambda x: float(x), aa_freqs)
        aa_freqs = map(lambda x: math.exp(-x), aa_freqs)  # aa_freq unit concersion (need to check possible numerical problem: sum(aa_freq) == 1 for each site)
        aa_freqs = map(lambda x: str(x), aa_freqs)
        seed_seq_master[i] = " ".join(aa_freqs)  #Construct the aa_freq string

seed_seq_master = filter(lambda x: x != "-", seed_seq_master)  #Discard deletion states "-"

# Mapping the aa_freqs to original alignment sites
while seed_seq_m_start <= seed_seq_m_end:
    for site_aa_freq in seed_seq_master:
        seed_to_ali_mapping = [site_aa_freq if x == seed_seq_m_start else x for x in seed_to_ali_mapping]
        seed_seq_m_start += 1

#Output the aa_freq for full & hmm_trimmed alignments
full_ali_obj = AlignIO.read(args.ali, "fasta")  # have to read the original alignment again!

# Executing Gblocks processes
# Make seq names shorter so Gblocks is happy (<= 15 characters, I think...)
i = 0
full_ali_seq_des = []  # getting sequence descriptions for restoring to the original descriptions after Gblocks trimming
for record in full_ali_obj:
    record.id = str(i)
    full_ali_seq_des.append(record.description)
    record.description = ""
    i += 1

AlignIO.write(full_ali_obj, args.ali+".temp_gb", "fasta")

print "\n---------------------------- Gblocks output ----------------------------"
try:
    subprocess.check_output([gblocks_exe, args.ali+".temp_gb", "-b4=5", "-b5=h"])
except subprocess.CalledProcessError as e:
    print e.output
print "-------------------------------------------------------------------------\n\n"

gbtrim_ali_obj = AlignIO.read(args.ali+".temp_gb-gb", "fasta")
i = 0
for record in gbtrim_ali_obj:
    record.description = full_ali_seq_des[0]
    i += 1

full_ali_raw_array = np.array([list(rec) for rec in full_ali_obj], np.character, order="F")  # Convert alignment object to np_array object
full_ali_col_array = np.transpose(full_ali_raw_array)  # align_columns.shape == (#col, #raw), i.e. (#site, #seq)
gbtrim_ali_raw_array = np.array([list(rec) for rec in gbtrim_ali_obj], np.character, order="F")  # Convert alignment object to np_array object
gbtrim_ali_col_array = np.transpose(gbtrim_ali_raw_array)  # align_columns.shape == (#col, #raw), i.e. (#site, #seq)
#print gbtrim_ali_col_array[0]

print "Generating aa_freq files..."  # some checking status...

aa_freqs_output_file = open(args.ali+".freq", 'w')  # for full original alignment
aa_freqs_output_file_hmmtrim = open(args.ali+".hmmtrim.freq", 'w')  # for trimmed alignment where hmm profile was mapped
aa_freqs_output_file_gbtrim = open(args.ali+".gbtrim.freq", 'w')  # for Gblocks trimmed alignment

hmmtrim_ali_obj = full_ali_obj[:,0:0]  # initiating an empty alignment for hmmtrim
hmmtrim_iqtree_site_pos = 1  # first element of the freq_file for hmmtrim in iqtree format
gbtrim_iqtree_site_pos = 1  # first element of the freq_file for gbtrim in iqtree format
for pos, aa_freq in enumerate(seed_to_ali_mapping):
    if type(aa_freq) is str:  # hmm_profile mapped
        if args.iqtree:  # freq_file in iqtree format
            aa_freqs_output_file.write(str(pos+1)+" "+aa_freq+"\n")
            aa_freqs_output_file_hmmtrim.write(str(hmmtrim_iqtree_site_pos)+" "+aa_freq+"\n")
            hmmtrim_iqtree_site_pos += 1
            try:
                if full_ali_col_array[pos].tolist() == gbtrim_ali_col_array[gbtrim_iqtree_site_pos-1].tolist():                    
                    aa_freqs_output_file_gbtrim.write(str(gbtrim_iqtree_site_pos)+" "+aa_freq+"\n")
                    gbtrim_iqtree_site_pos += 1

            except IndexError:  # capture the IndexError when gbtrim positions are running out, so the rest of processes can continue
                continue

        else:
            aa_freqs_output_file.write(aa_freq+"\n")
            aa_freqs_output_file_hmmtrim.write(aa_freq+"\n")
            try:
                if full_ali_col_array[pos].tolist() == gbtrim_ali_col_array[gbtrim_iqtree_site_pos-1].tolist():
                    while (gbtrim_iqtree_site_pos <= len(gbtrim_ali_col_array)):
                        aa_freqs_output_file_gbtrim.write(aa_freq+"\n")
                        gbtrim_iqtree_site_pos += 1
            
            except IndexError:  # capture the IndexError when gbtrim positions are running out, so the rest of processes can continue
                continue

        hmmtrim_ali_obj = hmmtrim_ali_obj[:,:] + full_ali_obj[:,pos-1:pos]
        
    else:
        if args.iqtree:  # freq_file in iqtree format
            aa_freqs_output_file.write(str(pos+1)+" "+"NA\n")
        else:
            aa_freqs_output_file.write("NA\n")

for record in hmmtrim_ali_obj:
    record.description = full_ali_seq_des[0]
    i += 1

AlignIO.write(hmmtrim_ali_obj, args.ali+".hmmtrim", "fasta")
AlignIO.write(gbtrim_ali_obj, args.ali+".gbtrim", "fasta")

aa_freqs_output_file.close()
aa_freqs_output_file_hmmtrim.close()
aa_freqs_output_file_gbtrim.close()
if args.iqtree:
    print '\nFull alignment:\t\t%s\nHMM_profile mapped alignment:\t%s\nGblock trimmed alignment:\t%s' % (args.ali+".freq", args.ali+".hmmtrim.freq", args.ali+".gbtrim.freq")
else:
    print '\nFull alignment:\t%s\nHMM_profile mapped alignment:\t%s' % (args.ali+".freq", args.ali+".hmmtrim.freq")
    
print "\nAll processes are completed."
