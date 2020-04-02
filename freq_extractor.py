#!/anaconda2/bin/python
import warnings, sys, os
import datetime
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
from Bio import SearchIO
import operator, argparse, subprocess, re, mmap, math
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from shutil import copyfile
import collections

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
#count time
starttime=datetime.datetime.now()
# Global variables


print __copyright__
###################################################################################
#################################   PARAMETERS   ##################################
###################################################################################

# to assure path exist

def assure_path_exists(path):
        dir = os.path.dirname(path)
        if not os.path.exists(dir):
                os.makedirs(dir)

def assure_slash_exists(path):
	
	if path[-1]!='/' :
		
		path+="/"
	return path	
# Commandline arguments
parser = argparse.ArgumentParser(description='Command-line options for freq_extractor')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('--ali', help='input alignment')
required.add_argument('--output_path', help='output files directory')
optional.add_argument('--database', help='hmm database file', default='example_files/database.hmm')
optional.add_argument('--file_name', help='output file name', default='dEfAuLt')
optional.add_argument('--hmmer_path', help='HMMER executable directory', default='bin/')
optional.add_argument('--gb_path', help='Gblocks executable directory', default='bin/')
optional.add_argument('--iq_path', help='IQ-TREE executable directory', default='bin/')
optional.add_argument('--seed_seq', help='--mode consensus (default) : a consensus sequence found by hmmemit with--symfra 0\
    --mode representative : a representative single sequence with highest pairwise similarity', default='consensus')
#optional.add_argument('--iqtree', action='store_true', help='flag to output frequency file in IQ-TREE supported format')
args = parser.parse_args()  # e.g. input alignment file: args.ali; output freq_file: args.out
# Demanding required input
if not (args.ali) or not (args.output_path) : parser.error(' Please provide necessary input, "python '+sys.argv[0]+' -h" to see more info')

# 3rd-party program directories
hmmscan_exe = str(args.hmmer_path)+'hmmscan'
hmmfetch_exe = str(args.hmmer_path)+'hmmfetch'
hmmbuild_exe = str(args.hmmer_path)+'hmmbuild'
hmmemit_exe = str(args.hmmer_path)+ 'hmmemit'
gblocks_exe = str(args.gb_path)+'Gblocks'
iqtree_exe=str(args.iq_path)+'iqtree'
if args.file_name=='dEfAuLt':                   #change default output file into input alignment name
  args.file_name=((args.ali).split('/'))[-1]

args.output_path=assure_slash_exists(args.output_path)
args.hmmer_path=assure_slash_exists(args.hmmer_path)
args.gb_path=assure_slash_exists(args.gb_path)
args.iq_path=assure_slash_exists(args.iq_path)
assure_path_exists(args.output_path)

print hmmscan_exe
print hmmfetch_exe
print gblocks_exe

# Empirical aa_freqs from LG model ordered as ARNDCQEGHILKMFPSTWYV
# Probably do not need it now...
#LG_aa_freqs = '2.34595 2.91338 3.49262 3.01205 4.45606 3.52632 2.83988 2.41146 3.78460 2.67538 2.34684 3.06526 3.47661 3.10027 3.17541 2.95790 2.93381 4.53976 3.46891 2.49357'

LG_MODEL="0.079066 0.055941 0.041977 0.053052 0.012937 0.040767 0.071586 0.057337 0.022355 0.062157 0.099081 0.064600 0.022951 0.042302 0.044040 0.061197 0.053287 0.012066 0.034155 0.069147"

###################################################################################
####################################   CLASS   ####################################
###################################################################################   
def HMM( hmm_file, seed_seq_id,hmm_database):
    # creating a HMM structure for the target hmm profile(s) -> {profile_name_1:[{hmm_aa_freq_1},profile_length_1], profile_name_2:[{hmm_aa_freq_2},profile_length_2],...}
    # self.hmm_database = {}
    match_freq_lines_pattern = re.compile(r"\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\d+\s+[A-Za-z].+")
    insert_freq_lines_pattern = re.compile(r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")

        
        # separate hmm_profile into blocks in the hmm_file, and parse profile information: NAME, LENG and aa_freq
    with open(hmm_file, 'r+') as hmmF:
        mmap_hmmF = mmap.mmap(hmmF.fileno(), 0)
        #target_hmm_profile_match = re.search('NAME\s+(.*)', mmap_hmmF)
        target_hmm_profile_match = re.search( seed_seq_id, mmap_hmmF)
        if target_hmm_profile_match:
            #profile_name = target_hmm_profile_match.group(1).strip('\r')
            #print 'Found target hmm profile: %s' % (profile_name)
            profile_name=seed_seq_id
            print 'Found target hmm profile: %s' % (profile_name)
            # running hmmfetch to retrieve the target hmm profile
             #target_hmm_lines = subprocess.check_output([hmmfetch_exe, args.database, profile_name])
            target_hmm_lines = subprocess.check_output([hmmfetch_exe, args.database, seed_seq_id])
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
                        
                
            hmm_database[profile_name] = [hmm_aa_freq, profile_len]

        else:
            print "No hmm_profile found in the database!"
            
            # checking...
#            for k, v in self.hmm_database.items():
#                print k
#               for aa, f in v[0].items():
#                    print '%s : match_poi_1 %s / insert_poi_1 %s' % (aa, f[0][0], f[1][0])
#                    print '%s : match_poi_-1 %s / insert_poi_-1 %s' % (aa, f[0][-1], f[1][-1])
    return hmm_database

###################################################################################
##################################   FUNCTIONS   ##################################
###################################################################################

#Parse the hmmscan output file(s) and return the information needed for alignment-hmm sites mapping
def ali_parser(alignment_file):
    

    # optional parameter  use score or hmmemit to choose seed sequence    #args.ali
    if args.seed_seq == 'consensus':
        
        subprocess.check_output([hmmbuild_exe, '--symfrac', '0', args.output_path+args.file_name+'.hmm', args.ali])    
        subprocess.check_output([hmmemit_exe, '-c', '-o', args.output_path+'temp_seed_seq.fasta', args.output_path+args.file_name+'.hmm'])
        seed_seq_record = AlignIO.read(args.output_path+'temp_seed_seq.fasta', 'fasta')
        seed_seq_record = str(seed_seq_record[0].seq)
        seed_seq_record = SeqRecord(Seq(seed_seq_record), id=args.file_name+'-consensus')

        seed_to_ali_mapping = [x+1 for x in range(len(seed_seq_record))]


    elif args.seed_seq == 'representative':

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
                if array_score[j][i] == '-': array_score[j][i] = 0.0
                else:
                    array_score[j][i] = float(col_scores[i])
            
                i += 1
            j += 1

        seq_score=np.zeros( j )
        #np_sum
        j = 0  # seq count
        while j < len(array_score):
            i = 0  # col count
            while i < len(col_scores):
                seq_score[j]=seq_score[j]+array_score[j][i]
            
                i += 1
            j += 1


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
       
        SeqIO.write(seed_seq_record, args.output_path+"temp_seed_seq.fasta", "fasta")
    
    # run hmmscan
    subprocess.check_output([hmmscan_exe, '-o', args.output_path+'temp_hmmscan.out', '-E', '1e-20' ,'--domE','1e-20', args.database, args.output_path+'temp_seed_seq.fasta'])
    return (seed_seq_record, seed_to_ali_mapping)
   
# parse the hmmscan output (TO-DO: need some sort of loop to handle multiple profile hits)
def hmmscan_parser(out_file, qurey_id):
    domain_info = {}  # store domain(s) info
    for QueryResult in SearchIO.parse(out_file,'hmmer3-text'):
       
        #if qurey_id in str(QueryResult.id):   
        for hit in QueryResult.hits:  # TO-DO: need to control for multiple hits, comparing e-values?
            domain_info[hit.id] = []

            for HSPFragments in hit.hsps:  # multiple domains

                domain_hit_score =HSPFragments.bitscore
                query_hit_start = HSPFragments.query_start+1 # parsed start position appears to be 1 site short
                query_hit_seq = HSPFragments.query.seq
                query_hit_end = HSPFragments.query_end
                domain_hit_start = HSPFragments.hit_start+1
                domain_hit_end = HSPFragments.hit_end
                domain_info[hit.id].append([domain_hit_score, query_hit_end, query_hit_seq, domain_hit_start, domain_hit_end,query_hit_start])
                    
	
	return domain_info

#add default frequency to the non-matching positions ( use LG model )
def add_default_model(alignment,default_model):   
    
    for pos, freq in enumerate(alignment):
        if type(freq) is str:
            continue

        else:
            freq=str(default_model) 
            alignment[pos]=str(default_model)
    return alignment

def check_domain_position(start,end,lists):
    flag=0
    if start==0 and end==0:
        return False
    while start<=end:
        if lists[start]==1:
            flag=1
            return False

        start+=1
    return True

def set_domain_position(start, end, lists):

    while start<=end:
        lists[start]=1
        start+=1
    return lists
###################################################################################
###############################   MAIN PROCESSES   ################################
###################################################################################

seed_seq, seed_to_ali_mapping = ali_parser(args.ali)[0], ali_parser(args.ali)[1]
domains = hmmscan_parser(args.output_path+'temp_hmmscan.out', seed_seq.id)
hmm_database = {}

for key, value in domains.iteritems():
#hmm_database = HMM(args.database, seed_seq.id).hmm_database  # Storing target hmm profile into a HMM Class which is structured as a dictionary
    
    hmm_database = HMM(args.database, key,hmm_database )

alimentObj = AlignIO.read(args.ali, "fasta")
aligmentObjLen=len(alimentObj[0])
local_seed_seq_master = [[]]*(aligmentObjLen+int(0.1*aligmentObjLen))
max_seed_seq_master = [[]]*aligmentObjLen

seed_seq_master = [[]]*(aligmentObjLen+int(0.1*aligmentObjLen)) #The master list storing aa_freqs for the entire seed_seq
domain_position_check=[0 for x in range(aligmentObjLen)]

##### LET HMM_database order as score
'''
domains_order = {}
dom_key=(domains.keys())


for i in dom_key[::-1]:
    domains_order[i]=domains[i]
    print i
#print domains_order
'''
higheststart=0
highestend=0

seq_start_end_list = []
#Getting the aa_freqs from respective hmms, storing structure see HMM.py
for hmm_name, hits in sorted(domains.items(), key=operator.itemgetter([1][0]),reverse=True):  #The aa_seqs should be recorded in from the N-terminal to C-terminal of the seed seq

       
    highestscore=0
    if hmm_database[hmm_name] is not None:
        for hit_info in hits:    
            seed_seq_m_start = int(hit_info[5])  #Seed seq match start
            score=float(hit_info[0])
            seq_count=seed_seq_m_start
            

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
                    pass
                    #continue  #Skipping the deletion state
                else:
                    print "Something is not right!!!"
                    break

                                
                if  local_seed_seq_master[seq_count]:
                    pass
                else:

                    local_seed_seq_master[seq_count]= hmmscan_seed_seq[i]
                seq_count+=1
            #seed_seq_master += hmmscan_seed_seq
            if score>highestscore:
                highestscore=score
                max_seed_seq_master=local_seed_seq_master
                higheststart= seed_seq_m_start
                highestend=seed_seq_m_end



    if check_domain_position(higheststart, highestend,domain_position_check):

        domain_position_check=set_domain_position(higheststart, highestend,domain_position_check)
        seq_start_end_list.append(higheststart)
        seq_start_end_list.append(highestend)
        while higheststart<=highestend:
            seed_seq_master[higheststart]=max_seed_seq_master[higheststart]
            higheststart+=1

        


    else:
        print "something is not right..." # Check

for i, aa_freqs in enumerate(seed_seq_master):
    #if isinstance(aa_freqs, list):
     if  seed_seq_master[i]=='-':
         continue
     elif seed_seq_master[i] :
         aa_freqs = map(lambda x: float(x), aa_freqs)
         aa_freqs = map(lambda x: math.exp(-x),aa_freqs)  # aa_freq unit concersion (need to check possible numerical problem: sum(aa_freq) == 1 for each site)
         aa_freqs = map(lambda x: str(x), aa_freqs)
         seed_seq_master[i] = " ".join(aa_freqs)  # Construct the aa_freq string
     else:

        pass
seed_seq_master = filter(lambda x: x != "-", seed_seq_master)  #Discard deletion states "-"
#print seed_seq_m_start, seed_seq_m_end
# Mapping the aa_freqs to original alignment sites
for i in range(0,len(seq_start_end_list),2):
    seed_seq_m_start=seq_start_end_list[i]

    seed_seq_m_end=seq_start_end_list[i+1]
    while seed_seq_m_start <= seed_seq_m_end:
        seed_to_ali_mapping = [seed_seq_master[seed_seq_m_start] if x == seed_seq_m_start else x for x in seed_to_ali_mapping]
        seed_seq_m_start += 1
  

'''
while seed_seq_m_start <= seed_seq_m_end:
    for site_aa_freq in seed_seq_master:
        seed_to_ali_mapping = [site_aa_freq if x == seed_seq_m_start else x for x in seed_to_ali_mapping]
        seed_seq_m_start += 1
'''
#print seed_to_ali_mapping
#Output the aa_freq for full & hmm_trimmed alignments


full_ali_obj = AlignIO.read(args.ali, "fasta")  # have to read the original alignment again!

# Executing Gblocks processes
# Make seq names shorter so Gblocks is happy (<= 15 characters, I think...)
i = 0
full_ali_seq_des = []  # getting sequence descriptions for restoring to the original descriptions after Gblocks trimming
alignment_name = []
for record in full_ali_obj:
    alignment_name.append(record.id)
    record.id = str(i)
    #full_ali_seq_des.append(record.description)
    record.description = ""
    i += 1

AlignIO.write(full_ali_obj, args.output_path+args.file_name+".temp_gb", "fasta")

print "\n---------------------------- Gblocks output ----------------------------"

try:
    subprocess.check_output([gblocks_exe, args.output_path+args.file_name+".temp_gb", "-b4=5", "-b5=h", "-p=t"])
except subprocess.CalledProcessError as e:
    print e.output
print "-------------------------------------------------------------------------\n\n"


gbtrim_ali_obj = AlignIO.read(args.output_path+args.file_name+".temp_gb-gb", "fasta")

i = 0
for record in gbtrim_ali_obj:
    record.id = alignment_name[i]
    #record.description = full_ali_seq_des[i]
    #gbtrim_len=len(record)
    #gbt_ali=record
    record.description = ""
    i += 1

for r in full_ali_obj:
    full_ali=r
full_ali_raw_array = np.array([list(rec) for rec in full_ali_obj], np.character, order="F")  # Convert alignment object to np_array object
full_ali_col_array = np.transpose(full_ali_raw_array)  # align_columns.shape == (#col, #raw), i.e. (#site, #seq)
gbtrim_ali_raw_array = np.array([list(rec) for rec in gbtrim_ali_obj], np.character, order="F")  # Convert alignment object to np_array object
gbtrim_ali_col_array = np.transpose(gbtrim_ali_raw_array)  # align_columns.shape == (#col, #raw), i.e. (#site, #seq)
#print gbtrim_ali_col_array[0]

print "Generating aa_freq files..."  # some checking status...

aa_freqs_output_file = open(args.output_path+args.file_name+".freq", 'w')  # for full original alignment
aa_freqs_output_file_hmmtrim = open(args.output_path+args.file_name+".hmmtrim.freq", 'w')  # for trimmed alignment where hmm profile was mapped
aa_freqs_output_file_gbtrim = open(args.output_path+args.file_name+".gbtrim.freq", 'w')  # for Gblocks trimmed alignment
#aa_test_full_array=open(args.output_path+args.file_name+"full.txt", 'w')
#aa_test_gb_array=open(args.output_path+args.file_name+"gb.txt", 'w')


hmmtrim_ali_obj = full_ali_obj[:,0:0]  # initiating an empty alignment for hmmtrim
hmmtrim_iqtree_site_pos = 1  # first element of the freq_file for hmmtrim in iqtree format
gbtrim_iqtree_site_pos = 1  # first element of the freq_file for gbtrim in iqtree format
match_count=0

#Write hmmtrim.freq
for pos, aa_freq in enumerate(seed_to_ali_mapping):
    if type(aa_freq) is str:
        aa_freqs_output_file_hmmtrim.write(str(hmmtrim_iqtree_site_pos)+" "+aa_freq+"\n")
        
        hmmtrim_iqtree_site_pos += 1
        hmmtrim_ali_obj = hmmtrim_ali_obj[:,:] + full_ali_obj[:,pos-1:pos]


#add LG_MODEL
seed_to_ali_mapping=add_default_model(seed_to_ali_mapping,LG_MODEL)


#Write .freq and gbtrim.freq

for pos, aa_freq in enumerate(seed_to_ali_mapping):
    if type(aa_freq) is str:  # hmm_profile mapped
            aa_freqs_output_file.write(str(pos+1)+" "+aa_freq+"\n")
            
            #aa_freqs_output_file_hmmtrim.write(str(hmmtrim_iqtree_site_pos)+" "+aa_freq+"\n")
            #hmmtrim_iqtree_site_pos += 1
            #hmmtrim_ali_obj = hmmtrim_ali_obj[:,:] + full_ali_obj[:,pos-1:pos]
            try:
                

                gbtrim_iqtree_site_pos=1


                for gaca in range(match_count,len(gbtrim_ali_col_array)):            
        
                    if full_ali_col_array[pos].tolist() == gbtrim_ali_col_array[gaca].tolist():
                        match_count+=1
                        #print "\n"+str(pos)+"\t"+str(gbtrim_iqtree_site_pos)
                        #print full_ali_col_array[pos].tolist()
                        aa_freqs_output_file_gbtrim.write(str(match_count)+" "+str(aa_freq)+"\n")
                        
                    
                        break
                    gbtrim_iqtree_site_pos=gbtrim_iqtree_site_pos+1    



            except IndexError:  # capture the IndexError when gbtrim positions are running out, so the rest of processes can continue
               continue

       
            
i=0
for record in hmmtrim_ali_obj:
    record.id = alignment_name[i]
    
    i += 1

AlignIO.write(hmmtrim_ali_obj, args.output_path+args.file_name+".hmmtrim", "fasta")
AlignIO.write(gbtrim_ali_obj, args.output_path+args.file_name+".gbtrim", "fasta")

aa_freqs_output_file.close()
aa_freqs_output_file_hmmtrim.close()
aa_freqs_output_file_gbtrim.close()

print '\nFull alignment:\t\t%s\nHMM_profile mapped alignment:\t%s\nGblock trimmed alignment:\t%s' % (args.output_path+args.file_name+".freq", args.output_path+args.file_name+".hmmtrim.freq", args.output_path+args.file_name+".gbtrim.freq")
























endtime=datetime.datetime.now()
print "executive time:",(endtime-starttime).seconds ,"seconds"
