import sys, re

class HMM:
# for hmm file the aa order is ACDEFGHIKLMNPQRSTVWY   

	def __init__(self, hmm_file):
				
		# (match_freq, insert_freq)
		self.hmm_aa_freq = {'A': [[],[]], 'R': [[],[]], 'N': [[],[]], 'D': [[],[]], 'C': [[],[]], 'Q': [[],[]], \
							'E': [[],[]], 'G': [[],[]], 'H': [[],[]], 'I': [[],[]], 'L': [[],[]], 'K': [[],[]], \
							'M': [[],[]], 'F': [[],[]], 'P': [[],[]], 'S': [[],[]], 'T': [[],[]], 'W': [[],[]], \
							'Y': [[],[]], 'V': [[],[]]}

		with open(hmm_file) as self.hmm_file:
			#Start parsing some meta info
			for line in self.hmm_file:
				if "NAME" in line: self.name = re.split('\s+', line)[1]
				if "LENG" in line: self.length = int(re.split('\s+', line)[1])
				if line.startswith("  COMPO "):
					break  # Pause reading lines --> next line is the aa insert-frequency line

			self.hmm_file.next()  # skip the next unwanted line			
			self.hmm_file.next()  # skip the next unwanted line
			
			#Start parsing the match/insert aa_frequency lines
			match_freq_lines_pattern = re.compile(r"^\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\d+\s+[A-Za-z].+$")
			insert_freq_lines_pattern = re.compile(r"^(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+$")

			for line in self.hmm_file:
				match_freq_line = re.match(match_freq_lines_pattern, line.lstrip())
				insert_freq_line = re.match(insert_freq_lines_pattern, line.lstrip())
				if match_freq_line is not None:
					#print "Match: "+line
					self.hmm_aa_freq['A'][0].append(match_freq_line.group(1))
					self.hmm_aa_freq['C'][0].append(match_freq_line.group(2))
					self.hmm_aa_freq['D'][0].append(match_freq_line.group(3))
					self.hmm_aa_freq['E'][0].append(match_freq_line.group(4))
					self.hmm_aa_freq['F'][0].append(match_freq_line.group(5))
					self.hmm_aa_freq['G'][0].append(match_freq_line.group(6))
					self.hmm_aa_freq['H'][0].append(match_freq_line.group(7))
					self.hmm_aa_freq['I'][0].append(match_freq_line.group(8))
					self.hmm_aa_freq['K'][0].append(match_freq_line.group(9))
					self.hmm_aa_freq['L'][0].append(match_freq_line.group(10))
					self.hmm_aa_freq['M'][0].append(match_freq_line.group(11))
					self.hmm_aa_freq['N'][0].append(match_freq_line.group(12))
					self.hmm_aa_freq['P'][0].append(match_freq_line.group(13))
					self.hmm_aa_freq['Q'][0].append(match_freq_line.group(14))
					self.hmm_aa_freq['R'][0].append(match_freq_line.group(15))
					self.hmm_aa_freq['S'][0].append(match_freq_line.group(16))
					self.hmm_aa_freq['T'][0].append(match_freq_line.group(17))
					self.hmm_aa_freq['V'][0].append(match_freq_line.group(18))
					self.hmm_aa_freq['W'][0].append(match_freq_line.group(19))
					self.hmm_aa_freq['Y'][0].append(match_freq_line.group(20))
				elif insert_freq_line is not None:
					#print "Insert: "+line
					self.hmm_aa_freq['A'][1].append(insert_freq_line.group(1))
					self.hmm_aa_freq['C'][1].append(insert_freq_line.group(2))
					self.hmm_aa_freq['D'][1].append(insert_freq_line.group(3))
					self.hmm_aa_freq['E'][1].append(insert_freq_line.group(4))
					self.hmm_aa_freq['F'][1].append(insert_freq_line.group(5))
					self.hmm_aa_freq['G'][1].append(insert_freq_line.group(6))
					self.hmm_aa_freq['H'][1].append(insert_freq_line.group(7))
					self.hmm_aa_freq['I'][1].append(insert_freq_line.group(8))
					self.hmm_aa_freq['K'][1].append(insert_freq_line.group(9))
					self.hmm_aa_freq['L'][1].append(insert_freq_line.group(10))
					self.hmm_aa_freq['M'][1].append(insert_freq_line.group(11))
					self.hmm_aa_freq['N'][1].append(insert_freq_line.group(12))
					self.hmm_aa_freq['P'][1].append(insert_freq_line.group(13))
					self.hmm_aa_freq['Q'][1].append(insert_freq_line.group(14))
					self.hmm_aa_freq['R'][1].append(insert_freq_line.group(15))
					self.hmm_aa_freq['S'][1].append(insert_freq_line.group(16))
					self.hmm_aa_freq['T'][1].append(insert_freq_line.group(17))
					self.hmm_aa_freq['V'][1].append(insert_freq_line.group(18))
					self.hmm_aa_freq['W'][1].append(insert_freq_line.group(19))
					self.hmm_aa_freq['Y'][1].append(insert_freq_line.group(20))

				