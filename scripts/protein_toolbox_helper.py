#!/usr/bin/python
import re, os, math, string, subprocess
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
from Bio import PDB

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

_hydrogen=re.compile("[123 ]*H.*") 

class ChainSelector(object): 
	"""
	Adapted from the extract() function of Biopython 
	"""

	def __init__(self, chain_id, start, end, model_id=0): 
		self.chain_id=chain_id 
		self.start=start 
		self.end=end 
		self.model_id=0 

	def accept_model(self, model): 
		# model - only keep model 0 
		if model.get_id()==self.model_id: 
			return 1 
		return 0 

	def accept_chain(self, chain): 
		if chain.get_id()==self.chain_id: 
			return 1 
		return 0 
		
	def accept_residue(self, residue): 
		# residue - between start and end 
		hetatm_flag, resseq, icode=residue.get_id() 
		if hetatm_flag!=" ": 
			# skip HETATMS 
			return 0 
		if self.start<=resseq<=self.end: 
			return 1 
			return 0 

	def accept_atom(self, atom): 
		# atoms - get rid of hydrogens 
		name=atom.get_id() 
		if _hydrogen.match(name): 
			return 0 
		else: 
			return 1 

def extract(structure, chain_id, start, end, filename): 
	""" 
	Write out selected residues from a particular given pdb structure to filename. 
	
	Args: 
		structure: The four letter pdb code for the protein (this pdb must currently exist!)
		chain_id: The chain in the pdb you are trying to extract
		start: The beginning residue that you want to start your extraction from
		end:
		filename: The filename (ex. my_protein.pdb) of the pdb that you will extract the selected residues to  	
	Returns:
		A pdb file that contains only the selected residues from the original pdb structure given as an argument
		
	Example of usage
		p=PDBParser()
		s=p.get_structure('X', '1RUZ.pdb')
		extract(s, 'H', 40, 60, 'extracted.pdb')
 
	""" 
	sel=ChainSelector(chain_id, start, end) 
	io=PDBIO() 
	io.set_structure(structure) 
	io.save(filename, sel) 


def get_pdb_info(data_file):
	pdbs = []
	chains = []
	pdb_file = open(data_file, "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list		
	return pdbs, chains
	
def get_sequences(file):
	"""
	This function takes a file with a bunch of sequences in fasta format and returns a list of the sequences. 
	
	Args: 
		file: The name of the fasta file of the given sequences
		
	Returns:
		Returns TWO lists
		It returns: all_sequences, headers
		all_sequences is a list of every sequence in the fasta file
		headers is a list of every header for each sequence in the fasta file
		The headers and the sequences are in corresponding elements in the two lists 
		(so the first element in the header list is the header for the first sequence in the sequence list.)				
	"""
		
	all_sequences = []  
	natural_sequences = []
	designed_sequences = []
	headers = []

	file_data = open(file, "r")
	seq_data = file_data.readlines()
	file_data.close()
	#print seq_data
	#This block of code basically just removes all of the headers from the array of sequence information.
	#It creates all_sequences, which is a list of all the sequences from the natural alignment file.
	string  = ''
	finished_sequence = ''
	for sequence in seq_data:
		if sequence[0] == '>': #If it is a header append the last sequence that was processed
			headers.append(sequence.rstrip())
			if(string != ''):
				finished_sequence = string #Just in case the sequence was translated from Stockholm format
				all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
				string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
		else:       
			string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
	all_sequences.append(string)
	num_sequences = len(all_sequences)
	return all_sequences, headers   

def get_sequence_array(seq):
	"""
	Turns a string representing a string into a list of characters.
	
	Args:
		seq: A string representing a sequence (ex. nucleotide, protein)
		
	Returns:
		seq_array: A list of single letter elements representing the sequence (ex. each amino acid in a protein is an element in the list)
	
	"""

	seq_array = []
	seq_length = len(seq)
	count = 0
	while count < seq_length:
		aa = seq[count]
		seq_array.append(aa)
		count = count + 1
	return seq_array

def array_to_seq(seq_array):
	"""
	Takes an array that represents a sequence and returns a string representation
	
	Args:
		seq_array: An array where each element is an character in a sequence
	
	Returns:
		seq: A string representation of the sequence
		
	"""

	seq = "" #Create an empty string
	for i in xrange(0, len(seq_array)): #For each element in the seq array
		seq = seq + seq_array[i]	 #Append the character to the empty string
	return seq

def get_info_from_pdb(pdb_file, chain):
	"""
	Extracts the sequences from of Protein Database File (PDB) file
	
	Args:
		pdb_file = A string that is the path to the file (just the name if in the same directory)
		chain = The chain in the pdb for the sequence you want to extract
			
	Returns:
		seq: A string representation of the sequence in one letter code
			
	"""
	p=PDB.PDBParser()
	s=p.get_structure('X', pdb_file)				
	ref_struct = s[0]
	ref_chain = ref_struct[chain]
	ref_residues = []
	ref_res_nums = []
	for res in ref_chain:
		#print res
		ref_residues.append( res.resname )
		ref_res_nums.append(res.id[1])
	seq_array = convert_to_one_letter_code(ref_residues)
	seq = array_to_seq(seq_array)
	return seq, ref_res_nums

def convert_to_one_letter_code(residue_list):
	"""
	Takes a list of amino acids in letter-code and changes it to one-letter code.
	
	Args:
		residue_list: A list with an amino acid sequence represented in three-letter code
	
	Returns:
		one_letter_residue_list: A list whose elements are the amino acids in one-letter code
	
	"""
	one_letter_residue_list = []
	
	for res in residue_list:
		abbrev = resdict[res]
		one_letter_residue_list.append(abbrev)
	return one_letter_residue_list	

def get_sequences(file):
	"""
	This function takes a file with a bunch of sequences in fasta format and returns a list of the sequences. 
	
	Args: 
		file: The name of the fasta file of the given sequences
		
	Returns:
		Returns TWO lists
		It returns: all_sequences, headers
		all_sequences is a list of every sequence in the fasta file
		headers is a list of every header for each sequence in the fasta file
		The headers and the sequences are in corresponding elements in the two lists 
		(so the first element in the header list is the header for the first sequence in the sequence list.)				
	"""
		
	all_sequences = []  
	natural_sequences = []
	designed_sequences = []
	headers = []

	file_data = open(file, "r")
	seq_data = file_data.readlines()
	file_data.close()
	
	#print seq_data
	#This block of code basically just removes all of the headers from the array of sequence information.
	#It creates all_sequences, which is a list of all the sequences from the natural alignment file.
	string  = ''
	finished_sequence = ''
	for sequence in seq_data:
		if sequence[0] == '>': #If it is a header append the last sequence that was processed
			headers.append(sequence.rstrip())
			if(string != ''):
				finished_sequence = string #Just in case the sequence was translated from Stockholm format
				all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
				string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
		else:       
			string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
	all_sequences.append(string)
	num_sequences = len(all_sequences)
	return all_sequences, headers   

def get_unaligned_seqs(fasta_file):
	"""
	This is a file that will open a fasta file and return the sequences sequences unaligned.
	
	Args:
		fasta_file: An alignment file in fasta format
		
	Returns:
		sequences: A list of sequences from the file
		headers: A list of the headers for the corresponding sequences	
	"""
	[sequences, headers] = get_sequences(fasta_file)
	new_seqs = []
	for seq in sequences:
		new_seq = ""
		seq_array = get_sequence_array(seq)
		for i in xrange(0, len(seq_array)):
			if(seq_array[i] != "-"):
				new_seq = new_seq + seq_array[i]
		new_seqs.append(new_seq)
		
	return new_seqs, headers

def write_seq_to_file(seqs, headers, out_fasta):
	"""
	Writes a group of sequences to a file
	
	Args:
		seqs: A list of sequences
		headers: A list of descriptors (often fasta headers > included)
		out_fasta: The name of the file that the sequences will be written to
		
	Returns:
		out_fasta: A string that is the filename that was given to the function
	
	"""
	file = open(out_fasta, "w")
	j = 0
	for seq in seqs:
		file.write(headers[j]+ "\n")
		file.write(seq + "\n")
		j = j+1
	file.close()
	return out_fasta

def align_seqs_mafft(seqs, headers, out_fasta):
	in_file = "temp_seqs_12345.txt"
	write_seq_to_file(seqs, headers, in_file)
	align_command = "mafft --auto " + in_file + " > " + out_fasta
	subprocess.call(align_command, shell = True)
	os.unlink(in_file)
