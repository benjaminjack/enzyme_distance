import protein_toolbox_helper as ph
import numpy as np
import re, os, sys

'''
Last Edited By: Eleisha Jackson on February 9, 2015
This is a script that extracts the amino acid rates from the individual Rate4Site file for each protein.
However, in this script that rates are not mapped back to the alignment
'''

def get_rates(rate_file):
	'''
	This is a function that extracts and maps the evolutionary rates to the pdb sequence
	
	Args:

		map_file: The file that maps the pdb sequence to the alignment	
		rate_file: The rate4site file
	Returns:
		sites: The sites in the pdb that have been mapped to the alignment	
		aligned_rates: The rates mapped back to sites in the pdb	
	'''
			
	align_sites = []
	rates = []
	amino_acids = []
	ref_seq_array  = []

	input = open(rate_file)
	data = input.readlines()
	data = data[13:]
	data.pop()
	data.pop()

	for line in data: #Get the evolutionary rates from the rate4site file
		site = line[0:4].strip()
		amino = line[5:10].strip()
		ref_seq_array.append(amino)
		rate = line[10:19].strip()
		align_sites.append(site)
		amino_acids.append(amino)
		rates.append(rate)	

	aligned_rates = []
	counter = 0
	# print ref_seq_array
	return ref_seq_array, rates

def extract_sites(map_file):
	'''
	This is a function extracts the site information from the mapped file
	
	Args:
		pdb: The pdb name
		map_file: The file that maps the pdb sequence to the alignment
		rate_file: the rate4site file
	Returns:
		sites: In this case just the pdb_aa 		
	'''

	site_data = np.genfromtxt(map_file, delimiter = "\t", dtype = "string", skip_header = 1, usecols = (0)) #Read in the data
	for s in site_data:
		sites.append(s.strip()) #String off the new lines and append the site to the sites list
	return sites

def main():
	if len( sys.argv ) != 3:
		print '''      Usage:      '''
		print "     ", sys.argv[0], "<input rate4site file>", "<output file>"
								
	else:
		rate_file = sys.argv[1] #rate4site file
		outfile = sys.argv[2] #outfile name
		
		out = open(outfile, "w")
		out.write('"pdb_aa","rate"\n')

		try:
			sites, rates = get_rates(rate_file) #Extract the rates
		except IOError:
			print "The File: " + pdb + " is not there."

		if (len(sites) != len(rates)):
			print "Sites not mapped correctly to Rates!"
			print "Length of Sites: ", len(sites)
			print "Length of Rates: ", len(rates)

		j = 0
		while ( j < len(sites)): #For each site
			out.write(str(sites[j]) + ',' + str(rates[j]) + "\n") #Write out the rate information for that site
			# print str(rates[j])
			j = j + 1
	
		out.close()		

if __name__ == "__main__":
	main()
