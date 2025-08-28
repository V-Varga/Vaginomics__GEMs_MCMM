#!/bin/python
# -*- coding: utf-8 -*-
"""

Title: assignFASTAheaders_GEMs.py
Date: 2025.08.11
Author: Vi Varga

Description:
	This program replaces adds a user-defined string to the start of a FASTA header.
		It is intended for use in GEM generation for MCMMs, where the NCBI 
		structure of the FASTA headers is important, but identical protein names
		in different organisms could cause issues.

List of functions:
	No functions are defined in this script.

List of standard and non-standard modules used:
	sys
	os
	re

Procedure:
	1. Loading required modules & assigning command line argument.
    2. Parsing the input FASTA file in order to extract headers and add the
		user-determined substring.
	3. Writing out the new FASTA file.

Known bugs and limitations:
	- There is no quality-checking integrated into the code.
	- There is no default substring - the user must provide one.

Usage
	./assignFASTAheaders_GEMs.py input_fasta user_string
	OR
	python assignFASTAheaders_GEMs.py input_fasta user_string
	
	Where user_string should be an identifying string for the FASTA file, without
		spaces or special characters. For eventual MCMM generation, it is 
		recommended to use a string identifying the species.

This script was written for Python 3.13.5, in Spyder 6.0.7. 

"""


# Part 1: Import modules & assign command line arguments

#import necessary modules
import sys #allows execution of script from command line
import os #allow access to computer files
import re #enables regex pattern matching


#load input and output files
input_fasta = sys.argv[1]
#input_fasta = "Data/gardnerella_vaginalis_GCF_001042655-1_prots.faa"
base = os.path.basename(input_fasta)
out_full = os.path.splitext(base)[0]
output_fasta = ".".join(input_fasta.split('.')[:-1]) + '_edit.fasta'

# determine the user-provided substring
user_string = str(sys.argv[2])


# Part 2: Assign the alphanumeric headers and write out results files

with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
	#open the input and output files
	for line in infile:
		#iterate through the input file line by line
		if line.startswith(">"):
			#identify the header lines and remove the end-line character
			header = line.strip()
			#remove the ">" character at the start of the line
			#this enables easier manipulation of the FASTA header
			header = re.sub(">", "", header)
			#create the new header
			new_header = user_string + "_" + header
			#now print the new header to the outfile
			outfile.write(">" + new_header + "\n")
		else:
			#sequence lines are copied to the outfile without changes
			outfile.write(line)
