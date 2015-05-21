#!/usr/bin/python

########################################################################
# 27 April 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools

def get_tss(g, g_list=None):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	if g == "mm10":
		genome="mmusculus_gene_ensembl"
	elif g == "hg19":
		genome = "hsapiens_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "description", "gene_biotype"]), 
		filters="ensembl_gene_id", mart=ensembl)
	gene = list(C1BM.rx(True,1))
	chr1 = list(C1BM.rx(True,2))
	tss = list(C1BM.rx(True,3))
	end = list(C1BM.rx(True,4))
	st = list(C1BM.rx(True,5))
	name = list(C1BM.rx(True,6))
	des = list(C1BM.rx(True,7))
	bio = list(C1BM.rx(True,8))
	data = {}
	for index, g in enumerate(gene):
		data[g] = (chr1[index], tss[index], end[index], st[index], name[index], des[index], bio[index])
	return data

def 
def main():
	parser = argparse.ArgumentParser(description='Takes BED files and intersect them with regions, uses TSS regions by default\n')
	parser.add_argument('-i', '--input', help='BED file', required=True)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include mm10/hg19', required=True)
	parser.add_argument('-o', '--output', help='Output name', required=True)
	parser.add_argument('-r', '--regions', help='Use specified regions instead of TSS in bed file format', required=False)
	parser.add_argument('-f', '--regions', help='Use specified regions instead of TSS in bed file format', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

