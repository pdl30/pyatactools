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
import HTSeq
from multiprocessing import Pool, Manager
import argparse
import pysam
import numpy
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
import pyatactools
import pkg_resources
import collections
from qcmodule import mystat

def sam_size(bam):
	results=  {}
	size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bam) ])
	return size

def read_tss(ifile, positions, halfwinwidth, count, return_dict):
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	constant = 1000000/float(sam_size(ifile))
	bamfile = HTSeq.BAM_Reader(ifile)
	for almnt in bamfile:	
		if almnt.aligned:
			s = set()
			for step_iv, step_set in positions[ almnt.iv ].steps():
				s |= step_set
			for p in s:
				if p.strand == "+":
					start_in_window = almnt.iv.start - p.start + halfwinwidth 
					end_in_window   = almnt.iv.end - p.start + halfwinwidth 
				else:
					start_in_window = p.start + halfwinwidth - almnt.iv.end
					end_in_window   = p.start + halfwinwidth - almnt.iv.start
				start_in_window = max( start_in_window, 0 )
				end_in_window = min( end_in_window, 2*halfwinwidth )
				profile[ start_in_window : end_in_window ] += constant
	profile = profile/float(count) #Average over the number of TSS regions
	return_dict[ifile] = profile

def read_tss_pysam(bam, position_dict, halfwinwidth, count, return_dict):
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	constant = 1000000/float(sam_size(bam))
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	for chrom, tss, strand in position_dict.values():
		coverage = {}
		chrom_start = tss-halfwinwidth
		if chrom_start <0: chrom_start=0
		chrom_end = tss+ halfwinwidth
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
		coverage = {}
		for i in range(1, 2*halfwinwidth):
			coverage[i] = 0.0
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos
			if strand == "+":
				ref_pos = pileupcolumn.pos - chrom_start
			else:
				ref_pos = chrom_end - pileupcolumn.pos 
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read += constant
		#	print chrom, tss, chrom_start, chrom_end, strand, ref_pos, cover_read
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
	#	if strand == '-':
	#		tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
	for key in aggreagated_cvg:
		profile[key] = aggreagated_cvg[key]
	return_dict[bam] = profile

def read_tss_function(args):
	return read_tss_pysam(*args)
def read_gene_function(args):
	return genebody_coverage(*args)

def genebody_percentile(anno, gene_filter, mRNA_len_cut = 100):
	'''
	return percentile points of gene body
	mRNA length < mRNA_len_cut will be skipped
	'''
	g_percentiles = {}
	g_filter = []
	if gene_filter:
		with open(gene_filter) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				g_filter.append(word[0])
	for line in open(anno,'r'):
		if line.startswith('Ensembl'):continue  
		# Parse fields from gene tabls
		fields = line.split()
		if fields[1] == "MT": chrom = "chrM"
		elif fields[1] == "X": chrom = "chrX"
		elif fields[1] == "Y": chrom = "chrY"
		elif fields[1].isdigit(): chrom = "chr" + fields[1]
		else: 
			continue
		tx_start  = int( fields[2] )
		tx_end    = int( fields[3] )
		geneName      = fields[0]
		if fields[4] == "1":
			strand = "+"
		else:
			strand = "-"
		geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
		gene_all_base=[]
		if g_filter:
			if geneName in g_filter:
				gene_all_base.extend(range(tx_start+1,tx_end+1))		#1-based coordinates on genome
				if len(gene_all_base) < mRNA_len_cut:
					continue
				g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
		else:
			gene_all_base.extend(range(tx_start+1,tx_end+1))		#1-based coordinates on genome
			if len(gene_all_base) < mRNA_len_cut:
				continue
			g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
	return g_percentiles

def genebody_coverage(bam, position_list, return_dict):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	constant = 1000000/float(sam_size(bam))
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	gene_finished = 0
	for chrom, strand, positions in position_list.values():
		coverage = {}
		for i in positions:
			coverage[i] = 0.0
		chrom_start = positions[0]-1
		if chrom_start <0: chrom_start=0
		chrom_end = positions[-1]
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions:
				continue
			if pileupcolumn.n == 0:
				coverage[ref_pos] = 0
				continue				
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=constant
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		if strand == '-':
			tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
		gene_finished += 1
		
		if gene_finished % 100 == 0:
			print >>sys.stderr, "\t%d transcripts finished\r" % (gene_finished),
	tmp2 = numpy.zeros( 100, dtype='i' ) 
	for key in aggreagated_cvg:
		tmp2[key] = aggreagated_cvg[key]
	return_dict[bam] = tmp2

def read_tss_anno(anno):
	positions = {}

def plot_tss_profile(conditions, anno, halfwinwidth, gene_filter, threads, outname):
	halfwinwidth = int(halfwinwidth)
	tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
	c = 0
	g_filter = []
	if gene_filter:
		with open(gene_filter) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				g_filter.append(word[0])
	positions = {}
	with open(anno) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			start = int(word[2]) - halfwinwidth
			if start < 0:
				start = 0
			if word[1] == "MT": chrom = "chrM"
			elif word[1] == "X": chrom = "chrX"
			elif word[1] == "Y": chrom = "chrY"
			elif word[1].isdigit(): chrom = "chr" + word[1]
			if word[4] == "1": 
				strand = "+"
			else: 
				strand = "-"
			if g_filter:
				if word[0] in g_filter:
					positions[word[0]] = (chrom, int(word[2]), strand)
					p = HTSeq.GenomicInterval( chrom, int(word[2]), int(word[2]), strand )
					window = HTSeq.GenomicInterval( chrom, start, int(word[2]) + halfwinwidth, "." )
					tsspos[ window ] += p
					c +=  1
			else:
				positions[word[0]] = (chrom, int(word[2]), strand)
				p = HTSeq.GenomicInterval( chrom, int(word[2]), int(word[2]), strand )
				window = HTSeq.GenomicInterval( chrom, start, int(word[2]) + halfwinwidth, "." )
				tsspos[ window ] += p
				c +=  1

	manager = Manager()
	return_dict = manager.dict()
#	read_tss_pysam(list(conditions)[0], positions, halfwinwidth, c, return_dict)
	pool = Pool(threads)
	pool.map(read_tss_function, itertools.izip(list(conditions.keys()), itertools.repeat(positions), itertools.repeat(halfwinwidth), 
		itertools.repeat(c), itertools.repeat(return_dict)))
	pool.close()
	pool.join()	
	for key in return_dict.keys():
		pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), return_dict[key], label=conditions[key])  
	pyplot.legend(prop={'size':8})
	pyplot.savefig(outname+".pdf")

def plot_genebody_profile(conditions, anno, gene_filter, threads, outname, filters):
	if isinstance(gene_filter, dict):
		for key1 in gene_filter:
			positions = genebody_percentile(anno, key1, mRNA_len_cut = 100)
			manager = Manager()
			return_dict = manager.dict()
			pool = Pool(threads)
			pool.map(read_gene_function, itertools.izip(list(conditions.keys()), itertools.repeat(positions), itertools.repeat(return_dict)))
			pool.close()
			pool.join()	
			for key in return_dict.keys():
				pyplot.plot( numpy.arange( 0, 100 ), return_dict[key], label="{}_{}".format(conditions[key], gene_filter[key1]))
	else
		positions = genebody_percentile(anno, gene_filter, mRNA_len_cut = 100)
		manager = Manager()
		return_dict = manager.dict()
		pool = Pool(threads)
		pool.map(read_gene_function, itertools.izip(list(conditions.keys()), itertools.repeat(positions), itertools.repeat(return_dict)))
		pool.close()
		pool.join()	
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( 0, 100 ), return_dict[key], label=conditions[key])
	pyplot.legend(prop={'size':8})
	pyplot.savefig(outname+".pdf")

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def main():
	parser = argparse.ArgumentParser(description='Takes BED files and intersect them with regions, uses TSS regions by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	tss_parser = subparsers.add_parser('tss', help='TSS plotter')
	tss_parser.add_argument('-c', '--config', help='BAM as keys, if using BEDS, include [size] section', required=True)
	tss_parser.add_argument('-g', '--genome', help='Options are mm10/hg19', required=True)
	tss_parser.add_argument('-f', '--filter', help='Gene name per line, filters TSS regions', required=False)
	tss_parser.add_argument('-a', action='store_true', help='Use filters from Config', required=False)
	tss_parser.add_argument('-o', '--output', help='Output name of pdf file', required=True)
	tss_parser.add_argument('-w', '--width', help='Width of region, default=1000', default=1000, required=False)
	tss_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	gene_parser = subparsers.add_parser('gene', help='Genebody plotter')
	gene_parser.add_argument('-c', '--config', help='BAM as keys, if using BEDS, include [size] section', required=False)
	gene_parser.add_argument('-g', '--genome', help='Options are mm10/hg19', required=False)
	gene_parser.add_argument('-f', '--filter', help='Gene name per line, filters TSS regions', required=False)
	gene_parser.add_argument('-a', action='store_true', help='Use filters from Config', required=False)
	gene_parser.add_argument('-o', '--output', help='Output name of pdf file', required=False)
	gene_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	rev_conds = reverse_dict(conditions)
	if args["a"]:
		filters = ConfigSectionMap("Filters", Config)
	elif args["filter"]
		filters = args["filter"]
	else:
		filters=None
	data = pkg_resources.resource_filename('pyatactools', 'data/{}_ensembl_80.txt'.format(args["genome"]))

	if args["subparser_name"] == "tss":
		plot_tss_profile(conditions, data, int(args["width"])/2.0, filters, int(args["threads"]), args["output"])#
	else:
		plot_genebody_profile(conditions, data, filters, int(args["threads"]), args["output"])
		
main()