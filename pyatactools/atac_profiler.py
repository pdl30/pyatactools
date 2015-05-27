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

def sam_size(ibam):
	results=  {}
	for bam in ibam:
		size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(bam) ])
		results[bam] = size
	return results

def read_bam_tss(bam, positions, halfwinwidth, return_dict):
	constant = 1000000/float(sam_size[bam])
	bamfile = HTSeq.BAM_Reader(bam)
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	c = 0
	for p in positions: #Problem if ensembl vs UCSC
		start = p.pos - halfwinwidth 
		end = p.pos + halfwinwidth 
		if start < 0:
			start = 0
		window = HTSeq.GenomicInterval( chrom, start, end, "." )
		for almnt in bamfile[ window ]:	
			if p.strand == "+":
				start_in_window = almnt.iv.start - p.pos + halfwinwidth 
				end_in_window   = almnt.iv.end - p.pos + halfwinwidth 
			else:
				start_in_window = p.pos + halfwinwidth - almnt.iv.end
				end_in_window   = p.pos + halfwinwidth - almnt.iv.start
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			if start_in_window >= 2*halfwinwidth or end_in_window < 0:
				continue
			profile[ start_in_window : end_in_window ] += constant
		c += 1
	profile = profile/float(c) #Average over the number of TSS regions
	return_dict[bam] = profile

def read_bed_tss(bed, positions, halfwinwidth, return_dict):
	constant = 1
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	c = 0
	for p in positions: #Problem if ensembl vs UCSC
		start = p.pos - halfwinwidth 
		end = p.pos + halfwinwidth 
		if start < 0:
			start = 0
		window = HTSeq.GenomicInterval( chrom, start, end, "." )
		with open(bed) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
			if p.strand == "+":
				start_in_window = int(word[1]) - p.pos + halfwinwidth 
				end_in_window   = int(word[2]) - p.pos + halfwinwidth 
			else:
				start_in_window = p.pos + halfwinwidth - int(word[2])
				end_in_window   = p.pos + halfwinwidth - int(word[1])
			start_in_window = max( start_in_window, 0 )
			end_in_window = min( end_in_window, 2*halfwinwidth )
			if start_in_window >= 2*halfwinwidth or end_in_window < 0:
				continue
			profile[ start_in_window : end_in_window ] += constant
		c += 1
	profile = profile/float(c) #Average over the number of TSS regions
	return_dict[bed] = profile

def read_bam_tss_function(args):
	return read_bam_tss(*args)

def read_bed_tss_function(args):
	return read_bed_tss(*args)

def plot_tss_profile(conditions, gff, halfwinwidth, gene_filter, threads, outname, bed=False):
	gtffile = HTSeq.GFF_Reader( gff )
	halfwinwidth = int(halfwinwidth)
	tsspos = set()
	#Add a filter for genes which are in the input file
	if gene_filter:
		g_filter = {}
		with open(gene_filter) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				g_filter[word[0]] = 1
		for feature in gtffile:
			if feature.name in g_filter:
				if feature.type == "exon" and feature.attr["exon_number"] == "1":
					tsspos.add( feature.iv.start_d_as_pos )
	else:
		for feature in gtffile:
			if feature.type == "exon" and feature.attr["exon_number"] == "1":
				tsspos.add( feature.iv.start_d_as_pos )
	manager = Manager()
	return_dict = manager.dict()
	pool = Pool(threads)
	if not bed:
		#pool.map(read_bam_tss_function, itertools.izip(list(conditions.keys()), itertools.repeat(tsspos), itertools.repeat(halfwinwidth), itertools.repeat(return_dict)))
		print list(conditions.keys())[0], tsspos, halfwinwidth, return_dict
		read_bam_tss(list(conditions.keys())[0], tsspos, halfwinwidth, return_dict)
	else:
		pool.map(read_bed_tss_function, itertools.izip(list(conditions.keys()), itertools.repeat(tsspos), itertools.repeat(halfwinwidth), itertools.repeat(return_dict)))
	pool.close()
	pool.join()	
	for key in return_dict.keys():
		pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), return_dict[key], label=conditions[key])  
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
	parser.add_argument('-c', '--config', help='BED/BAM as keys, if using BEDS, include [size] section', required=True)
	parser.add_argument('-g', '--gtf', help='GTF in ucsc format', required=True)
	parser.add_argument('-o', '--output', help='Output name of pdf file', required=True)
	parser.add_argument('-w', '--width', help='Width of region, default=1000', default=1000, required=False)
	parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	parser.add_argument('-b', action='store_true', help='Use if beds in config', required=False) 
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	rev_conds = reverse_dict(conditions)
	plot_tss_profile(conditions, args["gtf"], int(args["width"])/2.0, None, int(args["threads"]), args["output"], args["b"])

main()