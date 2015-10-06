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
from itertools import islice
import tempfile
#def sam_size(bam):
#	results=  {}
#	size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.flagstat(bam) ])
#	return size

def sam_size(bam):
	command = "samtools view -c -F 4 {}".format(bam)
	proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
	out = proc.communicate()[0]
	return int(out.upper())

def read_tss_pysam(bam, position_dict, halfwinwidth, norm, return_dict):
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	if norm:
		constant = 1e-6/float(norm[bam])
	else:
		constant = 1e-6/float(sam_size(bam))
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	count =0 
	for chrom, tss, strand in position_dict.values():
		count += 1
		coverage = {}
		chrom_start = tss-halfwinwidth
		if chrom_start <0: chrom_start=0
		chrom_end = tss+ halfwinwidth
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
		for i in range(1, 2*halfwinwidth):
			coverage[i] = 0.0
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			#ref_pos = pileupcolumn.pos
			if strand == "+":
				ref_pos = pileupcolumn.pos - chrom_start
			elif strand == "-":
				ref_pos = chrom_end - pileupcolumn.pos 
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read += constant
			coverage[ref_pos] = cover_read
		tmp1 = [coverage[k] for k in sorted(coverage)]

		for i in range(0,len(tmp1)):
			aggreagated_cvg[i] += tmp1[i]*10000

	for key in aggreagated_cvg:
		profile[key] = aggreagated_cvg[key]/float(count)
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
		tx_end = int( fields[3] )
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

def genebody_coverage(bam, position_list, norm, return_dict):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	if norm:
		constant = 1e-6/float(norm[bam])
	else:
		constant = 1e-6/float(sam_size(bam))
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	gene_finished = 0
	for chrom, strand, positions in position_list.values():
		before = positions[0]-1011
		if before ==0:
			before = 0
		after = positions[-1] + 1000
		before_list = range(before, positions[0]-1-10, 10)
		after_list = range(positions[-1]+10, after+10, 10)
		new_positions = before_list + positions + after_list
		if len(new_positions) != 300:
			continue
		coverage = {}
		for i in new_positions:
			coverage[i] = 0.0
		chrom_start = new_positions[0]
		if chrom_start <0: chrom_start=0
		chrom_end = new_positions[-1]
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in new_positions:
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
		
	tmp3 = numpy.zeros( 300, dtype='f' ) 
	for key in aggreagated_cvg:
		tmp3[key] = aggreagated_cvg[key]/float(len(position_list.keys()))
	return_dict[bam] = tmp3

def read_tss_anno(anno, gene_filter):
	positions = {}
	g_filter = []
	if gene_filter:
		with open(gene_filter) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				g_filter.append(word[0])
	with open(anno) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
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
					if strand == "+":
						positions[word[0]] = (chrom, int(word[2]), strand)
					elif strand == "-":
						positions[word[0]] = (chrom, int(word[3]), strand)
			else:
				if strand == "+":
					positions[word[0]] = (chrom, int(word[2]), strand)
				elif strand == "-":
					positions[word[0]] = (chrom, int(word[3]), strand)
	return positions

def plot_tss_profile(conditions, anno, halfwinwidth, gene_filter, threads, comb, outname, norm):
	halfwinwidth = int(halfwinwidth)
	positions = read_tss_anno(anno, gene_filter)
	fname = None
	manager = Manager()
	return_dict = manager.dict()
	pool = Pool(threads)
	pool.map(read_tss_function, itertools.izip(list(conditions.keys()), itertools.repeat(positions), itertools.repeat(halfwinwidth),itertools.repeat(norm), itertools.repeat(return_dict)))
	pool.close()
	pool.join()	
	if comb:
		combined_profiles = {}
		rev_conds = reverse_dict(conditions)
		for key in rev_conds:
			c = 0
			for bam in rev_conds[key]:
				if key not in combined_profiles:
					combined_profiles[key] = return_dict[bam]
				else:
					combined_profiles[key] += return_dict[bam]
				c+= 1
			combined_profiles[key] = combined_profiles[key]/float(c)
		pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		for key in combined_profiles.keys():
			pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), combined_profiles[key], label=key)
		pyplot.legend(prop={'size':8})
		pyplot.savefig(outname+".pdf")
		pyplot.close()
	else:
		if len(list(conditions.keys())) < 11:
			pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		else:
			colormap = pyplot.cm.gist_ncar
			pyplot.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, len(list(conditions.keys())))])
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), return_dict[key], label=conditions[key])
		pyplot.legend(prop={'size':8})
		pyplot.savefig(outname+".pdf")
		pyplot.close()

def plot_genebody_profile(conditions, anno, gene_filter, threads, comb, outname, norm):
	positions = genebody_percentile(anno, gene_filter, mRNA_len_cut = 100)
	fname = None
	manager = Manager()
	return_dict = manager.dict()
	pool = Pool(threads)
	pool.map(read_gene_function, itertools.izip(list(conditions.keys()), itertools.repeat(positions), itertools.repeat(norm), itertools.repeat(return_dict)))
	pool.close()
	pool.join()	
	if comb:
		pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		combined_profiles = {}
		rev_conds = reverse_dict(conditions)
		for key in rev_conds:
			c = 0
			for bam in rev_conds[key]:
				if key not in combined_profiles:
					combined_profiles[key] = return_dict[bam]
				else:
					combined_profiles[key] += return_dict[bam]
				c+= 1
			combined_profiles[key] = combined_profiles[key]/float(c)
		for key in combined_profiles.keys():
			pyplot.plot( numpy.arange( -100, 200 ), combined_profiles[key], label=key)
	else:
		if len(list(conditions.keys())) < 11:
			pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		else:
			colormap = pyplot.cm.gist_ncar
			pyplot.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, len(list(conditions.keys())))])
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( -100, 200 ), return_dict[key], label=conditions[key])
	pyplot.xticks([-100, 0, 100, 200], ['-1000', 'TSS', 'TES', "+1000"])
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

def get_insert_fun(args):
	return get_insert(*args)

def get_insert(sam, return_dict):
	cvg = collections.defaultdict(int)
	test2 = tempfile.NamedTemporaryFile(delete=False)
	command = "samtools view -h {} -o {}".format(sam, test2.name)
	subprocess.call(command, shell=True)
	with open(test2.name) as f:
		lines = [next(f) for x in xrange(5000000)]
		for line in lines:
			if line.startswith("@"):
				pass
			else:
				word = line.rstrip().split("\t")
				if len(word) < 9: #Presumes it removes unaligned reads
					pass
				else:
					if word[2] == "chrM" or word[2] == "M": #Filter because of not relevant
						pass
					else:
						if int(word[8]) < 1:
							pass
						elif int(word[8]) > 649:
							pass
						else:
							cvg[int(word[8])] += 1
	profile = numpy.zeros( 650, dtype='i' ) 
	for key in cvg:
		profile[key] = cvg[key]
	return_dict[sam] = profile

def plot_insert_bed2(sam, bedfile, return_dict):
	cvg = collections.defaultdict(int)
	test = tempfile.NamedTemporaryFile(delete = False)
	test2 = tempfile.NamedTemporaryFile(delete=False)
	command = "bedtools intersect -abam {} -b {} > {}".format(sam, bedfile, test.name)
	subprocess.call(command, shell=True)
	command = "samtools view -h {} -o {}".format(test.name, test2.name)
	subprocess.call(command, shell=True)
	with open(test2.name) as f:
		for line in f:
			if line.startswith("@"):
				pass
			else:
				word = line.rstrip().split("\t")
				if len(word) < 9: #Presumes it removes unaligned reads
					pass
				else:
					if word[2] == "chrM" or word[2] == "M": #Filter because of not relevant
						pass
					else:
						if abs(int(word[8])) < 1:
							pass
						elif abs(int(word[8])) > 100: 
							pass
						else:
							cvg[abs(int(word[8]))] += 1
	profile = numpy.zeros( 101, dtype='i' )
	for i in range(0, 101):
		if cvg[i]:
			profile[i] = cvg[i]
		else:
			profile[i] = 0
	return_dict[sam] = profile

def plot_insert_bed_function(args):
	return plot_insert_bed2(*args)

def plot_inserts(conditions, threads, output, bedfile, comb):
	if bedfile:
		positions = set()
		manager = Manager()
		return_dict = manager.dict()		
		pool = Pool(int(threads))
		pool.map(plot_insert_bed_function, itertools.izip(list(conditions.keys()), itertools.repeat(bedfile), itertools.repeat(return_dict)))
		pool.close()
		pool.join()	
		if comb:
			pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
			combined_profiles = {}
			rev_conds = reverse_dict(conditions)
			for key in rev_conds:
				c = 0
				for bam in rev_conds[key]:
					if key not in combined_profiles:
						combined_profiles[key] = return_dict[bam]
					else:
						combined_profiles[key] += return_dict[bam]
					c+= 1
				combined_profiles[key] = combined_profiles[key]/float(c)
			for key in combined_profiles.keys():
				pyplot.plot( numpy.arange( 0, 101), combined_profiles[key], label=key)
			pyplot.legend(prop={'size':8})
			pyplot.savefig(output+".pdf")
		else:
			for key in return_dict.keys():
				pyplot.plot( numpy.arange( 0, 101 ), return_dict[key], label=conditions[key])
			pyplot.legend(prop={'size':8})
			pyplot.savefig(output+".pdf")
	else:
		colormap = pyplot.cm.gist_ncar
		pyplot.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, len(list(conditions.keys())))])
		manager = Manager()
		return_dict = manager.dict()
		pool = Pool(int(threads))
		pool.map(get_insert_fun, itertools.izip(list(conditions.keys()), itertools.repeat(return_dict)))
		pool.close()
		pool.join()	
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( 0, 650 ), return_dict[key], label=conditions[key])
		pyplot.axvline(x=147.,color='k',ls='dashed')
		pyplot.axvline(x=294.,color='k',ls='dashed')
		pyplot.axvline(x=441.,color='k',ls='dashed')
		pyplot.legend(prop={'size':8})
		pyplot.savefig(output+".pdf")

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def read_peaks(bed):
	positions = {}
	c = 0
	with open(bed) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			med = int(word[2]) + int(word[1])
			med = round(med/float(2))
			positions[c] = (word[0], med)
			c +=  1
	return positions

def read_peak_pysam(bam, halfwinwidth, position_dict, norm, return_dict):
	profile = numpy.zeros( 2*halfwinwidth, dtype="f" )
	if norm:
		constant = 1e-6/float(norm[bam])
	else:
		constant = 1e-6/float(sam_size(bam))
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	for chrom, summit in position_dict.values():
		coverage = {}
		chrom_start = summit - halfwinwidth
		if chrom_start <0: chrom_start=0
		chrom_end = summit + halfwinwidth
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
		coverage = {}
		for i in range(1, 2*int(halfwinwidth)):
			coverage[i] = 0.0
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			#ref_pos = pileupcolumn.pos
			ref_pos = pileupcolumn.pos - chrom_start
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read += constant
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
	for key in aggreagated_cvg:
		profile[key] = aggreagated_cvg[key]/float(len(position_dict.keys()))
	return_dict[bam] = profile

def read_peak_function(args):
	return read_peak_pysam(*args)

def plot_peak_profile(conditions, bed, halfwinwidth, threads, comb, outname, norm):
	positions = read_peaks(bed)
	fname = None
	manager = Manager()
	return_dict = manager.dict()
	pool = Pool(threads)
	pool.map(read_peak_function, itertools.izip(list(conditions.keys()), itertools.repeat(halfwinwidth), itertools.repeat(positions), itertools.repeat(norm), itertools.repeat(return_dict)))
	pool.close()
	pool.join()	
	if comb:
		pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		combined_profiles = {}
		rev_conds = reverse_dict(conditions)
		for key in rev_conds:
			c = 0
			for bam in rev_conds[key]:
				if key not in combined_profiles:
					combined_profiles[key] = return_dict[bam]
				else:
					combined_profiles[key] += return_dict[bam]
				c+= 1
			combined_profiles[key] = combined_profiles[key]/float(c)
		for key in combined_profiles.keys():
			pyplot.plot( numpy.arange( 0, halfwinwidth*2 ), combined_profiles[key], label=key)
	else:
		if len(list(conditions.keys())) < 11:
			pyplot.rc('axes', color_cycle=['b','r', 'c', 'm', 'y', 'k', 'gray', "green", "darkred", "skyblue"])
		else:
			colormap = pyplot.cm.gist_ncar
			pyplot.gca().set_color_cycle([colormap(i) for i in numpy.linspace(0, 0.9, len(list(conditions.keys())))])
		for key in return_dict.keys():
			pyplot.plot( numpy.arange( 0, halfwinwidth*2 ), return_dict[key], label=conditions[key])
	pyplot.legend(prop={'size':8})
	pyplot.savefig(outname+".pdf")	

def main():
	parser = argparse.ArgumentParser(description='Takes BED files and intersect them with regions, uses TSS regions by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	tss_parser = subparsers.add_parser('tss', help='TSS plotter')
	tss_parser.add_argument('-c', '--config', help='BAM as keys', required=True)
	tss_parser.add_argument('-f', '--filter', help='Gene name per line, filters TSS regions', required=False)
	tss_parser.add_argument('-o', '--output', help='Output name of pdf file', required=True)
	tss_parser.add_argument('-d', action='store_true', help='Use combinations for plotting', required=False)
	tss_parser.add_argument('-w', '--width', help='Width of region, default=1000', default=1000, required=False)
	tss_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	tss_parser.add_argument('-n', action='store_true', help='Use [Norm] as constant from config', required=False)
	gene_parser = subparsers.add_parser('gene', help='Genebody plotter')
	gene_parser.add_argument('-c', '--config', help='BAM as keys', required=False)
	gene_parser.add_argument('-f', '--filter', help='Gene name per line, filters TSS regions', required=False)
	gene_parser.add_argument('-d', action='store_true', help='Use combinations for plotting', required=False)
	gene_parser.add_argument('-o', '--output', help='Output name of pdf file', required=False)
	gene_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	gene_parser.add_argument('-n', action='store_true', help='Use [Norm] as constant from config', required=False)
	insert_parser = subparsers.add_parser('insert', help='Insert histogram plotter')
	insert_parser.add_argument('-c', '--config', help='BAM as keys, must be bams!', required=False)
	insert_parser.add_argument('-o', '--output', help='Output name of pdf file', required=False)
	insert_parser.add_argument('-r', '--bedfile', help='Plot insert size of bed file regions', required=False)
	insert_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	insert_parser.add_argument('-d', action='store_true', help='Use combinations for plotting. Only works if supplying bed file', required=False)
	peak_parser = subparsers.add_parser('peak', help='Peak profiles plotter')
	peak_parser.add_argument('-c', '--config', help='BAM as keys', required=False)
	peak_parser.add_argument('-b', '--bed', help='Peak file, should be of standard width',required=False)
	peak_parser.add_argument('-w', '--width', help='Width of region, default=1000', default=1000, required=False)
	peak_parser.add_argument('-o', '--output', help='Output name of pdf file', required=False)
	peak_parser.add_argument('-t', '--threads', help='Threads, default=8', default=8, required=False)
	peak_parser.add_argument('-d', action='store_true', help='Use combinations for plotting', required=False)
	peak_parser.add_argument('-n', action='store_true', help='Use [Norm] as constant from config', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)

	if args["subparser_name"] == "tss":
		if args["filter"]:
			filters = args["filter"]
		else:
			filters=None
		if args["n"]:
			norm = ConfigSectionMap("Norm", Config)
		else:
			norm = None
		data = pkg_resources.resource_filename('pyatactools', 'data/mm10_ensembl_80.txt')
		plot_tss_profile(conditions, data, int(args["width"])/2.0, filters, int(args["threads"]), args["d"], args["output"], norm)
	elif args["subparser_name"] == "gene":
		if args["filter"]:
			filters = args["filter"]
		else:
			filters=None
		if args["n"]:
			norm = ConfigSectionMap("Norm", Config)
		else:
			norm = None
		data = pkg_resources.resource_filename('pyatactools', 'data/mm10_ensembl_80.txt')
		plot_genebody_profile(conditions, data, filters, int(args["threads"]), args["d"], args["output"], norm)
	elif args["subparser_name"] == "insert":
		plot_inserts(conditions, int(args["threads"]), args["output"], args["bedfile"], args["d"])
	elif args["subparser_name"] == "peak":
		if args["n"]:
			norm = ConfigSectionMap("Norm", Config)
		else:
			norm = None
		plot_peak_profile(conditions, args["bed"], int(args["width"])/2.0, int(args["threads"]), args["d"], args["output"], norm)