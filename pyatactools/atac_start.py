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
import argparse
import tempfile
import pkg_resources
from multiprocessing import Pool, Manager
import pybedtools


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

def ddup(conditions, sam=True):
	path = os.getcwd()
	ddup_sam = []
	#Convert sam to bam, rmdup bam, convert bam to sam and remove intermediates, could really do with a pipe!
	if sam:
		for sam in conditions:
			name = re.sub(".sam$", "", sam)
			command = "samtools view -bS {0}.sam | samtools sort - {0}_sort".format(name)
			subprocess.call(command, shell=True)
			command = "samtools rmdup {0}_sort.bam - | samtools view -h - > {0}_ddup.sam".format(name)
			subprocess.call(command, shell=True)
			os.remove("{0}_sort.bam".format(name))
			ddup_sam.append("{}_ddup.sam".format(name)) #Might find pipe useful here???
	else:
		for bam in conditions:
			name = re.sub(".bam$", "", bam)
			command = "samtools rmdup {0}.bam - | samtools view -h - > {0}_ddup.sam".format(name)
			subprocess.call(command, shell=True)
			ddup_sam.append("{}_ddup.sam".format(name)) #Might find pipe useful here???
	return ddup_sam

def transdense(sam, transdense_dir, return_dict):
	fh = tempfile.NamedTemporaryFile(delete = False)
	filename = os.path.basename(sam)
	name = re.sub(".sam$", "", filename)
	with open(sam) as f:
		for line in f:
			if line.startswith("@"):
				pass
			else:
				word = line.rstrip().split("\t")
				if len(word) > 8: #Presumes it removes unaligned reads
					if word[2] == "chrM" or word[2] == "M": #Filter because of not relevant
						pass
					else:
						if int(word[8]) > 0:		
							start = int(word[3]) - 11 
							end = int(word[3]) + 17
							if start > 0:
								fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], start, end)),
						elif int(word[8]) < 0:
							start = int(word[3]) - int(word[8])
							start -= 19
							end = int(word[3]) - int(word[8])
							end += 9

							if start > 0:
								fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], start, end)),
	fh.close()
	command = "sort -k1,1 -k2,2n {} | uniq > {}/{}_mytransDense.bed".format(fh.name, transdense_dir, name)
	return_dict["{}/{}_mytransDense.bed".format(transdense_dir, name)] = 1
	subprocess.call(command, shell=True)
	os.remove(fh.name)

def get_nfree(sam, nfree_dir, return_dict):
	filename = os.path.basename(sam)
	name = re.sub(".sam", "", filename)
	fh = tempfile.NamedTemporaryFile(delete = False)
	with open(sam) as f:
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
						if int(word[8]) == 0:
							pass
						else:
							if int(word[8]) < 100 and int(word[8]) > -100:
								if int(word[8]) > 0:		
									start = int(word[3]) - 11 
									end = int(word[3]) + 17
									if start > 0:
										fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], start, end)),
								elif int(word[8]) < 0:
									start = int(word[3]) - int(word[8])
									start -= 19
									end = int(word[3]) - int(word[8])
									end += 9
									if start > 0:
										fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], start, end)),
	fh.close()
	command = "sort -k1,1 -k2,2n {} | uniq > {}/{}_mynfree.bed".format(fh.name, nfree_dir, name)
	return_dict["{}/{}_mynfree.bed".format(nfree_dir, name)] = 1
	subprocess.call(command, shell=True)
	os.remove(fh.name)

def get_npres(sam, npres_dir, return_dict):
	fh = tempfile.NamedTemporaryFile(delete = False)
	filename = os.path.basename(sam)
	name = re.sub(".sam", "", filename)
	with open(sam) as f:
		for line in f:
			if line.startswith("@"):
				pass
			else:
				word = line.rstrip().split("\t")
				if len(word) > 9:
					if word[2] == "chrM" or word[2] == "M":
						pass
					else:
						if int(word[8]) > 180 and int(word[8]) < 247:
							mid = int(word[3])+ (float(word[8])/2)
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(mid-75-1)), int(round(mid+75-1)))), 
						elif int(word[8]) > 315 and int(word[8]) < 473:
							mid = float(word[8])/3
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(int(word[3])+mid-75-1)), int(round(int(word[3])+mid+75-1)))), 
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(int(word[3])+(mid*2)-75-1)), int(round(int(word[3])+(mid*2)+75-1)))), 
						elif int(word[8]) > 558 and int(word[8]) < 615:
							mid = int(word[3])+ (float(word[8])/4)
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(int(word[3])+mid-75-1)), int(round(int(word[3])+mid+75-1)))), 
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(int(word[3])+(mid*2)-75-1)), int(round(int(word[3])+(mid*2)+75-1)))), 
							fh.write("{}\t{}\t{}\tT\t0\t+\n".format(word[2], int(round(int(word[3])+(mid*3)-75-1)), int(round(int(word[3])+(mid*3)+75-1)))), 
	fh.close()
	command = "sort -k1,1 -k2,2n {} | uniq > {}/{}_mynPres.bed".format(fh.name, npres_dir, name)
	return_dict["{}/{}_mynPres.bed".format(npres_dir, name)] = 1
	subprocess.call(command, shell=True)                                             
	os.remove(fh.name) #Always do 

def convert_bed_bw(bed, chrom):
	name = re.sub(".bed", "", bed)
	#inbed = pybedtools.BedTool(bed)
	#outcov = inbed.genome_coverage(bg=True, genome=chrom)
	command = "bedtools genomecov -i {} -bg -g {} > {}.bedGraph".format(bed, chrom, name)
	#outcov.saveas(name+".bedGraph")
	subprocess.call(command, shell=True)
	command = ["bedGraphToBigWig", name+".bedGraph", chrom, name+".bw"]
	subprocess.call(command)
	os.remove(name+".bedGraph")

def function1(args):
	return transdense(*args)

def function2(args):
	return get_nfree(*args)

def function3(args):
	return get_npres(*args)

def function4(args):
	return convert_bed_bw(*args)

def main():
	parser = argparse.ArgumentParser(description='Takes bam files and preprocess\'s for analysis\n')
	parser.add_argument('-c', '--config', help='Conditions containing Sam/Bam files, values are naming', required=True)
	parser.add_argument('-g', '--genome', help='Genome the samples are aligned to, options include	 mm10/mm9/hg19', required=True)
	parser.add_argument('-o', '--outdir', help='Output directory, will create transdense, nfree and npres directories', required=True)
	parser.add_argument('-t', '--threads', help='threads, default=1', default=1, required=False)
	parser.add_argument('-b', action='store_true', help='Use if Config contains bam files', required=False) 
	parser.add_argument('-d', action='store_true', help='De-duplicate sam/bam files', required=False) 
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)

	chrom = pkg_resources.resource_filename('pyatactools', 'data/{}.chrom.sizes'.format(args["genome"]))
	if not os.path.isfile(chrom):
		raise Exception("Unsupported Genome!")
	#ddup(conditions)
	transdense_dir = os.path.join(args["outdir"], "transdense")
	nfree_dir = os.path.join(args["outdir"], "nfree")
	npres_dir = os.path.join(args["outdir"], "npres")
	if not os.path.isdir(transdense_dir):
		os.makedirs(transdense_dir)
		os.makedirs(nfree_dir)
		os.makedirs(npres_dir)
	if args["d"]:
		if args["b"]:
			ddup_bams = ddup(conditions, False)
		else:
			ddup_bams = ddup(conditions, True)
	else:
		ddup_bams = conditions

	pool = Pool(int(args["threads"]))
	manager = Manager()
	return_dict = manager.dict()
	pool.map(function1, itertools.izip(ddup_bams, itertools.repeat(transdense_dir), itertools.repeat(return_dict)))
	pool.map(function4, itertools.izip(list(return_dict.keys()), itertools.repeat(chrom)))
	return_dict = manager.dict()
	pool.map(function2, itertools.izip(ddup_bams, itertools.repeat(nfree_dir), itertools.repeat(return_dict)))
	pool.map(function4, itertools.izip(list(return_dict.keys()), itertools.repeat(chrom)))
	convert_bed_bw(args["genome"], chrom, nfree_list)
	return_dict = manager.dict()
	pool.map(function3, itertools.izip(ddup_bams, itertools.repeat(npres_dir), itertools.repeat(return_dict)))
	pool.map(function4, itertools.izip(list(return_dict.keys()), itertools.repeat(chrom)))

