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
from pychiptools.utilities import fastqc
import argparse

def align(fq1, threads, fasta, outdir, outname, fq2=None):
	if fq2:
		name = re.sub("_1.fastq", "", fq1)
		name = re.sub("_1.fq", "", fq1)
		command = "bwa aln -t {}, {} {} > {}/{}_1.bwa".format(threads, fasta, fq1, outdir, name)
		subprocess.call(command, shell=True)
		command = "bwa aln -t {}, {} {} > {}/{}_2.bwa".format(threads, fasta, fq2, outdir, name)
		subprocess.call(command, shell=True)
		command = "bwa sampe -a 2000 -n 1 {0} {4}/{1}_1.bwa {4}/{1}_2.bwa {2} {3}  > {4}/{1}.sam".format(fasta, name, fq1, fq2, outdir)
		subprocess.call(command, shell=True)
	else:
		name = re.sub("_1.fastq", "", fq1)
		name = re.sub("_1.fq", "", fq1)
		command = "bwa aln -t {} {} {} > {}/{}_1.bwa".format(threads, fasta, fq1, outdir, name)
		subprocess.call(command, shell=True)
		command = "bwa samse -a 2000 -n 1 {0} {3}/{1}_1.bwa {2} > {3}/{4}.sam".format(fasta, name, fq1, outdir, outname)
		subprocess.call(command, shell=True)

def main():
	parser = argparse.ArgumentParser(description='ChIP-seq fastqc, trimmer and bowtie wrapper\n ')
	parser.add_argument('-f','--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p','--pair', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-r','--ref', help='Path to reference fasta', required=True)
	parser.add_argument('-t','--threads', help='For bowtie2, how many threads to use (i.e. -p option on bowtie2', default=1, required=False)
	parser.add_argument('-n','--samname', help='Name of output sam file', required=True)
	parser.add_argument('-o','--outdir', help='Name of results directory', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	path = os.getcwd()
	if os.path.isdir(args["outdir"]):
		print("Results directory already exists!")
	else:
		subprocess.call(["mkdir", args["outdir"]])
	if args["pair"]:
		fq1 = args["pair"][0]
		fq2 = args["pair"][1]
		print "==> Running FastQC...\n"
		fastqc.run_fastqc(fq1, args["outdir"])
		fastqc.run_fastqc(fq2, args["outdir"])
		fwd_adapt = fastqc.find_adapters(args["outdir"], fq1)
		rev_adapt = fastqc.find_adapters(args["outdir"], fq2)
		if fwd_adapt or rev_adapt:
			print "==> Removing adapters...\n"
			fastqc.paired_cut_adapters(fwd_adapt, fq1, args["outdir"], rev_adapt, fq2)
			fq1 = args["outdir"]+"/trimmed_1.fastq" 
			fq2 = args["outdir"]+"/trimmed_2.fastq"
		print "==> Aligning fastq's...\n"
		align(fq1, int(args["threads"]), args["ref"], args["outdir"], args["samname"], fq2=fq2)
	elif args["fastq"]:
		fq1 = args["fastq"]
		print "==> Running FastQC...\n"
		fastqc.run_fastqc(fq1, args["outdir"])
		adapt = fastqc.find_adapters(args["outdir"], fq1)
		if adapt:
			print "==> Removing adapters...\n"
			fastqc.single_cut_adapters(adapt, fq1, args["outdir"])
			fq1 = args["outdir"]+"/trimmed.fastq" 
		print "==> Aligning fastq's...\n"
		align(fq1, int(args["threads"]), args["ref"], args["outdir"], args["samname"], fq2=None)