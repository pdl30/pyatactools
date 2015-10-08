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

def run_fastqc(fq1, outdir):
	command = "fastqc -q {} -o {}".format(fq1, outdir) 
	subprocess.call(command.split())

def find_adapters(outdir, fq):
	adapters = []
	name = re.sub(".fastq", "", fq)
	name = os.path.basename(name)
	command = "unzip -o -q {}/{}_fastqc.zip -d {}".format(outdir, name, outdir)
	subprocess.call(command.split())
	report = "{}/{}_fastqc/fastqc_data.txt".format(outdir, name)
	flist = open(report).readlines()
	parsing = False
	for line in flist:
		if line.startswith(">>Overrepresented sequences\tfail"):
			parsing = True
		elif line.startswith(">>END_MODULE"):
			parsing = False
		if parsing:
			if line.startswith(">>"):
				continue
			if line.startswith("#"):
				continue
			else:
				line = line.rstrip()
				word = line.split("\t")
				if word[3] != "No Hit":
					adapters.append(word[0])
	return adapters

def single_cut_adapters(adapters, fq1, outdir):
	trim = open('{}/trim_report.txt'.format(outdir), 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1
	command1 = "cutadapt -q 20 {0} --minimum-length=10 -o {1}/trimmed.fastq {2}".format(adapt1, outdir, fq1)
	p = subprocess.Popen(command1.split(), stdout=trim)
	p.communicate()

def paired_cut_adapters(adapters, fq1, outdir, rev_adapters, fq2):
	devnull = open('/dev/null', 'w')
	trim = open('{}/trim_report.txt'.format(outdir), 'w')
	adapt1 = ""
	for i in adapters:
		adapters = "-a {} ".format(i)
		adapt1 = adapters+adapt1

	adapt2 = ""
	for i in rev_adapters:
		adapters = "-a {} ".format(i)
		adapt2 = adapters+adapt2

	command1 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/tmp.2.fastq -o {1}/tmp.1.fastq {2} {3}".format(adapt1, outdir, fq1, fq2)
	p = subprocess.Popen(command1.split(), stdout=trim)
	p.communicate()
	command2 = "cutadapt -q 20 {0} --minimum-length=10 --paired-output {1}/trimmed_1.fastq -o {1}/trimmed_2.fastq {1}/tmp.2.fastq {1}/tmp.1.fastq".format(adapt2, outdir)
	p = subprocess.Popen(command2.split(), stdout=trim)
	p.communicate()
	cleanup = ["rm", "{0}/tmp.2.fastq".format(outdir), "{0}/tmp.1.fastq".format(outdir)]
	subprocess.call(cleanup, stdout=devnull)


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
		run_fastqc(fq1, args["outdir"])
		run_fastqc(fq2, args["outdir"])
		fwd_adapt = find_adapters(args["outdir"], fq1)
		rev_adapt = find_adapters(args["outdir"], fq2)
		if fwd_adapt or rev_adapt:
			print "==> Removing adapters...\n"
			paired_cut_adapters(fwd_adapt, fq1, args["outdir"], rev_adapt, fq2)
			fq1 = args["outdir"]+"/trimmed_1.fastq" 
			fq2 = args["outdir"]+"/trimmed_2.fastq"
		print "==> Aligning fastq's...\n"
		align(fq1, int(args["threads"]), args["ref"], args["outdir"], args["samname"], fq2=fq2)
	elif args["fastq"]:
		fq1 = args["fastq"]
		print "==> Running FastQC...\n"
		run_fastqc(fq1, args["outdir"])
		adapt = find_adapters(args["outdir"], fq1)
		if adapt:
			print "==> Removing adapters...\n"
			single_cut_adapters(adapt, fq1, args["outdir"])
			fq1 = args["outdir"]+"/trimmed.fastq" 
		print "==> Aligning fastq's...\n"
		align(fq1, int(args["threads"]), args["ref"], args["outdir"], args["samname"], fq2=None)