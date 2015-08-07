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
import tempfile
import argparse
import pkg_resources

def run_rcode(rscript):
	rcode = tempfile.NamedTemporaryFile(delete = False)
	rcode.write(rscript)
	rcode.close()
	print rcode.name
	try:
		subprocess.call(['Rscript', rcode.name])
	except:
		error("Error in running R script\n")
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)

def create_design_for_R(idict):
	output = tempfile.NamedTemporaryFile(delete = False)
	output.write("sampleName\tfileName\tcondition\n"),
	for key in sorted(idict.keys()):
		bam_name = os.path.basename(key)
		name = re.sub(".bam$", "", bam_name)
		output.write("{}\t{}\t{}\n".format(name, key, idict[key]))
	output.close()
	return output.name

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

def write_deseq(ifile, sample_dict, cond1, cond2, padj, f, outdir, gc_file, tmpdesign, cqn):
	print "==> Running differental expression analysis...\n"
	rscript =  "suppressMessages(library(cqn)); suppressMessages(library(scales))\n"
	rscript =  "suppressMessages(library(DESeq2))\n"
	rscript +=  "suppressMessages(library(RColorBrewer)); suppressMessages(library(ggplot2)); suppressMessages(library(gplots))\n"
	rscript += "pdata <- read.table('{}', header=T)\n".format(tmpdesign)
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(ifile)
	if f:
		rscript += "counts <- counts[,6:dim(counts)[2]]\n"
	else:
		pass
	if cqn:
		rscript += "x <- read.table({}, sep=\"\\t\", header=T)\n".format(gc_file)
		rscript += "g_length <- x[,4]-x[,3]+1; gc <- x[,2]/100; y <- cbind(g_length, gc); rownames(y) <- x[,1]\n" 
		rscript += "id <- match(rownames(counts), rownames(y)); z <- y[id,]; z2 <- z[complete.cases(z),]; id <- match(rownames(z2), rownames(counts)); counts2 <- counts[id,]\n"
		rscript += "cqn.subset <- cqn(counts2, lengths = z2[,1], x = z2[,2], sizeFactors = colsums(counts))\n"
		rscript += "cqnOffset <- cqn.subset$glm.offset; cqnNormFactors <- exp(cqnOffset)\n"
		rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts2, colData = data.frame(pdata), design = ~ condition)\n"
		rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
		rscript += "normalizationFactors(rnaseq_dds) <- cqnNormFactors"
	else:
		rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
		rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,3]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"	
	rscript += "rnaseq_res <- results(rnaseq_dds, contrast=c('condition','{0}','{1}'))\n".format(cond1, cond2)
	rscript += "rnaseq_sig <- rnaseq_res[which(rnaseq_res$padj <= {}),]\n".format(padj)
	rscript += "write.table(rnaseq_sig, file='{2}/{0}_vs_{1}_deseq2_significant.tsv', sep='\\t', quote=F)\n".format(cond1, cond2, outdir)
	rscript += "write.table(rnaseq_res, file='{2}/{0}_vs_{1}_deseq2_analysis.tsv', sep='\\t', quote=F)\n".format(cond1, cond2, outdir)
	rscript += "rld<-rlog(rnaseq_dds); colnames(rld) <- pdata$sampleName; rlogMat<-assay(rld); distsRL<-dist(t(assay(rld)))\n"
	rscript += "pdf('{0}/Sample-RLD-plots.pdf'); hmcol<-colorRampPalette(brewer.pal(9,'GnBu'))(100); mat<-as.matrix(distsRL)\n".format(outdir)
	rscript += "hc<-hclust(distsRL); par(cex.main=1);\n";
	rscript += "heatmap.2(mat,Rowv=as.dendrogram(hc),symm=TRUE,trace='none',col=rev(hmcol),margin=c(13,13),main='Sample-to-sample distances',cexRow=1,cexCol=1)\n"
	rscript += "data<-plotPCA(rld,intgroup=c('condition'),returnData=TRUE, ntop = nrow(rld)); percentVar<-round(100*attr(data,'percentVar'))\n"
	rscript += "ggplot(data,aes(PC1, PC2,color=condition,label=names))+geom_point(size=3)+xlab(paste0('PC1: ',percentVar[1],'% variance'))+ylab(paste0('PC2: ',percentVar[2],'% variance')) +geom_text(colour = 'black', alpha = 0.8, size = 2)\n"
	rscript += "dev.off()\n"
	return rscript

def get_cqn(ifile, sample_dict, f, outfile, gc_file, tmpdesign):
	print "==> Running differental expression analysis...\n"
	rscript =  "suppressMessages(library(cqn)); suppressMessages(library(scales))\n"
	rscript += "pdata <- read.table('{}', header=T)\n".format(tmpdesign)
	#Make sure they match!
	rscript += "counts <- read.table('{}', sep='\\t', header=T, row.names=1)\n".format(ifile)
	if f:
		rscript += "counts <- counts[,6:dim(counts)[2]]\n"
	else:
		pass
	rscript += "x <- read.table({}, sep=\"\\t\", header=T)\n".format(gc_file)
	rscript += "g_length <- x[,4]-x[,3]+1; gc <- x[,2]/100; y <- cbind(g_length, gc); rownames(y) <- x[,1]\n" 
	rscript += "id <- match(rownames(counts), rownames(y)); z <- y[id,]; z2 <- z[complete.cases(z),]; id <- match(rownames(z2), rownames(counts)); counts2 <- counts[id,]\n"
	rscript += "cqn.subset <- cqn(counts2, lengths = z2[,1], x = z2[,2], sizeFactors = colsums(counts))\n"
	rscript += "cqnOffset <- cqn.subset$glm.offset; cqnNormFactors <- exp(cqnOffset)\n"
	return rscript

def main():
	parser = argparse.ArgumentParser(description='Differential expression for ATAC-seq experiments. Runs DESEQ2\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	deseq2_parser = subparsers.add_parser('deseq', help="Runs DESEQ")
	deseq2_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	deseq2_parser.add_argument('-i','--input', help='Combined counts file from HTSeq or pyrna_count.py',required=True)
	deseq2_parser.add_argument('-p','--padj', help='Option for DESEQ2, default=0.05', default=0.05, required=False)
	deseq2_parser.add_argument('-n', action='store_true', help='Use CQN normalisation instead of DESEQ2', required=False)
	deseq2_parser.add_argument('-f', action='store_true', help='Use if featureCounts used as input', required=False)
	deseq2_parser.add_argument('-o','--output', help='Output counts file directory, default is current directory', required=False)

	norm_parser = subparsers.add_parser('norm', help="Outputs normalisation factors for samples in counts matrix")
	norm_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for usage!', required=False)
	norm_parser.add_argument('-i','--input', help='Combined counts file from HTSeq or pyrna_count.py',required=True)
	norm_parser.add_argument('-f', action='store_true', help='Use if featureCounts used as input', required=False)
	norm_parser.add_argument('-o','--output', help='Output counts file', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	gc_values= pkg_resources.resource_filename('pyatactools', 'data/gc_values.txt')
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)
	if args["output"]:
		output = args["output"]
	else:
		output = os.getcwd()
	if args["subparser_name"] == "deseq":
		design = create_design_for_R(conditions)
		comparisons = ConfigSectionMap("Comparisons", Config)
		for comp in comparisons:
			c = comparisons[comp].split(",")
			comps = [x.strip(' ') for x in c]
			rscript = write_deseq(args["input"], conditions, comps[0], comps[1], args["padj"], args["f"], output, gc_values, design, args["n"]) ##Needs changing!!!
			run_rcode(rscript)

	elif args["subparser_name"] == "norm":
		design = create_design_for_R(conditions)
		rscript = get_cqn(args["input"], conditions, args["f"], output, gc_values, design) ##Needs changing!!!
		run_rcode(rscript)
