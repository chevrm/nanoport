import sys, h5py, os
from glob import glob
from optparse import OptionParser
from Bio import SeqIO

## Print Fasta from Fastq
def getfasta (fastq, type, hd5name):
	fastqpart = fastq.splitlines()
	fastqpart[0] = fastqpart[0].replace("@", '')
	return ">%s_%s\n%s\n" % (fastqpart[0], type, fastqpart[1])

## Get Fastq from Fast5
def getfastq (hdfile, loc):
	if loc in hdfile:
		return hdfile[loc][()]

## Median function
def median(list):
        list.sort()
        if len(list) % 2:
                return list[int(len(list)/2)]
        else:
                return (list[int(len(list)/2)+1] + list[int(len(list)/2)]) / 2

## Parse Command Line Args
parser = OptionParser()
parser.add_option("--fasta", action="store_true", dest="fasta",
	default=False, help="Output to fasta (default fastq)")
parser.add_option("--with2d", action="store_true", dest="with2d",
	default=False, help="Output only reads with 2D basecalls")
parser.add_option("--template", action="store_true", dest="template",
	default=True, help="Output from Template basecalls")
parser.add_option("--complement", action="store_true", dest="complement",
	default=False, help="Output from Complement basecalls")
parser.add_option("--2d", action="store_true", dest="flag2d",
	default=False, help="Output from 2D basecalls")
(option, arg) = parser.parse_args()

## Read in Fast5, Print to File
fast5dir = arg[0]
fast5dir += "/*.fast5"
format = "fastq"
if (option.fasta):
        format = "fasta"
location = '/Analyses/Basecall_2D_000/BaseCalled_template/Fastq'
type = '1T'
if (option.complement):
	location = '/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'
        type = '1C'
elif (option.flag2d):
	location = '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
	type = '2D'
suffix = ".".join([type, format])
for fast5 in glob(fast5dir):
	hd5 = h5py.File(fast5, 'r')
	fastq = getfastq(hd5, location)
	if(fastq is not None):
		out = open(fast5.replace("fast5", suffix),'w')
		if (option.fasta):
			out.write(getfasta(fastq, type, fast5))
		else:
			out.write("%s" % fastq)
		out.close()
	hd5.close()

## Calculate summary stats
long, totseqs, totbases = 0, 0, 0
lenlist = []
over20dict = {}
fast5dir = fast5dir.replace("fast5", suffix )
combo = arg[0].replace("/", '')
catcommand = "cat %s > %s.%s 2>/dev/null" % (fast5dir, combo, suffix)
os.system(catcommand)
for seqfile in glob(fast5dir):
	seqhandle = open(seqfile, "rU")
	for record in SeqIO.parse(seqhandle, format):
		lenlist.append(len(record.seq))
		if len(record.seq) >= 20000:
			over20dict[record.id] = record.seq
		totseqs += 1
		totbases += len(record.seq)
		if len(record.seq) > long:
			long = len(record.seq)
if (totseqs == 0):
	sys.exit("No %s sequences found.\n" % type)
mean = totbases / totseqs
med = median(lenlist)
sys.stderr.write("\nTotal bases:\t%s\nTotal reads:\t%s\nMean read length:\t%s\nMedian read length:\t%s\nLongest read:\t%s\n" % (totbases, totseqs, mean, med, long) )
sys.stderr.write("Reads over 20kb:\n")
counter = 1
for id in over20dict:
	sys.stderr.write("\t%s\t%s\tlen=%s\n" % (counter, id, len(over20dict[id]) ) )
	counter += 1
