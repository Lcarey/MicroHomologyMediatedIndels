#!/home/lucas_pkuhpc/anaconda3/bin/python
import sys
import os
import argparse
import subprocess




# parse cmd line args
ap = argparse.ArgumentParser()
ap.add_argument("-bn", "--base_name", required=True, help="basename for all output files")
ap.add_argument("-g", "--genome_fasta_file", required=True, help="FASTA file with the genome / amplified sequence")
ap.add_argument("-r1", "--read1", required=True, help="FASTQ file first read")
ap.add_argument("-r2", "--read2", required=True, help="FASTQ file second read")
ap.add_argument("-fmh", "--find_mh", required=False, help="path to find_mh binary" , default='/lustre2/lucas_pkuhpc/bin/find_mh')
ap.add_argument("-bwa", "--bwa", required=False, help="path to bwa binary" , default='/lustre2/lucas_pkuhpc/bin/bwa')
ap.add_argument("-t", "--nthreads", required=False, help="number of threads" , default=1)
args = ap.parse_args()


# make sure we can find all the commands & files
die_with_error_flag = False
if not os.path.isfile( args.find_mh ):
	print( "find_mh not found, exiting\n   looked in:\t" ,  args.find_mh  , "\n" )
	die_with_error_flag = True
if not os.path.isfile( args.bwa ):
	print( "bwa not found, exiting\n   looked in:\t" ,  args.bwa  , "\n" )
	die_with_error_flag = True

if not os.path.isfile( args.read1 ):
	print( "read1 not found, exiting\n   looked in:\t" ,  args.read1  , "\n" )
	die_with_error_flag = True
if not os.path.isfile( args.read2 ):
	print( "read2 not found, exiting\n   looked in:\t" ,  args.read2  , "\n" )
	die_with_error_flag = True
if not os.path.isfile( args.genome_fasta_file ):
	print( "genome_fasta_file not found, exiting\n   looked in:\t" ,  args.genome_fasta_file  , "\n" )
	die_with_error_flag = True

if die_with_error_flag:
	sys.exit(1)



# start by indexing the genome and making the .cram file
print('--- indexing genome --- ')
cmd = args.bwa + ' index ' + args.genome_fasta_file 
print( cmd  )
subprocess.run(  cmd  ,  shell=True , check=True)
print(' ---- ---- ')

print(' ---- aligning reads to the genome ---- ')
cmd = args.bwa + ' mem -Y ' + ' -t ' + str(args.nthreads) + ' ' + args.genome_fasta_file + ' ' + args.read1 + ' ' + args.read2 \
	+ ' | samtools view -@ 3 -F 2048 -Sb ' \
	+ ' | samtools sort -O CRAM  --reference ' + args.genome_fasta_file + ' -@ 4 -m5G -o ' + args.base_name + '.cram'
print( cmd  )
subprocess.run(  cmd  ,  shell=True , check=True)
subprocess.run(  'samtools index '  + args.base_name + '.cram'  ,  shell=True , check=True)
print(' ---- ---- ')

