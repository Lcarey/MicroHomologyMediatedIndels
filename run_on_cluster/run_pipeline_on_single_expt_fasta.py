#!/home/lucas_pkuhpc/anaconda3/bin/python
import sys
import os
import argparse
import subprocess
import time




# parse cmd line args
ap = argparse.ArgumentParser()
ap.add_argument("-bn", "--base_name", required=True, help="basename for all output files")
ap.add_argument("-g", "--genome_fasta_file", required=True, help="FASTA file with the genome / amplified sequence")
ap.add_argument("-r1", "--read1", required=True, help="FASTQ file first read")
ap.add_argument("-r2", "--read2", required=True, help="FASTQ file second read")
ap.add_argument("-fmh", "--find_mh", required=False, help="path to find_mh binary" , default='/lustre2/lucas_pkuhpc/bin/find_mh')
ap.add_argument("-cs", "--catch_signatures", required=False, help="cmd for catch_signatures.awk" , default='/usr/bin/gawk -v SZ=25 -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/catch_signatures.awk')
ap.add_argument("-gs", "--generate_signatures", required=False, help="cmd for generate_signatures.awk" , default='/usr/bin/gawk -v SZ=25 -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/generate_signatures.awk')
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


# file names to create
bamcramfile = args.base_name + '.bam'
probesdir = 'probes___' + os.path.basename(args.genome_fasta_file)
signatures_file = os.path.join( probesdir ,  os.path.basename(args.genome_fasta_file) + '.signatures' )

# ## ## ## start by indexing the genome and making the .cram file
print("Start time: ", time.asctime(time.localtime()))
start_time = time.perf_counter()

print('--- indexing genome --- ')
cmd = args.bwa + ' index ' + args.genome_fasta_file 
print( cmd  )
subprocess.run(  cmd  ,  shell=True , check=True)
print(' ---- ---- ')

print(' ---- aligning reads to the genome ---- ')
cmd = args.bwa + ' mem -Y ' + ' -t ' + str(args.nthreads) + ' ' + args.genome_fasta_file + ' ' + args.read1 + ' ' + args.read2 \
	+ ' | samtools view -@ 3 -F 2048 -Sb ' \
	+ ' | samtools sort -O BAM  --reference ' + args.genome_fasta_file + ' -@ 4 -m5G -o ' + bamcramfile
print( cmd  )
if os.path.isfile( bamcramfile) & (os.path.getsize(bamcramfile) > 1000) :
	print( bamcramfile + " already exists. Skipping bwa mem\n")
else:
	subprocess.run(  cmd  ,  shell=True , check=True)
	subprocess.run(  'samtools index '  + bamcramfile ,  shell=True , check=True)
end_time_1 = time.perf_counter()
print(' ---- bwa mem took ' + str( end_time_1 - start_time)   + ' seconds.      ---- \n')

# ## ##  run find_mh & generate_signatures.awk
print(' ---- find_mh & extract signatures ------ ')
subprocess.run(  'mkdir -p ' + probesdir ,  shell=True , check=True)
probefiles = os.path.join( probesdir , os.path.basename(args.genome_fasta_file) + '.mh')
cmd = args.find_mh + ' ' + probefiles + ' < ' + args.genome_fasta_file 
print( cmd  )
subprocess.run(  cmd  ,  shell=True , check=True)
list_of_mh_files = subprocess.run( 'ls ' + probefiles  +  '*' + '  | sort -t. -k1n,1 ' , shell=True , capture_output=True , text=True )
list_of_mh_files = list_of_mh_files.stdout.replace( '\n' , ' ')
cmd = args.generate_signatures + ' ' + args.genome_fasta_file + ' ' + list_of_mh_files + ' > ' + signatures_file
print( cmd )
subprocess.run(  cmd  ,  shell=True , check=True)
end_time_2 = time.perf_counter()
print( ' ---------  find_mh & extract signatures took : ' +  str( end_time_2 - end_time_1) + ' seconds -------\n')


# ## ## catch_signatures.awk
cmd = 'samtools view -T ' + args.genome_fasta_file + ' ' + bamcramfile + ' | ' + args.catch_signatures + ' ' + signatures_file + ' - > ' + args.base_name + '.sign.count.tsv'
print( cmd ) 
subprocess.run(  cmd  ,  shell=True , check=True)
print( ' --------- catch_signatures took: ' +  str( time.perf_counter() - end_time_2) + ' seconds ------\n')
print( '  TOTAL time = ' +  str( time.perf_counter() - start_time) + ' seconds ------\n')




