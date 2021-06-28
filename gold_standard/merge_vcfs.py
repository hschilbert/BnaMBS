### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python merge_vcfs.py
					--in <FULL_PATH_TO_INPUT_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_VCF>
					--fasta <REFERENCE_FASTA_FILE>
					
					optional:
					--sort_script <FULL_PATH_TO_SORT_SCRIPT>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys, glob

# --- end of imports --- #

def load_all_variants( vcf_file, variants ):
	"""! @brief load all variants from VCF file (GATK default expected) """
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					variants[ parts[0] + '_%_' + parts[1].zfill( 9 ) ]
				except KeyError:
					variants.update( { parts[0] + '_%_' + parts[1].zfill( 9 ): line } )
			line = f.readline()
	return variants


def main( arguments ):
	"""! @brief runs everything """
	
	output_file = arguments[ arguments.index( '--out' )+1 ]
	input_folder = arguments[ arguments.index( '--in' )+1 ]
	sort_script = arguments[ arguments.index( '--sort_script' )+1 ]
	fasta_file = arguments[ arguments.index( '--fasta' )+1 ]
	
	if '--sort_script' in arguments:
		sort_script = arguments[ arguments.index( '--sort_script' )+1 ]
	else:
		sort_script = "sort_vcf_by_fasta.py"
	
	vcfs = glob.glob( input_folder + "*.vcf" )
	print "number of VCFs: " + str( len( vcfs ) )
	variants = {}
	for vcf in vcfs:
		variants = load_all_variants( vcf, variants )
	
	tmp_file = output_file + ".tmp"
	counter = 0
	with open( tmp_file, "w" ) as out:
		out.write( "##fileformat=VCFv4.2\n" )
		out.write( "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DATA_SET\n" )
		for key in sorted( variants.keys() ):
			out.write( variants[ key ] )
			counter += 1
	
	os.popen( "python " + sort_script + " --vcf " + tmp_file + " --fasta " + fasta_file + " --output " + output_file )
	print "number of variants: " + str( counter )


if '--in' in sys.argv and '--out' in sys.argv and '--fasta' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

