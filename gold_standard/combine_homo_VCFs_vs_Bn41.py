### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.3 ###

__usage__ = """
	python combine_single_VCFs.py
	--P1_VCF <FULL_PATH_TO_PARENT1_VCF_FILE>
	--P2_VCF <FULL_PATH_TO_PARENT2_VCF_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: hschilbe@cebitec.uni-bielefeld.de
					"""

import glob, sys

# --- end of imports --- #


def load_variants( vcf ):
	"""! @brief load all variants from given VCF """
	
	variants = {}
	with open( vcf, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts[3] ) == 1 and len( parts[4] ) == 1:
				variants.update( { parts[0]+'_%_'+(parts[1].zfill(8)): { 'chr':parts[0], 'pos':parts[1], 'id':parts[2], 'ref': parts[3], 'alt': parts[4], 'qual': parts[5], 'filter': parts[6], 'info': parts[7], 'format': parts[8], 'dataset': parts[9] } } )
			line = f.readline()
	return variants, header


def main( arguments ):
	"""! @brief run everything """
	
	vcf1 = arguments[ arguments.index('--P1_VCF')+1 ]
	vcf2 = arguments[ arguments.index('--P2_VCF')+1 ]

	output_file = arguments[ arguments.index('--out')+1 ]
	
	# --- load variants from all VCF files into separate dictionaries --- #

	variants1, header = load_variants( vcf1 ) 
	variants2, header = load_variants( vcf2 )
	
	# --- merge variant positions of all VCF files --- #
	
	lenvar1 = len(variants1)
	lenvar2 = len(variants2)
	
	# SNP amount before deleting all same SNVs between parents
	
	print str(vcf1.split('/')[-1].split('.')[0]) + ' is ' + str(len(variants1)) 
	print str(vcf2.split('/')[-1].split('.')[0]) + ' is ' + str(len(variants2))
	
	counter = 0
	counter1 = 0
	for key in sorted(variants2.keys()): 
		try:
			if variants1[key]['ref'] == variants2[key]['ref'] and variants1[key]['alt'] == variants2[key]['alt']: 
				counter1 = counter1 + 1
				del variants1[key] #getting unique SNPs in P1
				del variants2[key] #getting unique SNPs in P2
				#these are variants which should be filtered out as they are the same within the parents - thus can not be homozygous
			if variants1[key]['ref'] == variants2[key]['ref'] and variants1[key]['alt'] != variants2[key]['alt']: #triallelic variants
				counter = counter + 1
					
		except: 
			pass
	
	print counter
	print counter1	

	print str(vcf1.split('/')[-1].split('.')[0]) + ' is ' + str(len(variants1))
	print str((len(variants1)/float(lenvar1))*100) + ' percentage are filtered out.'
	
	print str(vcf2.split('/')[-1].split('.')[0]) + ' is ' + str(len(variants2))
	print str((len(variants2)/float(lenvar2))*100) + ' percentage are filtered out.'
	
	output_file1 = str(output_file) + str(vcf1.split('/')[-1].split('.')[0]) + '_unique.vcf'
	output_file2 = str(output_file) + str(vcf2.split('/')[-1].split('.')[0]) + '_unique.vcf'
	
	with open( output_file1 , "w" ) as out1:
		out1.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )	
		for key1 in sorted(variants1.keys()):
			try:
				new_line1 = [ key1.split('_%_')[0], str( int( key1.split('_%_')[-1] ) ), variants1[key1]['id'], variants1[key1]['ref'], variants1[key1]['alt'], variants1[key1]['qual'], variants1[key1]['filter'], variants1[key1]['info'], variants1[key1]['format'], variants1[key1]['dataset']]
				out1.write( "\t".join( new_line1 ) + '\n' )
			except:
				pass

	with open( output_file2 , "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )	
		for key2 in sorted(variants2.keys()):
			try:
				new_line2 = [ key2.split('_%_')[0], str( int( key2.split('_%_')[-1] ) ), variants2[key2]['id'], variants2[key2]['ref'], variants2[key2]['alt'], variants2[key2]['qual'], variants2[key2]['filter'], variants2[key2]['info'], variants2[key2]['format'], variants2[key2]['dataset']]
				out.write( "\t".join( new_line2) + '\n' )
			except:
				pass

if __name__ == '__main__':
	
	if '--P1_VCF' in sys.argv and '--P2_VCF' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
