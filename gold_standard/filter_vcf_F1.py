### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python combine_single_VCFs.py
	--F1_VCF <FULL_PATH_TO_F1_VCF_FILE>
	--P1_VCF <FULL_PATH_TO_PARENT1_VCF_FILE>
	--P2_VCF <FULL_PATH_TO_PARENT2_VCF_FILE>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug reports and feature requests: hschilbe@cebitec.uni-bielefeld.de
					"""

import sys, os

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
	
	vcf_f1 = arguments[ arguments.index( '--F1_VCF' ) + 1 ]
	vcf1 = arguments[ arguments.index( '--P1_VCF' ) + 1 ]
	vcf2 = arguments[ arguments.index( '--P2_VCF' ) + 1 ]
	prefix = arguments[ arguments.index('--out')+1 ]
	
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	# --- load variants from all VCF files into separate dictionaries --- #

	variants, header = load_variants( vcf_f1 ) 
	variants1, header = load_variants( vcf1 ) 
	variants2, header = load_variants( vcf2 ) 
	
	lenvar = len(variants)
	lenvar1 = len(variants1)
	lenvar2 = len(variants2)
	
	print str(vcf_f1.split('/')[-1].split('.')[0]) + ' raw vars nr is ' + str(len(variants)) 
	print str(vcf1.split('/')[-1].split('.')[0]) + ' raw vars nr is ' + str(len(variants1))  
	print str(vcf2.split('/')[-1].split('.')[0]) + ' raw vars nr is ' + str(len(variants2)) 
	
	# --- filter F1 vcf file for heterozygous SNVs --- #
	
	clean_variants = {}
	for key in sorted(variants.keys()):
		if variants[key]['dataset'].split(':')[0] == '0/1': #filter for heterozygous SNVs
			try:
				if 0.2 <= float(variants[key]['dataset'].split(':')[1].split(',')[1])/float((float(variants[key]['dataset'].split(':')[1].split(',')[0]) + float(variants[key]['dataset'].split(':')[1].split(',')[1]))) <= 0.8: #alt/(ref+alt)
					clean_variants.update( { key : variants[key] }) 
			except:
				pass
				#vars with no coverage but DP
	
	print str(vcf_f1.split('/')[-1].split('.')[0]) + ' is ' + str(len(clean_variants))

	# --- merge variant positions of all VCF files --- #
	
	counter = 0
	counter2 = 0
	

	with open( prefix + str(vcf2.split('/')[-1].split('.')[0]) + '_unique_f1_filtered.vcf' , "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )
		for key in sorted(variants2.keys()):
			try:
				if clean_variants[key]['pos'] == variants2[key]['pos']: 
					counter2 = counter2 + 1	
					new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variants2[key]['id'], variants2[key]['ref'], variants2[key]['alt'], variants2[key]['qual'], variants2[key]['filter'], variants2[key]['info'], variants2[key]['format'], variants2[key]['dataset']]
					out.write( "\t".join( new_line1 ) + '\n' )
			except: 
				counter = counter + 1
				pass
				
	print 'SNVs left in parent2 ' +str(counter2) 
	print 'SNVs filtered out ' + str(counter) 

	counter3 = 0
	counter4 = 0
	

	with open( prefix + str(vcf1.split('/')[-1].split('.')[0]) + '_unique_f1_filtered.vcf' , "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )
		for key in sorted(variants1.keys()):
			try:
				if clean_variants[key]['pos'] == variants1[key]['pos']: 
					counter3 = counter3 + 1	
					new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variants1[key]['id'], variants1[key]['ref'], variants1[key]['alt'], variants1[key]['qual'], variants1[key]['filter'], variants1[key]['info'], variants1[key]['format'], variants1[key]['dataset']]
					out.write( "\t".join( new_line1 ) + '\n' )
			except: 
				counter4 = counter4 + 1
				pass
				
	print 'SNVs left in parent1 ' +str(counter3) 
	print 'SNVs filtered out ' + str(counter4) 

if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv and '--F1_VCF' in sys.argv and 'P1_VCF' in sys.argv and 'P2_VCF' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
