### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python filter_pools_vcfs_for_gold_standard.py
	--out <FULL_PATH_TO_OUTPUT_FOLDER>
	
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
	
	prefix = arguments[ arguments.index('--out')+1 ]
	
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	# --- load variants from all VCF files into separate dictionaries --- #
	#vcf_hp = './vs_Bn41/NPZ/NPZ_high.vcf' 	
	#vcf_lp = './vs_Bn41/NPZ/NPZ_low.vcf' 
	
	#vcf_hp = './vs_Bn41/DSV/DSV_high.vcf' 
	#vcf_lp = './vs_Bn41/DSV/DSV_low.vcf' 
	
	vcf_hp = './G_high.vcf' 
	vcf_lp = './G_low.vcf' 
	
	vcfJ = './20200302_Janetzki_vs_Bn41_homo_variants_10_60_10_60_unique_unique_f1_filtered.vcf'
	vcfL = './20200302_Lorenz_vs_Bn41_homo_variants_10_60_10_60_unique_unique_f1_filtered.vcf'
	
	variantsHP, header = load_variants( vcf_hp ) 
	variantsLP, header = load_variants( vcf_lp ) 
	variantsJ, header = load_variants( vcfJ ) #JS
	variantsL, header = load_variants( vcfL ) #Lorenz

	raw_vars_HP = len(variantsHP) 
	raw_vars_LP = len(variantsLP) 
	raw_vars_J = len(variantsJ) 
	raw_vars_L = len(variantsL)	

	print str(vcf_hp.split('/')[-1].split('.')[0]) + ' raw SNPs nr is ' + str(len(variantsHP)) 
	print str(vcf_lp.split('/')[-1].split('.')[0]) + ' raw SNPs nr is ' + str(len(variantsLP)) 
	print str(vcfJ.split('/')[-1].split('.')[0]) + ' raw SNPs nr is ' + str(len(variantsJ)) 
	print str(vcfL.split('/')[-1].split('.')[0]) + ' raw SNPs nr is ' + str(len(variantsL))	

	# --- merge variant positions of all VCF files --- #
	
	counter = 0
	counter2 = 0
	
	with open( prefix + str(vcf_hp.split('/')[-1].split('.')[0]) + '_unique_f1_filtered.vcf' , "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )
		for key in sorted(variantsHP.keys()):
			try:
				if variantsHP[key]['pos'] == variantsJ[key]['pos']:
					variantsHP[key]['dataset'] = variantsHP[key]['dataset']+'_J'
					new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variantsHP[key]['id'], variantsHP[key]['ref'], variantsHP[key]['alt'], variantsHP[key]['qual'], variantsHP[key]['filter'], variantsHP[key]['info'], variantsHP[key]['format'], variantsHP[key]['dataset']]
					out.write( "\t".join( new_line1 ) + '\n' )
			except KeyError: 
				counter = counter + 1
				
				try:	
					if variantsHP[key]['pos'] == variantsL[key]['pos']: 
						variantsHP[key]['dataset'] = variantsHP[key]['dataset']+'_L'
						new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variantsHP[key]['id'], variantsHP[key]['ref'], variantsHP[key]['alt'], variantsHP[key]['qual'], variantsHP[key]['filter'], variantsHP[key]['info'], variantsHP[key]['format'], variantsHP[key]['dataset']]
						out.write( "\t".join( new_line1 ) + '\n' )
				except:
					counter2 = counter2 + 1
					pass
				
	print counter 
	print counter2 
	
	#nr of SNV which are found in HP from JS
	snv_in_HP_from_J = raw_vars_HP-counter
	print str(snv_in_HP_from_J)
	#percentage of SNVs of J which are also found in the pools
	per_snv_left_J = snv_in_HP_from_J/float(raw_vars_J)
	print str(per_snv_left_J)
	#nr of SNV which are found in LP from JS
	snv_in_LP_from_L = counter-counter2
	print str(snv_in_LP_from_L)
	#percentage of SNVs of LO which are also found in the pools
	per_snv_left_L = snv_in_LP_from_L/float(raw_vars_L)
	print str(per_snv_left_L)

	counter3 = 0
	counter4 = 0
	
	with open( prefix + str(vcf_lp.split('/')[-1].split('.')[0]) + '_unique_f1_filtered.vcf' , "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\n' )
		for key in sorted(variantsLP.keys()):
			try:
				if variantsLP[key]['pos'] == variantsJ[key]['pos']:
					variantsLP[key]['dataset'] = variantsLP[key]['dataset']+'_J'
					new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variantsLP[key]['id'], variantsLP[key]['ref'], variantsLP[key]['alt'], variantsLP[key]['qual'], variantsLP[key]['filter'], variantsLP[key]['info'], variantsLP[key]['format'], variantsLP[key]['dataset']]
					out.write( "\t".join( new_line1 ) + '\n' )
			except KeyError: 
				counter3 = counter3 + 1
				
				try:	
					if variantsLP[key]['pos'] == variantsL[key]['pos']: 
						variantsLP[key]['dataset'] = variantsLP[key]['dataset']+'_L'
						new_line1 = [ key.split('_%_')[0], str( int( key.split('_%_')[-1] ) ), variantsLP[key]['id'], variantsLP[key]['ref'], variantsLP[key]['alt'], variantsLP[key]['qual'], variantsLP[key]['filter'], variantsLP[key]['info'], variantsLP[key]['format'], variantsLP[key]['dataset']]
						out.write( "\t".join( new_line1 ) + '\n' )
				except:
					counter4 = counter4 + 1
					pass
				
	print counter3 
	print counter4 
	
	#nr of SNV which are found in HP from JS
	snv_in_LP_from_J = raw_vars_LP-counter3
	print str(snv_in_LP_from_J)
	#percentage of SNVs of J which are also found in the pools
	per_snv_left_J = snv_in_LP_from_J/float(raw_vars_J)
	print str(per_snv_left_J)
	#nr of SNV which are found in LP from JS
	snv_in_LP_from_L = counter3-counter4
	print str(snv_in_LP_from_L)
	#percentage of SNVs of LO which are also found in the pools
	per_snv_left_L = snv_in_LP_from_L/float(raw_vars_L)
	print str(per_snv_left_L)

if __name__ == '__main__':
	
	if '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
