### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python map_mean_exp_to_cand_genes_in_reg.py \
					--PAVs	<PAVs_FILE>
					--out <PATH_TO_OUTPUT_FILE>
					
					requires one of these two files:
					--anno_vars_w_reg_w_mean_exp	<FILE_ANNOTATED_VARS_W_REGION>
					--genes_in_regions <GENES_IN_REGIONS_FILE>
							
					"""

import sys, os, re

# --- end of imports --- #

def load_genes_in_regs( genes_in_regions ):
	
	genes_in_regs = {}
	counter = 0
	
	with open( genes_in_regions, "r" ) as f:
		header = f.readline() 
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes_in_regs.update( { parts[6]+'_'+ str(counter).zfill(5): line.strip() })
			counter = counter + 1
			line = f.readline()

	return header, genes_in_regs


def load_annotated_vars_w_region ( anno_vars_w_reg ):
	
	anno_vars = {}
	counter = 0
	print anno_vars_w_reg
	with open( anno_vars_w_reg, "r" ) as f:
		header = f.readline() 
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno_vars.update( { parts[4].split('_')[1]+'_'+ str(counter).zfill(5): line.strip() })
			counter = counter + 1
			line = f.readline()

	return header, anno_vars


def load_PAVs( PAVs ):
	
	PAV = {}
	counter1 = 0
	
	with open( PAVs, "r" ) as f:
		header = f.readline() 
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			PAV.update( { parts[0]+'_'+ str(counter1).zfill(5): { 'AveCovLowPool':parts[1], 'AveCovHighPool':parts[2] } } )
			counter1 = counter1 + 1
			line = f.readline()
	
	print counter1
	
	return PAV


def main( arguments ):
	"""! @brief run everything """
	
	PAVs = arguments[ arguments.index('--PAVs')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	PAV = load_PAVs( PAVs )
	
	if '--genes_in_regions' in arguments:
		genes_in_regions = arguments[ arguments.index( '--genes_in_regions' )+1 ]
		
		header, genes_in_regions = load_genes_in_regs( genes_in_regions )
		
		# --- generate output file --- #
		with open( output_file, "w" ) as out:
			out.write( header.strip() + '\tAverange_coverage_low_pool\tAverange_coverage_high_pool\n' )
			for key in sorted( genes_in_regions ):
				for gene_id in PAV:
					if key.split('_')[0] == gene_id.split('_')[0]:
						counter3 = counter3 + 1
						out.write( genes_in_regions[key] + '\t' + str(PAV[gene_id]['AveCovLowPool']) + '\t' + str(PAV[gene_id]['AveCovHighPool']) + '\n')
			
	if '--anno_vars_w_reg_w_mean_exp' in arguments:
		anno_vars_w_reg = arguments[ arguments.index( '--anno_vars_w_reg_w_mean_exp' )+1 ]
		
		header, anno_vars_w_reg = load_annotated_vars_w_region( anno_vars_w_reg )
		
		# --- generate output file --- #
		with open( output_file, "w" ) as out:
			out.write( header.strip() + '\tAverage_coverage_low_pool\tAverage_coverage_high_pool\n' )
			for key in sorted( anno_vars_w_reg ):
				for gene_id in PAV:
					if key.split('_')[0] == gene_id.split('_')[0]:
						counter2 = counter2 + 1
						out.write( anno_vars_w_reg[key] + '\t' + str(PAV[gene_id]['AveCovLowPool']) + '\t' + str(PAV[gene_id]['AveCovHighPool']) + '\n')
		

if __name__ == '__main__':
	
	if '--PAVs' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
