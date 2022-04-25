### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python map_mean_exp_to_cand_genes_in_reg.py \
					--mean_exp_table	<MEAN_TPM_OR_FPKM_FILE>
					--anno_vars_w_reg	<FILE_ANNOTATED_VARS_W_REGION>
					--out <PATH_TO_OUTPUT_FILE>
					"""


import sys, os, re

# --- end of imports --- #

def load_mean_exp ( mean_exp_tab ):
	
	mean_exp = {}
	
	with open( mean_exp_tab, "r" ) as f:
		header = f.readline() #header Gene_ID	leaf_28DAF_JS	leaf_35DAF_Exp	leaf_35DAF_SGDH14	seeds_23DAF_Exp	seeds_23DAF_SGDH14	seeds_28DAF_JS	seeds_35DAF_Exp	seeds_35DAF_SGDH14
		header_parts = header.strip().split('\t')
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mean_exp.update( { parts[0]: {'leaf_28DAF_JS':parts[1], 'seeds_28DAF_JS':parts[6] } })
			line = f.readline()
	
	return mean_exp

def load_annotated_vars_w_region ( anno_vars_w_reg ):
	
	anno_vars = {}
	counter = 0
	
	with open( anno_vars_w_reg, "r" ) as f:
		header = f.readline() #Chromosome	Position	ReferenceAllel	AlternativeAllel	GeneID	EffectType	Pool1_coverage_AF_info	Pool2_coverage_AF_info	dAF	Annotation	Region	Number of sig SNPs in Region	Number of HIGH SNPs in Region
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno_vars.update( { parts[4].split('_')[1]+'_'+ str(counter).zfill(5): line.strip() })
			counter = counter + 1
			line = f.readline()

	return header, anno_vars

def load_genes_in_regs( genes_in_regions ):
	
	genes_in_regs = {}
	counter = 0
	
	with open( genes_in_regions, "r" ) as f:
		header = f.readline() #Chromosome	Position	ReferenceAllel	AlternativeAllel	GeneID	EffectType	Pool1_coverage_AF_info	Pool2_coverage_AF_info	dAF	Annotation	Region	Number of sig SNPs in Region	Number of HIGH SNPs in Region
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes_in_regs.update( { parts[6]+'_'+ str(counter).zfill(5): line.strip() })
			counter = counter + 1
			line = f.readline()

	return header, genes_in_regs
	
def main( arguments ):
	"""! @brief run everything """
	
	mean_exp_tab = arguments[ arguments.index('--mean_exp_table')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	mean_exp_values = load_mean_exp( mean_exp_tab )

	
	if '--anno_vars_w_reg' in arguments:
		anno_vars_w_reg = arguments[ arguments.index( '--anno_vars_w_reg' )+1 ]
		
		header, anno_vars_w_reg = load_annotated_vars_w_region( anno_vars_w_reg )

		counter2 = 0
		
		# --- generate output --- #
		with open( output_file, "w" ) as out:
			out.write( header.strip() + '\tmean_expression_seeds_28DAF_JS\tmean_expression_leaves_28DAF_JS\n' )
			for key in sorted( anno_vars_w_reg ):
				for gene_id in mean_exp_values:
					if key.split('_')[0] == gene_id:
						counter2 = counter2+ 1
						out.write( anno_vars_w_reg[key] + '\t' + str(mean_exp_values[gene_id]['seeds_28DAF_JS']) + '\t' + str(mean_exp_values[gene_id]['leaf_28DAF_JS']) + '\n')
		
		#print counter2		
		
	
	if '--genes_in_regions' in arguments:
		genes_in_regions = arguments[ arguments.index( '--genes_in_regions' )+1 ]
		
		header, genes_in_regions = load_genes_in_regs( genes_in_regions )

		counter3 = 0
		
		# --- generate output --- #
		with open( output_file, "w" ) as out:
			out.write( header.strip() + '\tmean_expression_seeds_28DAF_JS\tmean_expression_leaves_28DAF_JS\n' )
			for key in sorted( genes_in_regions ):
				for gene_id in mean_exp_values:
					if key.split('_')[0] == gene_id:
						counter3 = counter3+ 1
						out.write( genes_in_regions[key] + '\t' + str(mean_exp_values[gene_id]['seeds_28DAF_JS']) + '\t' + str(mean_exp_values[gene_id]['leaf_28DAF_JS']) + '\n')
		
		#print counter3		
	

if __name__ == '__main__':
	
	if '--mean_exp_table' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
