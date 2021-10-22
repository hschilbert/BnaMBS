### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python generate_figure.py
					--in <FULL_PATH_TO_TPM_FILE>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					--genes <FULL_PATH_TO_GENE_FILE>
					--samples <FULL_PATH_TO_SAMPLE_FILE>
					"""

import matplotlib.pyplot as plt
import sys, os
import numpy as np 

# --- end of imports --- #

def load_TPMs( input_file ):
	"""! @brief load all TPMs from given file """
	
	TPMs = {}
	with open( input_file, "r" ) as f:
		samples = []
		headers = f.readline().strip().split('\t')
		for each in headers[1:]:
			if "_" in each:
				samples.append( each.split('_')[0] )
			else:
				samples.append( each )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			for idx, part in enumerate( parts[1:] ):
				try:
					TPMs[ parts[0] ].update( { samples[ idx ]: float( part ) } )
				except KeyError:
					TPMs.update( { parts[0]: { samples[ idx ]: float( part ) } } )
			line = f.readline()
	#print TPMs'chrC09-g83340': {'EsB35': 11.2935084982, 'ES35': 12.9809061688, 'ES23': 18.6847459001, 'SsS23': 13.2976871459, 'SSxSsS23': 18.7094005185, 'EB35': 5.09400008581, 'SsS': 19.1555441024, 'EBxEsB35': 4.76967975286, 'SS23': 22.9702966779, 'JSS228': 12.866219812, 'JSB3': 10.4127014609, 'JSB2': 11.5141679466, 'JSB1': 10.5356160501, 'SSxSsS': 15.6138427692, 'SsBxSB3': 5.65340997765, 'EsS35': 16.4265161585, 'EsS23': 17.5966963365, 'SS': 10.9409854187, 'JSS1': 17.2676691425, 'JSS2': 14.321625699, 'ESxEsS35': 16.7839698318, 'SB35': 5.55324034432, 'ESxEsS23': 17.7536731817, 'SsB35': 4.08367373229}
	
	return TPMs


def load_sample_groups( sample_group_file ):
	"""! @brief define sample groups based on input """
	
	sample_groups = {}
	with open( sample_group_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			#print parts
			try:
				sample_groups[ parts[1] ].append( parts[0] )
			except KeyError:
				sample_groups.update( { parts[1]: [ parts[0] ] } )
			line = f.readline()
	return sample_groups


def load_candidate_genes( candidate_gene_file ):
	"""! @brief load all candidate genes from given file """
	
	with open( candidate_gene_file, "r" ) as f:
		content = f.read().strip()
	return content.split('\n')


def generate_figure( data, samples, genes, figfile, key, means):
	"""! @brief generate boxplot for given groups of samples and genes """
	
	#fig, ax = plt.subplots()
	#new_xlabel = []
	values = []
	
	#print genes
	#print data
	
	for gene in genes:
		vals = []
		for sample in samples:
			#print samples
			try:
				vals.append( data[ gene ][ sample ] )
				#print gene
				#print sample
				#print vals
			except KeyError:
				pass
		values.append( vals )
		means.append( { 'gene': str(gene), 'mean': np.mean(vals) }  )
		# in order to shorten the gene IDs for x label #
		#new_xlabel.append(gene.split('_')[0])
	#print values  #[0.0, 0.20922154417, 0.307863361312]] R454660_chrA09p000010.1_BnaEXP SSxSsS23 [0.0, 0.20922154417, 0.307863361312]
	
	# ax.boxplot( values, showmeans=True )
	# ax.set_title( "n="+str( len( values[0] ) ) )
	# #print genes
	# ax.set_xticklabels( new_xlabel, rotation=90 )
	# ax.set_ylabel( "FPKMs" )
	# plt.subplots_adjust( left=0.1, bottom=0.3, right=0.99, top=0.9 )
	
	# fig.savefig( figfile, dpi=900 )
	# plt.close("all")

	return means

def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	candidate_gene_file = arguments[ arguments.index( '--genes' )+1 ]
	
	sample_group_file = arguments[ arguments.index( '--samples' )+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- load data --- #
	TPMs = load_TPMs( input_file )
	sample_groups = load_sample_groups( sample_group_file )
	candidate_genes = load_candidate_genes( candidate_gene_file )

	# --- generate figures --- #
	means = []
	mean_calc = {}
	for key in sample_groups.keys():
		#print key #tissues
		group = sample_groups[ key ]
		figfile = output_folder + key + ".png"
		means = generate_figure( TPMs, group, candidate_genes, figfile, key, means )
		mean_calc.update( {str(key): means } )
		#print mean_calc
		means = []
	
	''' out write mean values per gene per tissue '''
	#print mean_calc
	
	print mean_calc.keys()
	means_gene = []
	
	out_file = output_folder + 'means.txt'
	
	#means.append( { 'gene': str(gene), 'mean': np.mean(vals) }  )
	#candidate_genes
	#{ ,leav35_DAF: [ { 'gene': str(gene), 'mean': np.mean(vals) } , { 'gene': str(gene), 'mean': np.mean(vals) } ], }
	
	with open ( out_file , 'w' ) as out:
		out.write ('Gene_ID\t'+str(('\t').join(sorted(mean_calc.keys()))) + '\n')
		for gene in candidate_genes: #genes
			values_per_gene = [gene]
			for key in sorted(mean_calc.keys()): #tissues
				gene_exp_in_tissue = mean_calc[key]
				for entry in gene_exp_in_tissue : #means and genes
					#print value[0]['gene']
					#print key
					#print gene
					if entry['gene'] == gene:
						values_per_gene.append( str(entry['mean']) ) #mean_calc[key][0]['mean']
	
			out.write( ('\t').join( values_per_gene ) + '\n')
			


if '--in' in sys.argv and '--out' in sys.argv and '--genes' in sys.argv and '--samples' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
