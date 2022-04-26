### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###


import re, os, sys
from operator import itemgetter
from Bio import SeqIO

# --- end of imports --- #

__usage__ = """
	python fetch_gene_info_from_gff3_file.py
	--in <REGIONS>
	--anno <functional_annotation_file from orthofinder results>
	--RBH_BBH_file <functional_annotation_file from (reciprocal) best blast hits>
	--gff <GFF3_file>
	--out <FULL_PATH_TO_OUTPUT_DIR>
	
	bug reports and feature requests: hschilbe@cebitec.uni-bielefeld.de

	"""


def load_gene_positions( gff_file ):
	"""! @brief load genes from GFF file to construct query file """
	
	genes = []
	counter = 0
	
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#': 
				parts = line.strip().split('\t')
				if parts[2] == 'mRNA':

					genes.append( { 'id': parts[8].split('=')[1].split(';')[0], 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] )} )
					counter = counter + 1
			line = f.readline()
			
	print 'There are ' + str(counter) + ' genes in gff3 file.'
	
	return sorted( genes, key=itemgetter( 'chr', 'start', 'end' ) )

def load_regions ( regions ):

	regions_list = []
	
	with open( regions, "r" ) as f:
		line = f.readline() 
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			regions_list.append ({ 'regionID': parts[0].split('_')[1], 'chr': parts[0], 'start': int( parts[2] ), 'end': int( parts[3] )})
			line = f.readline()
	return sorted(regions_list) 

def load_RBH_BBH_file ( RBH_BBH ):
	
	RBH_BBHs = {}
	
	with open( RBH_BBH, "r" ) as f:
		header = f.readline() 
		header_parts = header.strip().split('\t')
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			RBH_BBHs.update( { parts[0].split('=')[1]: {'status':parts[2], 'annotation':parts[4] } })
			line = f.readline()
	
	return RBH_BBHs

def load_annotation ( annotation_file ):
	"""! @brief load functional gene annotation from given file """
	
	annotation = {}
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation.update( { parts[0]: "_%_".join( parts[1:] ) } )
			line = f.readline()
	return annotation
	
def get_genes_in_regions ( genes, regions, output_file_genes_in_regions, output_file_gene_IDs, orthofinder_annotation, RBH_BBH_file ):
	
	counter1 = 0
	genes_in_region = []
	
	with open( output_file_genes_in_regions, "w" ) as out:
		with open ( output_file_gene_IDs , "w" ) as out2:
			out.write ( ('\t').join( [ 'Region_ID', 'Chromosome', 'Start_region', 'End_region', 'Start_gene', 'End_gene', 'Gene_ID', 'Annotation', '\n' ] ) )
			for region in regions:
				for gene in genes:	
					if region['chr'].split('_')[-1] == gene['chr'] or ('_').join(region['chr'].split('_')[-2:]) == gene['chr']: 
						if region['start'] <= gene['start'] <= region['end'] and region['start'] <= gene['end'] <= region['end']:
							genes_in_region.append( gene['id'] )
							try:
								out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], orthofinder_annotation[gene['id']], '\n' ] ) )
							except KeyError:
								try:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], RBH_BBH_file[gene['id']]['annotation'], '\n' ] ) )
								except KeyError:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], 'n/a', '\n' ] ) )

							counter1 = counter1 + 1 
							
						if region['start'] <= gene['start'] <= region['end'] and region['start'] <= gene['end'] > region['end']:
							genes_in_region.append( gene['id'] )
							try:
								out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], orthofinder_annotation[gene['id']], '\n' ] ) )
							except KeyError:
								try:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], RBH_BBH_file[gene['id']]['annotation'], '\n' ] ) )
								except KeyError:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], 'n/a', '\n' ] ) )

	
							counter1 = counter1 + 1 
							
						if region['start'] > gene['start'] <= region['end'] and region['start'] <= gene['end'] <= region['end']:
							genes_in_region.append( gene['id'] )
							try:
								out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], orthofinder_annotation[gene['id']], '\n' ] ) )
							except KeyError:
								try:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], RBH_BBH_file[gene['id']]['annotation'], '\n' ] ) )
								except KeyError:
									out.write ( ('\t').join( [ region['regionID'], region['chr'], str(region['start']), str(region['end']), str(gene['start']), str(gene['end']), gene['id'], 'n/a', '\n' ] ) )

							counter1 = counter1 + 1 
			
			for gene_id in genes_in_region:
				out2.write( str(gene_id) + '\n') 
	
	print 'There are ' + str(counter1) + ' genes in all analysed regions.'
	
	return genes_in_region
	

def main( arguments ):
	"""! @brief run all parts of this script """
	
	regions = arguments[ arguments.index('--in')+1 ]
	annotation_file = arguments[ arguments.index('--anno')+1 ]
	RBH_BBH_file = arguments[ arguments.index('--RBH_BBH_file')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	
	prefix = arguments[ arguments.index('--out')+1 ]
	
	genes = load_gene_positions ( gff_file )	
	
	orthofinder_annotation = load_annotation ( annotation_file )

	RBH_BBH_file = load_RBH_BBH_file ( RBH_BBH_file )

	regions = load_regions ( regions )
	print 'regions loaded.'
	print 'now checking which genes are located inside the regions'
	
	output_file_genes_in_regions = prefix + 'genes_in_regions.txt'
	output_file_gene_IDs = prefix + 'genes_IDs.txt'

	genes_in_region = get_genes_in_regions ( genes, regions, output_file_genes_in_regions, output_file_gene_IDs, orthofinder_annotation, RBH_BBH_file )


if __name__ == "__main__":
	
	if '--in' in sys.argv and '--out' in sys.argv and '--RBH_BBH_file' in sys.argv and '--anno' in sys.argv and '--gff' in sys.argv :
		main( sys.argv )
	else:
		sys.exit( __usage__ )

	
	
