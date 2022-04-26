### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.4 ###

__usage__ = """
	python combine_single_VCFs_for_SnpEff.py
	--in <FULL_PATH_TO_VCF_FILE_FOLDER>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, sys

# --- end of imports --- #


def load_variants( vcf, variants ):
	"""! @brief load all variants from given VCF """
	
	ID = vcf.split('/')[-1].split('.')[0] 
	
	with open( vcf, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[6] == 'PASS':
					if parts[-1] !='./.':
						if parts[-1].split(':')[1] !='0,0':
							if len(parts[4]) == 1:
								try:
									variants[ parts[0]+'_%_'+(parts[1].zfill(8)) + '_%_' + parts[3] +  '_%_' + parts[4] ][ 'pools' ].append( str(ID) + '_' + parts[-1].split(':')[1].split(',')[0] + ':' +parts[-1].split(':')[1].split(',')[1] + ':' + str( int(parts[-1].split(':')[1].split(',')[0]) + int(parts[-1].split(':')[1].split(',')[1] ) ) + '_' + str(round((float(parts[-1].split(':')[1].split(',')[1])/float(int(parts[-1].split(':')[1].split(',')[0]) + int(parts[-1].split(':')[1].split(',')[1]) )) ,3 ))  )
								except KeyError:
									variants.update( { parts[0]+'_%_'+(parts[1].zfill(8)) + '_%_' + parts[3] +  '_%_' + parts[4] : { 'ref': parts[3], 'alt': parts[4], 'info': parts[-1], 'pools': [ str(ID) + '_%_' + parts[-1].split(':')[1].split(',')[0] + ':' +parts[-1].split(':')[1].split(',')[1] + ':' + str( int(parts[-1].split(':')[1].split(',')[0]) + int(parts[-1].split(':')[1].split(',')[1] ) ) + '_' + str( round((float(parts[-1].split(':')[1].split(',')[1])/float(int(parts[-1].split(':')[1].split(',')[0]) + int(parts[-1].split(':')[1].split(',')[1]) )), 3)) ], 'chr': parts[0], 'pos': parts[1], 'status': parts[6] } } )
									counter =+ 1		
			line = f.readline()
			
	return variants


def main( arguments ):
	"""! @brief run everything """
	
	input_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	# --- load variants from all VCF files into separate dictionaries --- #
	vcfs = sorted( glob.glob( input_folder + "*.vcf" ) )
	print vcfs
	global_variants = {}
	for vcf in vcfs:
		global_variants = load_variants( vcf, global_variants )
		
	print "number of variants: " + str( len( global_variants.keys() ) )
	
	# --- generate output file --- #
	with open( output_file, "w" ) as out:
		out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n" )
		for key in sorted( global_variants.keys() ):
			info = global_variants[ key ]
			try:
				info['pools'].append(str(float(info['pools'][0].split('_')[-1]) - float(info['pools'][1].split('_')[-1]))) 
				out.write( "\t".join( [ info['chr'], str( info['pos'] ), ".", info['ref'], info['alt'], ".", info['status'], ",".join( info['pools'] ), "GT:AD:DP:GQ:PL" ] ) + '\n' )
			except:
				out.write( "\t".join( [ info['chr'], str( info['pos'] ), ".", info['ref'], info['alt'], ".", info['status'], ",".join( info['pools'] ), "GT:AD:DP:GQ:PL" ] ) + '\n' )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
