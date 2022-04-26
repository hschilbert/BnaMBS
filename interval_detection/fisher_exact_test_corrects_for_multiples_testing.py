### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###

### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

### v0.3 ###

__usage__ = """
					python fisher_exact_test_corrects_for_multiples_testing.py \
					--out <FULL_PATH>
					--in <FULL_PATH_TO_VCF>
					--sig <INTEGER_EXAMPLE_0.05>
					--pool1 <COMMA_SEPARATED_LIST_OF_SAMPLES>
					--pool2 <COMMA_SEPARATED_LIST_OF_SAMPLES>
					
					but reports and feature requests: hschilbe@cebitec.uni-bielefeld.de
					"""

import scipy.stats as stats
import sys, os

# --- end of imports --- #


def load_header_and_samples( in_vcf ):
	"""! @brief get header row of given VCF """
	
	lines = []
	counter = 0
	with open( in_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == "#":
				lines.append( line )
			else:
				counter += 1
			line = f.readline()
	return lines[-1].strip().split('\t'), counter


def load_variants( vcf_file, sig_cutoff, pool1_indices, pool2_indices ):
	"""! @brief load all variants from VCF file, run Fisher's exact test and save significant variants in dictionary """
	
	variants = {}
	counter_00 = 0		# positions without sufficient coverage
	counter = 0			# all variants
	counter_empty = 0	# positions without sufficient information to compute test
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				
				# --- prepare variables to save information per row (variant) --- #
				
				status_pool1 = False
				status_pool2 = False
				ref1 = 0
				alt1 = 0
				ref2 = 0
				alt2 = 0
				
				# --- collecting values of pool1 --- #
				
				for idx in pool1_indices:
					if parts[ idx ] not in [ './.', "./._L", "./._J" ]:
						status_pool1 = True
						x, y = map( int, parts[ idx ].split(':')[1].split(',') )
						ref1 += x
						alt1 += y
				
				# --- collect values of pool2 --- #
				
				for idx in pool2_indices:
					if parts[ idx ] not in [ './.', "./._L", "./._J" ]:
						status_pool2 = True
						x, y = map( int, parts[ idx ].split(':')[1].split(',') )
						ref2 += x
						alt2 += y
				
				if status_pool1 + status_pool2 == 2:
					if ref1+alt1 > 0 and ref2+alt2 > 0:
						oddsratio, pvalue = stats.fisher_exact( ( ( ref1, alt1 ), ( ref2, alt2 ) ) )
						if pvalue <= sig_cutoff:
							variants.update( { parts[0]+"_%_"+parts[1].zfill( 9 ): parts[0:7] + [ str( pvalue ) ] + parts[8:] } )
					else:
						counter_00 += 1
				else:
					counter_empty += 1
				counter += 1
			line = f.readline()

		counter = counter - counter_00 - counter_empty # total nr of SNPs without the no coverage cases
		
		print 'There are ' + str(counter_00) + ' SNPs with no coverage for at least one pool, which will not be included in the analysis.\n'
		print 'There are ' + str(counter_empty) + ' SNPs with empty line for at least one pool, which will not be included in the analysis.\n'
		print 'There are ' + str( len( variants.keys() ) ) + ' SNPs which are significant - adjusted alpha:(' + str( sig_cutoff ) + ')\n'
		print str( 1 - float( len( variants.keys() ) ) / float( counter ) ) + ' percent of SNPs ('+ str( counter - len( variants.keys() ) ) +') are filtered out, because they are not significant.\n'
	
	return variants


def main( parameters ):
	"""! @brief run all functions """
	
	prefix = parameters[ parameters.index( '--out' )+1 ]
	if prefix[-1] != '/':
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if '--sig' in parameters:
		significance_level = float( parameters[ parameters.index( '--sig' )+1 ] )
	else:
		significance_level = 0.05
	
	pool1 = parameters[ parameters.index( '--pool1' )+1 ]
	if "," in pool1:
		pool1 = pool1.split(',')
	else:
		pool1 = [ pool1 ]
	pool2 = parameters[ parameters.index( '--pool2' )+1 ]
	if "," in pool2:
		pool2 = pool2.split(',')
	else:
		pool2 = [ pool2 ]
	
	in_vcf = parameters[ parameters.index( '--in' )+1 ]	
	output_vcf = prefix + in_vcf.split('/')[-1].lower().split('.vcf')[0] + "_sig_SNPs_" + str(significance_level) + "corrected_for_mul_testing.vcf"
	
	header, sample_size = load_header_and_samples( in_vcf )
	pool1_indices = []
	pool2_indices = []
	for idx, each in enumerate( header ):
		if each in pool1:
			pool1_indices.append( idx )
		elif each in pool2:
			pool2_indices.append( idx )
	
	
	adjusted_cutoff = significance_level / sample_size
	print "detected sample size: " + str( sample_size )
	print "adjusted cutoff: " + str( adjusted_cutoff )
	variants = load_variants( in_vcf , adjusted_cutoff, pool1_indices, pool2_indices )	
	
	
	# --- prepare output file --- #
	
	with open( output_vcf, "w" ) as out:
		out.write( "\t".join( header[:7] + [ "p-value" ] + header[8:] ) + '\n' )
		for key in sorted( variants.keys() ):
			out.write( '\t'.join( variants[ key ] ) + '\n')
	print 'all done.'


if __name__ == '__main__':
	if '--out' in sys.argv and '--in' in sys.argv and '--pool1' in sys.argv and '--pool2' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
