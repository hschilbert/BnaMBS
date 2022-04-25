### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v3 ###

# this version can now also include the total amount of vars present in the original/unfiltered pool vcfs file
# this is to visualize regions with no vars present
# this version will now also include a coverage plot under the x axis

__usage__ = """
					python dAF_with_unique_and_sig_SNP_v3.py\n
					--input_vcf <FILENAME>
					--input_vcf_sig_SNP <FILENAME>
					--in_merged_ori_vcf <FILENAME>
					--reference_file <FILENAME>
					--output_dir <DIRECTORY_NAME>[will be generated if required]
					--high_pool <sample name in VCF; multiple samples names can be provided comma-seperated>
					--low_pool <sample name in VCF; multiple samples names can be provided comma-seperated>
					
					
					bug reports and feature requests: hschilbe@cebitec.uni-bielefeld.de
					"""
					
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import sys, os
import numpy as np
from matplotlib.lines import Line2D

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


def load_variants( vcf_file, pool1_indices, pool2_indices ):
	"""! @brief load all variants from VCF file """
	
	variants = []
	counter_00 = 0		#positions without sufficient coverage
	counter = 0			#all variants
	counter_empty = 0	#positions without sufficient information to compute test

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
						# oddsratio, pvalue = stats.fisher_exact( ( ( ref1, alt1 ), ( ref2, alt2 ) ) )
						# we will not calc fisher exact test here, because it is time consuming and we only want to get the SNPs which were used for this test 
						# thus instead of a p value a zero is appended
						
						variants.append( parts[0:7] + [ str( 0 ) ] + parts[8:] )

					else:
						counter_00 += 1
				else:
					counter_empty += 1
				counter += 1
			line = f.readline()

		counter = counter - counter_00 - counter_empty # total nr of SNPs without the no coverage cases
		
		print 'There are ' + str(counter_00) + ' SNPs with no coverage for at least one pool, which will not be included in the analysis.\n'
		print 'There are ' + str(counter_empty) + ' SNPs with empty line for at least one pool, which will not be included in the analysis.\n'
		
		return variants


def get_coverage( input_vcf, high_name, low_name ):
	"""! @brief get average coverage for both samples of interest """
	
	coverage_high = []
	coverage_low = []	
	print input_vcf
	with open( input_vcf, "r" ) as f:
		line = f.readline()
		print line
		#break
		while line:
			if line[0] != '#':

				parts = line.strip().split('\t')
				if parts[6] == "PASS":
					sample_data = parts[9:]
					
					tmp_high = []
					for index in high_sample_idx:
						hd = sample_data[ index ]	#variant data for high sample
						if hd[:3] != "./.":
							hd_parts = hd.split(':')
							tmp_high.append( int( hd_parts[2] ) )
					if tmp_high > 0:
						coverage_high.append( sum( tmp_high ) )
					
					tmp_low = []
					for index in low_sample_idx:
						ld = sample_data[ index ]	#variant data for low sample
						if ld[:3] != "./.":
							ld_parts = ld.split(':')
							tmp_low.append( int( ld_parts[2] ) )
					if tmp_low > 0:
						coverage_low.append( sum( tmp_low ) )
			else:
				try:
					samples = line.strip().split('\t')[9:]
					high_sample_idx = []
					for each in high_name:
						try:
							high_sample_idx.append( samples.index( each ) )
						except ValueError:
							print "ERROR: sample not detected - " + each
					low_sample_idx = []
					for each in low_name:
						try:
							low_sample_idx.append( samples.index( each ) )
						except ValueError:
							print "ERROR: sample not detected - " + each
				except IndexError:
					pass
			line = f.readline()
	
	#plt.hist( coverage_high, bins = 20000 )
	#plt.title( "coverage high" )
	#plt.xlim( [ 0, 200 ] )
	#plt.show()
	
	#plt.hist( coverage_low, bins = 20000 )
	#plt.title( "coverage low" )
	#plt.xlim( [ 0, 200 ] )
	#plt.show()
	
	return sorted( coverage_high )[ len( coverage_high ) / 2 ], sorted( coverage_low )[ len( coverage_low ) / 2 ]


def get_delta_allel_frequencies( input_vcf, output_file, high_name, low_name, min_high_cov, max_high_cov, min_low_cov, max_low_cov):
	"""! @brief calculate all possible delta allele frequencies
	@note reference allele frequency is devided by alternative allele frequency
	@note triallelic variants (and higher) are skipped
	"""
	
	print 'The min_high_cov for pool1 filtering is  ' + str(min_high_cov)
	print 'The max_high_cov for pool1 filtering is  ' + str(max_high_cov)
	print 'The min_high_cov for pool2 filtering is  ' + str(min_high_cov)
	print 'The max_high_cov for pool2 filtering is  ' + str(max_high_cov)

	
	with open( output_file, "w" ) as out:
		with open( input_vcf, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[6] == "PASS" and not "," in parts[4]:	#avoid triallelic variants:
						sample_data = parts[9:]
						
						# --- combine information of all pool1 samples --- #
						h_cov = 0
						h_allele1 = 0
						h_allele2 = 0
						for each in high_sample_idx:
							hd = sample_data[ each ]	#variant data for high sample
							hd_parts = hd.split(':')
							try:
								try:
									h_cov += int( hd_parts[2] )
								except IndexError:
									pass
							except ValueError:
								pass
							if hd[:3] != "./.":
								h_allele1 += int( hd_parts[1].split(',')[0] )
								h_allele2 += int( hd_parts[1].split(',')[1] )
						
						# --- combine information of all pool2 samples --- #
						l_cov = 0
						l_allele1 = 0
						l_allele2 = 0
						for each in low_sample_idx:
							ld = sample_data[ each ]	#variant data for low sample
							ld_parts = ld.split(':')
							try:
								try:
									l_cov += int( ld_parts[2] )
								except IndexError:
									pass
							except ValueError:
								pass
							if ld[:3] != "./.":	
								l_allele1 += int( ld_parts[1].split(',')[0] )
								l_allele2 += int( ld_parts[1].split(',')[1] )
						
						# --- check if coverage of locus is within expected range ---- #
						if h_cov >= min_high_cov and h_cov <= max_high_cov and l_cov >= min_low_cov and l_cov <= max_low_cov:
							af_high_status = False
							try:
								af_high =  ( float( h_allele2 ) / float( h_allele1+h_allele2 ) )
							except ZeroDivisionError:
								af_high = 1
								af_high_status = True
							
							af_low_status = False
							try:
								af_low =  ( float( l_allele2 ) / float( l_allele1+l_allele2 ) )
							except ZeroDivisionError:
								af_low = 1
								af_low_status = True
							
							if af_high_status + af_low_status < 2:
								af_delta = af_high - af_low
								out.write( "\t".join( parts[:7] +  map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, af_delta ] ) ) + '\n' )
						elif h_cov >= min_high_cov and h_cov <= max_high_cov and l_cov == 0:
							out.write( "\t".join( parts[:7] + map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, 1 ] ) ) + '\n' )
						elif l_cov >= min_low_cov and l_cov <= max_low_cov and h_cov == 0:
							out.write( "\t".join( parts[:7] + map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, -1 ] ) ) + '\n' )
				else:
					try:
						try:
							samples = line.strip().split('\t')[9:]
							high_sample_idx = []
							for each in high_name:
								high_sample_idx.append( samples.index( each ) )
							low_sample_idx = []
							for each in low_name:
								low_sample_idx.append( samples.index( each ) )
							out.write( "\t".join( line.strip().split('\t')[:9] + [ "Pool1Coverage\tPool1RefCov\tPool1AltCov\tPool2Coverage\tPool2RefCov\tPool2AltCov\tdelta_AF" ] ) + '\n' )
						except ValueError:
							pass
					except IndexError:
						pass
				line = f.readline()


def plot_genome_wide_delta_allele_frequencies( coverage_threshold_to_plot, input_mapping_cov_high, input_mapping_cov_low, af_frequency_vcf, af_frequency_vcf_sig, af_frequency_vcf_merged_ori, chr_lengths, window_size, step_size, window_size_density, step_size_density, output_dir ):
	"""! @brief show genome wide distribution of AF frequencies """
	
	print window_size
	print chr_lengths
	
	# --- load information about chromosomes --- #
	chr_names = sorted( chr_lengths.keys() )
	raw_data_x = [ ]
	raw_data_y = [ ]
	af_data_y_pool1 = []	#same X values can be used for all
	af_data_y_pool2 = []
	
	raw_data_x_sig = [ ]
	raw_data_y_sig = [ ]
	af_data_y_pool1_sig = []	#same X values can be used for all
	af_data_y_pool2_sig = []
	
	for  k in chr_names:
		raw_data_x.append( [] )
		raw_data_y.append( [] )
		af_data_y_pool1.append( [] )
		af_data_y_pool2.append( [] )
		
	for  k in chr_names:
		raw_data_x_sig.append( [] )
		raw_data_y_sig.append( [] )
		af_data_y_pool1_sig.append( [] )
		af_data_y_pool2_sig.append( [] )
	
	# --- load normal variant information (variants detected in both pools) --- #
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline() #header
		line = f.readline()
		print line
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')

				raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] ) )				#position
				raw_data_y[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )	#dAF
				try:
					af_data_y_pool1[ chr_names.index( parts[0] ) ].append( float( parts[10] ) / float( parts[9] ) )
				except ZeroDivisionError:
					af_data_y_pool1[ chr_names.index( parts[0] ) ].append( 1 )
				try:
					af_data_y_pool2[ chr_names.index( parts[0] ) ].append( float( parts[13] ) / float( parts[12] ) )
				except ZeroDivisionError:
					af_data_y_pool2[ chr_names.index( parts[0] ) ].append( 1 )
			
			line = f.readline()
	print raw_data_x[0][0:10]

	# --- load sig variant information (variants detected in both pools) --- #
	with open( af_frequency_vcf_sig, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				# position for x axis
				raw_data_x_sig[ chr_names.index( parts[0] ) ].append( int( parts[1] ) )
				# position for y axis
				raw_data_y_sig[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
				try:
					af_data_y_pool1_sig[ chr_names.index( parts[0] ) ].append( float( parts[10] ) / float( parts[9] ) )
				except ZeroDivisionError: #if AF is 0 in one pool
					af_data_y_pool1_sig[ chr_names.index( parts[0] ) ].append( 1 )
				try:
					af_data_y_pool2_sig[ chr_names.index( parts[0] ) ].append( float( parts[13] ) / float( parts[12] ) )
				except ZeroDivisionError: #if AF is 0 in one pool
					af_data_y_pool2_sig[ chr_names.index( parts[0] ) ].append( 1 )
			line = f.readline()
	
	# --- load merged original vcfs vars variant information (variants detected in both pools) --- #
	raw_data_x_merged_ori = [ ]
	raw_data_y_merged_ori = [ ]
	#af_data_y_pool1_merged_ori = []	#same X values can be used for all
	#af_data_y_pool2_merged_ori = []
	
	for  k in chr_names:
		raw_data_x_merged_ori.append( [] )
		raw_data_y_merged_ori.append( [] )
		#af_data_y_pool1_merged_ori.append( [] )
		#af_data_y_pool2_merged_ori.append( [] )
	
	with open( af_frequency_vcf_merged_ori, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				raw_data_x_merged_ori[ chr_names.index( parts[0] ) ].append( int( parts[1] ) ) #position
			line = f.readline()


	# --- load coverage files information --- #
	raw_data_y_cov_high = [ ]
	raw_data_x_cov_high = [ ] 
	raw_data_y_cov_low = [ ] 
	raw_data_x_cov_low = [ ] 
	
	for  k in chr_names:
		raw_data_y_cov_high.append( [] )
		raw_data_y_cov_low.append( [] )
		raw_data_x_cov_low.append( [] )
		raw_data_x_cov_high.append( [] )

	with open( input_mapping_cov_high, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split(',')
				if float( parts[3] ) >= coverage_threshold_to_plot:
					raw_data_y_cov_high[ chr_names.index( parts[0] ) ].append( (float( coverage_threshold_to_plot )+1) ) #average coverage in window
					raw_data_x_cov_high[ chr_names.index( parts[0] ) ].append( (float( parts[1] ) + (window_size_density/2.0) )/ float( 1000000  ) ) #position in window
				else:
					raw_data_y_cov_high[ chr_names.index( parts[0] ) ].append( float( parts[3] ) ) #average coverage in window
					raw_data_x_cov_high[ chr_names.index( parts[0] ) ].append( (float( parts[1] ) + (window_size_density/2.0) )/ float( 1000000  ) ) #position in window
			line = f.readline()
			
	with open( input_mapping_cov_low, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split(',')
				if float( parts[3] ) >= coverage_threshold_to_plot:
					raw_data_y_cov_low[ chr_names.index( parts[0] ) ].append( (float( coverage_threshold_to_plot )+1)) #average coverage in window
					raw_data_x_cov_low[ chr_names.index( parts[0] ) ].append( (float( parts[1] ) + (window_size_density/2.0) )/ float( 1000000  ) ) #position in window

				else:
					raw_data_y_cov_low[ chr_names.index( parts[0] ) ].append( float( parts[3] ) ) #average coverage in window
					raw_data_x_cov_low[ chr_names.index( parts[0] ) ].append( (float( parts[1] ) + (window_size_density/2.0) ) / float( 1000000 ) ) #position in window
			line = f.readline()
	
	
	# --- prepare normal data (variants detected in both pools) for construction of plot --- #
	data_to_plot_x = [ ]
	data_to_plot_y = [ ]
	af_pool1 = []
	af_pool2 = []
	intervalls = []
	
	for k in chr_names:
		data_to_plot_x.append( [] )
		data_to_plot_y.append( [] )
		intervalls.append( [] )
		af_pool1.append( [] )
		af_pool2.append( [] )
	
	for idx, chr_data in enumerate( raw_data_x ): #iterate over positions of SNVs per chr
		
		start = 0
		end = 0 + window_size 
		
		while end < len( chr_data ):	
			intervalls[ idx ].append( ( chr_data[ start ], chr_data[ end ] ) )	
			
			if start > len( chr_data )-window_size:	#end analysis of this chromosome
				try:
					tmp_x = chr_data[ start: ]
					tmp_y = raw_data_y[ idx ][ start: ] 
					tmp_y_pool1 = af_data_y_pool1[ idx ][ start: ] 
					tmp_y_pool2 = af_data_y_pool2[ idx ][ start: ] 
					data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
					data_to_plot_y[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
					af_pool1[ idx ].append( np.median( tmp_y_pool1 ) )
					af_pool2[ idx ].append( np.median( tmp_y_pool2 ) )
				except ZeroDivisionError:
					data_to_plot_x[ idx ].append( 0 )
					data_to_plot_y[ idx ].append( 0 )
				break
				
			else:
				tmp_x = chr_data[ start:end ]
				tmp_y = raw_data_y[ idx ][ start:end ] 
				tmp_y_pool1 = af_data_y_pool1[ idx ][ start:end ] 
				tmp_y_pool2 = af_data_y_pool2[ idx ][ start:end ] 
				data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
				data_to_plot_y[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
				af_pool1[ idx ].append( np.median( tmp_y_pool1 ) )
				af_pool2[ idx ].append( np.median( tmp_y_pool2 ) )
				start += step_size
				end += step_size
	
	# --- prepare sig SNVs data for construction of plot --- #

	data_to_plot_x_sig = [ ]
	data_to_plot_y_sig = [ ]
	af_pool1_sig = []
	af_pool2_sig = []
	intervalls_sig = []
	
	for k in chr_names:
		data_to_plot_x_sig.append( [] )
		data_to_plot_y_sig.append( [] )
		intervalls_sig.append( [] )
		af_pool1_sig.append( [] )
		af_pool2_sig.append( [] )
	
	for idx, chr_data in enumerate( raw_data_x_sig ):
		start = 0
		end = 0 + window_size
		while end < len( chr_data ):	
			intervalls_sig[ idx ].append( ( chr_data[ start ], chr_data[ end ] ) )	
			if start > len( chr_data )-window_size:	#end analysis of this chromosome
				try:
					tmp_x = chr_data[ start: ]
					tmp_y = raw_data_y_sig[ idx ][ start: ] 
					tmp_y_pool1 = af_data_y_pool1_sig[ idx ][ start: ] 
					tmp_y_pool2 = af_data_y_pool2_sig[ idx ][ start: ] 
					data_to_plot_x_sig[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
					data_to_plot_y_sig[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
					af_pool1_sig[ idx ].append( np.median( tmp_y_pool1 ) )
					af_pool2_sig[ idx ].append( np.median( tmp_y_pool2 ) )
				except ZeroDivisionError:
					data_to_plot_x_sig[ idx ].append( 0 )
					data_to_plot_y_sig[ idx ].append( 0 )
				break
				
			else:
				tmp_x = chr_data[ start:end ]
				tmp_y = raw_data_y_sig[ idx ][ start:end ] 
				tmp_y_pool1 = af_data_y_pool1_sig[ idx ][ start:end ] 
				tmp_y_pool2 = af_data_y_pool2_sig[ idx ][ start:end ] 
				data_to_plot_x_sig[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
				data_to_plot_y_sig[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
				af_pool1_sig[ idx ].append( np.median( tmp_y_pool1 ) )
				af_pool2_sig[ idx ].append( np.median( tmp_y_pool2 ) )
				start += step_size
				end += step_size

	# --- for density plot of sig SNVs --- #
	
	raw_data_y_density = []		
	data_to_plot_x_density = [ ]
	data_to_plot_y_density = [ ]
	intervalls_density = []
	
	for k in chr_names:
		raw_data_y_density.append( [] )
		data_to_plot_x_density.append( [] )
		data_to_plot_y_density.append( [] )
		intervalls_density.append( [] )

	tmp_x = []

	for idx, chr_data in enumerate( raw_data_x_sig ): #iterate over positions of dARCs per chr
		start = 0
		end = 0 + window_size_density 
		if len(chr_data)>0:
			while start < chr_lengths[chr_names[idx]]: #max(chr_data)	is the highest position value in the list, alias the end position of the chromosome
				# get position of first and last SNVs in this window #
				intervalls_density[ idx ].append( ( start, end ) )	
				
				if  start > chr_lengths[chr_names[idx]] - window_size_density :	#end analysis of this chromosome
					try:
						for ele in chr_data:
							if  start <= ele :
								tmp_x.append( ele )
						data_to_plot_x_density[ idx ].append( (start + (window_size_density/2.0)  ) /1000000.0 )
						data_to_plot_y_density[ idx ].append( len( tmp_x ) )

					except ZeroDivisionError:
						data_to_plot_x_density[ idx ].append( 0 )
						data_to_plot_y_density[ idx ].append( 0 )
					break
					
				else:
					# geht all positions of SNVs between start and end position of current window
					for ele in chr_data:
						if  int(start) <= int(ele) <= int(end) :
							tmp_x.append( int(ele ))
					# here * 1000000.0 is calculated because Mbp will be plotted
					try:
						data_to_plot_x_density[ idx ].append( (start + (window_size_density/2.0) ) /1000000.0)
						# number of sig SNVs in that window is len(tmp_x)
						data_to_plot_y_density[ idx ].append( len( tmp_x ) )		#sum( tmp_y ) / float( len( tmp_y ) )
					
					except ZeroDivisionError:
						data_to_plot_x_density[ idx ].append( 0 )
						data_to_plot_y_density[ idx ].append( 0 )
					
				tmp_x = []
				start += step_size_density
				end += step_size_density
				
		else:
			continue
	
	# --- for density plot of merged ori SNVs --- #
	
	raw_data_y_density_merged_ori = []		
	#data_to_plot_x_density_merged_ori = [ ] is not nessesary since data_to_plot_x_density can be used as well
	data_to_plot_y_density_merged_ori = [ ]
	intervalls_density_merged_ori = []
	
	for k in chr_names:
		raw_data_y_density_merged_ori.append( [] )
		data_to_plot_y_density_merged_ori.append( [] )
		intervalls_density_merged_ori.append( [] )

	tmp_x_merged_ori = []

	for idx, chr_data in enumerate( raw_data_x_merged_ori ): #iterate over positions of SNVs per chr
		start = 0
		end = 0 + window_size_density 
		if len(chr_data)>0:
			while start < chr_lengths[chr_names[idx]]: #end < max( chr_data ): #max(chr_data)	
				intervalls_density_merged_ori[ idx ].append( ( start, end ) )	
				if start > chr_lengths[chr_names[idx]] - window_size_density: #start > max( chr_data )-window_size_density :	#end analysis of this chromosome
					try:
						for ele in chr_data:
							if  start <= ele :
								tmp_x_merged_ori.append( ele )
						data_to_plot_y_density_merged_ori[ idx ].append( len( tmp_x_merged_ori ) )

					except ZeroDivisionError:
						data_to_plot_y_density_merged_ori[ idx ].append( 0 )
					break
					
				else:
					# geht all positions of SNVs between start and end position of current window
					for ele in chr_data:
						if  int(start) <= int(ele) <= int(end) :
							tmp_x_merged_ori.append( int(ele ))
					try:
						data_to_plot_y_density_merged_ori[ idx ].append( len( tmp_x_merged_ori ) )		#sum( tmp_y ) / float( len( tmp_y ) )
					
					except ZeroDivisionError:
						data_to_plot_x_density_merged_ori[ idx ].append( 0 )
						data_to_plot_y_density_merged_ori[ idx ].append( 0 )
					
					tmp_x_merged_ori = []
					start += step_size_density
					end += step_size_density
		else:
			continue

	# --- construct single chromosome plots --- #
	
	output_density_dARCs = output_dir + "dARCs_density.txt"
	with open( output_density_dARCs, "w" ) as out:
		out.write( "Chr\tposition\tnorm_dARC_value\tunnorm_dARCs\n" )

		for idx, each in enumerate( chr_names ):
			fig_output_file = af_frequency_vcf + "." + str( window_size ) + "." + each + ".genome_wide.png"
			fig, ax = plt.subplots(  figsize=(20, 5) )
			
			# --- set whether top x spine and right spine should b visible and create a box --- #
			
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			
			# --- adding chromosomes --- #
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0, 0 ] , color="black", linewidth=.5 ) 		#plotting x axis
					
			# --- add helper lines --- #
			
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.25, 0.25 ] , color="black", linewidth=.1 )
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.5, 0.5 ] , color="black", linewidth=.1 )
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.75, 0.75 ] , color="black", linewidth=.1 )
			
			# --- adding more axis --- #
			
			ax2 = ax.twinx() #for raw desity of sig SNVs
			ax3 = ax.twinx() #for normalized desity of sig SNVs
			#ax4 = ax.twinx() #for coverage of corresponding mapping files (from high and low pool .bam files)
		
			# --- adding (variant) information --- #
			
			ax_dAF = ax.scatter( data_to_plot_x[ idx ], data_to_plot_y[ idx ], c=map( abs, data_to_plot_y[ idx ] ), s=1, cmap="cool", label="dAF of pools" )		
			###ax_raw_sig = ax.scatter( data_to_plot_x_sig[ idx ], data_to_plot_y_sig[ idx ], c=map( abs, data_to_plot_y_sig[ idx ] ), s=1, cmap="copper", label="dAF of dARCs" )
			###ax_pool1 = ax.scatter( data_to_plot_x[ idx ], af_pool1[ idx ], s=1, color="red", label="AF pool1" )
			###ax_pool2 = ax.scatter( data_to_plot_x[ idx ], af_pool2[ idx ], s=1, color="lime", label="AF pool2")
			
			###ax_raw_den_sig = ax2.plot( data_to_plot_x_density[ idx ], data_to_plot_y_density[ idx ], marker='x', color = 'black', label="density of dARCs")
			###ax2.scatter( data_to_plot_x_density[ idx ], data_to_plot_y_density[ idx ], linestyle='-', marker='X', s=1, color="black" )
			
			data_to_plot_y_normalized_sig_SNV_density = []
			
			for idx1, ele1 in enumerate ( data_to_plot_y_density[ idx ] ):
				for idx2, ele2 in enumerate( data_to_plot_y_density_merged_ori[ idx ] ):
					if idx1 == idx2:
						if ele2 != 0:
							norm_dARCs = float(ele1)/float(ele2)
							unnorm_dARCs = ele1
							data_to_plot_y_normalized_sig_SNV_density.append( norm_dARCs )
							out.write( "\t".join( [str(each) , str(idx1) , str(norm_dARCs) , str(unnorm_dARCs) , '\n' ]))
						else:
							data_to_plot_y_normalized_sig_SNV_density.append( 0 )
										
			ax_norm_den_sig = ax3.plot( data_to_plot_x_density[ idx ], data_to_plot_y_normalized_sig_SNV_density , marker='+', color = 'orange', label="normalized density of dARCs")

			#ax_cov_high_pool = ax4.scatter( raw_data_x_cov_high[ idx ], raw_data_y_cov_high[ idx ], s=1, color="black", label="coverage high pool" )
			#ax_cov_low_pool = ax4.scatter( raw_data_x_cov_low[ idx ], raw_data_y_cov_low[ idx ], s=1, color="blue", label="coverage low pool" )

			#--- add legend ---#
			
			#lns = [ax_dAF, ax_raw_sig, ax_cov_high_pool, ax_cov_low_pool, ax_raw_den_sig, ax_norm_den_sig, ax_pool1, ax_pool2]
			
			# --- only legend of ax plots are set here --- #

			legend_elements = [	Line2D([0], [0], color='orange', lw=4, label='normalized density of dARCs'),
			Line2D([0], [0], marker='o', color='w', label="dAF of pools", markerfacecolor='deepskyblue', markersize=10),
			Line2D([0], [0], marker='o', color='w', label='normalized coverage of high pool', markerfacecolor='red', markersize=10),
			Line2D([0], [0], marker='o', color='w', label='normalized coverage of low pool', markerfacecolor='lime', markersize=10)]


			#legend_elements = [Line2D([0], [0], color='black', lw=4, label='density of dARCs'), 
			#Line2D([0], [0], color='orange', lw=4, label='normalized density of dARCs'),
			#Line2D([0], [0], marker='o', color='w', label="dAF of pools", markerfacecolor='deepskyblue', markersize=10),
			#Line2D([0], [0], marker='o', color='w', label="dAF of dARCs", markerfacecolor='goldenrod', markersize=10),
			#Line2D([0], [0], marker='o', color='w', label='AF pool1', markerfacecolor='red', markersize=10),
			#Line2D([0], [0], marker='o', color='w', label='AF pool2', markerfacecolor='lime', markersize=10)]
			
			
			# legend_elements = [ Line2D([0], [0], marker='o', color='w', label='coverage high pool', markerfacecolor='black', markersize=10),
			# Line2D([0], [0], marker='o', color='w', label='coverage low pool', markerfacecolor='blue', markersize=10)]
			#Patch(facecolor='orange', edgecolor='r',label='Color Patch')
			###ax.legend(handles=legend_elements, bbox_to_anchor=( 0.8, 1.18 ), ncol=6, fancybox=True)
			
			ax.legend(handles=legend_elements, bbox_to_anchor=( 0.8, 1.18 ), ncol=4, fancybox=True)
		   
			# --- plot secons y axis behind the first one --- #
			
			ax3.spines['right'].set_position(('axes', 1.03) )      
			#ax4.spines['right'].set_position(('axes', 1.1) )  	    
			ax2.set_yticks([])
			ax3.set_yticks([0, 0.25, 0.5, 0.75, 1])

			# --- scale axis --- #
			if len(data_to_plot_y_density[ idx ]) != 0:
				print max(data_to_plot_y_density[ idx ])+1
				ax.set_ylim( -1 , 1 ) #dAF
				ax2.set_ylim( - ( max(data_to_plot_y_density[ idx ])+1) , max(data_to_plot_y_density[ idx ])+1) #raw density of sig SNVs
				ax3.set_ylim( -1 , 1 ) #normalized sig SNVs number
			if len(data_to_plot_y_density[ idx ]) == 0:
				
				ax.set_ylim( -1 , 1 ) #dAF
				ax2.set_ylim( 0 , 0) #raw density of sig SNVs
				ax3.set_ylim( -1 , 1 ) #normalized sig SNVs number

						
			bp_x_high = raw_data_x_cov_high[ idx ]	#position for point in Mbp (no modification required)
			bp_y_high = raw_data_y_cov_high[ idx ]	#coverage values (average per window)
			bp_x_low = raw_data_x_cov_low[ idx ]	#position for point in Mbp (no modification required)
			bp_y_low = raw_data_y_cov_low[ idx ]	#coverage values (average per window)
			
			bpnorm_y_high = []
			bpmax_value_high = float( max( bp_y_high ) )
			for bpvalue_high in bp_y_high:
				bpnorm_y_high.append( -1 * ( bpvalue_high / bpmax_value_high ) )
			
			ax.plot( bp_x_high, bpnorm_y_high, color="red" )
			
			bpnorm_y_low = []
			bpmax_value_low = float( max( bp_y_low ) )
			for bpvalue_low in bp_y_low:
				bpnorm_y_low.append( -1 * ( bpvalue_low / bpmax_value_low ) )
			
			ax.plot( bp_x_low, bpnorm_y_low, color="lime" )
						
			
			ax3.set_ylabel( "                                  normalized density of dARCs" , color="orange" )
			#ax2.set_ylabel( "                      density of dARCs" , color="black"  )
			ax.set_xlabel( "chromosome position [Mbp]" )
			ax.set_ylabel( "normalized coverage                    dAF               " )
			
			ax.set_xlim(  0, chr_lengths[ each ]/1000000.0 )
			
			start, end = ax.get_xlim()
			ax.xaxis.set_ticks( np.arange( start, end, 1 ) )
			ax.set_yticklabels(['1', '0.75', '0.5', '0.25', '0', '0.25', '0.5','0.75', '1'])

			plt.subplots_adjust( left=0.04, right=0.93, top=0.89, bottom=0.2 )
			fig.savefig( fig_output_file, dpi=300 )
			plt.close('all')
		
	return data_to_plot_x, data_to_plot_y, intervalls	


def plot_genome_wide_single_pos_dAF( af_frequency_vcf, af_frequency_vcf_sig, chr_lengths ):
	"""! @brief show genome wide distribution of AF frequencies """
	
	# --- load information about chromosomes --- #
	
	chr_names = sorted( chr_lengths.keys() )
	raw_data_x = [ ]
	raw_data_y = [ ]
	
	raw_data_x_sig = [ ]
	raw_data_y_sig = [ ]
	
	for  k in chr_names:
		raw_data_x.append( [] )
		raw_data_y.append( [] )
		
	for  k in chr_names:
		raw_data_x_sig.append( [] )
		raw_data_y_sig.append( [] )
	
	# --- load normal variant information (variants detected in both pools) --- #
	
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] )/1000000.0 )
				raw_data_y[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
			line = f.readline()
			
	with open( af_frequency_vcf_sig, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				raw_data_x_sig[ chr_names.index( parts[0] ) ].append( int( parts[1] )/1000000.0 )
				raw_data_y_sig[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
			line = f.readline()
		
	# --- construct single plots --- #
	
	for idx, each in enumerate( chr_names ):
		fig_output_file = af_frequency_vcf + "." + each + ".genome_wide.single_variants.png"
	
		fig, ax = plt.subplots(  figsize=(20, 3) )	
			
		# --- adding chromosomes --- #
		
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0, 0 ] , color="black", linewidth=.5 )
		for i in range(10):
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.1*i, 0.1*i ] , color="grey", linewidth=.1 )
		
		# --- adding variant information --- #
		
		ax.scatter( raw_data_x[ idx ], raw_data_y[ idx ], c=map( abs, raw_data_y[ idx ] ), s=1, cmap="cool" )	#hot, binary
		
		ax.scatter( raw_data_x_sig[ idx ], raw_data_y_sig[ idx ], c=map( abs, raw_data_y_sig[ idx ] ), s=1, cmap="copper" )	#hot, binary

		ax.set_xlabel( "chromosome position [Mbp]" )
		ax.set_ylabel( "delta Allele Frequency" )
		ax.set_title( each )
		
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)

		ax.set_ylim( 0, 1 )
		ax.set_xlim(  0, chr_lengths[ each ]/1000000.0 )
		
		start, end = ax.get_xlim()
		ax.xaxis.set_ticks( np.arange( start, end, 1 ) )
		
		plt.subplots_adjust( left=0.03, right=0.98, top=0.98, bottom=0.2 )
		fig.savefig( fig_output_file, dpi=300 )
		plt.close('all')


def load_seq_lengths( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	seq_lens = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split('.')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					seq_lens.update( { header: len( seq ) } )
					header = line.strip()[1:].split('.')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		seq_lens.update( { header: len( seq ) } )
	return seq_lens


def construct_overlapping_delta_AF_frequency_hist( af_frequency_vcf, af_frequency_vcf_sig, af_frequency_vcf_merged_ori ):
	"""! @brief construct histogram of delta AF distribution """

	fig_output_file = af_frequency_vcf + ".overlapping_AF_hist.png"
	print 'HI'
	delta_AFs = []
	delta_AFs_sig = []
	delta_AFs_merged_ori = []

	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				delta_AFs.append( float( line.strip().split('\t')[-1] ) )
			line = f.readline()

	with open( af_frequency_vcf_sig, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				delta_AFs_sig.append( float( line.strip().split('\t')[-1] ) )
			line = f.readline()

	with open( af_frequency_vcf_merged_ori, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				delta_AFs_merged_ori.append( float( line.strip().split('\t')[-1] ) )
			line = f.readline()
	
	fig, ax = plt.subplots()
	ax.hist( delta_AFs_merged_ori, bins=10000, color = 'black' )
	ax.hist( delta_AFs, bins=10000 )
	ax.hist( delta_AFs_sig, bins=10000, color = 'red' )

	ax.set_ylim( 0, 4900 ) # needs to be set according to expected value
	ax.set_ylabel( "number of variants" )
	ax.set_xlabel( "delta allele frequency (high pool - low pool)" )
	
	fig.savefig( fig_output_file, dpi=300 )
	print "total number of calculated delta allele frequencies: " + str( len( delta_AFs ) )


def construct_delta_AF_frequency_hist( af_frequency_vcf ):
	"""! @brief construct histogram of delta AF distribution """
	
	delta_AFs = []
	
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				delta_AFs.append( float( line.strip().split('\t')[-1] ) )
			line = f.readline()
	
	fig_output_file = af_frequency_vcf + ".AF_hist_for_sig.png"
	
	fig, ax = plt.subplots()
	ax.hist( delta_AFs, bins=10000 )

	ax.set_ylim( 0, 50 ) # needs to be adjusted depending on the input file, these figure margins might be to high/low
	ax.set_ylabel( "number of variants" )
	ax.set_xlabel( "delta allele frequency (pool1 - pool2)" )
	
	fig.savefig( fig_output_file, dpi=300 )



def main( arguments ):
	"""! @brief run all parts of this script """
	
	# ---- collecting inputs --- #
	
	input_vcf = arguments[ arguments.index( '--input_vcf' ) + 1 ]
	input_vcf_sig = arguments[ arguments.index( '--input_vcf_sig_SNP' ) + 1 ]
	fasta_ref_file = arguments[ arguments.index( '--reference_file' ) + 1 ]
	input_mapping_cov_high = arguments[ arguments.index( '--input_mapping_cov_high' ) + 1 ]
	input_mapping_cov_low = arguments[ arguments.index( '--input_mapping_cov_low' ) + 1 ]
	coverage_threshold_to_plot = float(arguments[ arguments.index( '--coverage_threshold_to_plot' ) + 1 ])
	
	output_dir = arguments[ arguments.index( '--output_dir' ) + 1 ]
	
	# --- if output directory does not already exist, generate it --- #
	
	if not output_dir[-1] == "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
		
	# --- get name of the pools ( must be equal to the last two column names in the vcfs )--- #
	
	raw_high_name = arguments[ arguments.index( '--high_pool' ) + 1 ]
	high_name = []
	for each in raw_high_name.split(','):
		if len( each ) > 0:
			high_name.append( each )
	
	raw_low_name = arguments[ arguments.index( '--low_pool' ) + 1 ]
	low_name = []
	for each in raw_low_name.split(','):
		if len( each ) > 0:
			low_name.append( each )

	# --- loading of vars from merged ori vcf used for sig snv calc --- #
	
	merged_ori_vcf = arguments[ arguments.index( '--in_merged_ori_vcf' )+1 ]	

	header, sample_size = load_header_and_samples( merged_ori_vcf )

	# --- get column (idx) for each pool --- #
	
	pool1_indices = []
	pool2_indices = []

	for idx, each in enumerate( header ):
		if each in high_name:
			pool1_indices.append( idx ) 	
		elif each in low_name:
			pool2_indices.append( idx )
	
	# --- load the vars which would be used for fisher exact test, but without filtering for an alpha, since we want to get all vars --- #
	
	variants = load_variants( merged_ori_vcf , pool1_indices, pool2_indices )	

	# --- prepare outputfile, since this vcf is now filtered as for the sig SNVs, but contains all vars which were tested. Thus, not filtered for an alpha --- #
	
	output_vcf = output_dir + merged_ori_vcf.split('/')[-1].lower().split('.vcf')[0] + "_merged_ori_vcf_filtered_as_for_sig_SNVs.vcf"

	with open( output_vcf, "w" ) as out:
		out.write( "\t".join( header[:7] + [ "p-value" ] +header[8:] ) + '\n' )
		for var in variants:
			#print var
			out.write(  ('\t').join(var) + '\n' )

	# --- calling all functions and running detection of loci --- #

	af_frequency_vcf = output_dir + "allele_frequencies.vcf"
	af_frequency_vcf_sig = output_dir + "allele_frequencies_sig.vcf"
	af_frequency_vcf_merged_ori = output_dir + "allele_frequencies_merged_ori.vcf"

	# --- setting window sizes and step sizes for dAF plots and plotting of amounts of sig SNVs --- #
	
	window_sizes = [100] 
	step_size = 5
	
	window_size_density = 100000
	step_size_density = 30000
	
	# --- calculate average coverage in vcfs --- #
	
	avg_cov_high, avg_cov_low = get_coverage( input_vcf, high_name, low_name )
	
	# coverage confidence intervall borders for high pool
	
	min_high_cov = 0.75*avg_cov_high	
	max_high_cov = 1.5*avg_cov_high	
	
	# coverage confidence intervall borders for low pool
	
	min_low_cov = 0.75 * avg_cov_low	
	max_low_cov = 1.5*avg_cov_low	
	
	get_delta_allel_frequencies( input_vcf, af_frequency_vcf, high_name, low_name, min_high_cov, max_high_cov, min_low_cov, max_low_cov )

	
	avg_cov_high_sig, avg_cov_low_sig = get_coverage( input_vcf_sig, high_name, low_name )
	
	avg_cov_merged_ori, avg_cov_merged_ori = get_coverage( output_vcf, high_name, low_name )
	
	# --- since for sig SNVs and merged ori SNVs we do not want to filter for an average coverage in these vcfs, we set the cut off here manually, so that every SNV will be included in the get_delta_allel_frequencies and plot_genome_wide_delta_allele_frequencies analysis --- #
	
	# coverage confidence intervall borders for pool1 (high)
	
	min_high_cov = 0 #0.75*avg_cov_high	#0.75
	max_high_cov = 1000 #1.5*avg_cov_high	#1.5
	
	# coverage confidence intervall borders for pool2 (low)
	
	min_low_cov = 0 #0.75 * avg_cov_low	#0.75
	max_low_cov = 1000 #1.5*avg_cov_low	#1.5

	get_delta_allel_frequencies( input_vcf_sig, af_frequency_vcf_sig, high_name, low_name, min_high_cov, max_high_cov, min_low_cov, max_low_cov )

	get_delta_allel_frequencies( output_vcf, af_frequency_vcf_merged_ori, high_name, low_name, min_high_cov, max_high_cov, min_low_cov, max_low_cov )

	# # --- construct dAF frequency histogram --- #
	
	construct_delta_AF_frequency_hist( af_frequency_vcf )
	construct_delta_AF_frequency_hist( af_frequency_vcf_sig )
	construct_delta_AF_frequency_hist( af_frequency_vcf_merged_ori )
	construct_overlapping_delta_AF_frequency_hist( af_frequency_vcf, af_frequency_vcf_sig, af_frequency_vcf_merged_ori )
	
	# --- load chromosome length as preparation for genome wide plots --- #
	chr_lengths = load_seq_lengths( fasta_ref_file )
	
	value_output_file = output_dir + "value_output_file.txt"
	
	print af_frequency_vcf_sig
	
	with open( value_output_file, "w" ) as out:
		plot_genome_wide_single_pos_dAF( af_frequency_vcf, af_frequency_vcf_sig, chr_lengths )
		for window_size in window_sizes:
			data_to_plot_x, data_to_plot_y, intervalls = plot_genome_wide_delta_allele_frequencies( coverage_threshold_to_plot, input_mapping_cov_high, input_mapping_cov_low, af_frequency_vcf, af_frequency_vcf_sig, af_frequency_vcf_merged_ori, chr_lengths, window_size, step_size , window_size_density, step_size_density, output_dir)	 
			chr_names = sorted( chr_lengths.keys() )
			for idx, k in enumerate( chr_names ):
				intervalls_of_interest = intervalls[ idx ]
				y_values_of_interest = data_to_plot_y[ idx ]
				for i, intervall in enumerate( intervalls_of_interest ):
					out.write( k + '\t' + str( y_values_of_interest[ i ] ) + '\t' + str( intervall[0] ) + '\t' + str( intervall[1] ) + '\n' )
					

if __name__ == '__main__':
	
	main( sys.argv )
	
	print "all done!"
