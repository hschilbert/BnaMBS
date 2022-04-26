### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###
### v0.4 ###

__usage__ = """
					python get_intervals_based_on_sig_snps.py
					--sig_snp_vcf <FULL_PATH_TO_SIG_SNPS_VCF_FILE>
					--snp_eff_res <FULL_PATH_TO_SNPEFF_RES>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					--anno <ANNOTATION_FILE_ORTHOFINDER>
					--ZCR <FULL_PATH_TO_ZCR_FILE_FROM_PAV_FINDER>
					
					--dis_out_reg <SET_DISTANCE_IN_BP_BETWEEN_STAND_ALONE_SNVS_AND_REGION>
					--min_nr_sig_snvs_in_reg <SET_MINIMUM_NUMBER_OF_SNVS_IN_REG>
					--dis_in_reg <SET_DISTANCE_IN_BP_OF_SNVS_IN_A_REGION>
					
					bug reports and feature requests: hschilbe@cebitec.uni-bielefeld.de
					"""

import numpy as np
import sys, os, re
import matplotlib.pyplot as plt

# --- end of imports --- #


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
	
def load_chr_lengths_and_names ( chr_len_file ):
	"""! @brief load chr length for plotting """
	
	chr_lengths = {}
	chr_names = []
	
	with open ( chr_len_file, 'r' ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			chr_lengths.update ({ parts[0]: parts[1] })
			chr_names.append (parts[0]) 
			line = f.readline()	
	return chr_lengths, chr_names
	
def find_high_impact_variants ( snpeff_res, output_file_vars, annotation ):
	"""! @brief load only first high impact effect per gene """
	
	ended_genes = []
	vars_w_anno = []
	effect_impact_positions = { }
		
	splice_variants_counter = 0
	premature_stop_counter = 0
	frame_shift_counter = 0
	lost_stop_counter = 0
	
	high_counter = 0
	moderate_counter = 0
	low_counter = 0
	modifier_counter = 0
	
	premature_stops = []
	
	data_for_extraction = []
	
	total_annotated_variant_counter = 0
	total_small_variant_counter = 0
	total_snp_counter = 0
		
	with open( snpeff_res, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				
				if len(parts[7].split(';')[0].split(',')) == 1:
					cov_pool1 = parts[7].split(';')[0].split(',')[0] 
					cov_pool2 = 'n/a'
					ratio_pools = 'n/a'
				
				if len(parts[7].split(';')[0].split(',')) > 1:
					cov_pool1 = parts[7].split(';')[0].split(',')[0]
					cov_pool2 = parts[7].split(';')[0].split(',')[1]
					ratio_pools = parts[7].split(';')[0].split(',')[2]
				total_annotated_variant_counter += 1
				
				try:
					x = abs( len( parts[3] ) - len( parts[4] ) )
				except IndexError:
					print line
					
				if abs( len( parts[3] ) - len( parts[4] ) ) < 100:	#filter out large InDels
					total_small_variant_counter += 1
					if len( parts[3] ) == len( parts[4] ):
						total_snp_counter += 1
					subparts = parts[7].split(',') 

					for subpart in subparts:
						
						# --- check all variant annotations and select only one (highest impact and first) --- #
						try:
							ID = re.findall( "Gene_Bna[A-Z]\d+g\d+[A-Z]", subpart )[0] 
							if ID not in ended_genes:

								# --- analyse high impact effects --- #
								if 'HIGH' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]							
									except:
										effect_annotation = "ERROR"
										if 'stop_gained' in subpart:
											premature_stop_counter += 1
											premature_stops.append( ID )
											ended_genes.append( ID )
											effect_annotation = 'stop_gained'
										elif 'stop_lost' in subpart:
											lost_stop_counter += 1
											effect_annotation = 'stop_lost'
										elif 'splice_region_variant' in subpart:
											splice_variants_counter += 1
											effect_annotation = 'splice_region_variant'
										elif 'frameshift' in subpart:
											frame_shift_counter += 1
											effect_annotation = 'frameshift'
										
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										high_counter += 1
										if effect_annotation != "ERROR":
											data_for_extraction.append( [ parts[0], parts[1], parts[3], parts[4], re.findall( "Gene_Bna[A-Z]\d+g\d+[A-Z]", ID )[0], effect_annotation, str(cov_pool1), str(cov_pool2), str(ratio_pools) ] )											
											#vars_w_anno.update ( { parts[0]+parts[1].zfill( 9 ) : { 'chr': parts[0], 'pos':parts[1], 'ref': parts[3] , 'alt': parts[4] , 'gene_ID': re.findall( "chr[A-Z]\d+_[A-Z]\d+.g\d+.t\d+", ID )[0], 'effect_annotation': effect_annotation } } )
								
								elif 'MODERATE' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										moderate_counter += 1
								elif 'LOW' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										low_counter += 1
								elif 'MODIFIER' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" })
										modifier_counter += 1
						except:
							pass 
			line = f.readline()
			
	print "RESULTS:"
	
	print "total annotated variants: " + str( total_annotated_variant_counter )
	print "total small variants: " + str( total_small_variant_counter )
	print "total number of SNPs: " + str( total_snp_counter )
	
	print "number of splice variants: " + str( splice_variants_counter )
	print "number of premature stops: " + str( premature_stop_counter )
	print "number of frame shifts: " + str( frame_shift_counter )
	print "number of lost stops: " + str( lost_stop_counter )
	
	print "HIGH: " + str( high_counter )
	print "MODERATE: " + str( moderate_counter )
	print "LOW: " + str( low_counter )
	print "MODIFIER: " + str( modifier_counter )
	
	print "premature stop check: " + str( len( premature_stops ) )
	print "premature stop check (unique): " + str( len( list( set( premature_stops ) ) ) )
	
	lengths_of_ref_alleles = []
	lengths_of_alt_alleles = []
	
	# --- write collected data into output file --- #
	
	with open( output_file_vars, "w" ) as out:
		header = [ "Chromosome", "Position", "ReferenceAllel", "AlternativeAllel", "GeneID", "EffectType",  "CoveragePool1","CoveragePool2", "RatiosInPools", "Annotation"]
		out.write( "\t".join( header ) + '\n' )
		for each in data_for_extraction:
			try:
				out.write( "\t".join( each + [ annotation[ each[4].split('_')[1] ] ] )  + '\n' )
				vars_w_anno.append( each + [ annotation[ each[4].split('_')[1] ] ] )
			except KeyError:
				out.write( "\t".join( each ) + '\tn/a\n' )
				vars_w_anno.append(   each + ['\tn/a\n'] )
				
			lengths_of_ref_alleles.append( len( each[2] ) )
			lengths_of_alt_alleles.append( len( each[3] ) )
	
	return vars_w_anno
			
def load_sig_variants ( sig_snps ):
	"""! @brief load all sig variant positions from given VCF """
		
	variants = []
	 
	with open( sig_snps, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				variants.append( [parts[0], parts[1], parts[3], parts[4] , parts[7]] ) 
			counter =+ 1
			line = f.readline()
	return variants
						
def def_intervals_on_sig_snps ( sig_snps, dis_out_reg, dis_in_reg ):
	""" define intervals/regions based on distance of sig snps to each other """

	counter = 0
	sig_snp_region = []
	sig_snp_w_no_region = []
	regions = {}
	
	for idx, element in enumerate(sorted(sig_snps)):
		if idx+1 < len(sig_snps) :
			if sig_snps[idx][0] == sig_snps[idx+1][0]: 
				# check how close is next variant, if it is less than X kb away assume it is still one region:
				
				if abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1])) < dis_out_reg:
					sig_snps[idx].append(abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1]))) 
					sig_snp_region.append(sig_snps[idx] )
				else: 
					if len(sig_snp_region) != 0: 
						sig_snps[idx].append(abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1])))
						sig_snp_region.append(sig_snps[idx]) 
						counter += 1
						regions.update( { 'REGION_' + str(counter).zfill(5): sig_snp_region } )	
						sig_snp_region = [] 
					else: 
						if abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1])) >= dis_out_reg:
							sig_snps[idx].append(abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1])))
							sig_snp_w_no_region.append(sig_snps[idx])
			else:
				if len(sig_snp_region) != 0: 
					sig_snps[idx].append(abs(int(sig_snps[idx][1]) - int(sig_snps[idx+1][1])))
					sig_snp_region.append(sig_snps[idx])
					counter += 1
					regions.update( { 'REGION_' + str(counter).zfill(5): sig_snp_region } )	
					sig_snp_region = []

	print 'region building based on sig SNPs is done'
	
	sig_snp_region.append(sig_snps[idx]) 
	counter += 0
	regions.update( { 'REGION_' + str(counter).zfill(5): sig_snp_region } )		
	
	# --- check for distance of the SNPs within a region --- #
	
	big_dis_snv = []
	counter1 = 0
	counter2 = 0
	
	print 'Before applying interal SNP distance filter there are ' + str( len(regions) ) +' regions.' 

	for region in sorted(regions):
		big_dis_snv = []
		for sig_snv in regions[region]:
			try:
				if sig_snv[5] >= dis_in_reg:
					big_dis_snv.append(sig_snv)
				else:
					continue
			except: 
				continue
		
		if len(big_dis_snv) >= 4:
			counter2 = counter2 + 1
		if len(big_dis_snv) < 4:
			del regions[region]
			counter1 = counter1 + 1
	
	print 'There are ' + str(counter1) + ' regions filtered out, because they had not at least 3 SNVs with a distance of ' + str(dis_in_reg) + ' bp between each other.\n'
	print 'There are ' + str(counter2) + ' regions kept in, because they had at least 3 SNVs with a distance of ' + str(dis_in_reg) +' bp between each other.\n'
	print 'Checking for internal SNPs distance is done.\n'
	
	return regions, sig_snp_w_no_region	

def load_my_ints ( my_int ):
	
	my_ints = {}
	
	with open( my_int, "r" ) as f:
		line = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len(parts[0].split('_')) < 4: 
				try:			
					my_ints[ parts[0].split('_')[2] ].append( { 'ID':parts[0].split('_')[2]+'_'+ parts[0].split('_')[1] , 'size':parts[1], 'start': parts[2] , 'end': parts[3], 'Nr. of sig SNPs': parts[4] } )
				except KeyError:
					my_ints.update( { parts[0].split('_')[2] : [ { 'ID':parts[0].split('_')[2]+'_'+ parts[0].split('_')[1] , 'size':parts[1], 'start': parts[2] , 'end': parts[3], 'Nr. of sig SNPs': parts[4] } ] } )
			else:
				try:			
					my_ints[ parts[0].split('_')[2]+'_'+ parts[0].split('_')[3]].append( { 'ID':parts[0].split('_')[2]+'_'+ parts[0].split('_')[3]+'_'+ parts[0].split('_')[1], 'size':parts[1], 'start': parts[2] , 'end': parts[3], 'Nr. of sig SNPs': parts[4] } )
				except KeyError:
					my_ints.update( { parts[0].split('_')[2]+'_'+ parts[0].split('_')[3] : [ { 'ID':parts[0].split('_')[2]+'_'+ parts[0].split('_')[3]+'_'+ parts[0].split('_')[1], 'size':parts[1], 'start': parts[2] , 'end': parts[3], 'Nr. of sig SNPs': parts[4] } ] } )
			
			line = f.readline()
	
	print 'The length of my int is the same as the number of input chromosomes: ' + str( len(my_ints) ) + '\n'
	
	return my_ints

def load_ZCR ( ZCR ):
	
	ZCRs = {}
	counter = 0
	
	with open( ZCR, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')	
			try:					
				ZCRs[ parts[0].split('_%_')[0] ].append( { 'ID':parts[0].split('_%_')[0], 'start': parts[0].split('_%_')[1] , 'end': parts[0].split('_%_')[2] } )
			except KeyError:
				ZCRs.update({ parts[0].split('_%_')[0] : [ { 'ID':parts[0].split('_%_')[0], 'start': parts[0].split('_%_')[1] , 'end': parts[0].split('_%_')[2] } ] } )
			counter = counter + 1
			line = f.readline()
			
	print 'The length of zcrs int is the same as the number of input chromosomes: ' + str( len(ZCRs) ) + '\n'
	
	return ZCRs


def main( arguments ):
	"""! @brief run everything """
	
	sig_snps = arguments[ arguments.index('--sig_snp_vcf')+1 ]
	snpeff_res = arguments[ arguments.index('--snp_eff_res')+1 ]
	ZCR = arguments[ arguments.index('--ZCR')+1 ]

	min_nr_sig_snvs_in_reg = int(arguments[ arguments.index('--min_nr_sig_snvs_in_reg')+1 ])
	dis_out_reg = int(arguments[ arguments.index('--dis_out_reg')+1 ])
	dis_in_reg = int(arguments[ arguments.index('--dis_in_reg')+1 ])

	prefix = arguments[ arguments.index('--out')+1 ]
	
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	
	""" load functional annotation in a dict """
	
	if '--anno' in arguments:
		annotation_file = arguments[ arguments.index('--anno')+1 ]
		annotation = load_annotation( annotation_file )
	else:
		annotation = {}
	print 'Functional annotation loaded.\n'
	
	
	""" load chr length and names in a dict """
	
	chr_len_file = '/vol/cluster-data/hschilbe/Rapeq/reference/Brassica_napus_v4.1.chromosomes_chr_lengths.txt'
	chr_lengths, chr_names = load_chr_lengths_and_names ( chr_len_file )
	print 'Chr lengths loaded.\n'


	""" load high impact variants in a dict """
	
	output_file_vars = prefix + 'annotated_vars.txt'
	vars_w_anno = find_high_impact_variants( snpeff_res, output_file_vars, annotation )
	print 'High impact vars loaded.\n'
	
	""" Interval detection """
		
	sig_snps = load_sig_variants ( sig_snps ) 

	regions, sig_snp_w_no_region = def_intervals_on_sig_snps ( sig_snps , dis_out_reg, dis_in_reg)
	print 'Regions based on sig SNPs loaded.\n'
	print 'There are ' + str(len(sig_snp_w_no_region)) + ' SNVs which are significant but not anchored in a region, because they stand alone.\n'
	
	""" outwrite regions with additional information (nr of sig SNVs) in a new file """
	
	counter_regions = 1
	counter_reg = 0
	regions_out = prefix + 'regions_corrected_for_dis_in_reg.txt'
	nr_of_sig_s_in_reg = []
	snp_positions = []

	with open ( regions_out , 'w' ) as out:
		header = "\t".join( ['Interval_ID', 'Size [bp]', 'Start','End','Number of sig SNPs'])
		out.write ( header + '\n')
		
		for idx, region in enumerate(sorted(regions)):
			if counter_regions <= len(regions) :				
				if abs(int(regions[region][0][1]) - int(regions[region][-1][1])) != 0 and len(regions[region]) >= min_nr_sig_snvs_in_reg : 
					out.write (  "\t".join( [ str(region) + '_' + str(regions[region][0][0]), str(abs(int(regions[region][0][1]) - int(regions[region][-1][1]))), str(regions[region][0][1]), str(regions[region][-1][1]), str(len(regions[region])) ] ) + '\n' )
					nr_of_sig_s_in_reg.append( str( len(regions[region]) ) )
					counter_reg += 1
			counter_regions += 1
		
	print 'There are ' + str(counter_reg) + ' regions in total. Now check for ZCRs between them and if we can combine some regions.'


	""" load ZCRs and combine intervals which distance was due to the ZCR longer than dis_out_reg """
	
	my_ints = load_my_ints (regions_out)
	ZCRs = load_ZCR (ZCR)

	counter = 0
	sum_zcrs = 0	

	new_combined = {}
	out_my_regions = []
	one_reg_on_chr = {}
	intervals_to_remove = []
	
	for key in sorted(my_ints.keys()):
		intervals = sorted(my_ints[key])
		if len(	intervals) > 1:
			for idx, my_interval in enumerate(intervals):
				if idx+1 < len(intervals):
					out_reg_size = int(intervals[idx+1]['start']) - int(my_interval['end'])
							
					for zcr in sorted(ZCRs[key]):
						if int(zcr['start']) > int(my_interval['end']) and int(zcr['end']) < int(intervals[idx+1]['start']):
							sum_zcrs = sum_zcrs + (int(zcr['end']) - int(zcr['start']))
							out_my_regions.append(zcr)
							
					if (int(out_reg_size) - int(sum_zcrs)) < dis_out_reg:
						intervals_to_remove.append(my_interval)
						intervals_to_remove.append(intervals[idx+1])
						
						if len(intervals[idx+1]['ID'].split('_')) < 3:
							try:					
								new_combined[ key ].append( {'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[1]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) } )
								intervals[intervals.index(intervals[idx+1])]={'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[1]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) }
							except KeyError:
								new_combined.update( { key : [ {'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[1]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) } ] } )
								intervals[intervals.index(intervals[idx+1])]={'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[1]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) }
						if len(intervals[idx+1]['ID'].split('_')) >= 3:
							try:					
								new_combined[ key ].append( {'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[2]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) } )
								intervals[intervals.index(intervals[idx+1])]={'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[2]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) }
							except KeyError:
								new_combined.update( { key : [ {'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[2]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) } ] } )
								intervals[intervals.index(intervals[idx+1])]={'ID': (my_interval['ID']+'%'+intervals[idx+1]['ID'].split('_')[2]), 'size': ( int( intervals[idx+1]['end'] )- int(my_interval['start'] ) ), 'start': my_interval['start'] , 'end': intervals[idx+1]['end'], 'Nr. of sig SNPs': ( int( intervals[idx+1]['Nr. of sig SNPs'] ) + int(my_interval['Nr. of sig SNPs'] ) ) }

					
					sum_zcrs = 0	
		else:
			one_reg_on_chr.update( {key : intervals} )

	print 'There are ' + str(len(one_reg_on_chr)) +' regions which are alone on one chromosome, thus can not be combined.'
	print 'These are: ' + str(one_reg_on_chr) + '\n'
	
	print 'Now I have to remove ' + str( len(intervals_to_remove) ) + ' regions/intervals, because they are now part of a combined region.'
	print 'These are: ' + str(intervals_to_remove) + '\n'


	""" write new (included combined) regions into a new file and remove thereby all just merged regions from the interval set """

	new_combined_final = {}
	
	output_file = prefix + 'Regions_corrected_for_ZCRs.txt'
	
	# write new regions into file #
	
	with open(  output_file , "w" ) as out:
		header = "\t".join( ['Interval_ID', 'Size [bp]', 'Start','End','Number of sig SNPs'])
		out.write ( header + '\n')
		for key in sorted(my_ints.keys()):
			intervals = sorted(my_ints[key])	
			for key2  in sorted(new_combined):
				if key == key2:
					combined_ints = sorted(new_combined[key2])
					for com_int in sorted(combined_ints):
						try:
							new_combined[key]
							intervals.append(com_int)
						except KeyError:
							continue
				
			for rem_interval in sorted(intervals_to_remove):
				try:
					intervals.remove(rem_interval)
				except ValueError:
					continue
			
			for interval in sorted(intervals):
				if len(interval['ID'].split('_')) < 3: 
					out.write ( 'REGION_' + str(interval['ID'].split('_')[1]) + '_' + str(interval['ID'].split('_')[0]) + '\t' + str(interval['size']) + '\t' + str(interval['start']) + '\t' + str(interval['end']) + '\t' + str(interval['Nr. of sig SNPs']) + '\n')
					if len(interval['ID'].split('%')) > 1: 
						try:
							new_combined_final[ key ].append(  interval  )
						except KeyError:
							new_combined_final.update( { key : [interval] } )
						counter = counter +1
						
				if len(interval['ID'].split('_')) >= 3: 
					out.write ('REGION_' +str(interval['ID'].split('_')[2]) + '_' + str(interval['ID'].split('_')[0]) + '_' + str(interval['ID'].split('_')[1]) + '\t' + str(interval['size']) + '\t' + str(interval['start']) + '\t' + str(interval['end']) + '\t' + str(interval['Nr. of sig SNPs']) + '\n')
					if len(interval['ID'].split('%')) > 1: 
						try:
							new_combined_final[ key ].append(  interval  )
						except KeyError:
							new_combined_final.update( { key : [interval] } )
						counter = counter +1
				
	
	print 'The to removed single regions (INCLUDING douple regions if the end result is a triple region): ' + str( len(intervals_to_remove) ) + ' are combined into ' + str(counter) + ' combined regions.\n'
	print 'The single regions are now removed.\n'
	print 'And ' + str(counter) +' regions are combined and added to the interval list: ' + str(new_combined_final) + '\n'	

	""" check in new combined regions whether there are more sig SNPs included now due to the fusion of the two regions """
		
	regions_to_combine = []
	snps_to_com = []
	regs_to_be_deleted = []
	snps_out_now_in_re_to_com =  []
	
	for key in sorted(new_combined_final.keys()): 
		intervals_new = sorted(new_combined_final[key])
		for interval_n in sorted(intervals_new): 		
			for combined_region_ids in interval_n['ID'].split('_')[-1].split('%'): 
				regions_to_combine.append('REGION_'+str(combined_region_ids)) 
			for regs_to_combine in sorted(regions_to_combine): # iterate over region ids, which should be combined		
				try: 
					for snp in sorted(regions[regs_to_combine]): #get list of all snps of one region which should be included in the new combined region
						snps_to_com.append(snp)					
				except KeyError:
					print 'Key does not exist in regions.'
					continue

				regs_to_be_deleted.append(regs_to_combine) #save the single region IDs, to delete them later because we only want to keep the combined regions			

			""" check whether sig SNPs outside regions are now within one of the newly combined regions, if so add these SNPs to the corresponding region """	

			for sig_snp_out_reg in sorted(sig_snp_w_no_region):
				if sig_snp_out_reg[0] == key:
					if interval_n['start'] <= int(sig_snp_out_reg[1]) <= interval_n['end']: # if there are new snp to attach to the region
						snps_out_now_in_re_to_com.append(sig_snp_out_reg) # only count the SNPs which were outside two regions and are now part of the combined region
						snps_to_com.append(sig_snp_out_reg) # adds the outside SNPs to their new combined region
			
			regions.update( { ('%').join(regions_to_combine) : snps_to_com } ) # the snps of the combined regions are now combined
			
			# before we get to the next new combined region, the SNPs and the region IDs need to be clean again, so that we can fill them in with the infos of the next new combined region 
			snps_to_com = []
			regions_to_combine = []		
	
	for reg_to_delete in set(regs_to_be_deleted): # the single region IDs (which are now combined) are removed from the region dict new entry are here removed
		del regions[reg_to_delete]

	print 'Checking in these new combined regions, whether there are more sig SNPs included now due to the fusion of the two regions. Done.\n'
	
	
	""" sort list after chr """
	""" count nr of as HIGH classified SNP in regions """

	chr_dic_reg = {}
	chr_list = []
	
	count_high_i_v = 0
	nr_h_i_v_per_region = {}
	
	for chr_name in chr_names:
		for idx, region in enumerate(sorted(regions.keys())):
			if str(regions[region][0][0]) == chr_name : 
				for snp in regions[region]: 
					chr_list.append (snp) 
			for var in  sorted(vars_w_anno): 
				if str(regions[region][0][0]) == str([var][0][0]): # check for same chr ID
				
					"""### THIS LINE IS IMPORTANT ###"""
					if int(regions[region][0][1])-5000 <= int([var][0][1]) <= int(regions[region][-1][1])+5000: #check if the high impact variant is located in a region +-5 kb before the start and end of this region
						count_high_i_v += 1
						
			nr_h_i_v_per_region.update( { str(region) : { 'sig_snps_in_reg': str(len(regions[region])), 'counts': int(count_high_i_v), 'chr':str(regions[region][0][0]) , 'start_reg': int(regions[region][0][1])-5000, 'end_reg': int(regions[region][-1][1])+5000 } } )
			count_high_i_v = 0	
			
		chr_dic_reg.update( { chr_name: chr_list } ) 
		chr_list = []
				
	""" briefly add region and nr of sig SNPs within this region into the annotated vars .txt file """
	
	output_file_w_reg = prefix + 'annotated_vars_w_reg.txt'
	
	with open (output_file_w_reg, 'w') as out:
		with open (output_file_vars, 'r') as f:
			header = "\t".join(['Chromosome',	'Position',	'ReferenceAllel'	,'AlternativeAllel',	'GeneID'	,'EffectType'	,'Pool1_coverage_AF_info' , 'Pool2_coverage_AF_info', 'dAF', 'Annotation', 'Region', 'Number of sig SNPs in Region', 'Number of HIGH SNPs in Region'])
			out.write ( header + '\n')
			line = f.readline() 
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				for reg in sorted(nr_h_i_v_per_region.keys()):
					if parts[0] == nr_h_i_v_per_region[reg]['chr'] :
						if int(nr_h_i_v_per_region[reg]['start_reg']) <= int(parts[1]) <= int(nr_h_i_v_per_region[reg]['end_reg']) and int(nr_h_i_v_per_region[reg]['sig_snps_in_reg']) >= int(min_nr_sig_snvs_in_reg) :
							out.write ( "\t".join(parts[0:]) + '\t' + str(reg) + '\t' + str(nr_h_i_v_per_region[reg]['sig_snps_in_reg']) + '\t' + str(nr_h_i_v_per_region[reg]['counts']) + '\n' ) 
				line = f.readline()
	
	total_regions = counter_reg - counter
	print 'There are ' + str(total_regions) + ' regions in total.'	

	""" PLOTTING hist (cumulative) for nr of sig SNPs in regions """
	
	nr_of_sig_s_in_reg
	fig_output_file = prefix + "hist.png"
	fig, ax = plt.subplots(figsize=(20, 3))	
	
	ax.hist ( map(int, nr_of_sig_s_in_reg), bins = 10000, cumulative=True)
	ax.set_xlabel('Amount of SNPs in each region')
	ax.set_ylabel('Frequency')
	ax.set_xlim( 0, max(map(int, nr_of_sig_s_in_reg)))
	
	plt.subplots_adjust( left=0.03, right=0.98, top=0.98, bottom=0.2 )
	fig.savefig( fig_output_file, dpi=300 )
	plt.close('all')
	
	""" PLOTTING significant snps along chr """ 
	
	for chr_name in chr_dic_reg:
		if len(chr_dic_reg[chr_name]) != 0:
			fig_output_file = prefix + str(chr_name) + ".png"
			fig, ax = plt.subplots(figsize=(20, 3))	
			
			# # --- adding chromosomes --- #
			ax.plot( [ 0,  float(chr_lengths[chr_name])/1000000.0], [ 0, 0 ] , color="black", linewidth=.5 )
			
			# # --- add helper lines --- #
			ax.plot( [ 0,  float(chr_lengths[chr_name])/1000000.0], [ 0.25, 0.25 ] , color="black", linewidth=.1 )
			ax.plot( [ 0,  float(chr_lengths[chr_name])/1000000.0], [ 0.5, 0.5 ] , color="black", linewidth=.1 )
			ax.plot( [ 0,  float(chr_lengths[chr_name])/1000000.0], [ 0.75, 0.75 ] , color="black", linewidth=.1 )

			# --- construct single plots and adding variant information --- #
			for snp in chr_dic_reg[chr_name][0:]:
				ax.scatter( int(snp[1])/1000000.0, 0.5, s=1, color='red' )
			
			ax.set_xlabel( "chromosome position [Mbp]" )			
			ax.spines["top"].set_visible(False)
			ax.spines["left"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.set_frame_on(False)
			ax.axes.get_yaxis().set_visible(False)
			
			ax.set_ylim( 0, 1 )
			ax.set_xlim( 0, float(chr_lengths[chr_name])/1000000.0 )
			
			start, end = ax.get_xlim()
			ax.xaxis.set_ticks( np.arange( start, end, 1 ) )
			
		plt.subplots_adjust( left=0.03, right=0.98, top=0.98, bottom=0.2 )
		fig.savefig( fig_output_file, dpi=300 )
		plt.close('all')


if __name__ == '__main__':
	
	if '--snp_eff_res' in sys.argv and '--sig_snp_vcf' in sys.argv and '--out' in sys.argv and '--min_nr_sig_snvs_in_reg' in sys.argv and '--dis_out_reg' in sys.argv and '--dis_in_reg' in sys.argv and '--ZCR' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
