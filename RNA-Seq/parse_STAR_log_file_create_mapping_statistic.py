### Hanna Schilbert ###
### hschilbe@cebitec.uni-bielefeld.de ###

import sys, glob

# --- end of imports --- #

__usage__ = """
	python parse_STAR_log_file_create_mapping_statistic.py
	--path_to_log_files <FULL_PATH_TO_LOG_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	"""
	
def load_log_file ( log_file ):
	"""! @brief load all infos from given log file """

	logs = {}
	
	with open( log_file, "r" ) as f:
		
		list_of_lines = f.readlines()
				
		total_mapped = (int(list_of_lines[8].strip().split('\t')[1]) + int(list_of_lines[23].strip().split('\t')[1]) + int(list_of_lines[25].strip().split('\t')[1]) )
		percentage_total_mapped = (float(total_mapped) / float(list_of_lines[5].strip().split('\t')[1]))*100

		total_unmapped = int(list_of_lines[28].strip().split('\t')[1]) + int(list_of_lines[30].strip().split('\t')[1]) + int(list_of_lines[32].strip().split('\t')[1])
		percentage_total_unmapped = (float(total_unmapped) / float(list_of_lines[5].strip().split('\t')[1]))*100

		total_multimapped = int(list_of_lines[23].strip().split('\t')[1]) + int(list_of_lines[25].strip().split('\t')[1])
		percentage_total_multimapped = (float(total_multimapped) / float(list_of_lines[5].strip().split('\t')[1]))*100

		logs.update( { 	
		"Number of input reads": int(list_of_lines[5].strip().split('\t')[1]),
		"Total mapped": total_mapped,
		"Total mapped percentage" : round(percentage_total_mapped, 2),
		"Uniquely mapped" : int(list_of_lines[8].strip().split('\t')[1]),
		"Uniquely mapped percentage" : float(list_of_lines[9].strip().split('\t')[1].split('%')[0]),
		"Multi mapped" : total_multimapped,
		"Multi mapped percentage" : round(percentage_total_multimapped,2),
		"Total unmapped" : total_unmapped,
		"Total unmapped percentage" : round(percentage_total_unmapped,2) } )

	return logs
	

def main( arguments ):
	"""! @brief run all parts of this script """
	
	path_to_log_files = arguments[ arguments.index('--path_to_log_files')+1 ]
	out = arguments[ arguments.index('--out')+1 ]
	
	log_files = sorted( glob.glob( path_to_log_files + "*Log.final.out" ) )

	IDs = []
	global_logs = {}
	
	with open ( out, "w" ) as out:
		
		header = "Sample name"+'\t'+"Number of input reads"+'\t'+"Total mapped"+'\t'+"Uniquely mapped"+'\t'+"Multi mapped"+'\t'+"Total unmapped"+'\n'
		out.write (header)
		
		for log_file in log_files:
			ID = log_file.split('/')[-1].split('.')[0].split('_')[0]
			IDs.append( ID )
			log_infos = load_log_file( log_file )
			
			out.write ( str(ID)+'\t'+str(log_infos["Number of input reads"])+'\t'+str(log_infos["Total mapped"])+str(" (")+str(log_infos["Total mapped percentage"])+str("%)")+'\t'+str(log_infos["Uniquely mapped"])+str(" (")+str(log_infos["Uniquely mapped percentage"])+str("%)")+'\t'+str(log_infos["Multi mapped"])+str(" (")+str(log_infos["Multi mapped percentage"])+str("%)")+'\t'+str(log_infos["Total unmapped"])+str(" (")+str(log_infos["Total unmapped percentage"])+str("%)")+'\n')

	
	
if __name__ == "__main__":
	
	if '--path_to_log_files' and '--out' in sys.argv :
		main( sys.argv )
	else:
		sys.exit( __usage__ )
