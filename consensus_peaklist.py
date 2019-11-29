import numpy as np
import pandas as pd
from os import path
import argparse
from interval import interval
from itertools import chain
from mir_parser import MirParser

# Pass tuples of peak ranges (.mir) and full dataframes(.h5)
def find_consensus_peaks(t, *args):
	# Notice current problems
	print("Chaining problems are not resolved yet!")
	print("Multiple overlaps are resolved by revoking any kind of merge, this is a temporary solution!")

	# Allow decimal and percent values
	if t > 1:
		t = t / 100

	# Parse Peak Ranges:
	parsed_mir = []
	for mir_path, _ in args:
		parsed_mir.append(MirParser(mir_path).parse_peaks())
		
	# Parse Dataframes:
	parsed_frames = []
	for _, frame_path in args:
		parsed_frames.append(pd.read_hdf(frame_path))

	# Dict for future peak lists
	peak_dict = dict((key,[]) for key in range(len(args)))

	# Dict to determine dataset membership
	peak_ranges_dict = dict((idx, peak_ranges) for idx, peak_ranges in enumerate(parsed_mir))

	joint_peak_ranges = list(chain(*parsed_mir))
	joint_peak_ranges.sort(key= lambda x: x[0])

	# Use reduces joint peak ranges list to reduce search time in find_max_overlap()
	reduced_joint_peak_ranges = joint_peak_ranges[:] # ?: not necessary because of index usage?
	
	# Initialize candidates as dict to ensure correkt if statement at the beginning
	# Candidates will change between dict and list a few times
	candidates = dict()
	
	# Calculate the overlap between ranges and join them if they exceed a defined threshold
	for idx, peak_ref in enumerate(joint_peak_ranges):
		if peak_ref in chain(*candidates.values()): # skip peaks that were already merged durin the maximum overlap detection
			continue # ?: should work because of ordered
		else:
			candidates = []
			ref = interval[peak_ref.min_mass, peak_ref.max_mass]
			for peak_compare in joint_peak_ranges[idx:]: # include self peak
				compare = interval[peak_compare.min_mass, peak_compare.max_mass]
				pct_ref, pct_compare = calc_percent_overlap(ref, compare)
				if pct_ref > t or pct_compare > t:
					candidates.append(peak_compare)
				else:
					# check if current candidate list can be extendet by spanning overlaps
					candidates = find_max_overlap(candidates, reduced_joint_peak_ranges[idx:], t) # ?: should work because of ordered # candidates is list
					candidates = check_for_split(candidates, t) # candidates is dict
					break
			else:
				# if the last iteration of the for loop is an overlap, candidates will remain as list but needs to be a dict
				# ?: Changing candidates list into simple dict may be enough, check if find_max_overlap() and check_for_split() is needed!
				if type(candidates) != dict:
					candidates = find_max_overlap(candidates, reduced_joint_peak_ranges[idx:], t) # ?: is this needed? # candidates is list
					#candidates = check_for_split(candidates, t) # candidates is dict
					candidates = check_for_revoke(candidates, t) # revokes any merge in case of chaining or multiple overlap, e.g. D1: x & y overlap with D2: z 
		# ?: check if double entries will occur -> should not because of starting if statement
		# Define new peak based on midpoint of all candidates
		# ?: weight by intensity could improve this
		for key, value in candidates.items():
			peak = define_peak(value)
			# add new peak to each data set that added at least one candidate
			for peak_ranges in value:
				for k, v in peak_ranges_dict.items():
					if peak_ranges in v:
						peak_dict[k].extend([peak])
	
	# ?: elemination of double new peaks necessary?
	'''
	#eleminate double entries
	for key, value in peak_dict:
		### test area ###
		test = list(set(value))[:]
		test.sort()
		if test != list(set(value)):
			print("not sorted!!!")
		### test end ###
		peak_dict[key] = list(set(value))
	'''

	for key, value in peak_dict.items(): # ?: should work because of ordered
		parsed_frames[key].columns = value



	for idx, frame in enumerate(parsed_frames):
		key = ""
		frame_path = args[idx][1]
		frame_path = frame_path.replace(path.sep, "/")
		with pd.HDFStore(frame_path) as f:
			key = f.keys()[0][1:]
		with pd.HDFStore(path.splitext(frame_path)[0] + "_aligned_" + str(t*100) + "pct.h5", complib='blosc', complevel=9) as store:
			store[key] = frame





# Calculate how many percent the overlap spans in each interval
def calc_percent_overlap(interval_a, interval_b):
	if interval_a & interval_b:
		span_a = interval_a[0].sup - interval_a[0].inf
		span_b = interval_b[0].sup - interval_b[0].inf
		overlap = interval_a & interval_b
		span_o = overlap[0].sup - overlap[0].inf
		pct_a = (span_o/span_a)
		pct_b = (span_o/span_b)
		return pct_a, pct_b
	else:
		return 0, 0

# Define peak based on a overlapping peak ranges
def define_peak(peak_ranges):
	min_peak = min([peak.min_mass for peak in peak_ranges])
	max_peak = max([peak.max_mass for peak in peak_ranges])
	return round(interval[min_peak, max_peak].midpoint[0][0], 3)

# check for each candidate in candidates if a candidate can extend the overlap by spanning across multiple disconnected peak ranges
def find_max_overlap(candidates, joint_peak_ranges, t):
	to_iter = candidates
	for candidate in to_iter:
		ref = interval[candidate.min_mass, candidate.max_mass]
		for peak_compare in joint_peak_ranges: # reduce compare by reducing to iter
			compare = interval[peak_compare.min_mass, peak_compare.max_mass]
			pct_ref, pct_compare = calc_percent_overlap(ref, compare)
			if pct_ref > t or pct_compare > t:
				if peak_compare not in candidates:
					candidates.append(peak_compare)
	if candidates != to_iter:
		candidates = find_max_overlap(candidates, joint_peak_ranges, t)
	return candidates

def check_for_split(candidates, t):
	split_dict = dict()
	key_counter = 0
	for idx, candidate in enumerate(candidates):
		ref = interval[candidate.min_mass, candidate.max_mass]
		for candidate_compare in candidates[idx:]: # ?: should work because of ordered
			compare = interval[candidate_compare.min_mass, candidate_compare.max_mass]
			pct_ref, pct_compare = calc_percent_overlap(ref, compare)
			if pct_ref < t and pct_compare < t:
				split_dict[key_counter] = [candidate]
				key_counter += 1
				split_dict[key_counter] = [candidate_compare]
				key_counter += 1
	if len(split_dict): # If no split happend, put the whole candidate array in one split dict entry
		for key, value in split_dict.items():
			ref = interval[value[0].min_mass, value[0].max_mass]
			for candidate_compare in candidates: # ?: should work because of ordered
				compare = interval[candidate_compare.min_mass, candidate_compare.max_mass]
				if compare != ref:
					pct_ref, pct_compare = calc_percent_overlap(ref, compare)
					if pct_ref > t or pct_compare > t:
						split_dict[key].extend([candidate_compare])
	else:
		split_dict[key_counter] = candidates
	return split_dict
	
def check_for_revoke(candidates, t):
	split_dict = dict()
	key_counter = 0
	for idx, candidate in enumerate(candidates):
		ref = interval[candidate.min_mass, candidate.max_mass]
		for candidate_compare in candidates[idx:]: # ?: should work because of ordered
			compare = interval[candidate_compare.min_mass, candidate_compare.max_mass]
			pct_ref, pct_compare = calc_percent_overlap(ref, compare)
			if pct_ref < t and pct_compare < t:
				split_dict[key_counter] = [candidate]
				key_counter += 1
				split_dict[key_counter] = [candidate_compare]
				key_counter += 1
	if len(split_dict): # If no split happend, put the whole candidate array in one split dict entry
		for key, value in split_dict.items(): # If split happend, do not merge anything and use mean for each individual value 
			split_dict[key] = [round(np.mean([value[0].min_mass, value[0].max_mass]),3)]
	else:
		split_dict[key_counter] = candidates
	return split_dict
				
				
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Align selected peak ranges of multiple data sets.")
	parser.add_argument("-t", type=float, action="store", dest="t", help="Percentual overlap between two selected ranges.")
	parser.add_argument("-p", nargs='+', action="store", dest="paths", help="Tuples of paths of form: file1.mir file1.h5 file2.mir file2.h5, ...")
	args = parser.parse_args()
	try:
		path_tuples = list(zip(args.paths[::2], args.paths[1::2]))
		find_consensus_peaks(args.t, *path_tuples)
	except:
		parser.print_help()


'''
print("#####")
mir_1 = "C:\\Users\\kwuellems\\Desktop\\test5\\20180221_healthy_gs_AF.mir"
mir_2 = "C:\\Users\\kwuellems\\Desktop\\test5\\20180221_PXE_gs_AF.mir"

frame_1 = "C:\\Users\\kwuellems\\Desktop\\test5\\skin_healthy_processed_picked\\20171115_DHB_RP_HS_healthy_3_oversampling_AF_processed.h5"
frame_2 = "C:\\Users\\kwuellems\\Desktop\\test5\\skin_pxe_processed_picked\\20171115_DHB_RP_HS_PXE_1_oversampling_AF_processed.h5"

print("begin")
find_consensus_peaks(0.5, (mir_1, frame_1), (mir_2, frame_2))
print("end")
print("###################")
'''