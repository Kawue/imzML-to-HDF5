from mir_parser import MirParser
import csv
from os import path, chmod
from sys import argv

def mir_to_csv(output_path, *args):
    mir_list = []
    for p in args:
        mir_list.append(MirParser(p).parse_peaks())

    peak_dict = dict()
    for idx, mir in enumerate(mir_list):
        single_peaks = []
        for peak in mir:
            #single_peaks.append([round((peak.min_mass + peak.max_mass)/2, 7)])
			single_peaks.append([(peak.min_mass + peak.max_mass)/2])
        peak_dict[idx] = single_peaks

    idx = 0
    for key, value in peak_dict.items():
        output = path.join(output_path, path.basename(args[idx]).split(".")[0]  + ".csv")
        idx += 1
        with open(output, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter =",", lineterminator='\n')
            writer.writerows(value)

mir_to_csv(argv[1], *argv[2:])
print("All csv created.")