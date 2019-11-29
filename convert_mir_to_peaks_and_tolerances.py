# -*- coding: UTF-8
#
import argparse
import numpy as np
import os.path
from mir_parser import MirParser


def parse_peaks_and_tolerances_from_mir(mir_path, out_path):
    peaks = sorted(MirParser(mir_path).parse_peaks())
    peak_list_name = os.path.basename(mir_path).split('.mir')[0]
    peak_file = os.path.join(out_path, '{}_peaks.txt').format(peak_list_name)
    np.savetxt(peak_file, np.around(np.mean(peaks, axis=1), decimals=3), delimiter=',')
    tolerance_file = os.path.join(out_path, '{}_tolerance.txt').format(peak_list_name)
    np.savetxt(tolerance_file, peaks, delimiter=',')
    print('Saved output to {} and {}.'.format(peak_file, tolerance_file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parse mir file (Bruker specfic) and peak and tolerance lists.')
    parser.add_argument('-i', action='store', dest='mir_path',
                        help='Path to mir file (Bruker specific).')
    parser.add_argument('-o', action='store', dest='out_path',
                        help='Path to the output directory.')

    args = parser.parse_args()
    try:
        parse_peaks_and_tolerances_from_mir(args.mir_path, args.out_path)
    except:
        parser.print_help()
