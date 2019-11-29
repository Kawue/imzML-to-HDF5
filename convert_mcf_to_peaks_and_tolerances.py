# -*- coding: UTF-8
#
import numpy as np
import os
from mcl_parser import MclParser

project_path = '../data/wuellems_paper_072017/17062016_N5203_MAPS_neg'
dataset_name = os.path.basename(project_path)
mcl_path = os.path.join(project_path, 'glio rn 100 3203.mcl')

mcl_peaks = MclParser(mcl_path).parse_peaks()
peaks = np.array([(p, (p - p * t/10**6), (p + p * t/10**6)) for p, t in mcl_peaks])

peak_file = os.path.join(project_path, '{}_peaks.txt').format(dataset_name)
np.savetxt(peak_file, peaks[:, 0], delimiter=',')
tolerance_file = os.path.join(project_path, '{}_tolerance.txt').format(dataset_name)
np.savetxt(tolerance_file, peaks[:, [1, 2]], delimiter=',')
