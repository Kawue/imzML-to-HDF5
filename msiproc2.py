# -*- coding: UTF-8
#
import argparse
import os
import os.path
import re
import cv2
import xml.etree.ElementTree as etree
import numpy as np
import pandas as pd
from collections import namedtuple
from pyimzml.ImzMLParser import ImzMLParser
from peakselection import select_peaks_from_msi_frame


def register_gridcoords_with_image(path_to_imzml_file, path_to_mis_file, out_path):
	"""
	Generates and saves images showing the coordinates of the measurements for each analysis from
	Bruker. Results are saved as png in the running directory with the name of the analysis and a
	fixed postfix.

	The results give information on the basic shape of the analyzed sectors and the resolution of
	the analyses (at least within the coordinate system in use). Each point in the pictures
	represents one m/z spectrum.
	"""

	mis_tree = etree.parse(path_to_mis_file)

	teach_points_yx_px = []
	teach_points_yx_um = []

	# Match two pairs of comma seperated coordinates.
	# The pairs themself are seperated by a semicolon.
	# The first coordinate is measured in px, the second in um.
	# For example: '5,2;800,600'
	coord_regex = re.compile(r'(\d+),(\d+);(-?\d+),(-?\d+)')
	for teach_point_xml in mis_tree.iter('TeachPoint'):
		coord_match = coord_regex.match(teach_point_xml.text)
		if coord_match:
			teach_point = coord_match.groups()
			teach_points_yx_px.append([int(teach_point[1]), int(teach_point[0])])
			teach_points_yx_um.append([int(teach_point[3]), int(teach_point[2])])

	for image_file_xml in mis_tree.iter('ImageFile'):
		image_file = os.path.join(os.path.dirname(path_to_mis_file), image_file_xml.text)

	if image_file:
		image = cv2.imread(image_file)
		image_shape_yx_px = np.array(image.shape[0:2], dtype=float)
		teach_points_yx_relative = teach_points_yx_px / image_shape_yx_px
		distances_yx_relative = np.abs(teach_points_yx_relative -
									   np.roll(teach_points_yx_relative, -1, axis=0))
		distances_yx_um = np.abs(teach_points_yx_um - np.roll(teach_points_yx_um, -1, axis=0))
		max_distance_rel_y, max_distance_rel_x = distances_yx_relative.max(axis=0)
		max_distance_um_y, max_distance_um_x = distances_yx_um.max(axis=0)
		um_height = max_distance_um_y / max_distance_rel_y
		um_width = max_distance_um_x / max_distance_rel_x
		image_shape_yx_um = np.array([um_height, um_width], dtype=float)

	coord_regex = re.compile(r'(\d+),(\d+)')
	for raster_xml in mis_tree.iter('Raster'):
		coord_match = coord_regex.match(raster_xml.text)
		raster_width_um, raster_height_um = coord_match.groups()
		raster_shape_yx_um = np.array([raster_height_um, raster_width_um], dtype=float)

	raster_shape_yx_px = raster_shape_yx_um / image_shape_yx_um * image_shape_yx_px
	raster_shape_yx_grid_coord = image_shape_yx_px / raster_shape_yx_px
	print('raster_shape_yx_px', raster_shape_yx_px)
	print('raster_shape_yx_grid_coord', raster_shape_yx_grid_coord)

	# using the imzML meta data
	p = ImzMLParser(path_to_imzml_file, parse_lib='ElementTree')
	Position = namedtuple('Position', ['x', 'y'])
	coord_pairs = []
	for coordinate in p.coordinates:
		coord_pairs.append(Position(x=coordinate[0], y=coordinate[1]))

	print('image.shape', image.shape)
	coord_mapping = np.zeros(shape=image.shape, dtype=int)
	coord_image = np.copy(image)
	coord_array = []
	for idx, coord in enumerate(coord_pairs):
		y = int(raster_shape_yx_px[0] + raster_shape_yx_px[0] * coord.y)
		x = int(raster_shape_yx_px[1] + raster_shape_yx_px[1] * coord.x)
		coord_image[y, x] = [0, 0, 255]
		coord_mapping[y, x] = [idx, coord.y, coord.x]
		coord_array.append([y, x])
	print(coord_mapping[:, :, 0].max(), coord_mapping[:, :, 1].max(), coord_mapping[:, :, 2].max())

	coord_array = np.array(coord_array)
	min_y, min_x = coord_array.min(axis=0)
	max_y, max_x = coord_array.max(axis=0)

	min_y -= int(raster_shape_yx_px[0] / 2)
	min_x -= int(raster_shape_yx_px[0] / 2)
	max_y += int(raster_shape_yx_px[1] / 2)
	max_x += int(raster_shape_yx_px[1] / 2)

	base_path = os.path.join(os.path.dirname(path_to_mis_file), out_path)
	mis_name, _ = os.path.splitext(os.path.basename(path_to_mis_file))
	raw_img_path = os.path.join(out_path, '{}_raw.png'.format(mis_name)) #os.path.join(base_path, '{}_raw.png'.format(mis_name))
	coord_img_path = os.path.join(out_path, '{}_coord.png'.format(mis_name)) #os.path.join(base_path, '{}_coord.png'.format(mis_name))

	trunc_raw_img_path = os.path.join(out_path, '{}_trunc_raw.png'.format(mis_name)) #os.path.join(base_path, '{}_trunc_raw.png'.format(mis_name))
	trunc_coord_img_path = os.path.join(out_path, '{}_trunc_coord.png'.format(mis_name)) #os.path.join(base_path, '{}_trunc_coord.png'.format(mis_name))
	cv2.imwrite(trunc_raw_img_path, image[min_y:max_y, min_x:max_x])
	cv2.imwrite(trunc_coord_img_path, coord_image[min_y:max_y, min_x:max_x])

	interpolation_methods = {
		# 'INTER_NEAREST': cv2.INTER_NEAREST,
		# 'INTER_LINEAR': cv2.INTER_LINEAR,
		'INTER_AREA': cv2.INTER_AREA,
		# 'INTER_CUBIC': cv2.INTER_CUBIC,
		# 'INTER_LANCZOS4': cv2.INTER_LANCZOS4
	}

	for i_name, i_method in interpolation_methods.items():
		i_name = str.lower(i_name)
		print('Resizing with', i_name)
		trunc_raw_img_path = os.path.join(out_path, '{}_trunc_raw_{}.png'.format(mis_name, i_name)) #os.path.join(base_path, '{}_trunc_raw_{}.png'.format(mis_name, i_name))
		trunc_coord_img_path = os.path.join(
			out_path, '{}_trunc_coord_{}.png'.format(mis_name, i_name))  #trunc_coord_img_path = os.path.join(base_path, '{}_trunc_coord_{}.png'.format(mis_name, i_name))
		fy = 1 / raster_shape_yx_px[0]
		fx = 1 / raster_shape_yx_px[1]
		print('by', ':', fx, '/', fy)
		trunc_raw_img = cv2.resize(
			image[min_y:max_y, min_x:max_x], (0, 0), fy=fy, fx=fx, interpolation=i_method)
		cv2.imwrite(trunc_raw_img_path, trunc_raw_img)
		trunc_coord_img = cv2.resize(
			coord_image[min_y:max_y, min_x:max_x], (0, 0), fy=fy, fx=fx, interpolation=i_method)
		cv2.imwrite(trunc_coord_img_path, trunc_coord_img)

	coord_map_path = os.path.join(out_path, '{}_coordmap'.format(mis_name)) #os.path.join(base_path, '{}_coordmap'.format(mis_name))
	print(coord_img_path)
	print()
	print()
	cv2.imwrite(raw_img_path, image)
	cv2.imwrite(coord_img_path, coord_image)
	np.savez_compressed(coord_map_path, coordmap=coord_mapping)


def intensities_generator(imzmlParser, mz_index, selection=slice(None)):
	for i in range(len(imzmlParser.coordinates)):
		yield imzmlParser.getspectrum(i)[1][selection]


def imzml_to_hdf5(imzml_file_path, out_path, mir_path):

	dataset_name, _ = os.path.splitext(os.path.basename(imzml_file_path))

	print()
	print('Loading', imzml_file_path)
	p = ImzMLParser(imzml_file_path, parse_lib='ElementTree')
	print()
	print('Loading done!')

	# check if all spectra have the same mz axis
	num_spectra = len(p.mzLengths)
	mz_index = np.array(p.getspectrum(0)[0])
	mz_index_length = len(mz_index)
	print()
	print('m/z consistency check ...')

	# '0' = mz values, '1' = intensities
	mz_index = np.unique(np.concatenate([p.getspectrum(i)[0] for i in range(num_spectra)]))
	
	if len(mz_index) != mz_index_length:
		print('WARNING: Not all spectra have the same mz values. Missing values are filled with zeros!')

	print()
	print('m/z consistency check done!')
	
	# DEV: use small range to test bigger datasets on little memory
	mz_selection = slice(None) # range(100)
	# load all intensities into a single data frame
	# resulting format:
	#   1 row = 1 spectrum
	#   1 column = all intensities for 1 mz, that is all values for a single intensity image
	print()
	print('DataFrame creation ...')
	msi_frame = pd.DataFrame(intensities_generator(p, mz_index, mz_selection), columns=mz_index[mz_selection])
	print('DataFrame creation done')
	print()
	print("DataFrame size equals: %i pixels, %i mz-values"%msi_frame.shape)
	print()

	if mir_path:
		print()
		print('Peak picking ...')
		msi_frame = select_peaks_from_msi_frame(msi_frame, mir_path)
		print()
		print('Peak picking done!')

	msi_frame = msi_frame.fillna(0)

	xycoordinates = np.asarray(p.coordinates)[:,[0,1]]
	multi_index = pd.MultiIndex.from_arrays(xycoordinates.T, names=("grid_x", "grid_y"))
	msi_frame.set_index(multi_index, inplace=True)

	msi_frame["dataset"] = [dataset_name]*msi_frame.shape[0]
	msi_frame = msi_frame.set_index("dataset", append=True)

	# For some data sets a small fraction of intensities (~0.1%) have been
	# negative, this might be a numerical issue in the imzml export by bruker.
	# DEV ad-hoc fix (couldn't figure out the cause or a more reasonable fix so far)
	msi_frame[msi_frame < 0] = 0

	print()
	print('Write DataFrame ...')
	h5_store_path = os.path.join(out_path, dataset_name + '.h5')
	save_name_frame = 'msi_frame_' + dataset_name
	with pd.HDFStore(h5_store_path, complib='blosc', complevel=9) as store:
		store[save_name_frame] = msi_frame
	print()
	print('done. Script completed!')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Convert imzml to hdf5.')
	parser.add_argument('-i', action='store', required=True, dest='imzml_path',
						help='Path to the input imzML file (required).')
	parser.add_argument('-o', action='store', required=True, dest='out_path',
						help='Path to the output directory (required).')
	parser.add_argument('-p', action='store', required=False, dest='mir_path',
						help='Path to the input mir file (peak selection) (optional).')
	parser.add_argument('-m', action='store', required=False, dest='mis_path',
						help='Path to the input mis file (Bruker specific) (optional).')

	args = parser.parse_args()

	try:
		# normalizing all paths
		args.imzml_path = os.path.normpath(args.imzml_path)
		args.out_path = os.path.normpath(args.out_path)

		if not os.path.exists(args.out_path):
			os.makedirs(args.out_path)

	except:
		parser.print_help()
	else:
		try:
			args.mir_path = os.path.normpath(args.mir_path)
		except:
			print("")
			print("")
			print("No peak list specified!")
			print("Conversion will be based on the full detection range.")
			print("")
			print("")
		
		if args.mis_path:
			args.mis_path = os.path.normpath(args.mis_path)
			register_gridcoords_with_image(args.imzml_path, args.mis_path, args.out_path)

		imzml_to_hdf5(args.imzml_path, args.out_path, args.mir_path)
