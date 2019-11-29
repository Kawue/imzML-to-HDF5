import numpy as np
import pandas as pd
from sys import argv
from os.path import join, normpath, exists
from os import makedirs
import argparse
from pyimzml.ImzMLWriter import ImzMLWriter

def convert(h5file, outpath, filename):
    data = pd.read_hdf(h5file)
    if data.index.names[0] != "grid_x" or data.index.names[1] != "grid_y":
        raise ValueError("Wrong index in h5, contact programmer!")

    with ImzMLWriter(join(outpath, filename)) as writer:
        for coords, series in data.iterrows():
            gx,gy,rx,ry = coords
            mzs = np.array(series.index)
            intensities = series.values
            writer.addSpectrum(mzs, intensities, (gx,gy))


#h5file = argv[1]
#outpath = argv[2]
#filename = argv[3]
#h5file = "C:\\Users\\kwuellems\\Desktop\\Orbi\\test\\20180705_AF_MV_DHB_RP_163_6_processed_norm.h5"
#outpath = "C:\\Users\\kwuellems\\Desktop\\Orbi\\test"
#filename = "output.imzML"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hdf5 to imzml.')
    parser.add_argument('-d', action='store', dest='h5file',
                        help='Path to the input hdf5 file.')
    parser.add_argument('-o', action='store', dest='outpath',
                        help='Path to the output directory.')
    parser.add_argument('-n', action='store', dest='filename',
                        help='Name of the output file.')
    args = parser.parse_args()
    try:
        # normalizing all paths
        args.h5file = normpath(args.h5file)
        args.outpath = normpath(args.outpath)
        if not exists(args.outpath):
            makedirs(args.outpath)
    except:
        parser.print_help()
    else:
        convert(args.h5file, args.outpath, args.filename)
