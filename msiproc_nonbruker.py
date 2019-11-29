import os
from sys import argv
import numpy as np
import pandas as pd
from pyimzml.ImzMLParser import ImzMLParser
import matplotlib.pyplot as plt

#Scheint soweit zu funktionieren, zumindest mit dem bladder datensatz.
#ggf. nochmal das binning angucken und mit anderem datensatz testen.
#ansonsten average spektrum exportieren und picken lassen mit mMass oder selber picken.

### This script is a very rough starting point. It consists only of very simple rebinning and averaging!
#TODO: Make rebinning optional.

def imzml_to_hdf5(imzml_path, out_path, rebinning=True, number_bins=20000):
    p = ImzMLParser(imzml_path)
    # +1 Because this number defines the edges
    bin_number = number_bins + 1
    print(p.__dict__.keys())
    xycoordinates = np.asarray(p.coordinates)[:,[0,1]]
    spectra = np.asarray([p.getspectrum(i) for i, _ in enumerate(p.coordinates)])
    mz_all = spectra[:, 0]
    mz_all = np.concatenate(mz_all)
    mz_all = np.unique(mz_all)

    bins, mz_repr, steps = create_mz_bins(mz_all, bin_number)

    print(xycoordinates.size//2)
    print(mz_all.size)
    # es muss gebinnt werden weil die anzahl aller mz zu groß werden kann
    data = np.empty((xycoordinates.size//2, mz_repr.size))
    data[:] = np.nan
    for i, _ in enumerate(p.coordinates):
        mzs, intensities = p.getspectrum(i)
        
        #plt.figure()
        #plt.title("Original")
        #plt.plot(mzs, intensities)
        #plt.show()

        binned_intensities = bin_mz_axis(mzs, intensities, bins)
        #spectrum = np.zeros(mz_all.size)
        #indices = np.searchsorted(mz_all, mzs)
        #spectrum[indices] = intensities
        
        data[i] = binned_intensities
        
        #print(binned_intensities.size)
        #print(mz_repr.size)
        
    print(data)
    print(data.shape)
    print(np.where(data == np.nan))
    print(xycoordinates.shape)
    #print(type(xycoordinates))
    mz_repr = np.around(mz_repr, 3)
    print(mz_repr)
    #print(type(mz_repr))
    pd_index = pd.MultiIndex.from_arrays(xycoordinates.T, names=("grid_x", "grid_y")) # order correct?
    #print(pd_index)
    dataframe = pd.DataFrame(data, index=pd_index, columns=mz_repr)

    dataset_name, _ = os.path.splitext(os.path.basename(imzml_path))
    out_path = os.path.join(out_path, dataset_name + ".h5")
    frame_name = "msi_frame_" + dataset_name

    with pd.HDFStore(out_path, complib="blosc", complevel=9) as store:
        store[frame_name] = dataframe

    '''
    for i, r in dataframe.iterrows():
        plt.figure()
        plt.title("Rebinned")
        plt.plot(r.index, r)
        plt.show()
    '''
    print(dataframe)

def bin_mz_axis(mz_values, intensity_values, bins):
    intensity_sum, _ = np.histogram(mz_values, bins, weights=intensity_values)
    mz_count, _ = np.histogram(mz_values, bins)
    binned_intensities = np.divide(intensity_sum, mz_count, out=np.zeros_like(intensity_sum), where=mz_count!=0)
    return binned_intensities


def create_mz_bins(mz_spectrum, number_bins):
    bins, steps = np.linspace(mz_spectrum.min(), mz_spectrum.max(), num=number_bins, endpoint=True, retstep=True)
    representatives = (bins[:-1] + np.roll(bins, -1)[:-1]) / 2
    return bins, representatives, steps

print("Start")
imzml_to_hdf5(argv[1], argv[2])

#imzml = "C:\\Users\\kwuellems\\Desktop\\Neuer Ordner (5)\\20180613_AF_MK_DHB_RP_139_1_processed.imzML"
#imzml = "D:\\Datasets\\ToF\\test\\20180613_AF_MK_DHB_RP_139_1_processed.imzML"
#imzml = "C:\\Users\\kwuellems\\Desktop\\Orbi\\20171017_mouseVibrissae_2_DHB_MatrixTestTime_24h_AF_msrange_10_1910_processed.imzML"
#outpath = "D:\\Datasets\\ToF\\test"
#outpath = "C:\\Users\\kwuellems\\Desktop\\Neuer Ordner (5)"
#imzml_to_hdf5(imzml, outpath)
print("Finished!")

#"C:\\Users\\kwuellems\\Desktop\\Neuer Ordner\\HR2MSI mouse urinary bladder S096.imzML"
#"C:\\Users\\kwuellems\\Desktop\\Neuer Ordner"