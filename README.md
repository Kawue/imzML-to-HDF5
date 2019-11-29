# Description
Parser to convert MSI imzML files into HDF5 files.

# Preconditions
Python is needed for execution. We recomment to install anaconda from https://www.anaconda.com/distribution/#download-section .
Afterwards, use the predefined anaconda environment. Therefore, open your console, navigate into the imzML-to-HDF5 package and call `create -f environment.yml`.
Activate the environment with `activate msiparse` (Windows users should use the CMD instead of powershell). Activation needs to be done every time before the parser is called.
Deactivate the environment with `deactivate`.

# Execution
1. Call `activate msiparse`.
2. Call `python msiproc.py` to convert from imzml + ibd to h5 data. This requires a lot of memory for bigger imzml/ibd files and should be executed on a cluster. For parameter details call `msiproc.py -h`

## Optional
  (a). Provide a mir file (peak selection file) to reduce data sets generated to the specified peak selection.

  (b). Call `python consensus_peaklist.py` to align multiple **peak selected** data sets. Peaks with a user specified percentual overlap will be merged. For parameter details call `python consensus_peaklist.py -h`. <br> **Important:** This method will only change the representative *m/z* value of peaks. The intensity will remain as the sum over the original peak range, i.e. a shifting will be simulated without changing the peaks in the original data. The new peak values should not be used in relation with the original data set!

  (c). Provide a mis file to generate position mapping images and scaled images.

Hints and issues:
 * The following files have to be named equally: imzml, ibd(, mis).
 * Computing consensus peaks can cause some issues with due to peak chaining effects and multiple peak overlapping, e.g. two peaks of D<sub>1</sub> overlap with one peak of D<sub>2</sub>. <br> None of these Problems should raise with enough spacing between selected peaks.
