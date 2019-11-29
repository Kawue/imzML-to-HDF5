# -*- coding: UTF-8
#
try:
    from lxml import etree
except ImportError:
    import xml.etree.ElementTree as etree
import collections
import os.path


Peak = collections.namedtuple('Peak', ('min_mass', 'max_mass'))

class MirParser:

  def __init__(self, mir_file=None):
    parser = etree.XMLParser()#recover=True encoding='utf-8'
    if os.path.isfile(mir_file):
      self._mir_root = etree.parse(mir_file, parser=parser).getroot()
    else:
      self._mir_root = etree.fromstring(mir_file, parser=parser)


  def parse_peaks(self):
    """Parse peaks from the .mir file.

    >>> example_mir = '''
    ... <ImagingResults flexImagingVersion="4.0.32.188" last_modified="2016-06-16T13:46:38">
    ...   <Result Type="PkFilter" Name="01" Color="#0000ff" Show="0" MinIntensity="0" IntensityThreshold="60.0022" AbsIntens="0" LogScale="0" MinMass="24.535144" MaxMass="24.735144" Integrate="0" UsePercentile="1" IntensityPercentile="0.95" FindMass="0" RelMass="0"></Result>
    ...   <Result Type="PkFilter" Name="02" Color="#00ff00" Show="0" MinIntensity="4.52205" IntensityThreshold="50.3247" AbsIntens="0" LogScale="0" MinMass="40.1662" MaxMass="40.3662" Integrate="0" UsePercentile="0" IntensityPercentile="0.95" FindMass="0" RelMass="0"></Result>
    ...   <Result Type="PkFilter" Name="03" Color="#ff0000" Show="0" MinIntensity="0" IntensityThreshold="10.693" AbsIntens="0" LogScale="0" MinMass="59.00335" MaxMass="59.20335" Integrate="0" UsePercentile="0" IntensityPercentile="0.95" FindMass="0" RelMass="0"></Result>
    ... </ImagingResults>'''
    >>> m = MirParser(example_mir)
    >>> m.parse_peaks()
    [Peak(min_mass='24.535144', max_mass='24.735144'), Peak(min_mass='40.1662', max_mass='40.3662'), Peak(min_mass='59.00335', max_mass='59.20335')]
    """
    peaks = []
    if self._mir_root.tag == 'ImagingResults':
      for peak in self._mir_root.findall("Result[@Type='PkFilter']"):
        peak_attributes  = dict(peak.items())
        peaks.append(Peak(float(peak_attributes['MinMass']), float(peak_attributes['MaxMass'])))
    return sorted(peaks, key=lambda x: x[0])

if __name__ == "__main__":
    import doctest
    doctest.testmod()

