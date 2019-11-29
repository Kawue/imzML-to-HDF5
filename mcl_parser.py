# -*- coding: UTF-8
#
try:
    from lxml import etree
except ImportError:
    import xml.etree.ElementTree as etree
import collections
import os.path


Peak = collections.namedtuple('Peak', ('mz', 'tolerance'))


class MclParser(object):

    def __init__(self, mcl_file=None):
        parser = etree.XMLParser()#recover=True encoding='utf-8'
        if os.path.isfile(mcl_file):
            self._mir_root = etree.parse(mcl_file, parser=parser).getroot()
        else:
            self._mir_root = etree.fromstring(mcl_file, parser=parser)


    def parse_peaks(self):
        """Parse peaks from the .mcl file.

        >>> example_mir = '''
        ... <method MethodName="glio rn 100 3203" Method_ID="1234567890">
        ...    <peak_list>
        ...        <peak background="false" calibrant="true" label="01" lock="false" m_z="147.023450" tolerance="680.164"/>
        ...        <peak background="false" calibrant="true" label="02" lock="false" m_z="89.058404" tolerance="1122.86"/>
        ...        <peak background="false" calibrant="true" label="03" lock="false" m_z="129.030900" tolerance="775.008"/>
        ...    </peak_list>
        ... </method>'''
        >>> m = MclParser(example_mir)
        >>> m.parse_peaks()
        [Peak(mz=89.058404, tolerance=1122.86), Peak(mz=129.0309, tolerance=775.008), Peak(mz=147.02345, tolerance=680.164)]
        """
        peaks = []
        if self._mir_root.tag == 'method':
            for peak in self._mir_root[0].findall("peak"):
                p = dict(peak.items())
                peaks.append(Peak(float(p['m_z']), float(p['tolerance'])))
        return sorted(peaks, key=lambda x: x[0])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
