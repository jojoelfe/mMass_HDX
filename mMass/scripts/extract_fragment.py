import sys
sys.path.append("../")
import mspy
import numpy as np
import pylab as pl
from matplotlib.collections import LineCollection

#Load mzml file
parser = mspy.parseMZML(
    "C:\Users\jojot_000\Documents\\011114\FIMC_digest_test04_1_100_pepsin.mzML")

parser.load()

ScanList = parser.scanlist()
ms1ScanList = []

for key, scan_info in ScanList.iteritems():
    if scan_info['msLevel'] == 1:
        ms1ScanList.append(key)
#Get fragment parameters
fragment = "SLSPHRPRHSRLQREPQVQWL"
charge = 4


sequence_obj = mspy.sequence(fragment)
pattern_obj = mspy.pattern(sequence_obj.formula(),
                           charge=charge, real=False)
#Perform calculations
time = []
basepeak = []
rmsd = []
profile = None
for scan_number in ms1ScanList:
    time.append(ScanList[scan_number]['retentionTime'])
    if mspy.checkpattern(signal=parser.scan(scan_number).profile,
                         pattern=pattern_obj) is not None:
        rmsd.append(mspy.checkpattern(signal=parser.scan(scan_number).profile,
                                      pattern=pattern_obj).rmsd)
        basepeak.append(mspy.checkpattern(
            signal=parser.scan(scan_number).profile,
            pattern=pattern_obj).basepeak)
        if rmsd[-1] < 0.15:
            if profile is None:
                profile = mspy.crop(parser.scan(scan_number).profile,
                                    pattern_obj[0][0] - 0.1, pattern_obj[-1][0] + 0.1)
            else:
                profile = mspy.combine(profile, mspy.crop(
                    parser.scan(scan_number).profile, pattern_obj[0][0] - 0.1,
                    pattern_obj[-1][0] + 0.1))
    else:
        rmsd.append(1)
        basepeak.append(0)
profile = mspy.reduce(profile)

pl.plot(time,rmsd)
pl.show()

pl.show()
pl.plot(*profile.T)
pl.show()
#Plot
