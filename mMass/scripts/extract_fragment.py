import sys
sys.path.append("../")
import mspy
import numpy as np
import pylab as pl

#Load mzml file
parser = mspy.parseMZML(
    "/Volumes/FAT/011114/IMC_digest_test04_1_100_pepsin.mzML")

parser.load()

ScanList = parser.scanlist()
ms1ScanList = []
profiles = {}

for key, scan_info in ScanList.iteritems():
    if scan_info['msLevel'] == 1:
        ms1ScanList.append(key)
        profiles[key] = parser.scan(key).profile
#Get fragment parameters


parser = None

fragment_list = ["SLSPHRPRHSRLQREPQVQWL"]
charge_min = 1
charge_max = 6


sequence_obj = mspy.sequence(fragment_list[0])
pattern_obj = mspy.pattern(sequence_obj.formula(),
                           charge=4, real=False)
#Perform calculations
time = []
basepeak = []
rmsd = []
profile = None
for scan_number in ms1ScanList:
    time.append(ScanList[scan_number]['retentionTime'])
    checkPatternResult = mspy.checkpattern(signal=profiles[scan_number],
                                           pattern=pattern_obj)
    if checkPatternResult is not None:
        rmsd.append(checkPatternResult.rmsd)
        basepeak.append(checkPatternResult.basepeak)
        if rmsd[-1] < 0.15:
            if profile is None:
                profile = mspy.crop(profiles[scan_number],
                                    pattern_obj[0][0] - 0.1,
                                    pattern_obj[-1][0] + 0.1)
            else:
                profile = mspy.combine(profile, mspy.crop(
                    profiles[scan_number], pattern_obj[0][0] - 0.1,
                    pattern_obj[-1][0] + 0.1))
    else:
        rmsd.append(1)
        basepeak.append(0)
profile = mspy.reduce(profile)

pl.plot(time, rmsd)
pl.show()
pl.plot(time, basepeak)
pl.show()
pl.plot(*profile.T)
pl.show()
#Plot
