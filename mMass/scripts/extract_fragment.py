import sys
sys.path.append("../")
import mspy
import numpy as np
import pylab as pl
import time as timer

#Load mzml files

t0 = timer.time()
path = "/Users/johannes/SHINDELAB/MassSpec/Johannes/011114_FIMC_test_bottom/"

files = ["FIMC_digest_test02_1_20_pepsin.mzML",
         "FIMC_digest_test03_1_50_pepsin.mzML",
         "FIMC_digest_test04_1_100_pepsin.mzML"]


ScanList = {}
ms1ScanList = {}
profiles = {}


for file_iter in files:

    parser = mspy.parseMZML(path + file_iter)

    parser.load()

    ScanList[file_iter] = parser.scanlist()
    ms1ScanList[file_iter] = []
    profiles[file_iter] = {}

    for key, scan_info in ScanList[file_iter].iteritems():
        if scan_info['msLevel'] == 1:
            ms1ScanList[file_iter].append(key)
            profiles[file_iter][key] = parser.scan(key).profile

t1 = timer.time() - t0

print 'Loaded files in %s ' % t1

#Get fragment parameters

t0 = timer.time()

fragment_list = ["YHFWHRGVT"]
charge_min = 1
charge_max = 6


sequence_obj = mspy.sequence(fragment_list[0])
pattern_obj = mspy.pattern(sequence_obj.formula(),
                           charge=3, real=False)
#Perform calculations
time = {}
basepeak = {}
rmsd = {}
profile = {}
for file_iter in files:
    time[file_iter] = []
    basepeak[file_iter] = []
    rmsd[file_iter] = []
    profile[file_iter] = None
    for scan_number in ms1ScanList[file_iter]:
        time[file_iter].append(ScanList[file_iter][scan_number]['retentionTime'])
        checkPatternResult = mspy.checkpattern(signal=profiles[file_iter][scan_number],
                                            pattern=pattern_obj)
        if checkPatternResult is not None:
            rmsd[file_iter].append(checkPatternResult.rmsd)
            basepeak[file_iter].append(checkPatternResult.basepeak)
            if rmsd[file_iter][-1] < 0.15:
                if profile[file_iter] is None:
                    profile[file_iter] = mspy.crop(profiles[file_iter][scan_number],
                                        pattern_obj[0][0] - 0.1,
                                        pattern_obj[-1][0] + 0.1)
                else:
                    profile[file_iter] = mspy.combine(profile[file_iter], mspy.crop(
                        profiles[file_iter][scan_number], pattern_obj[0][0] - 0.1,
                        pattern_obj[-1][0] + 0.1))
        else:
            rmsd[file_iter].append(1)
            basepeak[file_iter].append(0)
    if profile[file_iter] is not None:
        profile[file_iter] = mspy.reduce(profile[file_iter])


t1 = timer.time() - t0

print 'Extracted data in %s ' % t1

for file_iter in files:
    pl.plot(time[file_iter], rmsd[file_iter], label = file_iter)
pl.legend()
pl.show()
for file_iter in files:
    pl.plot(time[file_iter], basepeak[file_iter], label = file_iter)
pl.legend()
pl.show()
for file_iter in files:
    if profile[file_iter] is not None:
        pl.plot(*profile[file_iter].T,label = file_iter)
pl.legend()
pl.show()
#Plot
