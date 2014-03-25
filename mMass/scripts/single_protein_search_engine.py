import sys
sys.path.append("../")
import mspy
import time as timer
import os
import pickle

path = "/Users/johannes/SHINDELAB/MassSpec/Johannes/"

filename = "031814_FIMC_digesttests/TXIII_60min"

charge_min = 1
charge_max = 8

mz_min = 250
mz_max = 1800

RmsdThreshold = 0.15

start_time = 300
stop_time = 900

sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRLQREPQVQWLEQQVAKRRTKR"

max_length = 50

name = "031814_TXIII_60_min"



def load_mzml_file(path, filename, start_time, stop_time):
    """Load MS1 scans from MZML files"""

    t0 = timer.time()

    ScanList = {}
    ms1ScanList = {}
    profiles = {}



    parser = mspy.parseMZML(path + filename + '.mzML')

    parser.load()

    ScanList = parser.scanlist()
    ms1ScanList = []
    profiles = {}

    for key, scan_info in ScanList.iteritems():
        if scan_info['msLevel'] == 1 \
          and scan_info['retentionTime'] >= start_time \
          and scan_info['retentionTime'] <= stop_time:
            ms1ScanList.append(key)
            profiles[key] = parser.scan(key).profile

    t1 = timer.time() - t0

    print 'Loaded MS1 spectra in %s ' % t1
    return [ScanList, ms1ScanList, profiles]

def generate_peptide_list(sequence,max_length,mz_min,mz_max):
    """Generates list of all possible peptides and their profiles"""

    t0 = timer.time()

    seq_obj = mspy.sequence(sequence)
    peptide_objects = mspy.mod_proteo.digest(seq_obj,
                                    'Non-Specific',miscleavage=max_length)
    mass_indexed_peptide_list = {}

    for peptide in peptide_objects:
        compound = mspy.obj_compound.compound(peptide.formula())
        for z in range(charge_min,charge_max+1):
            pattern = compound.pattern(charge=z,real=False)
            if pattern[0][0] > mz_min and pattern[-1][0] < mz_max:
                highest_intensity_peak = max(pattern, key=lambda p: p[1])
                mass_indexed_peptide_list[highest_intensity_peak[0]] = (peptide,z,pattern)

    t1 = timer.time() - t0
    print 'Produced peptide isotopic distributions in %s ' % t1

    return mass_indexed_peptide_list

def match_peptide_patterns_to_ms1profile(mass_indexed_peptide_list,profile):
    for mass in mass_indexed_peptide_list.keys():
        if mspy.calculations.signal_intensity(profile, mass) > 5000:
            result = mspy.mod_pattern.checkpattern_fast(profile,mass_indexed_peptide_list[mass][2])
            if result.rmsd < RmsdThreshold:
                print mass_indexed_peptide_list[mass][0].format()
                print mass_indexed_peptide_list[mass][1]
                print mass


#mass_indexed_peptide_list =  generate_peptide_list(sequence,max_length,mz_min,mz_max)
reportDataFile = file(os.path.join(path,name+"_peppro.json"), 'r')
#pickle.dump(mass_indexed_peptide_list, reportDataFile, pickle.HIGHEST_PROTOCOL)
mass_indexed_peptide_list = pickle.load(reportDataFile)
reportDataFile.close()
ScanList, ms1ScanList, profiles = load_mzml_file(path, filename, start_time, stop_time)
print(ms1ScanList[157])
match_peptide_patterns_to_ms1profile(mass_indexed_peptide_list,profiles[ms1ScanList[157]])
