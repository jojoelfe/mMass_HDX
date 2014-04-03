library = "/Users/johannes/Documents/mMass_HDX/mMass/"

import sys
sys.path.append(library)
import mspy
import time as timer
import os
import pickle
from collections import defaultdict
import re
import json
import argparse



parser = argparse.ArgumentParser(description='Find proteolytic fragments in MS Data')
parser.add_argument('filename', metavar='filename', type=str, 
                   help='Filename')


args = parser.parse_args()



filename = args.filename.split(".")[0]

charge_min = 1
charge_max = 8

mz_min = 250
mz_max = 1800

RmsdThreshold = 0.15

start_time = 300
stop_time = 900

sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRLQREPQVQWLEQQVAKRRTKR"
sequence_name = "FIMC"

max_length = 50





def load_mzml_file(filename, start_time, stop_time):
    """Load MS1 scans from MZML files"""

    t0 = timer.time()

    ScanList = {}
    ms1ScanList = {}
    profiles = {}



    parser = mspy.parseMZML(filename + '.mzML')

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
    Match_result = defaultdict(dict)
    for mass in mass_indexed_peptide_list.keys():
        if mspy.calculations.signal_intensity(profile, mass) > 5000:
            result = mspy.mod_pattern.checkpattern_fast(profile,mass_indexed_peptide_list[mass][2])
            if result.rmsd < RmsdThreshold:
                Match_result[mass_indexed_peptide_list[mass][0].format()][mass_indexed_peptide_list[mass][1]] = result
    return Match_result

def generate_peptide_positions(sequence, peptides):
    #Generates position of peptides for logo

    position_list = []

    for peptide in peptides:
        pos = sequence.find(peptide)
        position_list.append([pos, pos + len(peptide),peptide])

    position_list = sorted(position_list, key=lambda x: x[1], reverse=True)
    position_list = sorted(position_list, key=lambda x: x[0])

    position_list_arranged = []

    end_values = [0]

    for position in position_list:
        for i, end in enumerate(end_values):
            if position[0] >= end:
                position_list_arranged.append([position[0],
                                               position[1],
                                               i,
                                               position[2]])
                end_values[i] = position[1]
                break
        else:
            position_list_arranged.append([position[0],
                                           position[1],
                                           len(end_values),
                                           position[2]])
            end_values.append(position[1])

    return str(position_list_arranged)

def generate_peptide_html(peptides):
    #Generate peptide html report
    buff = ""
    for peptide in peptides:

        buff_h = []
        buff_basepeak = []
        for z in peptides[peptide].keys():
            buff_basepeak.append("<td><div id='basepeak_{0}_{1}' data-peptide='{0}' \
                           class='basepeak' data-charge=\
                          '{1}' data-type='Basepeak' style='width:900px;\
                          height:220px'></div></td>\
                          <td><div id='spectrum_{0}_{1}' data-peptide='{0}' \
                           class='spectrum' data-charge=\
                          '{1}' data-type='Spectrum' style='width:300px;\
                          height:220px'></div></td>\
                          </tr>".format(peptide, z))
            buff_h.append("<tr><th>{0}</th>".format(z))
        buff += "<div id=\"{0}\" ><h3>{0}</h3>".format(peptide)
        buff += "<table>"
        buff += "".join([a+b for a,b in zip(buff_h,buff_basepeak)])
        buff += "</table></div>"
    return buff


try:
    os.stat(sequence_name+"_peppro.pickle")
    reportDataFile = file(os.path.join(sequence_name+"_peppro.pickle"), 'r')
    mass_indexed_peptide_list = pickle.load(reportDataFile)
    reportDataFile.close()
except:
    mass_indexed_peptide_list =  generate_peptide_list(sequence,max_length,mz_min,mz_max)
    reportDataFile = file(os.path.join(sequence_name+"_peppro.pickle"), 'wb')
    pickle.dump(mass_indexed_peptide_list, reportDataFile, pickle.HIGHEST_PROTOCOL)
    reportDataFile.close()

ScanList, ms1ScanList, profiles = load_mzml_file(filename, start_time, stop_time)
Retention_time = []
Match_results = []
peptides = defaultdict(dict)
t0 = timer.time()
for scan_number in ms1ScanList:
    Retention_time.append(ScanList[scan_number]['retentionTime'])
    Match_results.append(match_peptide_patterns_to_ms1profile(mass_indexed_peptide_list,profiles[scan_number]))
    for peptide in Match_results[-1].keys():
        for z in Match_results[-1][peptide].keys():
            peptides[peptide][z] = 1
t1 = timer.time() - t0
print 'Matched MS1 scans in %s ' % t1

Peptide_Intensities = {}
Peptide_IntIntensities = {}
for peptide in peptides.keys():
    Peptide_Intensities[peptide] = {}
    Peptide_IntIntensities[peptide] = 0.0
    for z in peptides[peptide].keys():
        Peptide_Intensities[peptide][z] = []


for i, time in enumerate(Retention_time):
    for peptide in peptides.keys():
        for z in peptides[peptide].keys():
            if peptide in Match_results[i].keys() and z in Match_results[i][peptide].keys():
                Peptide_Intensities[peptide][z].append(Match_results[i][peptide][z].basepeak)
                Peptide_IntIntensities[peptide] += Match_results[i][peptide][z].basepeak
            else:
                Peptide_Intensities[peptide][z].append(0.0)

Data_Json = {}
Data_Json['RetentionTime'] = Retention_time
Data_Json['Intensities'] = Peptide_Intensities
Data_Json['IntIntensities'] = Peptide_IntIntensities
Data_Json['peptides'] = peptides.keys()
Data_Json['pepcoords'] = generate_peptide_positions(sequence, peptides.keys())

regexp = re.compile(re.compile("\{\{([A-Z]+)\}\}"))

buff = ""
rep_strings = {}
rep_strings['name'] = filename
rep_strings['sequence'] = sequence

rep_strings['pepcoords'] = Data_Json['pepcoords']
rep_strings['pepplots'] = generate_peptide_html(peptides)

with open(os.path.join(library, 'scripts/single_protein_search_engine_template.html'), 'r') as template:
        for template_line in template:
            buff += re.sub(
                regexp, lambda x: rep_strings[x.group(1).lower()],
                template_line)

reportDir = filename + '_spse_results/'

try:
    os.stat(reportDir)
except:
    os.mkdir(reportDir)

reportPath = os.path.join(reportDir, 'report.html')
reportDataPath = os.path.join(reportDir, 'data.js')
reportFile = file(reportPath, 'w')

reportFile.write(buff.encode("utf-8"))
reportFile.close()

reportDataFile = file(reportDataPath, 'w')

reportDataFile.write("data=".encode("utf-8"))
reportDataFile.write(json.dumps(Data_Json))
reportDataFile.close()
