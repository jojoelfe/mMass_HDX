import sys
sys.path.append("../")
import mspy
import time as timer
import os
import re
import json
import numpy
#PARAMETER DECLARATION
path = "/Users/johannes/SHINDELAB/MassSpec/Johannes/011114_FIMC_test_bottom/"


files = ["FIMC_digest_test02_1_20_pepsin.mzML",
         "FIMC_digest_test03_1_50_pepsin.mzML",
         "FIMC_digest_test04_1_100_pepsin.mzML"]

peptides = []
with open("pepsin_fragments.txt", 'r') as f:
    for line in f:
        peptides.append(line.rstrip())
charge_min = 1
charge_max = 7
RmsdThreshold = 0.15

mz_min = 400
mz_max = 1500

sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRL" \
    "QREPQVQWLEQQVAKRRTKR"

name = "FIMC_pepsin_011114"

# FUNCTIONS


class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray) and obj.ndim == 2:
            return [y for y in [x for x in obj]]
        if isinstance(obj, numpy.ndarray) and obj.ndim == 1:
            return [x for x in obj]

        return json.JSONEncoder.default(self, obj)


def load_mzml_file(files):
    #Load MS1 scans from MZML files
    t0 = timer.time()

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

    print 'Extracted data in %s ' % t1
    return [ScanList, ms1ScanList, profiles]


def generate_peptide_positions(sequence, peptides):
    #Generates position of peptides for logo

    position_list = []

    for peptide in peptides:
        pos = sequence.find(peptide)
        position_list.append([pos, pos + len(peptide)])

    position_list = sorted(position_list, key=lambda x: x[1], reverse=True)
    position_list = sorted(position_list, key=lambda x: x[0])

    position_list_arranged = []

    end_values = [0]

    for position in position_list:
        for i, end in enumerate(end_values):
            if position[0] >= end:
                position_list_arranged.append([position[0],
                                               position[1],
                                               i])
                end_values[i] = position[1]
                break
        else:
            position_list_arranged.append([position[0],
                                           position[1],
                                           len(end_values)])
            end_values.append(position[1])

    return str(position_list_arranged)


def match_peptides_in_scans(peptides, PatternObjects, ScanList, ms1ScanList,
                            profiles, charge_min, charge_max, mz_min, mz_max,
                            files):
    #Matches expected isotopic distributions of peptides in MS1 spectra

    t0 = timer.time()
    RetentionTimeData = {}
    BasepeakData = {}
    RmsdData = {}
    ProfileData = {}

    # Iterate through peptide and initiate data structures
    for peptide in peptides:
        RetentionTimeData[peptide] = {}
        BasepeakData[peptide] = {}
        RmsdData[peptide] = {}
        ProfileData[peptide] = {}

        #Iterate through charge states and initiate data structures
        for z in range(charge_min, charge_max + 1):
            RetentionTimeData[peptide][z] = {}
            BasepeakData[peptide][z] = {}
            RmsdData[peptide][z] = {}
            ProfileData[peptide][z] = {}

            #Iterate through files and initiate data structures
            for file_iter in files:
                RetentionTimeData[peptide][z][file_iter] = []
                BasepeakData[peptide][z][file_iter] = []
                RmsdData[peptide][z][file_iter] = []
                ProfileData[peptide][z][file_iter] = None
                if (PatternObjects[peptide][z][0][0] < mz_min or
                        PatternObjects[peptide][z][-1][0] > mz_max):
                    continue
                #Iterate through scans
                for scan_number in ms1ScanList[file_iter]:
                    RetentionTimeData[peptide][z][file_iter].append(
                        ScanList[file_iter][scan_number]['retentionTime'])
                    checkPatternResult = mspy.checkpattern_fast(
                        signal=profiles[file_iter][scan_number],
                        pattern=PatternObjects[peptide][z])
                    if checkPatternResult is not None:
                        RmsdData[peptide][z][file_iter].append(
                            checkPatternResult.rmsd)
                        BasepeakData[peptide][z][file_iter].append(
                            checkPatternResult.basepeak)
                        if checkPatternResult.rmsd < RmsdThreshold:
                            if ProfileData[peptide][z][file_iter] is None:
                                ProfileData[peptide][z][file_iter] = mspy.crop(
                                    profiles[file_iter][scan_number],
                                    PatternObjects[peptide][z][0][0] - 0.1,
                                    PatternObjects[peptide][z][-1][0] + 0.1)
                            else:
                                ProfileData[peptide][z][file_iter] = mspy.combine(
                                    ProfileData[peptide][z][file_iter], mspy.crop(
                                        profiles[file_iter][scan_number],
                                        PatternObjects[peptide][z][0][0] - 0.1,
                                        PatternObjects[peptide][z][-1][0] + 0.1))
                    else:
                        RmsdData[peptide][z][file_iter].append(1)
                        BasepeakData[peptide][z][file_iter].append(0)
                if ProfileData[peptide][z][file_iter] is not None:
                    ProfileData[peptide][z][file_iter] = mspy.reduce(
                        ProfileData[peptide][z][file_iter])
    t1 = timer.time() - t0

    print 'Loaded files in %s ' % t1

    MatchData = {}
    MatchData['RetentionTime'] = RetentionTimeData
    MatchData['Basepeak'] = BasepeakData
    MatchData['Rmsd'] = RmsdData
    MatchData['Profile'] = ProfileData
    return MatchData


def generate_pattern_objects(peptides, charge_min, charge_max):
    #Calculates theoretical patterns

    SequenceObjects = {}
    PatternObjects = {}

    for peptide in peptides:
        SequenceObjects[peptide] = mspy.sequence(peptide)
        PatternObjects[peptide] = {}
        for z in range(charge_min, charge_max + 1):
            PatternObjects[peptide][z] = mspy.pattern(
                SequenceObjects[peptide].formula(), charge=z, real=False)

    return [SequenceObjects, PatternObjects]


def generate_peptide_html(peptides, MatchData, files, charge_min, charge_max):
    #Generate peptide html report
    buff = ""
    for peptide in peptides:

        buff_rmsd = []
        buff_h = []
        buff_basepeak = []
        for z in range(charge_min, charge_max + 1):
            buff_rmsd.append("<td><div id='rmsd_{0}_{1}' data-peptide='{0}' \
                           class='rmsd' data-charge=\
                          '{1}' data-type='Rmsd' style='width:200px;\
                          height:120px'></div></td>".format(peptide, z))
            buff_basepeak.append("<td><div id='basepeak_{0}_{1}' data-peptide='{0}' \
                           class='basepeak' data-charge=\
                          '{1}' data-type='Basepeak' style='width:200px;\
                          height:120px'></div></td>".format(peptide, z))
            buff_h.append("<th>{0}</th>".format(z))
        buff += "<h3>" + peptide + "</h3>"
        buff += "<table><tr>"
        buff += "".join(buff_h)
        buff += "</tr><tr>"
        buff += "".join(buff_rmsd)
        buff += "</tr><tr>"
        buff += "".join(buff_basepeak)
        buff += "</tr></table>"
    return buff



#Perform calculations

#ScanList, ms1ScanList, profiles = load_mzml_file(files)

#SequenceObjects, PatternObjects = generate_pattern_objects(peptides,
#                                                           charge_min,
#                                                           charge_max)
MatchData = {}
#MatchData = match_peptides_in_scans(peptides, PatternObjects, ScanList,
#                                    ms1ScanList, profiles, charge_min,
#                                    charge_max, mz_min, mz_max, files)

#Generating HTML report
buff = ""
rep_strings = {}
rep_strings['name'] = name
rep_strings['sequence'] = sequence

rep_strings['pepcoords'] = generate_peptide_positions(sequence, peptides)
rep_strings['pepplots'] = generate_peptide_html(peptides, MatchData, files,
                                                charge_min, charge_max)
confdir = "../configs/"
regexp = re.compile(re.compile("\{\{([A-Z]+)\}\}"))

with open(os.path.join(confdir, 'bottom_template.html'), 'r') as template:
        for template_line in template:
            buff += re.sub(
                regexp, lambda x: rep_strings[x.group(1).lower()],
                template_line)


reportDir = name + '_report/'

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
reportDataFile.write(json.dumps(MatchData, cls=NumpyAwareJSONEncoder))
reportDataFile.close()
