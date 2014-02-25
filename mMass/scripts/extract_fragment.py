import sys
sys.path.append("../")
import mspy
import time as timer
import os
import re
import json
import numpy
import pickle
import xml.etree.cElementTree as ET
from collections import defaultdict

#PARAMETER DECLARATION
path = "/Users/johannes/Downloads/"
#path = "F:/Documents/011114/"
#path = "/Volumes/Johannes-pc/Documents/011114/"
files = ["FIMC_digest_test_02231409_Re_TXIII_30"]


pepxml_file = "FIMC_digest_test_02231409_Re_TXIII_30.pep.xml"
peptides = []
with open("pepsin_fragments.txt", 'r') as f:
    for line in f:
        peptides.append(line.rstrip())
charge_min = 1
charge_max = 7
RmsdThreshold = 0.15

mz_min = 300
mz_max = 1500

start_time = 300
stop_time = 800

sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRL" \
    "QREPQVQWLEQQVAKRRTKR"

name = "FIMC_pepsin_022314_TXIII_30"

# FUNCTIONS


class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray) and obj.ndim == 2:
            return [y for y in [x for x in obj]]
        if isinstance(obj, numpy.ndarray) and obj.ndim == 1:
            return [x for x in obj]

        return json.JSONEncoder.default(self, obj)


def load_pepxml_file(peptides, pepxml_file):
    Ident_list = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for event, elem in ET.iterparse(pepxml_file):
        if elem.tag == '{http://regis-web.systemsbiology.net/pepXML}spectrum_query':
            condition = elem.attrib['spectrum'].split(".")[0]
            retention_time = float(elem.attrib['retention_time_sec']) * 60
            charge_state = int(elem.attrib['assumed_charge'])
            for peptide in elem[0]:
                if peptide.attrib['hit_rank'] == "1" and sequence.find(peptide.attrib['peptide']) != -1:
                    Ident_list[peptide.attrib['peptide']][charge_state][condition].append(retention_time)
    return Ident_list

def load_mzml_file(files, start_time, stop_time):
    """Load MS1 scans from MZML files"""

    t0 = timer.time()

    ScanList = {}
    ms1ScanList = {}
    profiles = {}

    for file_iter in files:

        parser = mspy.parseMZML(path + file_iter + '.mzML')

        parser.load()

        ScanList[file_iter] = parser.scanlist()
        ms1ScanList[file_iter] = []
        profiles[file_iter] = {}

        for key, scan_info in ScanList[file_iter].iteritems():
            if scan_info['msLevel'] == 1 \
                    and scan_info['retentionTime'] >= start_time \
                    and scan_info['retentionTime'] <= stop_time:
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

        buff_h = []
        buff_basepeak = []
        for z in range(charge_min, charge_max + 1):
            i = 0
            for data_i in MatchData['RetentionTime'][peptide][z].values():
                i += len(data_i)
            if i == 0:
                continue
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
        buff += "<h3 id=\"{0}\" >{0}</h3>".format(peptide)
        buff += "<table>"
        buff += "".join([a+b for a,b in zip(buff_h,buff_basepeak)])
        buff += "</table>"
    return buff


#Perform calculations

ScanList, ms1ScanList, profiles = load_mzml_file(files, start_time, stop_time)

MatchData = {}
IdentData = load_pepxml_file(peptides,path + pepxml_file)
peptides = IdentData.keys()
SequenceObjects, PatternObjects = generate_pattern_objects(peptides,
                                                           charge_min,
                                                           charge_max)
MatchData = match_peptides_in_scans(peptides, PatternObjects, ScanList,
                                    ms1ScanList, profiles, charge_min,
                                    charge_max, mz_min, mz_max, files)
MatchData["Ident"] = IdentData
MatchData["Pattern"] = PatternObjects
#with open("data.pickle", "wb") as f:
#    pickle.dump(MatchData,f,pickle.HIGHEST_PROTOCOL)
#with open("data.pickle", "r") as f:
#    MatchData = pickle.load(f)


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
