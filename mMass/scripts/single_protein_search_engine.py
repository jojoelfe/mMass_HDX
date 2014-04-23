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
import numpy
from bson import BSON

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray) and obj.ndim == 2:
            return [y for y in [x for x in obj]]
        if isinstance(obj, numpy.ndarray) and obj.ndim == 1:
            return [x for x in obj]

        return json.JSONEncoder.default(self, obj)

class SPSE:
    def __init__(self):
        self.charge_min = 1
        self.charge_max = 8

        self.mz_min = 250
        self.mz_max = 1800

        self.rmsd_threshold = 0.15

        self.start_time = 300
        self.stop_time = 900

        self.sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRLQREPQVQWLEQQVAKRRTKR"
        self.sequence_name = "FIMC"

        self.max_length = 50
        # Parse first argument as filename
        parser = argparse.ArgumentParser(description='Find proteolytic fragments in MS Data')
        parser.add_argument('filename', metavar='filename', type=str,
                           help='Filename')

        args = parser.parse_args()

        self.filename = args.filename.split(".")[0]


    def load_mzml_file(self):
        """Load MS1 scans from MZML files"""

        t0 = timer.time()

        self.parser = mspy.parseMZML(self.filename + '.mzML')

        self.parser.load()

        self.scan_list  = self.parser.scanlist()
        self.ms1_indices = []
        self.profiles = {}

        for key, scan_info in self.scan_list.iteritems():
            if scan_info['msLevel'] == 1 \
              and scan_info['retentionTime'] >= self.start_time \
              and scan_info['retentionTime'] <= self.stop_time:
                self.ms1_indices.append(key)
                self.profiles[key] = self.parser.scan(key).profile

        t1 = timer.time() - t0

        print 'Loaded MS1 spectra in %s ' % t1

    def load_peptide_list(self):
        try:
            os.stat(self.sequence_name+"_peppro.pickle")
            reportDataFile = file(os.path.join(self.sequence_name+"_peppro.pickle"), 'r')
            (self.peptide_indexed_iso_maxpeak, self.peptide_indexed_iso_dist) = pickle.load(reportDataFile)
            reportDataFile.close()
        except:
            self.generate_peptide_list()
            reportDataFile = file(os.path.join(self.sequence_name+"_peppro.pickle"), 'wb')
            pickle.dump((self.peptide_indexed_iso_maxpeak,self.peptide_indexed_iso_dist), reportDataFile, pickle.HIGHEST_PROTOCOL)
            reportDataFile.close()


    def generate_peptide_list(self):
        """Generates list of all possible peptides and their profiles"""

        t0 = timer.time()

        seq_obj = mspy.sequence(self.sequence)
        peptide_objects = mspy.mod_proteo.digest(seq_obj,
                                        'Non-Specific',miscleavage=self.max_length)
        self.peptide_indexed_iso_maxpeak = {}
        self.peptide_indexed_iso_dist = {}
        for peptide in peptide_objects:
            compound = mspy.obj_compound.compound(peptide.formula())
            self.peptide_indexed_iso_dist[peptide.format()] = {}
            self.peptide_indexed_iso_maxpeak[peptide.format()] = {}
            for z in range(self.charge_min,self.charge_max+1):
                pattern = compound.pattern(charge=z,real=False)
                if pattern[0][0] > self.mz_min and pattern[-1][0] < self.mz_max:
                    highest_intensity_peak = max(pattern, key=lambda p: p[1])
                    self.peptide_indexed_iso_maxpeak[peptide.format()][str(z)] = highest_intensity_peak[0]
                    self.peptide_indexed_iso_dist[peptide.format()][str(z)] = pattern

        t1 = timer.time() - t0
        print 'Produced peptide isotopic distributions in %s ' % t1

    def match_peptide_patterns_to_ms1profile(self,profile):
        match_result = defaultdict(dict)
        for peptide in self.peptide_indexed_iso_maxpeak.keys():
            for z in self.peptide_indexed_iso_maxpeak[peptide].keys():
                if mspy.calculations.signal_intensity(profile, self.peptide_indexed_iso_maxpeak[peptide][z]) > 5000:
                    result = mspy.mod_pattern.checkpattern_fast(profile,self.peptide_indexed_iso_dist[peptide][z])
                    if result.rmsd < self.rmsd_threshold:
                        match_result[peptide][z] = result
        return match_result

    def generate_peptide_positions(self):
        #Generates position of peptides for logo

        position_list = []

        for peptide in self.matched_peptides.keys():
            pos = self.sequence.find(peptide)
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

    def generate_peptide_html(self):
        #Generate peptide html report
        buff = ""
        for peptide in self.matched_peptides.keys():

            buff_h = []
            buff_basepeak = []
            for z in self.matched_peptides[peptide].keys():
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

    def perform_peptide_match(self):
        self.matched_peptides = {}
        self.peptide_ms1_profiles = {}
        self.retention_times = []
        self.match_data = []
        t0 = timer.time()
        for scan_number in self.ms1_indices:
            self.retention_times.append(self.scan_list[scan_number]['retentionTime'])
            self.match_data.append(self.match_peptide_patterns_to_ms1profile(self.profiles[scan_number]))
            for peptide in self.match_data[-1].keys():
                if peptide in self.matched_peptides:
                    for z in self.match_data[-1][peptide].keys():
                        self.matched_peptides[peptide][z] = 1
                        if z in self.peptide_ms1_profiles[peptide]:
                            self.peptide_ms1_profiles[peptide][z] = mspy.combine(self.peptide_ms1_profiles[peptide][z] ,mspy.crop(
                                self.profiles[scan_number],
                                self.peptide_indexed_iso_dist[peptide][z][0][0] - 0.1,
                                self.peptide_indexed_iso_dist[peptide][z][-1][0] + 0.1))
                        else:
                            self.peptide_ms1_profiles[peptide][z] = mspy.crop(
                                self.profiles[scan_number],
                                self.peptide_indexed_iso_dist[peptide][z][0][0] - 0.1,
                                self.peptide_indexed_iso_dist[peptide][z][-1][0] + 0.1)
                else:
                    self.matched_peptides[peptide] = {}
                    self.peptide_ms1_profiles[peptide] = {}
                    for z in self.match_data[-1][peptide].keys():
                        self.matched_peptides[peptide][z] = 1
                        self.peptide_ms1_profiles[peptide][z] = mspy.crop(
                            self.profiles[scan_number],
                            self.peptide_indexed_iso_dist[peptide][z][0][0] - 0.1,
                            self.peptide_indexed_iso_dist[peptide][z][-1][0] + 0.1)
        t1 = timer.time() - t0
        print 'Matched MS1 scans in %s ' % t1

    def integrate_match_data(self):
        # Creates peptide inensities dictionary with timeline of intensities
        # and sum of intensities for all charge states for each peptide
        self.peptide_intensities = {}
        self.peptide_intensities_int = {}
        # Preparing data structures
        for peptide in self.matched_peptides.keys():
            self.peptide_intensities[peptide] = {}
            self.peptide_intensities_int[peptide] = 0.0
            for z in self.matched_peptides[peptide].keys():
                self.peptide_intensities[peptide][z] = []
                self.peptide_ms1_profiles[peptide][z] = mspy.reduce(
                    self.peptide_ms1_profiles[peptide][z])
            # Summing over each data point
        for i, time in enumerate(self.retention_times):
            for peptide in self.matched_peptides.keys():
                for z in self.matched_peptides[peptide].keys():
                    if peptide in self.match_data[i].keys() and z in self.match_data[i][peptide].keys():
                        self.peptide_intensities[peptide][z].append(self.match_data[i][peptide][z].basepeak)
                        self.peptide_intensities_int[peptide] += self.match_data[i][peptide][z].basepeak
                    else:
                        self.peptide_intensities[peptide][z].append(0.0)

    def extract_relevant_ms2scans(self):

        # 1.key: Peptide, 2,key charge state, 3 key: activation method
        self.ms2_spectra = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        for key, scan_info in self.scan_list.iteritems():
            if scan_info['msLevel'] == 2 \
              and scan_info['retentionTime'] >= self.start_time \
              and scan_info['retentionTime'] <= self.stop_time:
                  parent_index = self.ms1_indices.index(scan_info['parentScanNumber'])
                  for peptide in self.match_data[parent_index].keys():
                      for z in self.match_data[parent_index][peptide].keys():
                          if abs(scan_info['precursorMZ'] - self.peptide_indexed_iso_maxpeak[peptide][z]) < 2:
                              scan = self.parser.scan(key)
                              scan.swap()
                              self.ms2_spectra[peptide][z][scan_info['activation']].append({
                              "retention_time" : scan_info['retentionTime'],
                              "peaklist" : scan.profile,
                              "scan_info" : scan_info
                              })

        self.ms2_spectra = dict(self.ms2_spectra)
        for peptide_key in self.ms2_spectra.keys():
            for charge_key in self.ms2_spectra[peptide_key].keys():
                self.ms2_spectra[peptide_key][charge_key] = dict(self.ms2_spectra[peptide_key][charge_key])
            self.ms2_spectra[peptide_key] = dict(self.ms2_spectra[peptide_key])



    def calculate_num_overlap_isodist(self,isodist1,isodist2):
        differences = [a[0] - b[0] for a in isodist1 for b in isodist2]
        num = 0
        for difference in differences:
            if abs(difference) < 0.03:
                num += 1
        return num

    def find_potential_clashes(self):
        self.overlaps = {}
        for peptide in self.matched_peptides.keys():
            self.overlaps[peptide] = {}
            for z in self.matched_peptides[peptide].keys():
                self.overlaps[peptide][z] = []
                for peptidei in self.matched_peptides.keys():
                    for zi in self.matched_peptides[peptidei].keys():
                        if abs(self.peptide_indexed_iso_dist[peptidei][zi][0][0] - self.peptide_indexed_iso_dist[peptide][z][0][0]) < 3:
                            numdiff = self.calculate_num_overlap_isodist(self.peptide_indexed_iso_dist[peptidei][zi],self.peptide_indexed_iso_dist[peptide][z])
                            if numdiff > 1 and peptide != peptidei:
                                self.overlaps[peptide][z].append((peptidei,zi))

    def write_json_and_html(self):
        reportDir = self.filename + '_spse_results/'

        try:
            os.stat(reportDir)
        except:
            os.mkdir(reportDir)

        Data_Json = {}
        Data_Json['RetentionTime'] = self.retention_times
        Data_Json['Intensities'] = self.peptide_intensities
        Data_Json['IntIntensities'] = self.peptide_intensities_int
        Data_Json['peptides'] = self.matched_peptides
        Data_Json['sequence'] = self.sequence
        Data_Json['ms1Profiles'] = self.peptide_ms1_profiles
        Data_Json['IsotopDistr'] = self.peptide_indexed_iso_dist
        Data_Json['overlaps'] = self.overlaps
        Data_Json['ms2scans'] = self.ms2_spectra


        reportDataPath = os.path.join(reportDir, 'data.js')
        reportDataFile = file(reportDataPath, 'w')

        reportDataFile.write("data=".encode("utf-8"))
        jsondump = json.dumps(Data_Json,cls=NumpyAwareJSONEncoder)
        reportDataFile.write(jsondump)
        reportDataFile.close()

        #Write as BSON (not as effective as I thought)
        #reportDataPath = os.path.join(reportDir, 'data.bson')
        #reportDataFile = file(reportDataPath, 'wb')


        #reportDataFile.write(BSON.encode(json.loads(jsondump)))
        #reportDataFile.close()

        regexp = re.compile(re.compile("\{\{([A-Z]+)\}\}"))

        buff = ""
        rep_strings = {}
        rep_strings['name'] = self.filename
        rep_strings['sequence'] = self.sequence

        rep_strings['pepcoords'] = ""
        rep_strings['pepplots'] = ""

        with open(os.path.join(library, 'scripts/single_protein_search_engine_template.html'), 'r') as template:
                for template_line in template:
                    buff += re.sub(
                        regexp, lambda x: rep_strings[x.group(1).lower()],
                        template_line)



        reportPath = os.path.join(reportDir, 'report.html')
        reportDataPath = os.path.join(reportDir, 'data.js')
        reportFile = file(reportPath, 'w')

        reportFile.write(buff.encode("utf-8"))
        reportFile.close()


spse_obj = SPSE()

spse_obj.load_peptide_list()
spse_obj.load_mzml_file()
spse_obj.perform_peptide_match()
spse_obj.integrate_match_data()
spse_obj.find_potential_clashes()
spse_obj.extract_relevant_ms2scans()
spse_obj.write_json_and_html()
