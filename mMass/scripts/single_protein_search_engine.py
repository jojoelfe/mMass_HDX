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
import numpy as np

class XCorr:
    def __init__(self,bin_width=1.0011413,bin_offset=0.68):
        self.bin_width=bin_width
        self.bin_offset=bin_offset

        self.NUM_REGIONS = 10
        self.MAX_XCORR_OFFSET = 75

    def __integerize(self,value):
        return int((value / self.bin_width + 1) - self.bin_offset)

    def __normalize_each_region(self,data, max_int_overall, max_int_per_regions, region_selector):

        region_idx = 0
        max_int = max_int_per_regions[region_idx]

        for bin_idx, intensity in enumerate(data):

            if bin_idx >= region_selector*(region_idx+1) and region_idx < self.NUM_REGIONS-1:
                region_idx += 1
                max_int = max_int_per_regions[region_idx]

            if max_int != 0.0 and intensity > 0.05 * max_int_overall:
                data[bin_idx] = (intensity / max_int) * 50
            else:
                data[bin_idx] = 0.0


    def process_ms2(self,ms2_profile,max_mass,max_bin,precursor_mz):

        max_int_overall = 0
        max_int_per_regions = np.zeros(self.NUM_REGIONS)

        region_selector = int(max_bin / self.NUM_REGIONS)

        region = 0

        observed = np.zeros(max_bin)

        for peak in ms2_profile:
            if peak[0] > max_mass:
                continue

            if abs(peak[0]-precursor_mz) < 15:
                continue

            mz = self.__integerize(peak[0]);
            if mz >= max_bin:
              continue
            region = int(mz/region_selector)

            if(region >= self.NUM_REGIONS):
                region = self.NUM_REGIONS - 1

            intensity = np.sqrt(peak[1])

            if intensity > max_int_overall:
                max_int_overall = intensity

            if observed[mz] < intensity:
                observed[mz] = intensity

            if max_int_per_regions[region] < intensity:
                max_int_per_regions[region] = intensity

        self.__normalize_each_region(observed,max_int_overall,max_int_per_regions,region_selector)

        new_observed = observed.copy()

        for ix in range(1,self.MAX_XCORR_OFFSET+1):
            new_observed -= np.pad(observed,(ix,0),'constant',constant_values=(0,0))[0:-ix] / (self.MAX_XCORR_OFFSET * 2)

        for ix in range(1,self.MAX_XCORR_OFFSET+1):
            new_observed -= np.pad(observed,(0,ix),'constant',constant_values=(0,0))[ix:] / (self.MAX_XCORR_OFFSET * 2)

        return new_observed

    def create_theoretical_spectrum(self,fragment_masses,fragment_loss_masses,max_bin):

        theoretical = np.zeros(max_bin)

        for mass in fragment_masses:
            mz = self.__integerize(mass)
            if mz >= max_bin:
                continue
            if theoretical[mz] < 50:
                theoretical[mz] = 50


        for mass in fragment_loss_masses:
            mz = self.__integerize(mass)
            if mz >= max_bin:
                continue
            if theoretical[mz] < 10:
                theoretical[mz] = 10

        return theoretical

    def get_xcorr(self,ms2_profile,precursor_mz,precursor_charge,fragment_masses,fragment_loss_masses):

        max_mass = min(ms2_profile[-1][0],precursor_mz*precursor_charge)
        max_bin = self.__integerize(max_mass)
        observed = self.process_ms2(ms2_profile,max_mass,max_bin,precursor_mz)
        theoretical = self.create_theoretical_spectrum(fragment_masses, fragment_loss_masses,max_bin)

        return np.dot(observed,theoretical)/10000


class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 2:
            return [y for y in [x for x in obj]]
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
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

        t0 = timer.time()
        # 1.key: Peptide, 2,key charge state, 3 key: activation method
        self.ms2_spectra = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        for key, scan_info in self.scan_list.iteritems():
            if scan_info['msLevel'] == 2 \
              and scan_info['retentionTime'] >= self.start_time \
              and scan_info['retentionTime'] <= self.stop_time:
                  try:
                      parent_index = self.ms1_indices.index(scan_info['parentScanNumber'])
                  except ValueError:
                      continue
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
        t1 = timer.time() - t0
        print 'Extracted MS2 scans in %s ' % t1


    def annotate_ms2scans(self):
        # Checks each peak against list of potential fragments
        t0 = timer.time()
        self.ms2_annotation = {}
        xcorr = XCorr()
        for peptide in self.ms2_spectra.keys():
            self.ms2_annotation[peptide] = {}
            for z in self.ms2_spectra[peptide].keys():
                self.ms2_annotation[peptide][z]={}
                for activ in self.ms2_spectra[peptide][z].keys():
                    fragments = mspy.mod_proteo.fragment(mspy.sequence(peptide),"by")
                    fragments_loss = mspy.mod_proteo.fragmentlosses(fragments,["H20","NH3"])
                    fragment_mzs = []
                    fragment_loss_mzs = []
                    for z2 in range(1,int(z)+1):
                        for fragment in fragments:
                            fragment_mzs.append([
                                fragment.mz(charge=z2)[1],
                                fragment.format('f ') +' '+ str(z2) + '+'
                            ])
                        for fragment in fragments_loss:
                            fragment_loss_mzs.append([
                                fragment.mz(charge=z2)[1],
                                fragment.format('f ') +' '+ str(z2) + '+'
                            ])
                    self.ms2_annotation[peptide][z][activ] = {}
                    self.ms2_annotation[peptide][z][activ] = fragment_mzs
                    self.ms2_annotation[peptide][z][activ] += fragment_loss_mzs
                    self.ms2_annotation[peptide][z][activ].append((
                    mspy.mod_proteo.fragment(mspy.sequence(peptide),"M")[0].mz(charge=int(z))[1],
                    mspy.mod_proteo.fragment(mspy.sequence(peptide),"M")[0].format('f '+' '+str(z)+'+')
                    ))
                    for scan_obj in self.ms2_spectra[peptide][z][activ]:
                        xcorr_score = xcorr.get_xcorr(scan_obj['peaklist'],scan_obj['scan_info']['precursorMZ'],int(z),map(lambda x: x[0],fragment_mzs),map(lambda x: x[0],fragment_loss_mzs))
                        scan_obj['scan_info']['Xcorr'] = xcorr_score


        t1 = timer.time() - t0
        print 'Calculated annotations for MS2 scans in %s ' % t1

    def calculate_num_overlap_isodist(self,isodist1,isodist2):
        differences = [a[0] - b[0] for a in isodist1 for b in isodist2]
        num = 0
        for difference in differences:
            if abs(difference) < 0.03:
                num += 1
        return num

    def find_potential_clashes(self):

        t0 = timer.time()

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
        t1 = timer.time() - t0
        print 'Calculated peptide clashes in %s ' % t1


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
        Data_Json['ms2annotation'] = self.ms2_annotation


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
spse_obj.annotate_ms2scans()
spse_obj.write_json_and_html()
