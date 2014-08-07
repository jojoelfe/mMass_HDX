import sys
sys.path.append('/Users/johannes/Documents/mMass_HDX/mMass/')
import numpy as np
import mspy
import json
import argparse

ms2_tolerance = 0.1

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 2:
            return [y for y in [x for x in obj]]
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
            return [x for x in obj]

        return json.JSONEncoder.default(self, obj)




def extract_intensities(files,peak_list):
    intensities_data = {}
    isotopic_distributions = {}
    for peak in peak_list:
        seq_obj = mspy.sequence(peak['sequence'])
        if 'fragment' not in peak or not peak['fragment']:
            isotopic_distributions[peak['name']] = mspy.mod_pattern.pattern(seq_obj.formula(),charge=peak['charge'])
            peak["isotopic_distribution"] = isotopic_distributions[peak['name']]
        else:
            if peak['fragmentserie'] == 'c':
                fragment = seq_obj[:peak['fragmentsite']]
                fragment.cTermFormula = mspy.blocks.fragments[peak['fragmentserie']].cTermFormula
                isotopic_distributions[peak['name']] = mspy.mod_pattern.pattern(fragment.formula(),charge=peak['charge'])
                peak["isotopic_distribution"] = isotopic_distributions[peak['name']]
            if peak['fragmentserie'] == 'z':
                fragment = seq_obj[-peak['fragmentsite']:]
                fragment.nTermFormula = mspy.blocks.fragments[peak['fragmentserie']].nTermFormula
                isotopic_distributions[peak['name']] = mspy.mod_pattern.pattern(fragment.formula(),charge=peak['charge'])
                peak["isotopic_distribution"] = isotopic_distributions[peak['name']]
    # Apply for each file
    for file in files:
        #Parse file
        print("Processing "+file+" ...")
        run = mspy.parseMZML(file)
        run.load()
        print("Finished parsing")
        scan_list = run.scanlist()
        intensities_data[file] = {}
        #Prepare data structure
        for peak in peak_list:
            intensities_data[file][peak['name']] = []
        for spec in scan_list.values():
            for peak in peak_list:
                if (peak['spectrum']['level'] == 1
                    and spec['msLevel'] == 1) or (peak['spectrum']['level'] == spec['msLevel']
                    and abs(peak['spectrum']['precursormz'] - spec['precursorMZ']) < ms2_tolerance ):
                    if spec['retentionTime'] >= peak["quantify_region"][0] and spec['retentionTime'] <= peak["quantify_region"][1]:
                        profile = run.scan(spec['scanNumber']).profile
                        subprofile = mspy.calculations.signal_crop(profile,isotopic_distributions[peak['name']][0][0]-0.1,isotopic_distributions[peak['name']][-1][0]+0.1)
                        intensities_data[file][peak['name']].append((spec['retentionTime'], mspy.calculations.signal_area(subprofile), spec['scanNumber']))

    return(intensities_data)

def create_usable_scan_list(intensity_data):
    scan_list = {}
    for file in intensity_data.keys():
        scan_list[file] = {}
        for peak in intensity_data[file].keys():
            scan_list[file][peak] = []
            maximum = max(intensity_data[file][peak],key=lambda x: x[1])
            for intensity in intensity_data[file][peak]:
                if intensity[1] > maximum[1]/2:
                    scan_list[file][peak].append(intensity[2])
    return scan_list

def average_profiles_from_usable_scans(peak_list, usable_scan_list):
    averaged_profiles = {}
    for condition in usable_scan_list.keys():
        averaged_profiles[condition] = {}
        for filename in usable_scan_list[condition].keys():
            print("Processing "+filename+" ...")
            run = mspy.parseMZML(filename)
            run.load()
            print("Finished parsing")
            averaged_profiles[condition][filename] = {}
            for peak in peak_list:
                for scanID in usable_scan_list[condition][filename][peak['name']]:
                    profile = run.scan(scanID).profile
                    subprofile = mspy.calculations.signal_crop(profile,peak['isotopic_distribution'][0][0]-0.1,peak['isotopic_distribution'][-1][0]+1)
                    if peak['name'] in averaged_profiles[condition][filename]:
                        averaged_profiles[condition][filename][peak['name']] = mspy.calculations.signal_combine(averaged_profiles[condition][filename][peak['name']],subprofile)
                    else:
                        averaged_profiles[condition][filename][peak['name']] = subprofile
                averaged_profiles[condition][filename][peak['name']] = mspy.calculations.signal_reduce(averaged_profiles[condition][filename][peak['name']],0.002)
    return averaged_profiles

def generate_seq_objects(peak_list):
    seq_objects = {}
    for peak in peak_list:
        seq_obj = mspy.sequence(peak['sequence'])
        if 'fragment' not in peak or not peak['fragment']:
            seq_objects[peak['name']] = seq_obj
        else:
            if peak['fragmentserie'] == 'c':
                fragment = seq_obj[:peak['fragmentsite']]
                fragment.cTermFormula = mspy.blocks.fragments[peak['fragmentserie']].cTermFormula
                seq_objects[peak['name']] = fragment
            if peak['fragmentserie'] == 'z':
                fragment = seq_obj[-peak['fragmentsite']:]
                fragment.nTermFormula = mspy.blocks.fragments[peak['fragmentserie']].nTermFormula
                seq_objects[peak['name']] = fragment
    return seq_objects

def quantify_uptake(files,peak_list):
    uptake_data = {}
    fitters = {}
    for peak in peak_list:
        seq_obj = mspy.sequence(peak['sequence'])
        if 'fragment' not in peak or not peak['fragment']:
            fitters[peak['name']] = exchange_quantifier_average_mass(seq_obj=seq_obj,charge=peak['charge'],scales=2)
        else:
            if peak['fragmentserie'] == 'c':
                fragment = seq_obj[:peak['fragmentsite']]
                fragment.cTermFormula = mspy.blocks.fragments[peak['fragmentserie']].cTermFormula
                fitters[peak['name']] = exchange_quantifier_average_mass(seq_obj=fragment,charge=peak['charge'],scales=2)
            if peak['fragmentserie'] == 'z':
                fragment = seq_obj[-peak['fragmentsite']:]
                fragment.nTermFormula = mspy.blocks.fragments[peak['fragmentserie']].nTermFormula
                fitters[peak['name']] = exchange_quantifier_average_mass(seq_obj=fragment,charge=peak['charge'],scales=2)

    # Apply for each file
    for file in files:
        #Parse file
        run = mspy.parseMZML(file)
        run.load()
        scan_list = run.scanlist()
        uptake_data[file] = {}
        #Prepare data structure
        for peak in peak_list:
            uptake_data[file][peak['name']] = []
        for spec in scan_list.values():
            for peak in peak_list:
                if (peak['spectrum']['level'] == 1
                    and spec['msLevel'] == 1) or (peak['spectrum']['level'] == spec['msLevel']
                    and abs(peak['spectrum']['precursormz'] - spec['precursorMZ']) < ms2_tolerance ):
                    if spec['retentionTime'] >= peak["quantify_region"][0] and spec['retentionTime'] <= peak["quantify_region"][1]:
                        profile = run.scan(spec['scanNumber']).profile
                        res = fitters[peak['name']].quantify_from_signal(profile)
                        uptake_data[file][peak['name']].append((spec['retentionTime'], res))

    return(uptake_data)

def quantify_uptake_in_file(file,peaklist):
    pass

def quantify_uptake_in_scan(spec,peaklist):
    pass

def pick_peaks_from_profiles(profiles, peak_list, seq_objects):
    result = {}
    for condition in profiles.keys():
        result[condition] = {}
        for filename in profiles[condition].keys():
            result[condition][filename] = {}
            for peak in peak_list:
                picker_label_height = peak_picker_label_height(seq_objects[peak["name"]],peak["charge"])
                picker_integral = peak_picker_integral(seq_objects[peak["name"]],peak["charge"])

                result[condition][filename][peak["name"]]["label_height"] = picker_label_height.pick_from_signal(profiles[condition][filename][peak["name"]])
                result[condition][filename][peak["name"]]["integral"] = picker_integral.pick_from_signal(profiles[condition][filename][peak["name"]])
    return result

class peak_picker( object ):
    """Abstract class that defines an interface for peak fitters. Initialized using a sequence object to generate isotopic distribution.
      Returns array of intensities for I, I+1, etc peaks.
    """

    def __init__(self,seq_obj=None,charge=None,scales=1):
        self.scales = scales

        if seq_obj is not None and charge is not None:
            self.seq_obj = seq_obj
            self.charge = charge
        else:
            raise StandardError("You have to provide either a sequence object or a no-exchange model")
        self._generate_pickmzs()

    def _generate_pickmzs(self):
        pattern = mspy.mod_pattern.pattern(self.seq_obj.formula(),charge=self.charge,real=False)
        self.mzs = []
        for i in range(self.scales+len(pattern)):
            av_list = []
            for idx in range(self.scales+1):
                if (i-idx) >= 0 and i-idx < len(pattern):
                    av_list.append(pattern[i-idx][0]+idx*1.003/self.charge)
            self.mzs.append(sum(av_list) / len(av_list))


class peak_picker_label_height( peak_picker ):

    def pick_from_signal(self, signal):
        result = []
        for idx, mz in enumerate(self.mzs):
            peak = mspy.labelpeak(signal,mz)
            if peak is not None:
                result.append([peak.mz,peak.ai])
            else:
                result.append([mz,0.0])

        return(result)

class peak_picker_integral( peak_picker ):
    def pick_from_signal(self, signal):
        result = []
        for idx, mz in enumerate(self.mzs):
            subprofile = mspy.calculations.signal_crop(signal, mz-0.1, mz+0.1)
            result.append([0,mspy.calculations.signal_area(subprofile)])
        return(result)





class exchange_quantifier( object ):
    """Abstract class that defines interface for hdx quantifiers and implements common methods.
    In general it quantifies the intensity of I, I+1 etc peaks assuming that dmass is unity"""

    def __init__(self,exchange_model=None,seq_obj=None,charge=None,scales=1):

        self.scales = scales
        if exchange_model is not None:
            self.no_exchange_model = exchange_model
        elif seq_obj is not None and charge is not None:
            self.no_exchange_model = mspy.mod_pattern.pattern(seq_obj.formula(),charge=charge,real=False)
        else:
            raise StandardError("You have to provide either a sequence object or a no-exchange model")
        self._generate_models()

    def quantify_from_peaklist(self):
        raise NotImplementedError( "Quantification from peaklist should be implemented by children" )

    def quantify_from_signal(self, signal):

        peaklist = np.zeros(len(self.mzs))

        for idx, mz in enumerate(self.mzs):
            peak = mspy.labelpeak(signal,mz)
            if peak is not None:
                peaklist[idx] = peak.ai

        return(self.quantify_from_peaklist(peaklist))


    def _generate_models(self):
        self.mzs = map(lambda x: x[0], self.no_exchange_model)
        for idx in range(self.scales):
            self.mzs.append(self.mzs[-1]+1.003)
        self.models = np.zeros((self.scales+1,len(self.mzs)))
        for idx in range(self.scales+1):
            for jdx, intensity in enumerate(map(lambda x: x[1], self.no_exchange_model)):
                self.models[idx][idx+jdx] = intensity

class exchange_quantifier_fit(exchange_quantifier):
    pass

class exchange_quantifier_peak_ratio(exchange_quantifier):
    pass

class exchange_quantifier_average_mass(exchange_quantifier):

    def quantify_from_peaklist(self,peaklist):
        rel_masses = np.arange(peaklist.shape[0])
        weighted_masses = np.dot(peaklist,rel_masses)
        av_mass = np.sum(weighted_masses) / np.sum(peaklist)
        return(av_mass)

# Main program starts

#Parsing arguments
parser = argparse.ArgumentParser(description="Analyze HDX exchange in proteolytic peptides")
parser.add_argument('conditions', metavar ='conditions', type=str, help="File with list of conditions in json")
parser.add_argument('peaks', metavar ='peaks', type=str, help="File with list of peaks to analyze")

args = parser.parse_args()


#Reading input files
with open(args.conditions) as conditions_file:
    conditions = json.load(conditions_file)
with open(args.peaks) as peaks_file:
    peaklist = json.load(peaks_file)


int_data = {}
usable_scan_list = {}

seq_objects = generate_seq_objects(peaklist)


for cond in conditions:
    int_data[cond["name"]] = extract_intensities(cond["files"],peaklist)
    usable_scan_list[cond["name"]] = create_usable_scan_list(int_data[cond["name"]])

profiles = average_profiles_from_usable_scans(peaklist,usable_scan_list)

picking = pick_peaks_from_profiles(profiles, peaklist, seq_objects)

#Output results
result = {}
result["input"] = {"conditions":conditions,"peaklist":peaklist}
result["int_data"] = int_data
result["scans_used_for_averaging"] = usable_scan_list
result["profiles"] = profiles
result["picking"] = picking

reportDataPath = 'data.js'
reportDataFile = file(reportDataPath, 'w')

reportDataFile.write("data=".encode("utf-8"))
jsondump = json.dumps(result,cls=NumpyAwareJSONEncoder)
reportDataFile.write(jsondump)
reportDataFile.close()

