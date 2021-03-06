import sys
sys.path.append('/Users/johannes/Documents/mMass_HDX/mMass/')
import numpy as np
import mspy
import json
import argparse
import pickle
import scipy.optimize
from numpy.linalg import solve as solveLinEq
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
            if peak['fragmentserie'] == 'b':
                fragment = seq_obj[:peak['fragmentsite']]
                fragment.cTermFormula = mspy.blocks.fragments[peak['fragmentserie']].cTermFormula
                isotopic_distributions[peak['name']] = mspy.mod_pattern.pattern(fragment.formula(),charge=peak['charge'])
                peak["isotopic_distribution"] = isotopic_distributions[peak['name']]
            if peak['fragmentserie'] == 'z':
                fragment = seq_obj[-peak['fragmentsite']:]
                fragment.nTermFormula = mspy.blocks.fragments[peak['fragmentserie']].nTermFormula
                isotopic_distributions[peak['name']] = mspy.mod_pattern.pattern(fragment.formula(),charge=peak['charge'])
                peak["isotopic_distribution"] = isotopic_distributions[peak['name']]
            if peak['fragmentserie'] == 'y':
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
                                and abs(peak['spectrum']['precursormz'] - spec['precursorMZ']) < ms2_tolerance and ((not "activation" in peak['spectrum']) or peak['spectrum']['activation'] == spec['activation']) ):
                            if spec['retentionTime'] >= peak["quantify_region"][0] and spec['retentionTime'] <= peak["quantify_region"][1]:
                                profile = run.scan(spec['scanNumber']).profile
                                print("Cropping: scan " + str(spec['scanNumber']) + "| "+str(isotopic_distributions[peak['name']][0][0]-0.1)+"-"+str(isotopic_distributions[peak['name']][-1][0]+0.1))
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
                averaged_profiles[condition][filename][peak['name']] = mspy.calculations.signal_reduce(averaged_profiles[condition][filename][peak['name']],0.0005)
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
            if peak['fragmentserie'] == 'b':
                fragment = seq_obj[:peak['fragmentsite']]
                fragment.cTermFormula = mspy.blocks.fragments[peak['fragmentserie']].cTermFormula
                seq_objects[peak['name']] = fragment
            if peak['fragmentserie'] == 'y':
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
                result[condition][filename][peak["name"]] = {}
                result[condition][filename][peak["name"]]["label_height"] = picker_label_height.pick_from_signal(profiles[condition][filename][peak["name"]])
                result[condition][filename][peak["name"]]["integral"] = picker_integral.pick_from_signal(profiles[condition][filename][peak["name"]])
    return result

def quantify_exchange_from_picking(picking, peak_list, seq_objects, conditions):
    peaks_control = {}

    # Extract peak pickings for controls (Water)
    for condition in conditions:
        if condition["control"]:
            for peak in peak_list:
                peaks_control[peak['name']] = {} 
                for filename in picking[condition['name']].keys():
                    for pick_alg in picking[condition['name']][filename][peak['name']].keys():
                        if pick_alg in peaks_control[peak['name']]:
                            peaks_control[peak['name']][pick_alg].append(picking[condition['name']][filename][peak['name']][pick_alg])
                        else:
                            peaks_control[peak['name']][pick_alg] = [picking[condition['name']][filename][peak['name']][pick_alg]]

    #Average control peak values
    avg_control_peaks = {}
    for peak in peaks_control.keys():
        avg_control_peaks[peak] = {}
        for pick_alg in peaks_control[peak].keys():
            avg_control_peaks[peak][pick_alg] = np.zeros(np.array(peaks_control[peak][pick_alg][0]).shape)
            for peak_sample in peaks_control[peak][pick_alg]:
                avg_control_peaks[peak][pick_alg] += np.array(peak_sample)


    
    result = {}

    for condition in picking.keys():
        result[condition] = {}
        for filename in picking[condition].keys():
            result[condition][filename] = {}
            for peak in peak_list:
                result[condition][filename][peak["name"]] = {}
                for pick_alg in picking[condition][filename][peak["name"]].keys():
                    result[condition][filename][peak["name"]][pick_alg] = {}
                    quantifier_av_mass = exchange_quantifier_average_mass(seq_obj= seq_objects[peak["name"]], charge = peak["charge"])
                    quantifier_av_mass_exp_cont = exchange_quantifier_average_mass(charge = peak["charge"],exchange_model = avg_control_peaks[peak['name']][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["av_mass"] = quantifier_av_mass.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["av_mass_exp_cont"] = quantifier_av_mass_exp_cont.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
                    quantifier_peak_ratio = exchange_quantifier_peak_ratio(seq_obj= seq_objects[peak["name"]], charge = peak["charge"])
                    quantifier_peak_ratio_exp_cont = exchange_quantifier_peak_ratio(charge = peak["charge"],exchange_model = avg_control_peaks[peak['name']][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["peak_ratio"] = quantifier_peak_ratio.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["peak_ratio_exp_cont"] = quantifier_peak_ratio_exp_cont.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
                    quantifier_fit = exchange_quantifier_fit(seq_obj= seq_objects[peak["name"]], charge = peak["charge"])
                    quantifier_fit_exp_cont = exchange_quantifier_fit(charge = peak["charge"],exchange_model = avg_control_peaks[peak['name']][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["fit"] = quantifier_fit.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
                    result[condition][filename][peak["name"]][pick_alg]["fit_exp_cont"] = quantifier_fit_exp_cont.quantify_from_peaklist(picking[condition][filename][peak["name"]][pick_alg])
    return result


def hendersson(x, kmax, pka):
    return kmax/(1.0 + np.power(10,pka-x))

def determine_pka(quantification,peak_list,conditions,histidines):
    result = {}
    pickalgdic = {}
    quantalgsdict = {}
    for histidine in histidines:
        result[histidine["name"]] = {}
        for fragment in histidine["fragments"]:
            result[histidine["name"]][fragment] = {}
            result[histidine["name"]][fragment]["pH"] = []
            result[histidine["name"]][fragment]["filename"] = []
            result[histidine["name"]][fragment]["exchange"] = {}
            result[histidine["name"]][fragment]["uptake"] = {}
            result[histidine["name"]][fragment]["pKa"] = {}
            result[histidine["name"]][fragment]["kmax"] = {}
            for condition in conditions:
                if condition["control"]:
                    continue
                for filename in condition["files"]:
                    result[histidine["name"]][fragment]["pH"].append(condition["ph"])
                    result[histidine["name"]][fragment]["filename"].append(filename)
                    for pickalg in quantification[condition["name"]][filename][fragment].keys():
                        pickalgdic[pickalg] = 1
                        if not pickalg in result[histidine["name"]][fragment]["exchange"]:
                            result[histidine["name"]][fragment]["exchange"][pickalg] = {}
                        if not pickalg in result[histidine["name"]][fragment]["uptake"]:
                            result[histidine["name"]][fragment]["uptake"][pickalg] = {}
                        for quant_alg in quantification[condition["name"]][filename][fragment][pickalg].keys():
                            quantalgsdict[quant_alg] = 1
                            if not quant_alg in result[histidine["name"]][fragment]["exchange"][pickalg]:
                                result[histidine["name"]][fragment]["exchange"][pickalg][quant_alg] = []
                            if not quant_alg in result[histidine["name"]][fragment]["uptake"][pickalg]:
                                result[histidine["name"]][fragment]["uptake"][pickalg][quant_alg] = []
                            exchange = quantification[condition["name"]][filename][fragment][pickalg][quant_alg]
                            rate = - np.log((1 - exchange)/1) / 72
                            result[histidine["name"]][fragment]["exchange"][pickalg][quant_alg].append(rate)
                            result[histidine["name"]][fragment]["uptake"][pickalg][quant_alg].append(exchange)
            for pickalg in pickalgdic.keys():
                result[histidine["name"]][fragment]["pKa"][pickalg] = {}
                result[histidine["name"]][fragment]["kmax"][pickalg] = {}
                for quantalg in quantalgsdict:  
                    try:    
                        fitpars, covmat = scipy.optimize.curve_fit(hendersson, result[histidine["name"]][fragment]["pH"],result[histidine["name"]][fragment]["exchange"][pickalg][quantalg], p0=[0.02,6.5])
                        result[histidine["name"]][fragment]["pKa"][pickalg][quantalg] = fitpars[1]
                        result[histidine["name"]][fragment]["kmax"][pickalg][quantalg] = fitpars[0]
                    except:
                        print "Upps"
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
        #self._generate_models()

    def quantify_from_peaklist(self):
        raise NotImplementedError( "Quantification from peaklist should be implemented by children" )

    def quantify_from_signal(self, signal):

        peaklist = np.zeros(len(self.mzs))

        for idx, mz in enumerate(self.mzs):
            peak = mspy.labelpeak(signal,mz)
            if peak is not None:
                peaklist[idx] = peak.ai

        return(self.quantify_from_peaklist(peaklist))


    def _generate_models(self,numpeaks):
        self.mzs = map(lambda x: x[0], self.no_exchange_model)
        self.mzs.append(self.mzs[-1]+1.003)
        self.models = np.zeros((self.scales+1,numpeaks))
        for idx in range(self.scales+1):
            for jdx, intensity in enumerate(map(lambda x: x[1], self.no_exchange_model)):
                if idx+jdx < numpeaks: 
                    self.models[idx][idx+jdx] = intensity

class exchange_quantifier_fit(exchange_quantifier):
    def quantify_from_peaklist(self, peaklist):
        self._generate_models(len(peaklist))
        try:
            result = self._leastSquare(np.array([a[1] for a in peaklist]),self.models,iterLimit=3000)
        except np.linalg.linalg.LinAlgError:
            result = np.array([1,0])
        result = result/sum(result)
        a=0
        for i,value in enumerate(result):
           a += i*value 
        return a 

    def _leastSquare(self, data, models, iterLimit=None, chiLimit=1e-3):
        """Least-square fitting. Adapted from the original code by Konrad Hinsen."""

        normf = 100./np.max(data)
        data *= normf

        params = [50.] * len(models)
        id = np.identity(len(params))
        chisq, alpha = self._chiSquare(data, models, params)
        l = 0.001

        niter = 0
        while True:


            niter += 1
            delta = solveLinEq(alpha+l*np.diagonal(alpha)*id,-0.5*np.array(chisq[1]))
            next_params = map(lambda a,b: a+b, params, delta)

            for x in range(len(next_params)):
                if next_params[x] < 0.:
                    next_params[x] = 0.

            next_chisq, next_alpha = self._chiSquare(data, models, next_params)
            if next_chisq[0] > chisq[0]:
                l = 5.*l
            elif chisq[0] - next_chisq[0] < chiLimit:
                break
            else:
                l = 0.5*l
                params = next_params
                chisq = next_chisq
                alpha = next_alpha

            if iterLimit and niter == iterLimit:
                break

        next_params /= normf

        return next_params
    
    def _chiSquare(self, data, models, params):
        """Calculate fitting chi-square for current parameter set."""

        # calculate differences and chi-square value between calculated and real data
        differences = np.sum(models * [[x] for x in params], axis=0) - data
        chisq_value = np.sum(differences**2)

        # calculate chi-square deriv and alpha
        cycles = len(models)
        chisq_deriv = cycles*[0]
        alpha = np.zeros((len(params), len(params)))
        for x in range(len(data)):

            deriv = cycles*[0]
            for i in range(cycles):
                p_deriv = cycles*[0]
                p_deriv[i] = models[i][x]
                deriv = map(lambda a,b: a+b, deriv, p_deriv)
            chisq_deriv = map(lambda a,b: a+b, chisq_deriv, map(lambda x,f=differences[x]*2:f*x, deriv))

            d = np.array(deriv)
            alpha = alpha + d[:,np.newaxis]*d

        return [chisq_value, chisq_deriv], alpha

class exchange_quantifier_peak_ratio(exchange_quantifier):
    def quantify_from_peaklist(self, peaklist):
        if peaklist[0][1] == 0:
            return(0)
        ratio = peaklist[1][1]/peaklist[0][1]
        orig_ratio = self.no_exchange_model[1][1]/self.no_exchange_model[0][1]
        inc = 1 - ( 1 / (1 + ratio - orig_ratio))
        return(inc)

class exchange_quantifier_average_mass(exchange_quantifier):

    def quantify_from_peaklist(self,peaklist):

        intensities_noex = np.array(map(lambda x: x[1], self.no_exchange_model))
        intensities = np.array(map(lambda x: x[1], peaklist))
        rel_masses = np.arange(intensities.shape[0])
        rel_masses_noex = np.arange(intensities_noex.shape[0])
        weighted_masses = np.dot(intensities,rel_masses)
        weighted_masses_noex = np.dot(intensities_noex,rel_masses_noex)
        av_mass = np.sum(weighted_masses) / np.sum(intensities) - np.sum(weighted_masses_noex) / np.sum(intensities_noex)
        return(av_mass)

# Main program starts

#Parsing arguments
parser = argparse.ArgumentParser(description="Analyze HDX exchange in proteolytic peptides")
parser.add_argument('conditions', metavar ='conditions', type=str, help="File with list of conditions in json")
parser.add_argument('peaks', metavar ='peaks', type=str, help="File with list of peaks to analyze")
parser.add_argument('histidines', metavar ='peaks', type=str, help="File with list histidines and corresponding ions ")
parser.add_argument('--write-parsing-pickle', action="store", dest="writepickle")
parser.add_argument('--read-parsing-pickle', action="store", dest="readpickle")
args = parser.parse_args()


#Reading input files
with open(args.conditions) as conditions_file:
    conditions = json.load(conditions_file)
with open(args.peaks) as peaks_file:
    peaklist = json.load(peaks_file)
with open(args.histidines) as histidines_file:
    histidines = json.load(histidines_file)


int_data = {}
usable_scan_list = {}

seq_objects = generate_seq_objects(peaklist)

if args.readpickle is None:
    for cond in conditions:
        int_data[cond["name"]] = extract_intensities(cond["files"],peaklist)
        usable_scan_list[cond["name"]] = create_usable_scan_list(int_data[cond["name"]])

    profiles = average_profiles_from_usable_scans(peaklist,usable_scan_list)

    picking = pick_peaks_from_profiles(profiles, peaklist, seq_objects)
    if args.writepickle is not None:
        with open(args.writepickle,'wb') as picklefile:
            pickle.dump((int_data,usable_scan_list,profiles,picking,peaklist),picklefile)
else:
    with open(args.readpickle,'rb') as picklefile:
        int_data,usable_scan_list,profiles,picking,peaklist = pickle.load(picklefile)
quantification = quantify_exchange_from_picking(picking, peaklist, seq_objects, conditions)
pka_determination = determine_pka(quantification, peaklist, conditions, histidines)

#Output results
result = {}
result["input"] = {"conditions":conditions,"peaklist":peaklist,"histidines":histidines}
result["int_data"] = int_data
result["scans_used_for_averaging"] = usable_scan_list
result["profiles"] = profiles
result["picking"] = picking
result["quantification"] = quantification
result["pka_determination"] = pka_determination

reportDataPath = 'data.js'
reportDataFile = file(reportDataPath, 'w')

reportDataFile.write("data=".encode("utf-8"))
jsondump = json.dumps(result,cls=NumpyAwareJSONEncoder)
reportDataFile.write(jsondump)
reportDataFile.close()

