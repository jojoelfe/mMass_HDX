import pymzml
import numpy as np
import mod_signal
import calculations
from mspy import *
import matplotlib.pyplot as plt

            

def extract_peaks(file,time_range):
    peak_data = {}
    
    run = pymzml.run.Reader(file+'.mzML')
    peak_data[file] = {}
    print(file)
    for spec in run:

        if spec.get('filter string')!= None and spec.get('total ion current') !=  None and spec.get('scan time') >= time_range[0] and spec.get('scan time') <= time_range[1]:

            if spec.get('filter string') in peak_data[file]:
                peak_data[file][spec.get('filter string')]['scan_count'] += 1
                peak_data[file][spec.get('filter string')]['scan'] = mod_signal.combine(peak_data[file][spec.get('filter string')]['scan'], np.array(spec.peaks))
                calculations.signal_filter(peak_data[file][spec.get('filter string')]['scan'],0.005)
            else:
                peak_data[file][spec.get('filter string')] = {}
                peak_data[file][spec.get('filter string')]['scan_count'] = 1
                peak_data[file][spec.get('filter string')]['scan'] = np.array(spec.peaks)

    
    for scan in peak_data[file].keys():
        peak_data[file][scan]['scan'] = mod_signal.multiply(peak_data[file][scan]['scan'],y=1.0/peak_data[file][scan]['scan_count'])
        with open(file+scan+'.np', 'wb') as f:
            np.save(f,peak_data[file][scan]['scan'])


def create_fragment_list(sequence_string,scanfiles):
    scan_data = {}
    for scan in scanfiles.keys():
        with open(scanfiles[scan],'rb') as f:
            scan_data[scan] = np.load(f)
    sequence_obj = sequence(sequence_string)
    series = []
    series.append('b')
    series.append('y')
    fragments = mod_proteo.fragment(
                sequence = sequence_obj,
                series = series,
                scrambling = False
            )
    for fragment in fragments:
        print(fragment.format('f'))
        for z in range(1,8):
            monomass = fragment.mz(z)[0]
            if (monomass> 350 and monomass < 1200):
                print(z)
                for scan in scanfiles.keys():
                    plt.plot(*mod_signal.crop(scan_data[scan],monomass-0.1,monomass+2).T)
                    plt.show()
                    plt.clf()
                
                
    
            
