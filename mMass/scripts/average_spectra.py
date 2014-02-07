#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      jojot_000
#
# Created:     03/11/2013
# Copyright:   (c) jojot_000 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import sys
sys.path.append('/Users/johannes/Documents/mMass_HDX/mMass/')
import pymzml
import numpy as np
import mspy




def extract_peaks(file,time_range):
    peak_data = {}

    run = pymzml.run.Reader(file)
    peak_data[file] = {}
    print(file)
    for spec in run:

        if spec.get('filter string')!= None and spec.get('total ion current') !=  None and spec.get('scan time') >= time_range[0] and spec.get('scan time') <= time_range[1]:

            if spec.get('filter string') in peak_data[file]:
                peak_data[file][spec.get('filter string')]['scan_count'] += 1
                peak_data[file][spec.get('filter string')]['scan'] = mspy.combine(peak_data[file][spec.get('filter string')]['scan'], np.array(spec.peaks))

            else:
                peak_data[file][spec.get('filter string')] = {}
                peak_data[file][spec.get('filter string')]['scan_count'] = 1
                peak_data[file][spec.get('filter string')]['scan'] = np.array(spec.peaks)


    for scan in peak_data[file].keys():
        peak_data[file][scan]['scan'] = mspy.reduce(peak_data[file][scan]['scan'])
        peak_data[file][scan]['scan'] = mspy.multiply(peak_data[file][scan]['scan'],y=1.0/peak_data[file][scan]['scan_count'])
        with open(file+scan+'.np', 'wb') as f:
            np.save(f,peak_data[file][scan]['scan'])


def main():
    file_list=['01_H20_0','02_H20_72','03_3_5','04_4_0','05_4_5','06_5_0','07_5_5','08_6_0','09_6_5','10_7_0','11_7_5','12_8_0','13_8_5','14_9_0',]
    elution_time = (3.5,6.5)

    #Average spectra during elution and create data structurre
    for file in file_list:
        extract_peaks('../fimc_hdx_pH_rep3' + file +'.mzML',elution_time)





    pass

if __name__ == '__main__':
    main()
