#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      jojot_000
#
# Created:     02/11/2013
# Copyright:   (c) jojot_000 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import pymzml
import matplotlib.pyplot as plt

def extract_ic(files,mass_list):
    ic_data = {}
    for file in files:
        run = pymzml.run.Reader('../fimc_hdx_pH_rep3'+file+'.mzML')
        ic_data[file] = {}
        for spec in run:
            if spec.get('filter string')!= None and spec.get('total ion current') !=  None:
                if spec.get('filter string') in ic_data[file]:
                    ic_data[file][spec.get('filter string')]['scan time'].append(spec.get('scan time'))
                    ic_data[file][spec.get('filter string')]['tic'].append(spec.get('total ion current'))
                    for mass in mass_list:
                        isum = 0
                        for mz,i in spec.peaks:
                            if mz >= mass[0] and mz <= mass[1]:
                                isum += i


                        ic_data[file][spec.get('filter string')][mass].append(isum)
                else:
                    ic_data[file][spec.get('filter string')] = {}
                    ic_data[file][spec.get('filter string')]['scan time'] = [spec.get('scan time')]
                    ic_data[file][spec.get('filter string')]['tic'] = [spec.get('total ion current')]
                    for mass in mass_list:
                        isum = 0
                        for mz,i in spec.peaks:
                            if mz >= mass[0] and mz <= mass[1]:
                                isum += i
                        ic_data[file][spec.get('filter string')][mass] = [isum]
    return(ic_data)


def main():
    mass_list=[(706.8,707.7),(708.4,708.5)]
    file_list=['01_H20_0','02_H20_72','03_3_5','04_4_0','05_4_5','06_5_0','07_5_5','08_6_0','09_6_5','10_7_0','11_7_5','12_8_0','13_8_5','14_9_0',]
    ic_data = extract_ic(file_list,mass_list)
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(3,1,1)
    for file in file_list:
        ax.plot(ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]']['scan time'],ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]']['tic'],label=file)
    ax.legend()
    ax.set_title('Total Ion Current')
    ax = fig.add_subplot(3,1,2)
    for file in file_list:
        ax.plot(ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]']['scan time'],ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]'][(706.8, 707.7)])
    ax.set_title('Ion current mz 706.8-707.7 M14')
    ax = fig.add_subplot(3,1,3)
    for file in file_list:
        ax.plot(ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]']['scan time'],ic_data[file]['FTMS + p ESI sid=35.00  Full ms [400.00-1300.00]'][(708.4,708.5)])
    ax.set_title('Ion current mz 708.4-708.5 PEG contamination')
    fig.tight_layout()
    fig.savefig('Chromatograms.pdf', bbox_inches=0,dpi=300)
    pass

if __name__ == '__main__':
    main()
