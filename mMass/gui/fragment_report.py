
# -------------------------------------------------------------------------
#     Copyright (C) 2013-2013 Johannes Elferich

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file LICENSE.TXT in the
#     main directory of the program.
# -------------------------------------------------------------------------

#TODO:
# Add more infotmation (baspeak intensity, error in x)
# Add plots (matplotlib or d3.js) ?
# Put panel before html generation

from mspy import * # NOQA
from collections import defaultdict
import re
import sys

# set config folder for MAC OS X
if sys.platform == 'darwin':
    confdir = 'configs'
    support = os.path.expanduser("~/Library/Application Support/")
    userconf = os.path.join(support, 'mMass')
    if os.path.exists(support) and not os.path.exists(userconf):
        try:
            os.mkdir(userconf)
        except:
            pass
    if os.path.exists(userconf):
        confdir = userconf

# set config folder for Linux
elif sys.platform.startswith('linux') or sys.platform.startswith('freebsd'):
    confdir = 'configs'
    home = os.path.expanduser("~")
    userconf = os.path.join(home, '.mmass')
    if os.path.exists(home) and not os.path.exists(userconf):
        try:
            os.mkdir(userconf)
        except:
            pass
    if os.path.exists(userconf):
        confdir = userconf

# set config folder for Windows
else:
    confdir = os.path.sep
    for folder in os.path.dirname(
            os.path.realpath(__file__)).split(os.path.sep)[:-1]:
        path = os.path.join(confdir, folder)
        if os.path.isdir(path):
            confdir = path
        if os.path.isfile(path):
            break
    confdir = os.path.join(confdir, 'configs')
    if not os.path.exists(confdir):
        try:
            os.mkdir(confdir)
        except:
            pass

if not os.path.exists(confdir):
    raise IOError, "Configuration folder cannot be found!"


class FragmentInformation:
    def __init__(self, fragment_name, fragment_object):
        self.fragment_name = fragment_name
        self.fragment_object = fragment_object
        self.charges = []
        self.checkpatternresults = defaultdict(dict)
        self.patternobjects = {}
        self.confident = False
        self.potential = False
        self.bestDocID = None
        self.bestCharge = None
        self.bestRMSD = None
        self.bestResult = None

    def updateStatus(self, result, charge, docID):
        if self.bestRMSD is None or result.rmsd < self.bestRMSD:
            self.bestDocID = docID
            self.bestCharge = charge
            self.bestRMSD = result.rmsd
            self.bestResult = result
            if self.bestRMSD < 0.1:
                self.confident = True
                self.potential = False
            elif self.bestRMSD < 0.2:
                self.potential = True


def runSignalMatch(fragments, documents):
    """Performs matching of framgents against signal"""

    # 1. Dimension: fragment by name
    # Content fragment_information
    fragment_dict = {}
    for item in fragments:
        # Item [0]: fragment name
        # Item [3]: charge
        # Item [6]: fragment_object

        #Initialize object and append charge state
        if item[0] in fragment_dict:
            fragment_dict[item[0]].charges.append(item[3])
        else:
            fragment_dict[item[0]] = FragmentInformation(item[0], item[6])
            fragment_dict[item[0]].charges.append(item[3])

        #Create pattern object
        pattern_obj = pattern(
            item[6].formula(), charge=int(item[3]), real=False)
        fragment_dict[item[0]].patternobjects[item[3]] = pattern_obj

        #Check in each document
        for i, document_obj in enumerate(documents):
            result = checkpattern(signal=document_obj.spectrum.profile,
                                  pattern=pattern_obj)
            fragment_dict[item[0]].checkpatternresults[item[3]][i] = result
            if result is not None:
                fragment_dict[item[0]].updateStatus(result, item[3], i)

    fragment_list_sorted = sorted(fragment_dict.keys(),
                                  key=lambda fragment_name: fragment_dict
                                  [fragment_name].fragment_object
                                  .history[-1][1])
    fragment_list_sorted = sorted(fragment_list_sorted,
                                  key=lambda fragment_name: fragment_dict
                                  [fragment_name].fragment_object
                                  .history[-1][2])

    return fragment_list_sorted, fragment_dict


def get_nfrag_list(fragment_list_sorted, fragment_dict, sequence):
    """Returns json list of fragments to show in logo"""

    nfrag_list = []
    for fragment_name in fragment_list_sorted:
        history = fragment_dict[fragment_name].fragment_object.history[-1]
        stroke = None
        if fragment_dict[fragment_name].confident:
            stroke = 4
        elif fragment_dict[fragment_name].potential:
            stroke = 1
        if stroke is not None:
            if history[1] == 0 and history[2] != len(sequence.format()):
                nfrag_list.append(
                    "{{ \"i\": {0}, \"s\": {1}, \"name\": \"{2}\" }}"
                    .format(history[2], stroke, fragment_name))
    return ",".join(nfrag_list)


def get_cfrag_list(fragment_list_sorted, fragment_dict, sequence):
    """Returns json list of fragments to show in logo"""

    cfrag_list = []
    for fragment_name in fragment_list_sorted:
        history = fragment_dict[fragment_name].fragment_object.history[-1]
        stroke = None
        if fragment_dict[fragment_name].confident:
            stroke = 4
        elif fragment_dict[fragment_name].potential:
            stroke = 1
        if stroke is not None:
            if history[1] != 0 and history[2] == len(sequence.format()):
                cfrag_list.append(
                    "{{ \"i\": {0}, \"s\": {1}, \"name\": \"{2}\" }}"
                    .format(history[2], stroke, fragment_name))
    return ",".join(cfrag_list)


def get_html_conffrags(fragment_list_sorted, fragment_dict, documents):
    """Returns html code to build table of fragments found with confidence"""

    buff = ""
    for fragment_name in fragment_list_sorted:
        if fragment_dict[fragment_name].confident:
            buff += "<tr>"
            buff += "<td>{0}</td>".format(fragment_name)
            buff += "<td>{}</td>".format(fragment_dict[fragment_name]
                                         .bestResult.format())
            buff += "<td>{0}</td>".format(fragment_dict[fragment_name]
                                          .bestCharge)
            buff += "<td>{0}</td>".format(fragment_dict[fragment_name]
                                          .fragment_obj.mass(massType=0)
                                          / fragment_dict[fragment_name]
                                          .bestCharge)
            buff += "<td>{0}</td>".format(documents[
                    fragment_dict[fragment_name].bestDocID].title)
            buff += "</tr>"
    return buff


def get_html_potfrags(fragment_list_sorted, fragment_dict, documents):
    """Returns html code to build table of fragments potentially found"""

    buff = ""
    for fragment_name in fragment_list_sorted:
        if fragment_dict[fragment_name].potential:
            buff += "<tr>"
            buff += "<td>{0}</td>".format(fragment_name)
            buff += "<td>{}</td>".format(fragment_dict[fragment_name]
                                         .bestResult.format())
            buff += "<td>{0}</td>".format(fragment_dict[fragment_name]
                                          .bestCharge)
            buff += "<td>{0}</td>".format(fragment_dict[fragment_name]
                                          .fragment_obj.mass(massType=0)
                                          / fragment_dict[fragment_name]
                                          .bestCharge)
            buff += "<td>{0}</td>".format(documents[
                    fragment_dict[fragment_name].bestDocID].title)
            buff += "</tr>"
    return buff


def get_html_fragtables(fragment_list_sorted, fragment_dict, documents):
    """Returns html code to buidl list of scores for all fragmens"""

    buff = ""
    for fragment_name in fragment_list_sorted:
        buff += "<h3>{0}</h3>".format(fragment_name)
        buff += "<p>{0}</p>".format(fragment_dict[fragment_name]
                                    .fragment_object.format())
        buff += "<table><tr><th>Scan</th>"
        for z in sorted(fragment_dict[fragment_name].charges):
            buff += "<th>{0}</th>".format(z)
        buff += "</tr>"
        for i, document in sorted(enumerate(documents),
                                  key=lambda x: x[1]):
            buff += "<tr><td>{0}</td>".format(document.title)
            for z in sorted(fragment_dict[fragment_name].charges):
                if fragment_dict[fragment_name]\
                        .checkpatternresults[z][i] is not None:
                    buff += "<td>{}</td>".format(
                        fragment_dict[fragment_name].checkpatternresults[z][i]
                        .format())
                else:
                    buff += "<td></td>"
            buff += "</tr>"
        buff += "</table>"
    return buff


def generate_fragment_report(fragments, documents, sequence, rms_cutoff=0.2):
    buff = ""
    fragment_list_sorted, fragment_dict = runSignalMatch(fragments, documents)
    rep_strings = {}
    rep_strings['sequence'] = sequence.format()
    rep_strings['nfraglist'] = get_nfrag_list(fragment_list_sorted,
                                              fragment_dict, sequence)
    rep_strings['cfraglist'] = get_cfrag_list(fragment_list_sorted,
                                              fragment_dict, sequence)
    rep_strings['conffrags'] = get_html_conffrags(fragment_list_sorted,
                                                  fragment_dict, documents)
    rep_strings['potfrags'] = get_html_potfrags(fragment_list_sorted,
                                                fragment_dict, documents)
    rep_strings['fragtables'] = get_html_fragtables(fragment_list_sorted,
                                                    fragment_dict, documents)
    regexp = re.compile(re.compile("\{\{([A-Z]+)\}\}"))
    with open(os.path.join(confdir, 'frag_template.html'),'r') as template:
        for template_line in template:
            buff += re.sub(
                regexp, lambda x: rep_strings[x.group(1).lower()],
                template_line)
    return buff
# ----
