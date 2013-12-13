
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


def generate_fragment_report(fragments, documents, sequence, rms_cutoff=0.2):
    buff = REPORT_HEADER
    buff += '<h1>Fragment report of mMass spectra</h1>'
    fragment_dict = defaultdict(dict)
    check_pattern_dict = defaultdict(dict)
    conf_dict = {}
    det_dict = {}
    for item in fragments:
        fragment_dict[item[0]][int(item[3])] = item[6]
        check_pattern_dict[item[0]][int(item[3])] = {}

        pattern_obj = pattern(item[6].formula(), charge=int(item[3]), real=False)
        for document_obj in documents:
            check_pattern = checkpattern(signal=document_obj.spectrum.profile,
                                         pattern=pattern_obj)
            check_pattern_dict[item[0]][
                int(item[3])][document_obj.title] = check_pattern
            if check_pattern is not None and check_pattern.rms < 0.2:
                if item[0] in det_dict:
                    if det_dict[item[0]] is not None \
                            and check_pattern.rms < det_dict[item[0]][0]:
                        det_dict[item[0]] = (check_pattern,
                                             document_obj.title, item)
                else:
                    det_dict[item[0]] = (check_pattern,
                                         document_obj.title, item)
            if check_pattern is not None and check_pattern.rms < 0.1:
                det_dict[item[0]] = None
                if item[0] in conf_dict:
                    if check_pattern.rms < conf_dict[item[0]][0]:
                        conf_dict[item[0]] = (check_pattern,
                                              document_obj.title, item)
                else:
                    conf_dict[item[0]] = (check_pattern,
                                          document_obj.title, item)

    fragment_list_sorted = sorted(fragment_dict.keys(),
                                  key=lambda fragment_name: fragment_dict
                                  [fragment_name][1].history[-1][1])
    fragment_list_sorted = sorted(fragment_list_sorted,
                                  key=lambda fragment_name: fragment_dict
                                  [fragment_name][1].history[-1][2])

    # Add sequence logo

    buff += "<h2>Sequence Logo</h2>"
    buff += D3_SCRIPT_HEADER
    buff += "var alphabet = \"{0}\".split(\"\");".format(sequence.format())
    buff += "var width = 960, height = 54 * {0} + 32;\n"\
            .format(int(len(sequence.format()) / 40) + 1)
    nfrag_list = []
    cfrag_list = []
    for fragment_name in fragment_list_sorted:
        history = fragment_dict[fragment_name][1].history[1]
        stroke = None
        if fragment_name in conf_dict:
            stroke = 4
        elif fragment_name in det_dict:
            stroke = 1
        if stroke is not None:
            if history[2] == len(sequence.format()) and history[1] != 0:
                cfrag_list.append(
                    "{{ \"i\": {0}, \"s\": {1}, \"name\": \"{2}\" }}"
                    .format(history[1], stroke, fragment_name))
            elif history[1] == 0 and history[2] != len(sequence.format()):
                nfrag_list.append(
                    "{{ \"i\": {0}, \"s\": {1}, \"name\": \"{2}\" }}"
                    .format(history[2], stroke, fragment_name))

    buff += "var nfrag_list = [{0}];".format(",".join(nfrag_list))
    buff += "var cfrag_list = [{0}];".format(",".join(cfrag_list))
    buff += TRANS_FUNCTION
    buff += "update(alphabet);"
    buff += "</script>"

    # Add list of certain and potential identifications
    buff += "<h2>List of fragments found with confidence \
        (rms &lt; 0.1)</h2><table><tr><th>\
        Fragment</th><th>RMS</th><th>z</th><th>m/z</th><th>Document</th></tr>"
    for fragment_name in fragment_list_sorted:
        if fragment_name in conf_dict:
            buff += "<tr>"
            buff += "<td>{0}</td>".format(fragment_name)
            buff += "<td>{}</td>".format(conf_dict[fragment_name][0].format())
            buff += "<td>{0}</td>".format(conf_dict[fragment_name][2][3])
            buff += "<td>{0}</td>".format(conf_dict[fragment_name][2][2])
            buff += "<td>{0}</td>".format(conf_dict[fragment_name][1])
            buff += "</tr>"
    buff += "</table>"

    # Add list of certain and potential identifications
    buff += "<h2>List of fragments potentially found \
        (0.1 &lt; rms &lt; 0.2)</h2><table><tr><th>\
        Fragment</th><th>RMS</th><th>z</th><th>m/z</th><th>Document</th></tr>"
    for fragment_name in fragment_list_sorted:
        if fragment_name in det_dict and det_dict[fragment_name] is not None:
            buff += "<tr>"
            buff += "<td>{0}</td>".format(fragment_name)
            buff += "<td>{}</td>".format(det_dict[fragment_name][0].format())
            buff += "<td>{0}</td>".format(det_dict[fragment_name][2][3])
            buff += "<td>{0}</td>".format(det_dict[fragment_name][2][2])
            buff += "<td>{0}</td>".format(det_dict[fragment_name][1])
            buff += "</tr>"
    buff += "</table>"

    # Add Fragment Tables
    for fragment_name in fragment_list_sorted:
        buff += "<h3>{0}</h3>".format(fragment_name)
        buff += "<p>{0}</p>".format(fragment_dict[fragment_name][1].format())
        buff += "<table><tr><th>Scan</th>"
        for z in sorted(fragment_dict[fragment_name].keys()):
            buff += "<th>{0}</th>".format(z)
        buff += "</tr>"
        for document_name in sorted(check_pattern_dict
                                    [fragment_name][1].keys()):
            buff += "<tr><td>{0}</td>".format(document_name)
            for z in sorted(fragment_dict[fragment_name].keys()):
                if check_pattern_dict[
                        fragment_name][z][document_name] is not None:
                    buff += "<td>{}</td>".format(
                        check_pattern_dict[fragment_name][z][document_name]
                        .format())
                else:
                    buff += "<td></td>"
            buff += "</tr>"
        buff += "</table>"

    buff += '  <p id="footer">Generated by mMass &bull; Open Source Mass \
        Spectrometry Tool &bull; <a href="http://www.mmass.org/" \
        title="mMass homepage">www.mmass.org</a></p>\n'
    buff += '</body>\n'
    buff += '</html>'
    return buff


TRANS_FUNCTION = """
var svg = d3.select("body").append("svg")
    .attr("width", width)
    .attr("height", height)
  .append("g")
    .attr("transform", "translate(32,32)");
"""
REPORT_HEADER = """<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="cs" lang="cs">
<head>
  <meta http-equiv="content-type" content="text/html; charset=utf-8" />
  <meta name="author" content="Created by mMass - Open Source Mass Spectrometry Tool; www.mmass.org" />
  <title>mMass Report</title>
  <style type="text/css">
  <!--
    body{margin: 5%; font-size: 8.5pt; font-family: Arial, Verdana, Geneva, Helvetica, sans-serif;}
    h1{font-size: 1.5em; text-align: center; margin: 1em 0; border-bottom: 3px double #000;}
    h1 span{font-style: italic;}
    h2{font-size: 1.2em; text-align: left; margin: 2em 0 1em 0; border-bottom: 1px solid #000;}
    h2 span{font-style: italic;}
    table{border-collapse: collapse; margin: 1.5em auto; width: 100%; background-color: #fff;}
    thead{display: table-header-group;}
    th,td{font-size: .75em; border: 1px solid #aaa; padding: .3em; vertical-align: top; text-align: left;}
    html>body th, html>body td{font-size: .9em;}
    th{text-align: center; color: #000; background-color: #ccc;}
    th a{text-align: center; color: #000; background-color: #ccc; text-decoration: none;}
    #tableMainInfo th{text-align: right; width: 15%;}
    #tableMainInfo td{text-align: left;}
    #spectrum{text-align: center;}
    #footer{font-size: .8em; font-style: italic; text-align: center; color: #aaa; margin: 2em 0 1em 0; padding-top: 0.5em; border-top: 1px solid #000;}
    .left{text-align: left;}
    .right{text-align: right;}
    .center{text-align: center;}
    .nowrap{white-space:nowrap;}
    .sequence{font-size: 1.1em; font-family: monospace;}
    .modified{color: #f00; font-weight: bold;}
    .matched{text-decoration: underline;}

 text {
  font: bold 24px monospace;
}

.enter {
  fill: black;
}

.fragment_label  {
   font-size: 7pt;
   font-family: Arial, Verdana, Geneva, Helvetica, sans-serif;
   text-anchor: middle;
   }

.nfrag {
  stroke: blue;
  fill: none;
}

.cfrag {
  stroke: red;
  fill: none;
}

.update {
  fill: #000;
}

  -->
  </style>
  <script type="text/javascript">
    // This script was adapted from the original script by Mike Hall (www.brainjar.com)
    //<![CDATA[

    // for IE
    if (document.ELEMENT_NODE == null) {
      document.ELEMENT_NODE = 1;
      document.TEXT_NODE = 3;
    }

    // sort table
    function sortTable(id, col) {

      // get table
      var tblEl = document.getElementById(id);

      // init sorter
      if (tblEl.reverseSort == null) {
        tblEl.reverseSort = new Array();
      }

      // reverse sorting
      if (col == tblEl.lastColumn) {
        tblEl.reverseSort[col] = !tblEl.reverseSort[col];
      }

      // remember current column
      tblEl.lastColumn = col;

      // sort table
      var tmpEl;
      var i, j;
      var minVal, minIdx;
      var testVal;
      var cmp;

      for (i = 0; i < tblEl.rows.length - 1; i++) {
        minIdx = i;
        minVal = getTextValue(tblEl.rows[i].cells[col]);

        // walk in other rows
        for (j = i + 1; j < tblEl.rows.length; j++) {
          testVal = getTextValue(tblEl.rows[j].cells[col]);
          cmp = compareValues(minVal, testVal);

          // reverse sorting
          if (tblEl.reverseSort[col]) {
            cmp = -cmp;
          }

          // set new minimum
          if (cmp > 0) {
            minIdx = j;
            minVal = testVal;
          }
        }

        // move row before
        if (minIdx > i) {
          tmpEl = tblEl.removeChild(tblEl.rows[minIdx]);
          tblEl.insertBefore(tmpEl, tblEl.rows[i]);
        }
      }

      return false;
    }

    // get node text
    function getTextValue(el) {
      var i;
      var s;

      // concatenate values of text nodes
      s = "";
      for (i = 0; i < el.childNodes.length; i++) {
        if (el.childNodes[i].nodeType == document.TEXT_NODE) {
          s += el.childNodes[i].nodeValue;
        } else if (el.childNodes[i].nodeType == document.ELEMENT_NODE && el.childNodes[i].tagName == "BR") {
          s += " ";
        } else {
          s += getTextValue(el.childNodes[i]);
        }
      }

      return s;
    }

    // compare values
    function compareValues(v1, v2) {
      var f1, f2;

      // lowercase values
      v1 = v1.toLowerCase()
      v2 = v2.toLowerCase()

      // try to convert values to floats
      f1 = parseFloat(v1);
      f2 = parseFloat(v2);
      if (!isNaN(f1) && !isNaN(f2)) {
        v1 = f1;
        v2 = f2;
      }

      // compare values
      if (v1 == v2) {
        return 0;
      } else if (v1 > v2) {
        return 1;
      } else {
        return -1;
      }
    }

    //]]>
  </script>
</head>

<body>
"""

D3_SCRIPT_HEADER = """
<script src="http://d3js.org/d3.v2.min.js?2.10.1"></script>
<script>
function update(data) {

  // DATA JOIN
  // Join new data with old elements, if any.
  var text = svg.selectAll("text.seq")
      .data(data);
  var nfrag = svg.selectAll("path").data(nfrag_list)
  var cfrag = svg.selectAll("path").data(cfrag_list)
  // UPDATE
  // Update old elements as needed.
  var line = d3.svg.line()
    .interpolate("linear")
    .x(function(d,i) { return d.x; })
    .y(function(d,i) { return d.y; });
  // ENTER
  // Create new elements as needed.
  text.enter().append("text")
      .attr("class", "enter")
      .attr("x", function(d, i) { var j = i % 40;return j * 22 + ~~(j/10) * 8; })
      .attr("y", function(d, i) { var j = ~~(i / 40);return j * 54; })
      .attr("dy", ".35em");
  nfrag.enter().append("path")
      .attr("class", "nfrag")
      .attr("d", function(d) {return line(nfrag_path);})
      .attr("transform", function(d) { var j = (d.i-1) % 40; dx = j * 22 + ~~((j)/10) * 8; dy = ~~((d.i-1)/40) * 54; return "translate(" + dx + "," + dy + ")"; })
      .attr("stroke-width", function(d) {return d.s})
  svg.selectAll("text.labeln").data(nfrag_list).enter().append("text")
     .attr("class", "fragment_label")
     .attr("x", function (d) { var j = (d.i-1) % 40; dx = j * 22 + ~~((j)/10) * 8 + 12; return dx; } )
     .attr("y",function (d) {  dy = ~~((d.i-1)/40) * 54 + 25; return dy })
     .text( function(d) { return d.name; })
  cfrag.enter().append("path")
      .attr("class", "cfrag")
      .attr("d", function(d) {return line(cfrag_path);})
      .attr("transform", function(d,i) { var j = (d.i) % 40; dx = j * 22 + ~~((j)/10) * 8 - 22; dy = ~~((d.i)/40) * 54; return "translate(" + dx + "," + dy + ")"; })
      .attr("stroke-width", function(d) {return d.s})
  svg.selectAll("text.labelc").data(cfrag_list).enter().append("text")
     .attr("class", "fragment_label")
     .attr("x", function (d) { var j = (d.i) % 40; dx = j * 22 + ~~((j)/10) * 8 +2  ; return dx; } )
     .attr("y",function (d) {  dy = ~~((d.i)/40) * 54 - 14 ; return dy })
     .text( function(d) { return d.name; })
  // ENTER + UPDATE
  // Appending to the enter selection expands the update selection to include
  // entering elements; so, operations on the update selection after appending to
  // the enter selection will apply to both entering and updating nodes.
  text.text(function(d) { return d; });

  // EXIT
  // Remove old elements as needed.
  text.exit().remove();
}
var cfrag_path = [ { "x": 30, "y": -12 }, { "x" : 18, "y": -12}, { "x" : 18, "y" : 0}]
var nfrag_path = [ { "x": 6, "y": 14 }, { "x" : 18, "y": 14}, { "x" : 18, "y" : 3}]
"""
