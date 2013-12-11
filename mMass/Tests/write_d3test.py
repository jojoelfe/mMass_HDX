import sys
print 'Importing libaries'
sys.path.append('../')
import mspy
import tempfile
import os


def export_html():
    tmpDir = ''
    reportPath = os.path.join(tmpDir, 'mmass_report.html')
    reportHTML=REPORT_HEADER
    with open('d3.v3.min.js','r') as d3file:
        reportHTML += d3file.read()
    reportHTML += '</script>\n<body><div id="graph" class="aGraph" style="position:absolute;top:0px;left:0; float:left;"></div>'
    mzfile = mspy.parseMZML('t0r3t1.mzML')
    mzfile.load()
    print 'Loaded MZML file'
    scanlist = mzfile.scanlist()
    relevant_scans_20 = []
    for i_scan in scanlist.keys():
        if scanlist[i_scan]['filterString'] == u'FTMS + p ESI sid=35.00  Full ms2 707.17@etd20.00 [500.00-1600.00]' and scanlist[i_scan]['retentionTime'] > 300 and scanlist[i_scan]['retentionTime'] < 480:
            relevant_scans_20.append(i_scan)
    av_scan_20 = None
    for i_scan in relevant_scans_20:
        if av_scan_20 == None:
            av_scan_20 = mzfile.scan(i_scan)
        else:
            av_scan_20 += mzfile.scan(i_scan)
    reportHTML += '<script>\ndata_x = ['
    reportHTML += ','.join(map(str,av_scan_20.profile.T[0]))
    reportHTML += '];\n'
    reportHTML += 'data_y = ['
    reportHTML += ','.join(map(str,av_scan_20.profile.T[1]))
    reportHTML += '];\n'
    reportHTML += GRAPH_SCRIPT
    reportHTML += '</script>'


    reportHTML += '</body>'
    reportFile = file(reportPath, 'w')
    reportFile.write(reportHTML.encode("utf-8"))
    reportFile.close()

GRAPH_SCRIPT = """var m = [80, 80, 80, 80]; // margins
		var w = 1000 - m[1] - m[3]; // width
		var h = 400 - m[0] - m[2]; // height
                var x = d3.scale.linear().domain([d3.min(data_x), d3.max(data_x)]).range([0, w]);
		// Y scale will fit values from 0-10 within pixels h-0 (Note the inverted domain for the y-scale: bigger is up!)
		var y = d3.scale.linear().domain([0, d3.max(data_y)]).range([h, 0]);
                // automatically determining max range can work something like this
                // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);
                data = d3.zip(data_x,data_y)
		// create a line function that can convert data[] into x and y points
		var line = d3.svg.line()
			// assign the X function to plot our line as we wish
			.x(function(d) { 
				// verbose logging to show what's actually being done
				// return the X coordinate where we want to plot this datapoint
				return x(d[0]); 
			})
			.y(function(d) { 
				// verbose logging to show what's actually being done
				// return the Y coordinate where we want to plot this datapoint
				return y(d[1]); 
			})

			// Add an SVG element with the desired dimensions and margin.
			var graph = d3.select("#graph").append("svg:svg")
			      .attr("width", w + m[1] + m[3])
			      .attr("height", h + m[0] + m[2])
			    .append("svg:g")
			      .attr("transform", "translate(" + m[3] + "," + m[0] + ")");

			// create yAxis
			var xAxis = d3.svg.axis().scale(x).tickSize(-h).tickSubdivide(true);
			// Add the x-axis.
			graph.append("svg:g")
			      .attr("class", "x axis")
			      .attr("transform", "translate(0," + h + ")")
			      .call(xAxis);


			// create left yAxis
			var yAxisLeft = d3.svg.axis().scale(y).ticks(4).orient("left");
			// Add the y-axis to the left
			graph.append("svg:g")
			      .attr("class", "y axis")
			      .attr("transform", "translate(-25,0)")
			      .call(yAxisLeft);
			
  			// Add the line by appending an svg:path element with the data line we created above
			// do this AFTER the axes above so that the line is above the tick-lines
  			graph.append("svg:path").attr("d", line(data));
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
			/* tell the SVG path to be a thin blue line without any area fill */
			path {
				stroke: steelblue;
				stroke-width: 1;
				fill: none;
			}
			
			.axis {
			  shape-rendering: crispEdges;
			}

			.x.axis line {
			  stroke: lightgrey;
			}

			.x.axis .minor {
			  stroke-opacity: .5;
			}

			.x.axis path {
			  display: none;
			}

			.y.axis line, .y.axis path {
			  fill: none;
			  stroke: #000;
			}
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
  -->
  </style>
  <script type="text/javascript">
"""

if __name__ == "__main__":
    export_html()
