function binarySearch(items, value){

    var startIndex  = 0,
        stopIndex   = items.length - 1,
        middle      = Math.floor((stopIndex + startIndex)/2);

    while(items[middle] != value && startIndex < stopIndex){

        //adjust search area
        if (value < items[middle]){
            stopIndex = middle - 1;
        } else if (value > items[middle]){
            startIndex = middle + 1;
        }

        //recalculate middle
        middle = Math.floor((stopIndex + startIndex)/2);
    }

    //make sure it's the right value
    return middle;
}

function initiate_draw_time_graph() {
    setTimeout(draw_time_graph.bind(this), 50);
}

function draw_time_graph() {
    var peptide = this.getAttribute("data-peptide");
    var z = this.getAttribute("data-charge");
    var m = [10, 50, 30, 60]; // margins
    var w = parseInt(d3.select(this)
        .style("width")
        .replace('px', '')) - m[1] -
        m[3]; // width
    var h = parseInt(d3.select(this)
        .style("height")
        .replace('px', '')) - m[0] -
        m[2]; // height
    var data_time = data["RetentionTime"];
    var data_basepeak = data["Intensities"][peptide][z];
    //var data_rmsd = data["Rmsd"][peptide][z];
    //var data_ident = {};
    //if (peptide in data["Ident"]) {
    //data_ident = data["Ident"][peptide][z];
    //}
    var cond_data = [{
        name: name,
        values: d3.zip(data_time, data_basepeak)
    }];
    var ymax = d3.max(cond_data[0].values.map(function (datapair) {
        return datapair[1]
    }));
    if (cond_data[0].values.length == 0) {
        return;
    }
    var minx = cond_data[0].values[0][0];
    var maxx = cond_data[0].values[cond_data[0].values.length - 1][0];
    var x = d3.scale.linear()
        .domain([
            minx,
            maxx
        ])
        .range([0, w]);
    // Y scale will fit values from 0-10 within pixels h-0 (Note the inverted domain for the y-scale: bigger is up!)
    var y_basepeak = d3.scale.linear()
        .domain([0, ymax])
        .range([h, 0]);
    var y_rmsd = d3.scale.linear()
        .domain([0, 0.45])
        .range([h, 0]);
    // automatically determining max range can work something like this
    // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);
    // create a line function that can convert data[] into x and y points
    var line_basepeak = d3.svg.line()
    // assign the X function to plot our line as we wish
    .x(function (d) {
        // verbose logging to show what's actually being done
        // return the X coordinate where we want to plot this datapoint
        return x(d[0]);
    })
        .y(function (d) {
            // verbose logging to show what's actually being done
            // return the Y coordinate where we want to plot this datapoint
            return y_basepeak(d[1]);
        });
    var line_rmsd = d3.svg.line()
    // assign the X function to plot our line as we wish
    .x(function (d) {
        // verbose logging to show what's actually being done
        // return the X coordinate where we want to plot this datapoint
        return x(d[0]);
    })
        .y(function (d) {
            // verbose logging to show what's actually being done
            // return the Y coordinate where we want to plot this datapoint
            return y_rmsd(0.5)
        });
    // Add an SVG element with the desired dimensions and margin.
    var graph = d3.select(this)
        .append("svg:svg")
        .attr("width", w + m[1] + m[3])
        .attr("height", h + m[0] + m[2])
        .append("svg:g")
        .attr("transform", "translate(" + m[3] + "," + m[0] + ")");
    // create yAxis
    var xAxis = d3.svg.axis()
        .scale(x)
        .tickSize(6, -h)
        .ticks(20);
    // Add the x-axis.
    graph.append("svg:g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + h + ")")
        .call(xAxis);
    // create left yAxis
    var yAxisLeft = d3.svg.axis()
        .scale(y_basepeak)
        .ticks(4)
        .orient("left")
        .tickSize(
            6, -w)
        .tickFormat(d3.format(",.2e"));
    // Add the y-axis to the left
    graph.append("svg:g")
        .attr("class", "y axis")
        .attr("transform", "translate(0,0)")
        .call(yAxisLeft);
    // create right yAxis
    var yAxisRight = d3.svg.axis()
        .scale(y_rmsd)
        .ticks(4)
        .orient("right")
        .tickSize(
            6, -w)
        .tickFormat(d3.format("1.1"));
    // Add the y-axis to the left
    graph.append("svg:g")
        .attr("class", "y axis")
        .attr("transform", "translate(" + w + ",0)")
        .call(yAxisRight);
    // Add the line by appending an svg:path element with the data line we created above
    // do this AFTER the axes above so that the line is above the tick-lines
    var condition = graph.selectAll(".condition")
        .data(cond_data)
        .enter()
        .append("g")
        .attr("class", "condition")
    // Create Clip-path to hide thick line where rmsd is bad
    condition.append("defs")
        .append("clipPath")
        .attr("id", function (d) {
            return "clippath" + peptide + z + d.name;
        })
        .append("svg:path")
        .attr("class", "rmsd line")
        .attr("d", function (d) {
            return line_rmsd(d.values) + "L" + w + "," + y_rmsd(0.5) + "Z";
        })
        .style("stroke", function (d) {
            return "#FF0000";
        });
    // Create thick line
    condition.append("svg:path")
        .attr("class", "basepeaks line")
        .attr("d",
            function (d) {
                return line_basepeak(d.values);
            })
        .style("stroke", function (d) {
            return "#FF0000";
        })
        .attr("clip-path", function (d) {
            return "url(#clippath" + peptide + z + d.name + ")";
        });
    // Create thin line
    condition.append("svg:path")
        .attr("class", "basepeaku line")
        .attr("d",
            function (d) {
                return line_basepeak(d.values);
            })
        .style("stroke", function (d) {
            return "#FF0000";
        });
    // Create dashed lines where MS2 identification was found
}

function initiate_draw_spectrum() {
    setTimeout(draw_spectrum.bind(this),50);
}

function draw_spectrum() {
    //Initialize variables
    var peptide = this.getAttribute("data-peptide");
    var z = this.getAttribute("data-charge");
    var m = [10, 50, 30, 60]; // margins
    var w = parseInt(d3.select(this)
        .style("width")
        .replace('px', '')) - m[1] -
        m[3]; // width
    var h = parseInt(d3.select(this)
        .style("height")
        .replace('px', '')) - m[0] -
        m[2]; // height
    var data_spectrum = data["ms1Profiles"][peptide][z];
    var cond_data = [{
        name: name,
        values: data_spectrum
    }];
    var minx = 0;
    var maxx = 0;
    var ymax = d3.max(cond_data.map(function (cond) {
        if (cond.values == null) {
            return null;
        }
        minx = cond.values[0][0];
        maxx = cond.values[cond.values.length - 1][0];
        return d3.max(cond.values.map(function (datapair) {
            return datapair[1]
        }));
    }));
    if (ymax == null) {
        return;
    }
    //Setting up scales
    var x = d3.scale.linear()
        .domain([
            minx,
            maxx
        ])
        .range([0, w]);
    // Y scale will fit values from 0-10 within pixels h-0 (Note the inverted domain for the y-scale: bigger is up!)
    var y = d3.scale.linear()
        .domain([0, ymax])
        .range([h, 0]);
    // automatically determining max range can work something like this
    // var y = d3.scale.linear().domain([0, d3.max(data)]).range([h, 0]);
    // create a line function that can convert data[] into x and y points
    //creating line helper
    var line = d3.svg.line()
    // assign the X function to plot our line as we wish
    .x(function (d) {
        // verbose logging to show what's actually being done
        // return the X coordinate where we want to plot this datapoint
        return x(d[0]);
    })
        .y(function (d) {
            // verbose logging to show what's actually being done
            // return the Y coordinate where we want to plot this datapoint
            return y(d[1]);
        })
    // Add an SVG element with the desired dimensions and margin.
    var graph = d3.select(this)
        .append("svg:svg")
        .attr("width", w + m[1] + m[3])
        .attr("height", h + m[0] + m[2])
        .append("svg:g")
        .attr("transform", "translate(" + m[3] + "," + m[0] + ")");
    // create yAxis
    var xAxis = d3.svg.axis()
        .scale(x)
        .tickSize(6, -h)
        .ticks(5);
    // Add the x-axis.
    graph.append("svg:g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + h + ")")
        .call(xAxis);
    // create left yAxis
    var yAxisLeft = d3.svg.axis()
        .scale(y)
        .ticks(4)
        .orient("left")
        .tickSize(6, -
            w)
        .tickFormat(d3.format(",.2e"));
    // Add the y-axis to the left
    graph.append("svg:g")
        .attr("class", "y axis")
        .attr("transform", "translate(0,0)")
        .call(yAxisLeft);
    // Add the line by appending an svg:path element with the data line we created above
    // do this AFTER the axes above so that the line is above the tick-lines
    var condition = graph.selectAll(".condition")
        .data(cond_data)
        .enter()
        .append("g")
        .attr("class", "condition")
    condition.append("svg:path")
        .attr("class", "spectrum line")
        .attr("d",
            function (d) {
                if (d.values == null) {
                    return "M0,0"
                };
                return line(d.values);
            })
        .style("stroke", function (d) {
            return "#FF0000";
        });
    //Create annotation for theoretical isotopic peaks
    graph.selectAll(".peptide_iso_peaks")
        .data(data['IsotopDistr'][peptide][z])
        .enter()
        .append("line")
        .attr("class", "iso line")
        .attr("x1", function(d) {return x(d[0])})
        .attr("x2", function(d) {return x(d[0])})
        .attr("y1", function(d) {return y(0) })
        .attr("y2", function(d) {return y(ymax)});
    graph.selectAll(".peptide_iso_height")
        .data(data['IsotopDistr'][peptide][z])
        .enter()
        .append("line")
        .attr("class", "iso height")
        .attr("x1", function(d) {return x(d[0])-5})
        .attr("x2", function(d) {return x(d[0])+5})
        .attr("y1", function(d) {return y(ymax * d[1]) })
        .attr("y2", function(d) {return y(ymax * d[1])});
    //Show overlapping isotopic peaks
    min_index = binarySearch(data['MassSortedIsotops']['Masses'],minx)-1;
    max_index = binarySearch(data['MassSortedIsotops']['Masses'],maxx);

    overlap_isotopes = data['MassSortedIsotops']['Masses'].slice(min_index,max_index);
    graph.selectAll(".peptide_iso_peaks_overlap")
        .data(overlap_isotopes)
        .enter()
        .append("line")
        .attr("class", "iso line")
        .attr("x1", function(d) {return x(d)})
        .attr("x2", function(d) {return x(d)})
        .attr("y1", function(d) {return y(0) })
        .attr("y2", function(d) {return y(ymax)})
        .style("stroke","grey");
}

function reset_letter_clicked() {
  residue_status = Array.apply(null, new Array(data['sequence'].length))
      .map(Number.prototype
          .valueOf, 0);
      var svg = d3.select("#logo").select("svg");
      var text = svg.selectAll("text.enter");
      text.style("fill", function (d, i) {
          if (residue_status[i] == 2) {
              return 'rgb(255,0,0)';
          } else if (residue_status[i] == 1) {
              return 'rgb(0,255,0)';
          } else {
              return 'rgb(0,0,0)';
          }
      });
      var lines = svg.selectAll(".peptide")
          .data(data['pepcoords'],function (d) { return d[3]});
      var no_ones = residue_status.reduce(function (a, b) {
          if (b == 1) {
              return a + b;
          } else {
              return a;
          }
      });
      lines.style("display", function (d, i) {
          var i_no_ones = 0;
          for (var k = d[0]; k < d[1]; k++) {
              if (residue_status[k] == 2) {

                  return "none";
              }
              if (residue_status[k] == 1) {
                  i_no_ones += 1;
              }
          }
          if (i_no_ones == no_ones) {

              return "block";
          } else {

              return "none";
          }
      });
      d3.select("#peptide_list").selectAll(".peptide-detail")
      .data(data['pepcoords'], function(d) { if (typeof(d) == "string") { return d;} else {return d[3];} })
      .style("display", function (d, i) {
          var i_no_ones = 0;
          for (var k = d[0]; k < d[1]; k++) {
              if (residue_status[k] == 2) {

                  return "none";
              }
              if (residue_status[k] == 1) {
                  i_no_ones += 1;
              }
          }
          if (i_no_ones == no_ones) {

              return "block";
          } else {

              return "none";
          }
      });
}

function letter_clicked(d, i) {
    residue_status[i] += 1;
    if (residue_status[i] == 3) {
        residue_status[i] = 0;
    }
    var svg = d3.select("#logo").select("svg");
    var text = svg.selectAll("text.enter");
    text.style("fill", function (d, i) {
        if (residue_status[i] == 2) {
            return 'rgb(255,0,0)';
        } else if (residue_status[i] == 1) {
            return 'rgb(0,255,0)';
        } else {
            return 'rgb(0,0,0)';
        }
    });
    var lines = svg.selectAll(".peptide")
        .data(data['pepcoords'],function (d) { return d[3]});
    var no_ones = residue_status.reduce(function (a, b) {
        if (b == 1) {
            return a + b;
        } else {
            return a;
        }
    });
    lines.style("display", function (d, i) {
        var i_no_ones = 0;
        for (var k = d[0]; k < d[1]; k++) {
            if (residue_status[k] == 2) {

                return "none";
            }
            if (residue_status[k] == 1) {
                i_no_ones += 1;
            }
        }
        if (i_no_ones == no_ones) {

            return "block";
        } else {

            return "none";
        }
    });
    d3.select("#peptide_list").selectAll(".peptide-detail")
    .data(data['pepcoords'], function(d) { if (typeof(d) == "string") { return d;} else {return d[3];} })
    .style("display", function (d, i) {
        var i_no_ones = 0;
        for (var k = d[0]; k < d[1]; k++) {
            if (residue_status[k] == 2) {

                return "none";
            }
            if (residue_status[k] == 1) {
                i_no_ones += 1;
            }
        }
        if (i_no_ones == no_ones) {

            return "block";
        } else {

            return "none";
        }
    });
}

function update_logo_sort(way) {
    if (way == 0) {
        sort_peptides_by_intensity();
    } else {
        sort_peptides_by_position();
    }
    height = d3.max(data['pepcoords'].map(function (currentValue) {
        return currentValue[2];
    })) * 4 + 36;
    var svg = d3.select("#logo")
        .select("svg")
        .attr("height", height);
    var lines = svg.selectAll("line")
        .data(data['pepcoords'], function (d) {
            return d[3];
        });
    lines.transition()
        .attr("x1", function (d) {
            return d[0] * 11;
        })
        .attr("x2", function (d) {
            return d[1] * 11 - 2;
        })
        .attr("y1", function (d) {
            return d[2] * 4 + 14;
        })
        .attr("y2", function (d) {
            return d[2] * 4 + 14;
        });
}

function render_logo() {
    // DATA JOIN
    // Join new data with old elements, if any.
    generate_peptide_positions();
    height = d3.max(data['pepcoords'].map(function (currentValue) {
        return currentValue[2];
    })) * 4 + 36;
    var svg = d3.select("#logo")
        .append("svg")
        .attr("width", "960px")
        .attr("height", height)
        .append("g")
        .attr("transform", "translate(9,9)");
    var text = svg.selectAll("text.seq")
        .data(data['sequence']);
    // ENTER
    var lines = svg.selectAll("line")
        .data(data['pepcoords'], function (d) {
            return d[3];
        });
    // Create new elements as needed.
    text.enter()
        .append("text")
        .attr("class", "enter")
        .attr("x", function (d, i) {
            var j = i;
            return j * 11;
        })
        .attr("y", function (d, i) {
            return 0;
        })
        .attr("dy", ".35em")
        .style("fill", function (d, i) {
            return 'rgb(0,0,0)';
        });
    text.text(function (d) {
        return d;
    })
        .on("click", letter_clicked);
    lines.enter()
        .append("line")
        .attr("class", "peptide")
        .attr("x1", function (d) {
            return d[0] * 11;
        })
        .attr("x2", function (d) {
            return d[1] * 11 - 2;
        })
        .attr("y1", function (d) {
            return d[2] * 4 + 14;
        })
        .attr("y2", function (d) {
            return d[2] * 4 + 14;
        })
        .attr("dy", ".45em")
        .on("click", function (d) {
            location.hash = "#" + d[3];
        });
    // EXIT
    // Remove old elements as needed.
    var colscale = d3.scale.log()
        .domain([1, 100, 10000, 100000, 1000000, 10000000, 100000000,
            1000000000
        ])
        .range(['rgb(255,255,204)', 'rgb(255,237,160)', 'rgb(254,217,118)',
            'rgb(254,178,76)', 'rgb(253,141,60)', 'rgb(252,78,42)',
            'rgb(227,26,28)', 'rgb(177,0,38)'
        ]);
    lines.style("stroke", function (d) {
        return colscale(data['IntIntensities'][d[3]])
    });
}
var residue_status;

function generate_peptide_positions() {
    var peptides = d3.keys(data['peptides']);
    data['pepcoords'] = peptides.map(function (peptide) {
        start = data['sequence'].indexOf(peptide);
        stop = start + peptide.length;
        return [start, stop, 0, peptide];
    });
    sort_peptides_by_position();
}

function sort_peptides_by_position() {
    var positions = data['pepcoords'].sort(function (a, b) {
        return b[1] - a[1];
    });
    positions = positions.sort(function (a, b) {
        return a[0] - b[0];
    });
    end_values = [];
    data['pepcoords'] = positions.map(function (pos) {
        for (var i = 0; i < end_values.length; i++) {
            if (pos[0] >= end_values[i]) {
                end_values[i] = pos[1];
                return ([pos[0], pos[1], i, pos[3]]);
            }
        }
        end_values.push(pos[1]);
        return ([pos[0], pos[1], end_values.length - 1, pos[3]]);
    });
}

function sort_peptides_by_intensity() {
    positions = data['pepcoords'].sort(function (a, b) {
        return data['IntIntensities'][b[3]] - data['IntIntensities'][a[3]];
    });
    end_values = [];
    data['pepcoords'] = positions.map(function (pos) {
        outer :
        for (var i = 0; i < end_values.length; i++)  {
          for (var j = 0; j < end_values[i].length; j++) {
            if (pos[0] < end_values[i][j][1] && end_values[i][j][0] < pos[1]) {
              continue outer;
            }
          }
          end_values[i].push([pos[0],pos[1]]);
          return [pos[0],pos[1],i,pos[3]]
        }
        end_values.push([]);
        end_values[end_values.length-1].push([pos[0],pos[1]]);
        return ([pos[0], pos[1], end_values.length - 1, pos[3]]);
    });
}

function make_plots() {
    clearInterval(intervalId);
    var peptide_divs = d3.select("#peptide_list").selectAll("div")
                                                 .data(d3.keys(data['peptides'],function(d) {
                                                   return d;
                                                 }))
                                                 .enter().append("div")
                                                 .attr("id",function (d) { return d;})
                                                 .attr("class", function (d) { return 'peptide-detail'});

    peptide_divs.insert("h3").text(function(d) {return d;});
    var peptide_table_rows = peptide_divs.insert("table").selectAll("tr")
                                         .data(function (d) {
                                           return d3.keys(data['peptides'][d])
                                                    .map( function (e) {
                                                       return([d,e])
                                                     })
                                           })
                                          .enter().append("tr");
    peptide_table_rows.insert("th").text(function (d) { return d[1];});
    peptide_table_rows.insert("td").insert("div")
      .attr("id", function (d) { return 'basepeak_' + d[0] + '_' + d[1];})
      .attr("data-peptide", function (d) { return d[0] })
      .attr("data-charge", function (d) { return d[1] })
      .attr("data-type", function (d) { return 'basepeak' })
      .attr("class", function (d) { return 'basepeak' })
      .style("width", function (d) { return '960px'})
      .style("height", function (d) { return '220px'})
      .each(initiate_draw_time_graph);
    peptide_table_rows.insert("td").insert("div")
      .attr("id", function (d) { return 'spectrum_' + d[0] + '_' + d[1];})
      .attr("data-peptide", function (d) { return d[0] })
      .attr("data-charge", function (d) { return d[1] })
      .attr("data-type", function (d) { return 'spectrum' })
      .attr("class", function (d) { return 'spectrum' })
      .style("width", function (d) { return '320px'})
      .style("height", function (d) { return '220px'})
      .each(initiate_draw_spectrum);

    residue_status = Array.apply(null, new Array(data['sequence'].length))
    .map(Number.prototype
        .valueOf, 0);
    //basepeakplots = d3.selectAll(".basepeak");
    //basepeakplots.each(initiate_draw_time_graph);
    //spectrumplots = d3.selectAll(".spectrum");
    //spectrumplots.each(initiate_draw_spectrum);
    render_logo();
}


var intervalId = setInterval(function () {
    if (typeof data != "undefined") { //once data is detected, process data
        make_plots();
    }
}, 1000)
