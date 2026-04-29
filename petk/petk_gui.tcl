# PETK GUI
package require Tk

package provide petk_gui 1.0

namespace eval ::PETK::gui {
    variable window     ".petk_main_window"
}

set ::PETK::gui::module_dir [file normalize [file dirname [info script]]]
foreach module {utils.tcl tab1_nanopore.tcl tab2_analyte.tcl tab3_sem.tcl} {
    set module_path [file join $::PETK::gui::module_dir $module]
    if {[file exists $module_path]} {
        source $module_path
    } else {
        puts "Warning: PETK module not found: $module_path"
    }
}

proc ::PETK::gui::resourcePath {args} {
    if {![info exists ::PETK::gui::module_dir]} {
        set ::PETK::gui::module_dir [file normalize [file dirname [info script]]]
    }
    if {[llength $args] == 0} {
        return $::PETK::gui::module_dir
    }
    return [file join $::PETK::gui::module_dir {*}$args]
}

namespace eval ::PETK::gui::plot {
    # Plot options - track selected series
    variable selectedSeries "current"
    variable plotCurrent 1
    variable plotBlockage 0
    
    # NEW: Open pore current reference line
    variable showOpenPoreLine 1
    variable openPoreCurrent ""  # Will be auto-calculated or user-set
    variable openPoreCurrentMode "auto"  # "auto" or "manual"

    variable fontsInitialized 0
    variable titleFontName "::PETK::PlotTitleFont"
    variable axisFontName "::PETK::PlotAxisFont"
    variable tickFontName "::PETK::PlotTickFont"
    variable legendFontName "::PETK::PlotLegendFont"
}

proc ::PETK::gui::plot::ensureFonts {} {
    variable fontsInitialized
    variable titleFontName
    variable axisFontName
    variable tickFontName
    variable legendFontName

    if {$fontsInitialized} {
        return
    }

    foreach {name opts} [list \
        $titleFontName  {-family Helvetica -size 18 -weight bold} \
        $axisFontName   {-family Helvetica -size 14 -weight normal} \
        $tickFontName   {-family Helvetica -size 12 -weight normal} \
        $legendFontName {-family Helvetica -size 12 -weight normal}] {
        catch {font delete $name}
        catch {font create $name {*}$opts}
    }
    set fontsInitialized 1
}

proc ::PETK::gui::plot::updateSeriesFlags {} {
    variable selectedSeries
    variable plotCurrent
    variable plotBlockage
    if {$selectedSeries eq "conductance"} {
        set plotCurrent 0
        set plotBlockage 1
    } else {
        set plotCurrent 1
        set plotBlockage 0
    }
}

proc ::PETK::gui::plot::setSeries {newval} {
    variable selectedSeries
    set selectedSeries $newval
    ::PETK::gui::plot::updateSeriesFlags
}

::PETK::gui::plot::updateSeriesFlags

# Core plotting functions (simplified version)
proc ::PETK::gui::plot::parseFile {filename open_current voltage_v} {
    if {![file exists $filename]} {
        return {}
    }
    
    set infile [open $filename "r"]
    set metadata {}
    
    while {[gets $infile line] >= 0} {
        set line [string trim $line]
        if {[string match "#*" $line]} {
            if {[regexp {# Position (\d+)/(\d+)} $line -> pos total]} {
                dict set metadata position $pos
                dict set metadata total_positions $total
            } elseif {[regexp {# Z_position: ([\d\.\-e\+]+)} $line -> zpos]} {
                dict set metadata z_position $zpos
            } elseif {[regexp {# Current: ([\d\.\-e\+]+)} $line -> current]} {
                dict set metadata current $current
            } elseif {[regexp {# Blockage: ([\d\.]+)%} $line -> blockage]} {
                dict set metadata blockage $blockage
            }
        } elseif {[llength [split $line]] >= 6} {
            set values [split $line]
            dict set metadata z_position [lindex $values 0]
            dict set metadata current [lindex $values 1]
            if {[llength $values] >= 3} {
                dict set metadata blockage [lindex $values 2]
            }
            break
        }
    }
    
    close $infile

    if {[dict exists $metadata current] && $open_current ne "" && $voltage_v ne ""} {
        if {[catch {expr {double($voltage_v)}}]} {
            # Ignore invalid voltage
        } else {
            set current_val [expr {double([dict get $metadata current])}]
            if {[expr {$voltage_v != 0.0}]} {
                # Currents are logged in nA; dividing the difference by volts yields nS
                set delta_g_ns [expr {(double($open_current) - $current_val) / double($voltage_v)}]
                dict set metadata delta_g $delta_g_ns
            }
        }
    }

    return $metadata
}

proc ::PETK::gui::plot::parseResultsSummary {results_file} {
    set summary [dict create]
    if {![file exists $results_file]} {
        return $summary
    }

    set fh [open $results_file r]
    set currents {}
    set z_positions {}
    set open_current ""

    while {[gets $fh line] >= 0} {
        set line [string trim $line]
        if {$line eq ""} {
            continue
        }
        if {[string match "#*" $line]} {
            if {[regexp {Open_pore_current:\s*([-0-9.eE\+]+)} $line -> val]} {
                if {![catch {expr {double($val)}} parsed]} {
                    set open_current $parsed
                }
            }
            continue
        }

        set parts [split $line]
        if {[llength $parts] >= 2} {
            if {![catch {expr {double([lindex $parts 0])}} zval]} {
                lappend z_positions $zval
            }
            if {![catch {expr {double([lindex $parts 1])}} cur]} {
                lappend currents $cur
            }
        }
    }

    close $fh

    if {$open_current ne ""} {
        dict set summary open_current $open_current
    }
    if {[llength $currents] > 0} {
        dict set summary currents $currents
    }
    if {[llength $z_positions] > 0} {
        dict set summary z_positions $z_positions
    }

    return $summary
}

proc ::PETK::gui::plot::parseResultsData {results_file} {
    set alldata {}
    if {$results_file eq "" || ![file exists $results_file]} {
        return $alldata
    }

    set fh [open $results_file r]
    set header_indices [dict create z 0 current 1 blockage -1]
    set header_seen 0

    while {[gets $fh line] >= 0} {
        set trimmed [string trim $line]
        if {$trimmed eq ""} {
            continue
        }
        if {[string match "#*" $trimmed]} {
            if {!$header_seen && [string match "*Z_position*" $trimmed]} {
                set header_seen 1
                set header_line [string trim [string trimleft $trimmed "#"]]
                set header_line [string trim [string trimleft $header_line "#"]]
                set headers [split $header_line]
                for {set idx 0} {$idx < [llength $headers]} {incr idx} {
                    set name [lindex $headers $idx]
                    if {[string match "*Z_position*" $name]} {
                        dict set header_indices z $idx
                    } elseif {[string match "*Blockage*" $name]} {
                        dict set header_indices blockage $idx
                    } elseif {[string match "*Current*" $name] && ![string match "*Normalized*" $name]} {
                        dict set header_indices current $idx
                    }
                }
            }
            continue
        }

        set fields [split $trimmed]
        set z_idx [dict get $header_indices z]
        set cur_idx [dict get $header_indices current]
        set block_idx [dict get $header_indices blockage]
        if {$z_idx >= [llength $fields] || $cur_idx >= [llength $fields]} {
            continue
        }

        set entry [dict create]
        dict set entry z_position [lindex $fields $z_idx]
        dict set entry current [lindex $fields $cur_idx]
        if {$block_idx >= 0 && $block_idx < [llength $fields]} {
            dict set entry blockage [lindex $fields $block_idx]
        }
        lappend alldata $entry
    }
    close $fh

    return $alldata
}

proc ::PETK::gui::plot::escapePythonString {text} {
    return [string map {"\\" "\\\\" "\"" "\\\""} $text]
}

proc ::PETK::gui::plot::formatPythonNumberList {values} {
    set formatted {}
    foreach val $values {
        if {[catch {expr {double($val)}} parsed]} {
            continue
        }
        lappend formatted [format "%.10g" $parsed]
    }
    return [join $formatted ", "]
}

proc ::PETK::gui::plot::renderSeriesWithMatplotlib {args} {
    array set opts {
        title ""
        xlabel "Z Position (Å)"
        ylabel ""
        dataDir ""
        seriesLabel ""
        lineColor "#C62828"
        openPore ""
        scriptName "plot_sem_series.py"
        outputName "sem_results_plot.png"
        zValues {}
        yValues {}
    }
    array set opts $args

    set z_list [::PETK::gui::plot::formatPythonNumberList $opts(zValues)]
    set y_list [::PETK::gui::plot::formatPythonNumberList $opts(yValues)]
    if {$z_list eq "" || $y_list eq ""} {
        tk_messageBox -icon warning -title "Plot Error" -message "No numeric data available to plot."
        return ""
    }

    set data_dir $opts(dataDir)
    if {$data_dir eq ""} {
        if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
            set data_dir $::PETK::gui::workdir
        } else {
            set data_dir [pwd]
        }
    }

    if {![file exists $data_dir]} {
        tk_messageBox -icon error -title "Plot Error" -message "Data directory $data_dir was not found."
        return ""
    }

    set python_exec "python3"
    if {[info exists ::PETK::gui::pythonExecutable] && $::PETK::gui::pythonExecutable ne ""} {
        set python_exec $::PETK::gui::pythonExecutable
    }

    set script_path [file join $data_dir $opts(scriptName)]
    set output_path [file join $data_dir $opts(outputName)]

    set open_pore_literal "None"
    if {$opts(openPore) ne ""} {
        if {![catch {expr {double($opts(openPore))}} parsed]} {
            set open_pore_literal [format "%.10g" $parsed]
        }
    }

    set python_code [format {import matplotlib.pyplot as plt

z_values = [%s]
y_values = [%s]
title = "%s"
xlabel = "%s"
ylabel = "%s"
series_label = "%s"
line_color = "%s"
output_path = r"%s"
open_pore = %s

if not z_values or not y_values:
    print("No data to plot.")
else:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_values, y_values, color=line_color, linewidth=2, marker='o', label=series_label)
    if open_pore is not None:
        ax.axhline(open_pore, color="gray", linestyle="--", linewidth=2,
                   label=f"Open Pore Current ({open_pore:.3f} nA)")
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.show()
    print(f"Plot saved to {output_path}")
} $z_list $y_list [::PETK::gui::plot::escapePythonString $opts(title)] [::PETK::gui::plot::escapePythonString $opts(xlabel)] [::PETK::gui::plot::escapePythonString $opts(ylabel)] [::PETK::gui::plot::escapePythonString $opts(seriesLabel)] [::PETK::gui::plot::escapePythonString $opts(lineColor)] [string map {"\\" "\\\\"} $output_path] $open_pore_literal]

    if {[catch {
        set fh [open $script_path w]
        puts $fh $python_code
        close $fh
    } err]} {
        tk_messageBox -icon error -title "Plot Error" -message "Failed to write Python plot script:\n$err"
        return ""
    }

    set plot_output ""
    if {[catch {set plot_output [exec $python_exec $script_path]} err]} {
        tk_messageBox -icon error -title "Plot Error" -message "Python plotting failed:\n$err"
    } else {
        puts $plot_output
        tk_messageBox -icon info -title "Plot Generated" -message "Matplotlib figure launched.\nImage saved to:\n$output_path"
    }

    catch {file delete -force $script_path}
    return $output_path
}

proc ::PETK::gui::plot::readVoltageFromConfig {config_path} {
    if {![file exists $config_path]} {
        return ""
    }
    set fh [open $config_path r]
    set content [read $fh]
    close $fh
    if {[regexp {"voltage"\s*:\s*([-0-9.eE\+]+)} $content -> val]} {
        if {[catch {expr {double($val)/1000.0}} result]} {
            return ""
        }
        return $result
    }
    return ""
}

proc ::PETK::gui::plot::readOutputPrefixFromConfig {config_path} {
    if {![file exists $config_path]} {
        return ""
    }
    set fh [open $config_path r]
    set content [read $fh]
    close $fh
    if {[regexp {"output_prefix"\s*:\s*"([^"]+)"} $content -> val]} {
        return $val
    }
    return ""
}

proc ::PETK::gui::plot::collectRunMetadata {directory prefix} {
    set metadata [dict create open_current "" voltage_v "" results_file "" config_path ""]

    set candidates {}
    if {$prefix ne ""} {
        set candidates [glob -nocomplain -directory $directory "${prefix}_results.txt"]
    }
    if {[llength $candidates] == 0} {
        set candidates [glob -nocomplain -directory $directory "*_results.txt"]
    }
    if {[llength $candidates] > 0} {
        set results_file [lindex $candidates 0]
        dict set metadata results_file $results_file
        set summary [::PETK::gui::plot::parseResultsSummary $results_file]
        if {[dict exists $summary open_current]} {
            dict set metadata open_current [dict get $summary open_current]
        }
    }

    set config_path [file join $directory "config.json"]
    if {![file exists $config_path]} {
        set parent_config [file join [file dirname $directory] "config.json"]
        if {[file exists $parent_config]} {
            set config_path $parent_config
        } else {
            set config_path ""
        }
    }
    if {$config_path ne ""} {
        dict set metadata config_path $config_path
        set voltage_v [::PETK::gui::plot::readVoltageFromConfig $config_path]
        if {$voltage_v ne ""} {
            dict set metadata voltage_v $voltage_v
        }
        if {![dict exists $metadata output_prefix]} {
            set out_prefix [::PETK::gui::plot::readOutputPrefixFromConfig $config_path]
            if {$out_prefix ne ""} {
                dict set metadata output_prefix $out_prefix
            }
        }
    }

    return $metadata
}

proc ::PETK::gui::plot::computeRotationDeltas {base_dir} {
    set entries {}
    set rotation_dirs [lsort -dictionary [glob -nocomplain -types d [file join $base_dir "rot_*"]]]

    foreach dir $rotation_dirs {
        set config_path [file join $dir "config.json"]
        set voltage_v [::PETK::gui::plot::readVoltageFromConfig $config_path]
        if {$voltage_v eq "" || $voltage_v == 0.0} {
            continue
        }

        set results_candidates [glob -nocomplain -directory $dir "*_results.txt"]
        if {[llength $results_candidates] == 0} {
            continue
        }
        set results_path [lindex $results_candidates 0]
        set summary [::PETK::gui::plot::parseResultsSummary $results_path]
        if {![dict exists $summary open_current]} {
            continue
        }
        set open_current [dict get $summary open_current]
        if {![dict exists $summary currents] || [llength [dict get $summary currents]] == 0} {
            continue
        }
        set currents [dict get $summary currents]
        set measured_current [lindex $currents end]
        set delta_ns [expr {(double($open_current) - double($measured_current)) / double($voltage_v) * 1.0e9}]
        set label [file tail $dir]
        set entry [dict create label $label delta $delta_ns results_file $results_path voltage_v $voltage_v open_current $open_current config_path $config_path]
        lappend entries $entry
    }

    return $entries
}

proc ::PETK::gui::plot::plotRotationHistogram {} {
    if {![info exists ::PETK::gui::workdir] || $::PETK::gui::workdir eq ""} {
        tk_messageBox -icon error -title "Plot Error" -message "No working directory set."
        return
    }

    set candidate_dirs {}
    set workdir [file normalize $::PETK::gui::workdir]

    set results_root [file join $workdir "results"]
    lappend candidate_dirs $results_root
    lappend candidate_dirs [file join $results_root "rotations"]
    lappend candidate_dirs [file join $workdir "rotation" "results"]
    lappend candidate_dirs [file join $workdir "rotation" "results" "rotations"]

    if {[info procs ::PETK::gui::resourcePath] ne ""} {
        set install_root [::PETK::gui::resourcePath rotation results]
        lappend candidate_dirs $install_root
        lappend candidate_dirs [file join $install_root "rotations"]
    }

    set norm_candidates {}
    foreach dir $candidate_dirs {
        if {$dir eq ""} {
            continue
        }
        if {[catch {set norm_dir [file normalize $dir]}]} {
            continue
        }
        if {[lsearch -exact $norm_candidates $norm_dir] == -1} {
            lappend norm_candidates $norm_dir
        }
    }

    set base_dir ""
    set entries {}
    foreach dir $norm_candidates {
        if {![file isdirectory $dir]} {
            continue
        }
        set dir_entries [::PETK::gui::plot::computeRotationDeltas $dir]
        if {[llength $dir_entries] > 0} {
            set base_dir $dir
            set entries $dir_entries
            break
        }
    }

    if {$base_dir eq ""} {
        if {[llength $norm_candidates] == 0} {
            set searched "none"
        } else {
            set searched [join $norm_candidates "\n  • "]
        }
        tk_messageBox -icon info -title "Plot Rotation Scan" \
            -message "No rotation scan data found in any of the expected locations:\n  • $searched\nRun a rotation scan first, or copy the results into one of these directories."
        return
    }

    set record_literals {}
    foreach entry $entries {
        if {![dict exists $entry results_file]} {
            continue
        }
        set label [dict get $entry label]
        set results_file [dict get $entry results_file]
        set config_path ""
        if {[dict exists $entry config_path]} {
            set config_path [dict get $entry config_path]
        }
        set label_literal [string map {"\\" "\\\\" "'" "\\'"} $label]
        set results_literal [string map {"\\" "\\\\" "\"" "\\\""} $results_file]
        set config_literal [string map {"\\" "\\\\" "\"" "\\\""} $config_path]
        lappend record_literals [format "    {'label': '%s', 'results': r\"%s\", 'config': r\"%s\"}" $label_literal $results_literal $config_literal]
    }

    if {[llength $record_literals] == 0} {
        tk_messageBox -icon info -title "Plot Rotation Scan" -message "No rotation scan files were found."
        return
    }

    set entries_py [join $record_literals ",\n"]
    set output_path [file join $base_dir "conductance_diff_histogram.png"]
    set output_path_py [string map {"\\" "\\\\"} $output_path]
    set plot_title [format "Rotation ΔG Histogram (%s)" [file tail $workdir]]
    set plot_title_escaped [string map {"\"" "\\\""} $plot_title]
    set python_script_path [file join $base_dir "plot_rotation_hist.py"]

    set python_code [format {import json
import matplotlib.pyplot as plt

entries = [
%s
]

output_path = r"%s"
plot_title = "%s"

def parse_voltage(config_path):
    if not config_path:
        return None
    try:
        with open(config_path, 'r') as fh:
            data = json.load(fh)
        voltage = None
        if isinstance(data, dict):
            for key in ('simulation', 'simulation_parameters'):
                if key in data and isinstance(data[key], dict):
                    block = data[key]
                    if 'voltage' in block:
                        voltage = block['voltage']
                        break
                    if 'voltage_mV' in block:
                        voltage = block['voltage_mV']
                        break
            if voltage is None and 'voltage' in data:
                voltage = data['voltage']
        if voltage is None:
            return None
        return float(voltage) / 1000.0
    except Exception as exc:
        print(f"Warning: could not read voltage from {config_path}: {exc}")
        return None

def parse_results(results_path):
    open_current = None
    measured = None
    try:
        with open(results_path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    if 'Open_pore_current' in line:
                        try:
                            open_current = float(line.split(':', 1)[1].split()[0])
                        except Exception:
                            pass
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        measured = float(parts[1])
                    except Exception:
                        pass
        return open_current, measured
    except Exception as exc:
        print(f"Warning: could not read {results_path}: {exc}")
        return None, None

deltas = []
labels = []
for entry in entries:
    label = entry['label']
    results_path = entry['results']
    config_path = entry.get('config') or ''
    open_pore, measured = parse_results(results_path)
    if open_pore is None or measured is None:
        print(f"{label}: Missing current data.")
        continue
    voltage_v = parse_voltage(config_path)
    if voltage_v is None or voltage_v == 0.0:
        print(f"{label}: Missing voltage information.")
        continue
    delta_g = abs((open_pore - measured) / voltage_v)
    deltas.append(delta_g)
    labels.append(label)
    print(f"{label}: I={measured:.6f} nA, ΔG={delta_g:.6f} nS (open pore {open_pore:.6f} nA, V={voltage_v:.6f} V)")

if not deltas:
    print("No data found.")
else:
    plt.figure(figsize=(8, 6))
    plt.hist(deltas, bins='auto', edgecolor='black')
    plt.xlabel(r"$\Delta G$ (nS)", fontsize=16, weight='bold')
    plt.ylabel('Count', fontsize=16, weight='bold')
    plt.title(plot_title, fontsize=16, weight='bold')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.show()
    print(f"Histogram saved as {output_path}")
} $entries_py $output_path_py $plot_title_escaped]

    if {[catch {
        set py_fh [open $python_script_path w]
        puts $py_fh $python_code
        close $py_fh
    } write_error]} {
        tk_messageBox -icon error -title "Plot Error" -message "Failed to write temporary Python script:\n$write_error"
        return
    }

    set python_exec "python3"
    if {[info exists ::PETK::gui::pythonExecutable] && $::PETK::gui::pythonExecutable ne ""} {
        set python_exec $::PETK::gui::pythonExecutable
    }

    set plot_output ""
    if {[catch {set plot_output [exec $python_exec $python_script_path]} err]} {
        tk_messageBox -icon error -title "Plot Error" -message "Python plotting failed:\n$err"
    } else {
        puts $plot_output
        tk_messageBox -icon info -title "Plot Rotation Scan" -message "Matplotlib histogram opened.\nSaved to:\n$output_path"
    }

    catch {file delete -force $python_script_path}
}

proc ::PETK::gui::plot::parseDirectory {directory prefix} {
    set alldata {}
    set pattern [file join $directory "${prefix}_*.dat"]
    set files [lsort -dictionary [glob -nocomplain $pattern]]

    set run_meta [::PETK::gui::plot::collectRunMetadata $directory $prefix]
    set open_current [dict get $run_meta open_current]
    set voltage_v [dict get $run_meta voltage_v]
    set results_file ""
    if {[dict exists $run_meta results_file]} {
        set results_file [dict get $run_meta results_file]
    }

    foreach filename $files {
        set filedata [::PETK::gui::plot::parseFile $filename $open_current $voltage_v]
        if {[llength $filedata] > 0} {
            if {$open_current ne ""} {
                dict set filedata open_current $open_current
            }
            if {$voltage_v ne ""} {
                dict set filedata voltage_v $voltage_v
            }
            lappend alldata $filedata
        }
    }
    if {[llength $alldata] == 0 && $results_file ne ""} {
        set table_data [::PETK::gui::plot::parseResultsData $results_file]
        foreach entry $table_data {
            if {$open_current ne ""} {
                dict set entry open_current $open_current
            }
            if {$voltage_v ne ""} {
                dict set entry voltage_v $voltage_v
            }
            if {$open_current ne "" && $voltage_v ne "" && $voltage_v != 0.0 && [dict exists $entry current]} {
                if {![catch {expr {double($voltage_v)}} safe_voltage]} {
                    set current_val [expr {double([dict get $entry current])}]
                    set delta_g_ns [expr {(double($open_current) - $current_val) / $safe_voltage}]
                    dict set entry delta_g $delta_g_ns
                }
            }
            lappend alldata $entry
        }
    }

    return $alldata
}

proc ::PETK::gui::plot::createPlot {alldata title {data_dir ""}} {
    variable plotCurrent
    variable plotBlockage
    variable showOpenPoreLine
    variable openPoreCurrent
    variable openPoreCurrentMode
    
    if {[llength $alldata] == 0} {
        tk_messageBox -icon info -title "Plot Data" -message "No SEM results were found to plot."
        return ""
    }
    
    set zpos_array {}
    set current_array {}
    set delta_g_array {}
    
    foreach entry $alldata {
        if {[dict exists $entry z_position]} {
            lappend zpos_array [dict get $entry z_position]
            lappend current_array [dict get $entry current]
            if {[dict exists $entry delta_g]} {
                lappend delta_g_array [dict get $entry delta_g]
            } elseif {[dict exists $entry blockage]} {
                # Backwards compatibility
                lappend delta_g_array [dict get $entry blockage]
            }
        }
    }
    
    if {[llength $zpos_array] == 0} {
        tk_messageBox -icon info -title "Plot Data" -message "Parsed files did not contain Z-position data."
        return ""
    }
    
    if {$openPoreCurrentMode eq "auto" && $plotCurrent && [llength $current_array] > 0} {
        set max_current [::tcl::mathfunc::max {*}$current_array]
        set openPoreCurrent $max_current
        puts "Auto-calculated open pore current: [format "%.6f" $openPoreCurrent] nA"
    }

    set ylabel ""
    set series_values {}
    set series_label ""
    set line_color "#C62828"
    set script_name "plot_sem_series.py"
    set output_name "sem_results_plot.png"

    if {$plotCurrent} {
        set ylabel "Current (nA)"
        set series_label "Current (nA)"
        set series_values $current_array
        set line_color "#C62828"
        set script_name "plot_sem_current.py"
        set output_name "sem_results_plot_current.png"
    } elseif {$plotBlockage} {
        set ylabel "Δ Conductance (nS)"
        set series_label "Δ Conductance (nS)"
        set series_values $delta_g_array
        set line_color "#1565C0"
        set script_name "plot_sem_conductance.py"
        set output_name "sem_results_plot_conductance.png"
    }

    if {[llength $series_values] == 0} {
        tk_messageBox -icon warning -title "Plot Data" -message "No data points found for the selected series."
        return ""
    }

    set open_line ""
    if {$plotCurrent && $showOpenPoreLine && $openPoreCurrent ne "" && [string is double -strict $openPoreCurrent]} {
        set open_line $openPoreCurrent
    }
    
    return [::PETK::gui::plot::renderSeriesWithMatplotlib \
        title $title \
        ylabel $ylabel \
        seriesLabel $series_label \
        lineColor $line_color \
        dataDir $data_dir \
        scriptName $script_name \
        outputName $output_name \
        openPore $open_line \
        zValues $zpos_array \
        yValues $series_values]
}

proc ::PETK::gui::petk_gui {} {
    variable window

    # If window already exists, bring it up
    if { [winfo exists $window] } {
        wm deiconify $window
        return
    }

    # Destroy and create window
    catch {destroy $window}
    toplevel $window

    grid columnconfigure $window 0 -weight 1
    grid rowconfigure    $window 0 -weight 1

    wm title      $window "Pore Explorer"
    wm resizable  $window 1 1
    wm geometry   $window 1000x800
    wm minsize    $window 800 600

    ### Menubar
    menu $window.menubar -tearoff 0
    $window configure -menu $window.menubar

    set menu $window.menubar
    menu $menu.file    -tearoff 0
    menu $menu.options -tearoff 0
    menu $menu.help    -tearoff 0

    $menu add cascade -menu $menu.file    -label "Pore Explorer"
    $menu add cascade -menu $menu.options -label "Options"
    $menu add cascade -menu $menu.help    -label "Help"

    if {![string eq "Darwin" $::tcl_platform(os)]} {
        $menu.file add separator
        $menu.file add command -label "Preferences..."      -command  ::PETK::gui::showPref -font TkMenuFont
    } else {
        $menu add cascade -menu $menu.apple
    }

    $menu.file add separator
    $menu.file add command -label "Close window"    -command { wm withdraw $::PETK::gui::window }
    $menu.file add command -label "Quit"            -command { wm withdraw $::PETK::gui::window }

    $menu.options add checkbutton -label "Show advanced settings" -variable ::PETK::gui::showadvanced

    ### MAIN frame holding notebook
    ttk::frame $window.hlf
    grid $window.hlf  -column 0 -row 0 -sticky nsew
    grid columnconfigure $window.hlf 0 -weight 1
    grid rowconfigure    $window.hlf 0 -weight 1  ; # Give weight to row 0 (notebook row)

    ### Notebook
    ttk::notebook $window.hlf.nb
    grid $window.hlf.nb -column 0 -row 0 -sticky nsew

    ### CRITICAL FIX: Configure notebook to expand its content
    grid columnconfigure $window.hlf.nb 0 -weight 1
    grid rowconfigure    $window.hlf.nb 0 -weight 1   ; # FIXED: Was 1, should be 0!

    $window.hlf.nb add [ttk::frame $window.hlf.nb.tab1] -text "Nanopore Setup"
    $window.hlf.nb add [ttk::frame $window.hlf.nb.tab2] -text "Analyte Setup"
    $window.hlf.nb add [ttk::frame $window.hlf.nb.tab3] -text "SEM Setup"
    ttk::notebook::enableTraversal $window.hlf.nb

    foreach kid [winfo children $window.hlf.nb] {
        set name [string range [file extension $kid] 1 end]
        set [set name] $kid
        grid columnconfigure $kid 0 -weight 1
        grid rowconfigure $kid 0 -weight 1   ; # ADDED: Give each tab weight too
        
        # Store as global variables for access from other procedures
        set ::PETK::gui::$name $kid
    }

    ### Navigation across tabs
    ttk::frame  $window.hlf.nav
    ttk::button $window.hlf.nav.back -text "Previous" -command [list ttk::notebook::CycleTab $window.hlf.nb -1]
    ttk::button $window.hlf.nav.next -text "Next"     -command [list ttk::notebook::CycleTab $window.hlf.nb  1]
    grid $window.hlf.nav -column 0 -row 1 -sticky se -padx "0 15" -pady "0 10"
    grid $window.hlf.nav.back $window.hlf.nav.next
    
    ####################################################
    # Tab 1: Nanopore Setup
    ####################################################
    buildTab1 $tab1

    ####################################################
    # Tab 2: Analyte Setup
    ####################################################
    buildTab2 $tab2

    ####################################################
    # Tab 3: SEM Setup
    ####################################################
    buildTab3 $tab3

    return $window
}

####################################################
# Tab 1: Nanopore Setup
####################################################
# Implementation sourced from petk/tab1_nanopore.tcl

####################################################
# Tab 2: Analyte Setup
####################################################
# Implementation sourced from petk/tab2_analyte.tcl

####################################################
# Tab 3 SEM Setup
####################################################
# Implementation sourced from petk/tab3_sem.tcl

####################################################
# Scroll Region Update Function (from Tab 2)
####################################################

proc ::PETK::gui::onCanvasConfigured {canvas canvas_window {container ""}} {
    # Get canvas dimensions
    set canvas_width [winfo width $canvas]
    set canvas_height [winfo height $canvas]
    
    # Don't process if canvas isn't properly sized yet
    if {$canvas_width <= 1 || $canvas_height <= 1} {
        return
    }
    
    # Resolve content container
    if {$container eq ""} {
        if {![info exists ::PETK::gui::contentContainer]} {
            return
        }
        set container $::PETK::gui::contentContainer
    }
    
    if {![winfo exists $container]} {
        return
    }
    
    # Update content width to match canvas
    set content_width [winfo reqwidth $container]
    if {$canvas_width > 1} {
        set content_width $canvas_width
    }
    
    # Configure the canvas window
    $canvas itemconfig $canvas_window -width $content_width
    
    # Update scroll region
    update idletasks
    set content_height [winfo reqheight $container]
    $canvas configure -scrollregion [list 0 0 $content_width $content_height]
}

proc ::PETK::gui::onContentConfigured {canvas canvas_window {container ""}} {
    if {$container eq ""} {
        if {![info exists ::PETK::gui::contentContainer]} {
            return
        }
        set container $::PETK::gui::contentContainer
    }
    
    if {![winfo exists $container]} {
        return
    }
    
    # Force update and recalculate
    update idletasks
    set content_width [winfo reqwidth $container]
    set content_height [winfo reqheight $container]
    
    # Update canvas scroll region
    $canvas configure -scrollregion [list 0 0 $content_width $content_height]
}

proc ::PETK::gui::resizeCanvasContent {canvas canvas_window} {
    # Get the content frame widget
    set content_frame [lindex [$canvas itemcget $canvas_window -window] 0]
    
    if {![winfo exists $content_frame]} {
        return
    }
    
    update idletasks
    
    # Get current dimensions
    set canvas_width [winfo width $canvas]
    set content_req_width [winfo reqwidth $content_frame]
    set content_req_height [winfo reqheight $content_frame]
    
    # Skip if canvas not sized yet
    if {$canvas_width <= 1} {
        return
    }
    
    # Make content frame fill the canvas width (minus padding for scrollbar)
    set target_width [expr {$canvas_width - 5}]
    if {$target_width < $content_req_width} {
        set target_width $content_req_width
    }
    
    # Configure the canvas window item to use the calculated width
    $canvas itemconfig $canvas_window -width $target_width
    
    # Update the scroll region
    $canvas configure -scrollregion [list 0 0 $target_width $content_req_height]
}

proc ::PETK::gui::forceInitialCanvasResize {canvas canvas_window {container ""}} {
    # Force a complete layout recalculation
    update idletasks
    
    if {$container eq ""} {
        if {![info exists ::PETK::gui::contentContainer]} {
            return
        }
        set container $::PETK::gui::contentContainer
    }
    
    if {![winfo exists $container]} {
        return
    }
    
    set canvas_width [winfo width $canvas]
    set canvas_height [winfo height $canvas]
    set content_width [winfo reqwidth $container]
    set content_height [winfo reqheight $container]
    
    # Ensure content width matches canvas width
    if {$canvas_width > 1} {
        set content_width $canvas_width
    }
    
    # Configure canvas window and scroll region
    $canvas itemconfig $canvas_window -width $content_width
    $canvas configure -scrollregion [list 0 0 $content_width $content_height]
}

# === NOTEBOOK HIERARCHY FIX ===
proc ::PETK::gui::fixNotebookExpansion {} {
    # Fix the entire widget hierarchy for proper expansion
    
    # Main window
    catch {grid rowconfigure .petk_main_window 0 -weight 1}
    catch {grid columnconfigure .petk_main_window 0 -weight 1}
    
    # HLF frame 
    catch {grid rowconfigure .petk_main_window.hlf 0 -weight 1}
    catch {grid columnconfigure .petk_main_window.hlf 0 -weight 1}
    
    # Notebook
    set nb .petk_main_window.hlf.nb
    catch {
        grid $nb -sticky nsew
        grid rowconfigure $nb 0 -weight 1
        grid columnconfigure $nb 0 -weight 1
    }
    
    puts "Applied notebook expansion fix"
}

# === CONVENIENCE FUNCTION ===
proc ::PETK::gui::updateScrollRegion {canvas canvas_window} {
    # This maintains compatibility with your existing code
    ::PETK::gui::resizeCanvasContent $canvas $canvas_window
}
