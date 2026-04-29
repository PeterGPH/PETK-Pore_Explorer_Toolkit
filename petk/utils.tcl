# utils.tcl
# Common utility procedures and variable initialization

namespace eval ::PETK::gui::utils {
    variable showadvanced 0
    
    proc create_label_entry {parent lbl_text var width args} {
        ttk::label $parent.lbl -text $lbl_text -width 18
        ttk::entry $parent.entry -textvariable $var -width $width -justify center
        grid $parent.lbl $parent.entry -sticky ew -pady 3 {*}$args
        return $parent.entry
    }
    
    proc initializeMembraneVariables {} {
        # Initialize pore type
        if {![info exists ::PETK::gui::membraneType]} {
            set ::PETK::gui::membraneType "cylindrical"
        }
       
        # Initialize cylindrical parameters
        if {![info exists ::PETK::gui::cylindricalDiameter]} {
            set ::PETK::gui::cylindricalDiameter "200.0"
        }
        if {![info exists ::PETK::gui::cornerRadius]} {
            set ::PETK::gui::cornerRadius "50.0"
        }
       
        # Initialize double cone parameters
        if {![info exists ::PETK::gui::innerDiameter]} {
            set ::PETK::gui::innerDiameter "100.0"
        }
        if {![info exists ::PETK::gui::outerDiameter]} {
            set ::PETK::gui::outerDiameter "300.0"
        }
       
        # Initialize common parameters
        if {![info exists ::PETK::gui::nanoporeThickness]} {
            set ::PETK::gui::nanoporeThickness "200.0"
        }
       
        # Initialize box dimensions
        if {![info exists ::PETK::gui::autoCalculateBoxDimensions]} {
            set ::PETK::gui::autoCalculateBoxDimensions 1
        }
        if {![info exists ::PETK::gui::boxSizeX]} {
            set ::PETK::gui::boxSizeX "300.0"
        }
        if {![info exists ::PETK::gui::boxSizeY]} {
            set ::PETK::gui::boxSizeY "300.0"
        }
        if {![info exists ::PETK::gui::boxSizeZ]} {
            set ::PETK::gui::boxSizeZ "300.0"
        }
        if {![info exists ::PETK::gui::sysPadding]} {
            set ::PETK::gui::sysPadding "50.0"
        }
       
        # Initialize movement range for auto-calculation
        if {![info exists ::PETK::gui::zStartRange]} {
            set ::PETK::gui::zStartRange "150.0"
        }
        if {![info exists ::PETK::gui::zEndRange]} {
            set ::PETK::gui::zEndRange "-150.0"
        }
        if {![info exists ::PETK::gui::zStep]} {
            set ::PETK::gui::zStep "10.0"
        }
        if {![info exists ::PETK::gui::semPreviewFrames]} {
            set ::PETK::gui::semPreviewFrames 5
        }
       
        # Initialize auto-calculated box display
        if {![info exists ::PETK::gui::autoBoxX]} {
            set ::PETK::gui::autoBoxX ""
        }
        if {![info exists ::PETK::gui::autoBoxY]} {
            set ::PETK::gui::autoBoxY ""
        }
        if {![info exists ::PETK::gui::autoBoxZ]} {
            set ::PETK::gui::autoBoxZ ""
        }
       
        # Initialize output settings
        if {![info exists ::PETK::gui::workdir]} {
            if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
                set ::PETK::gui::workdir $::PETK::gui::workdir
            } else {
                set ::PETK::gui::workdir [pwd]
            }
        }
        if {![info exists ::PETK::gui::outputPrefix]} {
            set ::PETK::gui::outputPrefix "vertical_movement"
        }
       
        # Initialize display variables
        if {![info exists ::PETK::gui::currentPoreType]} {
            set ::PETK::gui::currentPoreType "Cylindrical"
        }
        if {![info exists ::PETK::gui::displayThickness]} {
            set ::PETK::gui::displayThickness ""
        }
        if {![info exists ::PETK::gui::displayParam1]} {
            set ::PETK::gui::displayParam1 ""
        }
        if {![info exists ::PETK::gui::displayParam2]} {
            set ::PETK::gui::displayParam2 ""
        }
        if {![info exists ::PETK::gui::calculatedVolume]} {
            set ::PETK::gui::calculatedVolume ""
        }
        if {![info exists ::PETK::gui::estimatedGridPoints]} {
            set ::PETK::gui::estimatedGridPoints ""
        }
       
        # Update display
        ::PETK::gui::utils::updateParameterDisplay
        ::PETK::gui::utils::calculateBoxDimensions
    }
    
    proc selectWorkDir {} {
        set tempdir [tk_chooseDirectory -title "Select project folder" -initialdir [pwd]]
        if {![string eq $tempdir ""]} {
            set ::PETK::gui::workdir $tempdir
        }
    }
    
    proc toggleBoxDimensionMode {} {
        set container $::PETK::gui::contentContainer
       
        # Hide both frames first
        catch {grid forget $container.boxdim.manual}
        catch {grid forget $container.boxdim.auto_display}
       
        if {$::PETK::gui::autoCalculateBoxDimensions} {
            # Show auto-calculated display
            grid $container.boxdim.auto_display -row 1 -column 0 -columnspan 6 -sticky ew
            ::PETK::gui::utils::calculateBoxDimensions
            ::PETK::gui::tab1::updateMembraneStatus "Switched to auto-calculated box dimensions"
        } else {
            # Show manual input fields
            grid $container.boxdim.manual -row 1 -column 0 -columnspan 6 -sticky ew
            ::PETK::gui::tab1::updateMembraneStatus "Switched to manual box dimensions"
        }
    }
    
    proc calculateBoxDimensions {} {
        if {!$::PETK::gui::autoCalculateBoxDimensions} {
            return
        }

        # Get pore radius based on type
        set pore_radius 0.0
        if {$::PETK::gui::membraneType eq "cylindrical"} {
            if {[string is double $::PETK::gui::cylindricalDiameter]} {
                set pore_radius [expr {$::PETK::gui::cylindricalDiameter / 2.0}]
            }
        } elseif {$::PETK::gui::membraneType eq "doublecone"} {
            if {[string is double $::PETK::gui::outerDiameter]} {
                set pore_radius [expr {$::PETK::gui::outerDiameter / 2.0}]
            }
        }
       
        # Calculate XY dimensions based on pore size and padding
        set padding $::PETK::gui::sysPadding
        set xy_size [expr {max(0.0, $pore_radius + 50.0 + $padding)}]
       
        # Calculate Z dimensions based on movement range and membrane
        set membrane_half_thickness 0.0
        if {[string is double $::PETK::gui::nanoporeThickness]} {
            set membrane_half_thickness [expr {$::PETK::gui::nanoporeThickness / 2.0}]
        }
       
        set z_start 0.0
        set z_end 0.0
        if {[string is double $::PETK::gui::zStartRange]} {
            set z_start $::PETK::gui::zStartRange
        }
        if {[string is double $::PETK::gui::zEndRange]} {
            set z_end $::PETK::gui::zEndRange
        }
       
        set z_min [expr {min(-$membrane_half_thickness, $z_end)}]
        set z_max [expr {max($membrane_half_thickness, $z_start)}]
       
        # Update display variables
        set ::PETK::gui::autoBoxX [format "%.1f to %.1f (%.1f)" [expr {-$xy_size}] $xy_size [expr {2*$xy_size}]]
        set ::PETK::gui::autoBoxY [format "%.1f to %.1f (%.1f)" [expr {-$xy_size}] $xy_size [expr {2*$xy_size}]]
        set ::PETK::gui::autoBoxZ [format "%.1f to %.1f (%.1f)" $z_min $z_max [expr {$z_max - $z_min + $padding}]]
       
        # Store actual values for calculations
        set ::PETK::gui::calculatedBoxSizeX [expr {2*$xy_size}]
        set ::PETK::gui::calculatedBoxSizeY [expr {2*$xy_size}]
        set ::PETK::gui::calculatedBoxSizeZ [expr {$z_max - $z_min}]
           
        puts "Auto-calculated box dimensions:"
        puts " X: -$xy_size to $xy_size Å"
        puts " Y: -$xy_size to $xy_size Å"
        puts " Z: $z_min to $z_max Å"
        puts " Pore radius: $pore_radius Å"
        puts " Movement range: $z_end to $z_start Å"
    }
    
    proc updateParameterDisplay {} {
        set ::PETK::gui::displayThickness "$::PETK::gui::nanoporeThickness Å"
       
        if {$::PETK::gui::membraneType eq "cylindrical"} {
            set ::PETK::gui::displayParam1 "Diameter: $::PETK::gui::cylindricalDiameter Å"
            set ::PETK::gui::displayParam2 "Corner R: $::PETK::gui::cornerRadius Å"
           
            # Calculate cylindrical volume
            if {[string is double $::PETK::gui::cylindricalDiameter] && [string is double $::PETK::gui::nanoporeThickness]} {
                set radius [expr {$::PETK::gui::cylindricalDiameter / 2.0}]
                set height $::PETK::gui::nanoporeThickness
                set volume [expr {3.14159 * $radius * $radius * $height}]
                set ::PETK::gui::calculatedVolume [format "%.1f Å³" $volume]
            } else {
                set ::PETK::gui::calculatedVolume "Invalid params"
            }
           
        } elseif {$::PETK::gui::membraneType eq "doublecone"} {
            set ::PETK::gui::displayParam1 "Inner: $::PETK::gui::innerDiameter Å"
            set ::PETK::gui::displayParam2 "Outer: $::PETK::gui::outerDiameter Å"
           
            # Calculate double cone volume (approximation)
            if {[string is double $::PETK::gui::innerDiameter] && [string is double $::PETK::gui::outerDiameter] && [string is double $::PETK::gui::nanoporeThickness]} {
                set r1 [expr {$::PETK::gui::innerDiameter / 2.0}]
                set r2 [expr {$::PETK::gui::outerDiameter / 2.0}]
                set h [expr {$::PETK::gui::nanoporeThickness / 2.0}]
                set volume [expr {3.14159/3.0 * $h * ($r1*$r1 + $r1*$r2 + $r2*$r2) * 2}]
                set ::PETK::gui::calculatedVolume [format "%.1f Å³" $volume]
            } else {
                set ::PETK::gui::calculatedVolume "Invalid params"
            }
        }
       
        ::PETK::gui::utils::calculateBoxDimensions
    }
    
    proc loadPoreImages {} {
        set ::PETK::gui::poreImages(cylindrical) ""
        set ::PETK::gui::poreImages(doublecone) ""
       
        set cyl_path [::PETK::gui::resourcePath shapes shape2.gif]
        if {[file exists $cyl_path]} {
            if {[catch {
                set ::PETK::gui::poreImages(cylindrical) [image create photo -file $cyl_path]
            } error]} {
                puts "Warning: Could not load cylindrical pore image: $error"
            }
        } else {
            puts "Warning: Cylindrical pore image not found at: $cyl_path"
        }
       
        set cone_path [::PETK::gui::resourcePath shapes shape1.gif]
        if {[file exists $cone_path]} {
            if {[catch {
                set ::PETK::gui::poreImages(doublecone) [image create photo -file $cone_path]
            } error]} {
                puts "Warning: Could not load double cone pore image: $error"
            }
        } else {
            puts "Warning: Double cone pore image not found at: $cone_path"
        }
    }
    
    proc drawPoreSchematic {canvas} {
        set width [winfo reqwidth $canvas]
        set height [winfo reqheight $canvas]
        if {$width <= 1} {set width 300}
        if {$height <= 1} {set height 250}
       
        set cx [expr {$width / 2}]
        set cy [expr {$height / 2}]
       
        if {$::PETK::gui::membraneType eq "cylindrical"} {
            # Draw cylindrical pore schematic
            # Membrane (gray rectangles)
            $canvas create rectangle 20 [expr {$cy-40}] [expr {$cx-30}] [expr {$cy+40}] -fill gray70 -outline black
            $canvas create rectangle [expr {$cx+30}] [expr {$cy-40}] [expr {$width-20}] [expr {$cy+40}] -fill gray70 -outline black
           
            # Pore (white rectangle)
            $canvas create rectangle [expr {$cx-30}] [expr {$cy-40}] [expr {$cx+30}] [expr {$cy+40}] -fill white -outline black
           
            # Labels
            $canvas create text $cx [expr {$cy-60}] -text "Cylindrical Pore" -font {TkDefaultFont 10 bold}
            $canvas create text 60 [expr {$cy-60}] -text "Membrane" -font {TkDefaultFont 9}
           
        } elseif {$::PETK::gui::membraneType eq "doublecone"} {
            # Draw double cone pore schematic
            # Membrane (gray rectangles)
            $canvas create rectangle 20 [expr {$cy-40}] [expr {$cx-40}] [expr {$cy+40}] -fill gray70 -outline black
            $canvas create rectangle [expr {$cx+40}] [expr {$cy-40}] [expr {$width-20}] [expr {$cy+40}] -fill gray70 -outline black
           
            # Double cone pore (trapezoid)
            set coords [list [expr {$cx-40}] [expr {$cy-40}] \
                             [expr {$cx-15}] $cy \
                             [expr {$cx-40}] [expr {$cy+40}] \
                             [expr {$cx+40}] [expr {$cy+40}] \
                             [expr {$cx+15}] $cy \
                             [expr {$cx+40}] [expr {$cy-40}]]
            $canvas create polygon $coords -fill white -outline black
           
            # Labels
            $canvas create text $cx [expr {$cy-60}] -text "Double Cone Pore" -font {TkDefaultFont 10 bold}
            $canvas create text 60 [expr {$cy-60}] -text "Membrane" -font {TkDefaultFont 9}
        }
    }
    
    proc addParameterLabels {canvas} {
        set width [winfo reqwidth $canvas]
        set height [winfo reqheight $canvas]
        if {$width <= 1} {set width 300}
        if {$height <= 1} {set height 250}
       
        # Add parameter annotations on the visualization
        if {$::PETK::gui::membraneType eq "cylindrical"} {
            $canvas create text [expr {$width-10}] 20 -text "D = $::PETK::gui::cylindricalDiameter Å" -anchor ne -font {TkDefaultFont 9} -fill blue
            $canvas create text [expr {$width-10}] 35 -text "R = $::PETK::gui::cornerRadius Å" -anchor ne -font {TkDefaultFont 9} -fill blue
            $canvas create text [expr {$width-10}] 50 -text "T = $::PETK::gui::nanoporeThickness Å" -anchor ne -font {TkDefaultFont 9} -fill blue
        } elseif {$::PETK::gui::membraneType eq "doublecone"} {
            $canvas create text [expr {$width-10}] 20 -text "Inner = $::PETK::gui::innerDiameter Å" -anchor ne -font {TkDefaultFont 9} -fill blue
            $canvas create text [expr {$width-10}] 35 -text "Outer = $::PETK::gui::outerDiameter Å" -anchor ne -font {TkDefaultFont 9} -fill blue
            $canvas create text [expr {$width-10}] 50 -text "T = $::PETK::gui::nanoporeThickness Å" -anchor ne -font {TkDefaultFont 9} -fill blue
        }
    }
    
    proc updateMembraneStatus {message} {
        if {[info exists ::PETK::gui::membraneStatusLabel] && [winfo exists ::PETK::gui::membraneStatusLabel]} {
            $::PETK::gui::membraneStatusLabel configure -text $message
            update
        }
        puts "PETK Membrane Status: $message"
    }
    
    proc onCanvasConfigured {canvas canvas_window} {
        set canvas_width [winfo width $canvas]
        set canvas_height [winfo height $canvas]
        if {$canvas_width <= 1 || $canvas_height <= 1} {
            return
        }
        if {[info exists ::PETK::gui::contentContainer]} {
            set container $::PETK::gui::contentContainer
            set content_width [winfo reqwidth $container]
            if {$content_width < $canvas_width} {
                set content_width $canvas_width
            }
            $canvas itemconfig $canvas_window -width $content_width
            update idletasks
            set content_height [winfo reqheight $container]
            $canvas configure -scrollregion [list 0 0 $content_width $content_height]
        }
    }
    
    proc onContentConfigured {canvas canvas_window} {
        if {[info exists ::PETK::gui::contentContainer]} {
            set container $::PETK::gui::contentContainer
            update idletasks
            set content_width [winfo reqwidth $container]
            set content_height [winfo reqheight $container]
            $canvas configure -scrollregion [list 0 0 $content_width $content_height]
        }
    }
    
    proc forceInitialCanvasResize {canvas canvas_window} {
        update idletasks
        if {[info exists ::PETK::gui::contentContainer]} {
            set container $::PETK::gui::contentContainer
            set canvas_width [winfo width $canvas]
            set canvas_height [winfo height $canvas]
            set content_width [winfo reqwidth $container]
            set content_height [winfo reqheight $container]
            if {$content_width < $canvas_width} {
                set content_width $canvas_width
            }
            $canvas itemconfig $canvas_window -width $content_width
            $canvas configure -scrollregion [list 0 0 $content_width $content_height]
        }
    }
    
    # Add other utility procs if needed for tab1, like transitionToPoreOption, but since they are tab1 specific, they go in tab1_nanopore.tcl
}
