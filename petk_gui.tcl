# PETK GUI
package require Tk

namespace eval ::PETK::gui {
    variable window     ".petk_main_window"
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
# Tab 2: Membrane and Nanopore Builder
####################################################

proc ::PETK::gui::buildTab1 {tab1} {
    ::PETK::gui::initializeMembraneVariables

    # === CRITICAL: Configure tab to expand ===
    grid columnconfigure $tab1 0 -weight 1
    grid rowconfigure $tab1 0 -weight 1

    # Create scrollable canvas container
    canvas $tab1.canvas -highlightthickness 0
    ttk::scrollbar $tab1.vscroll -orient vertical -command [list $tab1.canvas yview]
    ttk::scrollbar $tab1.hscroll -orient horizontal -command [list $tab1.canvas xview]
    
    # Configure canvas scrolling
    $tab1.canvas configure -yscrollcommand [list $tab1.vscroll set]
    $tab1.canvas configure -xscrollcommand [list $tab1.hscroll set]
    
    # Create the actual content frame inside the canvas
    ttk::frame $tab1.canvas.content
    set canvas_window [$tab1.canvas create window 0 0 -anchor nw -window $tab1.canvas.content]
    
    # Grid the canvas and scrollbars with proper expansion
    grid $tab1.canvas -row 0 -column 0 -sticky nsew
    grid $tab1.vscroll -row 0 -column 1 -sticky ns
    grid $tab1.hscroll -row 1 -column 0 -sticky ew
    
    # Configure grid weights - CRITICAL for expansion
    grid rowconfigure $tab1 0 -weight 1
    grid columnconfigure $tab1 0 -weight 1

    # Now use the content frame as your container
    set container $tab1.canvas.content
    grid columnconfigure $container 0 -weight 1
    set row 0
    
    # === OUTPUT SETTINGS SECTION ===
    ttk::labelframe $container.output -text "Output Settings" -padding 10
    grid $container.output -row $row -column 0 -sticky ew -padx 10 -pady 5
    grid columnconfigure $container.output {1 2} -weight 1
    incr row

    # Output directory
    ttk::label $container.output.dirlbl -text "Output directory:" -width 18
    ttk::entry $container.output.dir -textvariable ::PETK::gui::workdir -width 35
    ttk::button $container.output.browsedir -text "Browse..." -command {::PETK::gui::selectWorkDir}

    # Output prefix
    ttk::label $container.output.prefixlbl -text "Output prefix:" -width 18
    ttk::entry $container.output.prefix -textvariable ::PETK::gui::outputPrefix -width 25

    grid $container.output.dirlbl $container.output.dir $container.output.browsedir -sticky ew -pady 3
    grid $container.output.prefixlbl $container.output.prefix - -sticky ew -pady 3
    
    # === BOX DIMENSIONS SECTION ===
    ttk::labelframe $container.boxdim -text "Box Dimensions" -padding 10
    grid $container.boxdim -row $row -column 0 -sticky ew -padx 10 -pady 5
    grid columnconfigure $container.boxdim {1 3 5} -weight 1
    incr row

    # Auto-calculate checkbox
    ttk::checkbutton $container.boxdim.auto -text "Auto-calculate box dimensions" \
        -variable ::PETK::gui::autoCalculateBoxDimensions -command {::PETK::gui::toggleBoxDimensionMode}
    ttk::label $container.boxdim.cutofflbl -text "Distance cutoff (Å):" -width 18
    ttk::entry $container.boxdim.cutoff -textvariable ::PETK::gui::sysCutoff -width 12 -justify center

    grid $container.boxdim.auto $container.boxdim.cutofflbl $container.boxdim.cutoff -columnspan 6 -sticky w -pady "0 5"

    # Manual box dimension inputs
    ttk::frame $container.boxdim.manual
    ttk::label $container.boxdim.manual.xlbl -text "X size (Å):" -width 12
    ttk::entry $container.boxdim.manual.x -textvariable ::PETK::gui::boxSizeX -width 12 -justify center
    ttk::label $container.boxdim.manual.ylbl -text "Y size (Å):" -width 12
    ttk::entry $container.boxdim.manual.y -textvariable ::PETK::gui::boxSizeY -width 12 -justify center
    ttk::label $container.boxdim.manual.zlbl -text "Z size (Å):" -width 12
    ttk::entry $container.boxdim.manual.z -textvariable ::PETK::gui::boxSizeZ -width 12 -justify center

    grid $container.boxdim.manual.xlbl $container.boxdim.manual.x $container.boxdim.manual.ylbl $container.boxdim.manual.y $container.boxdim.manual.zlbl $container.boxdim.manual.z -sticky ew -pady 2

    # Auto-calculated dimension display
    ttk::frame $container.boxdim.auto_display
    ttk::label $container.boxdim.auto_display.title -text "Auto-calculated dimensions:" -font {TkDefaultFont 10 bold}
    grid $container.boxdim.auto_display.title -columnspan 6 -sticky w -pady "0 5"

    ttk::label $container.boxdim.auto_display.xlbl -text "X range:" -width 12
    ttk::label $container.boxdim.auto_display.xval -textvariable ::PETK::gui::autoBoxX -width 25 -anchor w -relief sunken -background white

    ttk::label $container.boxdim.auto_display.ylbl -text "Y range:" -width 12
    ttk::label $container.boxdim.auto_display.yval -textvariable ::PETK::gui::autoBoxY -width 25 -anchor w -relief sunken -background white

    ttk::label $container.boxdim.auto_display.zlbl -text "Z range:" -width 12
    ttk::label $container.boxdim.auto_display.zval -textvariable ::PETK::gui::autoBoxZ -width 25 -anchor w -relief sunken -background white

    grid $container.boxdim.auto_display.xlbl $container.boxdim.auto_display.xval - - - - -sticky ew -pady 2
    grid $container.boxdim.auto_display.ylbl $container.boxdim.auto_display.yval - - - - -sticky ew -pady 2
    grid $container.boxdim.auto_display.zlbl $container.boxdim.auto_display.zval - - - - -sticky ew -pady 2

    ttk::button $container.boxdim.auto_display.calc -text "Recalculate Box" -command {::PETK::gui::calculateBoxDimensions}
    grid $container.boxdim.auto_display.calc -row 3 -column 10 -columnspan 2 -sticky w -pady "5 0"

    # Grid the appropriate frame based on mode
    if {$::PETK::gui::autoCalculateBoxDimensions} {
        grid $container.boxdim.auto_display -row 1 -column 0 -columnspan 6 -sticky ew
    } else {
        grid $container.boxdim.manual -row 1 -column 0 -columnspan 6 -sticky ew
    }

    # === PORE TYPE SELECTION ===
    ttk::labelframe $container.poretype -text "Nanopore Type Selection" -padding 10
    grid $container.poretype -row $row -column 0 -sticky ew -padx 10 -pady "10 5"
    grid columnconfigure $container.poretype {0 1} -weight 1
    incr row

    ttk::radiobutton $container.poretype.cylindrical -text "Cylindrical Pore" -value "cylindrical" \
        -variable ::PETK::gui::membraneType -command {::PETK::gui::updateMembraneTypeDisplay}
    ttk::radiobutton $container.poretype.doublecone -text "Double Cone Pore" -value "doublecone" \
        -variable ::PETK::gui::membraneType -command {::PETK::gui::updateMembraneTypeDisplay}

    grid $container.poretype.cylindrical $container.poretype.doublecone -sticky w -pady 3 -padx 20

    # === PORE PARAMETERS SECTION ===
    ttk::labelframe $container.params -text "Pore Parameters" -padding 10
    grid $container.params -row $row -column 0 -sticky ew -padx 10 -pady 5
    grid columnconfigure $container.params {1 3 5} -weight 1
    incr row

    # Common parameter: Nanopore thickness
    ttk::label $container.params.thicklbl -text "Nanopore thickness (Å):" -width 20
    ttk::entry $container.params.thick -textvariable ::PETK::gui::nanoporeThickness -width 12 -justify center

    grid $container.params.thicklbl $container.params.thick - - - -sticky ew -pady 5

    # Cylindrical pore parameters
    ttk::frame $container.params.cyl
    ttk::label $container.params.cyl.diamlbl -text "Pore diameter (Å):" -width 20
    ttk::entry $container.params.cyl.diam -textvariable ::PETK::gui::cylindricalDiameter -width 12 -justify center
    ttk::label $container.params.cyl.radiuslbl -text "Corner radius (Å):" -width 15
    ttk::entry $container.params.cyl.radius -textvariable ::PETK::gui::cornerRadius -width 12 -justify center

    grid $container.params.cyl.diamlbl $container.params.cyl.diam $container.params.cyl.radiuslbl $container.params.cyl.radius -sticky ew -pady 3 -padx 5

    # Double cone pore parameters
    ttk::frame $container.params.cone
    ttk::label $container.params.cone.innerlbl -text "Inner diameter (Å):" -width 20
    ttk::entry $container.params.cone.inner -textvariable ::PETK::gui::innerDiameter -width 12 -justify center
    ttk::label $container.params.cone.outerlbl -text "Outer diameter (Å):" -width 15
    ttk::entry $container.params.cone.outer -textvariable ::PETK::gui::outerDiameter -width 12 -justify center

    grid $container.params.cone.innerlbl $container.params.cone.inner $container.params.cone.outerlbl $container.params.cone.outer -sticky ew -pady 3 -padx 5

    # === VISUALIZATION SECTION ===
    ttk::labelframe $container.visualization -text "Pore Geometry Visualization" -padding 10
    grid $container.visualization -row $row -column 0 -sticky nsew -padx 10 -pady 5
    grid columnconfigure $container.visualization {0 1} -weight 1
    grid rowconfigure $container.visualization 1 -weight 1
    incr row

    # Image display frame
    ttk::frame $container.visualization.image
    grid $container.visualization.image -column 0 -row 0 -rowspan 2 -sticky nsew -padx "0 10"
    grid rowconfigure $container.visualization.image 0 -weight 1
    grid columnconfigure $container.visualization.image 0 -weight 1

    # Create canvas for image display
    canvas $container.visualization.image.canvas -width 300 -height 250 -bg white -relief sunken -borderwidth 2
    grid $container.visualization.image.canvas -sticky nsew

    # Store canvas reference for image updates
    set ::PETK::gui::poreImageCanvas $container.visualization.image.canvas
    
    # Parameter summary frame
    ttk::frame $container.visualization.summary
    grid $container.visualization.summary -column 1 -row 0 -sticky nsew -padx "10 0"
    grid rowconfigure $container.visualization.summary 6 -weight 1

    ttk::label $container.visualization.summary.title -text "Current Configuration" -font {TkDefaultFont 11 bold}
    grid $container.visualization.summary.title -sticky w -pady "0 10"

    # Parameter display labels
    ttk::label $container.visualization.summary.typelbl -text "Pore type:" -width 15
    ttk::label $container.visualization.summary.type -textvariable ::PETK::gui::currentPoreType -width 15 -anchor w -relief sunken -background white
    grid $container.visualization.summary.typelbl $container.visualization.summary.type -sticky ew -pady 2

    ttk::label $container.visualization.summary.thicklbl -text "Thickness:" -width 15
    ttk::label $container.visualization.summary.thickval -textvariable ::PETK::gui::displayThickness -width 15 -anchor w -relief sunken -background white
    grid $container.visualization.summary.thicklbl $container.visualization.summary.thickval -sticky ew -pady 2

    ttk::label $container.visualization.summary.param1lbl -text "Parameter 1:" -width 15
    ttk::label $container.visualization.summary.param1val -textvariable ::PETK::gui::displayParam1 -width 15 -anchor w -relief sunken -background white
    grid $container.visualization.summary.param1lbl $container.visualization.summary.param1val -sticky ew -pady 2

    ttk::label $container.visualization.summary.param2lbl -text "Parameter 2:" -width 15
    ttk::label $container.visualization.summary.param2val -textvariable ::PETK::gui::displayParam2 -width 15 -anchor w -relief sunken -background white
    grid $container.visualization.summary.param2lbl $container.visualization.summary.param2val -sticky ew -pady 2

    # Volume calculation
    ttk::label $container.visualization.summary.vollbl -text "Pore volume:" -width 15
    ttk::label $container.visualization.summary.volval -textvariable ::PETK::gui::calculatedVolume -width 15 -anchor w -relief sunken -background white
    grid $container.visualization.summary.vollbl $container.visualization.summary.volval -sticky ew -pady 2

    # Pore validity status
    ttk::label $container.visualization.summary.validitylbl -text "Status:" -width 15
    ttk::label $container.visualization.summary.validityval -textvariable ::PETK::gui::poreValidityStatus -width 15 -anchor w -relief sunken
    grid $container.visualization.summary.validitylbl $container.visualization.summary.validityval -sticky ew -pady 2

    # Update button
    ttk::button $container.visualization.summary.update -text "Update Preview" -command {::PETK::gui::updatePoreVisualization}
    grid $container.visualization.summary.update -sticky ew -pady "10 0"

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab1.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab1.canvas $canvas_window]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab1.canvas $canvas_window]
    
    # Store references for later use
    set ::PETK::gui::mainCanvas $tab1.canvas
    set ::PETK::gui::canvasWindow $canvas_window
    set ::PETK::gui::contentContainer $container
    
    # Initialize variables and display
    ::PETK::gui::updateMembraneTypeDisplay
    ::PETK::gui::loadPoreImages
    
    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab1.canvas $canvas_window]
}

####################################################
# Variable Initialization
####################################################

proc ::PETK::gui::initializeMembraneVariables {} {
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
    if {![info exists ::PETK::gui::sysCutoff]} {
        set ::PETK::gui::sysCutoff "50.0"
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
    ::PETK::gui::updateParameterDisplay
    ::PETK::gui::calculateBoxDimensions
}

proc ::PETK::gui::selectWorkDir {} {
    set tempdir [tk_chooseDirectory -title "Select project folder" -initialdir [pwd]]
    if {![string eq $tempdir ""]} { 
        set ::PETK::gui::workdir $tempdir 
    }
}

####################################################
# Box Dimension Management Functions
####################################################
proc ::PETK::gui::toggleBoxDimensionMode {} {
    set container $::PETK::gui::window.hlf.nb.tab1.canvas.content
    
    # Hide both frames first
    catch {grid forget $container.boxdim.manual}
    catch {grid forget $container.boxdim.auto_display}
    
    if {$::PETK::gui::autoCalculateBoxDimensions} {
        # Show auto-calculated display
        grid $container.boxdim.auto_display -row 1 -column 0 -columnspan 6 -sticky ew
        ::PETK::gui::calculateBoxDimensions
        ::PETK::gui::updateMembraneStatus "Switched to auto-calculated box dimensions"
    } else {
        # Show manual input fields
        grid $container.boxdim.manual -row 1 -column 0 -columnspan 6 -sticky ew
        ::PETK::gui::updateMembraneStatus "Switched to manual box dimensions"
    }
}

proc ::PETK::gui::calculateBoxDimensions {} {
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
    set padding $::PETK::gui::sysCutoff
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
    puts "  X: -$xy_size to $xy_size Å"
    puts "  Y: -$xy_size to $xy_size Å"
    puts "  Z: $z_min to $z_max Å"
    puts "  Pore radius: $pore_radius Å"
    puts "  Movement range: $z_end to $z_start Å"
}

####################################################
# Pore Type Management Functions
####################################################

proc ::PETK::gui::updateMembraneTypeDisplay {} {
    set container $::PETK::gui::window.hlf.nb.tab1.canvas.content
    
    # Hide all parameter frames first
    catch {grid forget $container.params.cyl}
    catch {grid forget $container.params.cone}
    
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        grid $container.params.cyl -row 1 -column 0 -columnspan 6 -sticky ew -pady 5
        set ::PETK::gui::currentPoreType "Cylindrical"
        ::PETK::gui::updateMembraneStatus "Switched to cylindrical pore mode"
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        grid $container.params.cone -row 1 -column 0 -columnspan 6 -sticky ew -pady 5
        set ::PETK::gui::currentPoreType "Double Cone"
        ::PETK::gui::updateMembraneStatus "Switched to double cone pore mode"
    }
    
    ::PETK::gui::updateParameterDisplay
    ::PETK::gui::updatePoreVisualization
    ::PETK::gui::calculateBoxDimensions
}

proc ::PETK::gui::updateParameterDisplay {} {
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
            # Volume of truncated cone: π/3 * h * (r1² + r1*r2 + r2²)
            set volume [expr {3.14159/3.0 * $h * ($r1*$r1 + $r1*$r2 + $r2*$r2) * 2}]
            set ::PETK::gui::calculatedVolume [format "%.1f Å³" $volume]
        } else {
            set ::PETK::gui::calculatedVolume "Invalid params"
        }
    }
    
    # Update box dimensions and grid points
    ::PETK::gui::calculateBoxDimensions
}

####################################################
# Image Loading and Visualization Functions
####################################################

proc ::PETK::gui::loadPoreImages {} {
    # Load pore shape images
    set ::PETK::gui::poreImages(cylindrical) ""
    set ::PETK::gui::poreImages(doublecone) ""
    
    # Try to load cylindrical pore image
    set cyl_path "./shapes/shape2.gif"
    if {[file exists $cyl_path]} {
        if {[catch {
            set ::PETK::gui::poreImages(cylindrical) [image create photo -file $cyl_path]
        } error]} {
            puts "Warning: Could not load cylindrical pore image: $error"
        }
    } else {
        puts "Warning: Cylindrical pore image not found at: $cyl_path"
    }
    
    # Try to load double cone pore image
    set cone_path "./shapes/shape1.gif"
    if {[file exists $cone_path]} {
        if {[catch {
            set ::PETK::gui::poreImages(doublecone) [image create photo -file $cone_path]
        } error]} {
            puts "Warning: Could not load double cone pore image: $error"
        }
    } else {
        puts "Warning: Double cone pore image not found at: $cone_path"
    }
    
    # Load initial image
    ::PETK::gui::updatePoreVisualization
}

proc ::PETK::gui::updatePoreVisualization {} {
    if {![info exists ::PETK::gui::poreImageCanvas] || ![winfo exists $::PETK::gui::poreImageCanvas]} {
        return
    }
    
    set canvas $::PETK::gui::poreImageCanvas
    
    # Clear previous image
    $canvas delete all
    
    # Get current pore type image
    set image_key $::PETK::gui::membraneType
    
    if {[info exists ::PETK::gui::poreImages($image_key)] && $::PETK::gui::poreImages($image_key) ne ""} {
        set image $::PETK::gui::poreImages($image_key)
        
        # Get canvas and image dimensions
        set canvas_width [winfo reqwidth $canvas]
        set canvas_height [winfo reqheight $canvas]
        if {$canvas_width <= 1} {set canvas_width 300}
        if {$canvas_height <= 1} {set canvas_height 250}
        
        set img_width [image width $image]
        set img_height [image height $image]
        
        # Calculate scaling to fit canvas while maintaining aspect ratio
        set scale_x [expr {double($canvas_width) / $img_width}]
        set scale_y [expr {double($canvas_height) / $img_height}]
        set scale [expr {min($scale_x, $scale_y) * 0.8}]  ; # 80% of canvas size
        
        # Create scaled image if needed
        if {$scale < 1.0} {
            set new_width [expr {int($img_width * $scale)}]
            set new_height [expr {int($img_height * $scale)}]
            
            set scaled_image [image create photo]
            $scaled_image copy $image -subsample [expr {int(1.0/$scale)}]
            
            # Center the image on canvas
            set x [expr {($canvas_width - $new_width) / 2}]
            set y [expr {($canvas_height - $new_height) / 2}]
            $canvas create image $x $y -anchor nw -image $scaled_image
        } else {
            # Center the original image
            set x [expr {($canvas_width - $img_width) / 2}]
            set y [expr {($canvas_height - $img_height) / 2}]
            $canvas create image $x $y -anchor nw -image $image
        }
        
    } else {
        # Draw a simple schematic if no image available
        ::PETK::gui::drawPoreSchematic $canvas
    }
    
    # Add parameter labels on the image
    # ::PETK::gui::addParameterLabels $canvas
    ::PETK::gui::updateParameterDisplay
    ::PETK::gui::updatePoreVmdVisualization 
}

proc ::PETK::gui::drawPoreSchematic {canvas} {
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

proc ::PETK::gui::addParameterLabels {canvas} {
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

proc ::PETK::gui::updateMembraneStatus {message} {
    if {[info exists ::PETK::gui::membraneStatusLabel] && [winfo exists $::PETK::gui::membraneStatusLabel]} {
        $::PETK::gui::membraneStatusLabel configure -text $message
        update
    }
    puts "PETK Membrane Status: $message"
}



proc ::PETK::gui::updatePoreVmdVisualization {} {
    # Clear any existing graphics
    draw delete all
    
    # Debug: Print current parameter values
    # puts "=== VMD Visualization Debug ==="
    # puts "autoCalculateBoxDimensions: $::PETK::gui::autoCalculateBoxDimensions"
    # if {$::PETK::gui::autoCalculateBoxDimensions} {
    #    puts "autoBoxX: '$::PETK::gui::autoBoxX'"
    #    puts "autoBoxY: '$::PETK::gui::autoBoxY'"
    #    puts "autoBoxZ: '$::PETK::gui::autoBoxZ'"
    #} else {
    #    puts "boxSizeX: $::PETK::gui::boxSizeX"
    #    puts "boxSizeY: $::PETK::gui::boxSizeY" 
    #    puts "boxSizeZ: $::PETK::gui::boxSizeZ"
    #}
    #puts "membraneType: $::PETK::gui::membraneType"
    #puts "nanoporeThickness: $::PETK::gui::nanoporeThickness"
    #if {$::PETK::gui::membraneType eq "cylindrical"} {
    #    puts "cylindricalDiameter: $::PETK::gui::cylindricalDiameter"
    #    if {[info exists ::PETK::gui::cornerRadius]} {
    #        puts "cornerRadius: $::PETK::gui::cornerRadius"
    #    }
    #} else {
    #    puts "innerDiameter: $::PETK::gui::innerDiameter"
    #    puts "outerDiameter: $::PETK::gui::outerDiameter"
    #}
    #puts "=============================="
    
    # Get box dimensions based on current mode
    if {$::PETK::gui::autoCalculateBoxDimensions} {
        # Try different parsing approaches for auto-calculated dimensions
        if {[catch {
            # Method 1: Extract size from parentheses - format: "min to max (size)"
            regexp {\(([0-9.-]+)\)} $::PETK::gui::autoBoxX xbox_match xbox
            regexp {\(([0-9.-]+)\)} $::PETK::gui::autoBoxY ybox_match ybox  
            regexp {\(([0-9.-]+)\)} $::PETK::gui::autoBoxZ zbox_match zbox
            puts "Parsed box dimensions: $xbox x $ybox x $zbox"
        } error]} {
            puts "Error parsing auto dimensions: $error"
            # Fallback to manual dimensions if auto parsing fails
            set xbox $::PETK::gui::boxSizeX
            set ybox $::PETK::gui::boxSizeY
            set zbox $::PETK::gui::boxSizeZ
            puts "Using fallback manual dimensions: $xbox x $ybox x $zbox"
        }
    } else {
        # Use manual dimensions
        set xbox $::PETK::gui::boxSizeX
        set ybox $::PETK::gui::boxSizeY
        set zbox $::PETK::gui::boxSizeZ
        puts "Using manual dimensions: $xbox x $ybox x $zbox"
    }
    
    # Validate dimensions
    if {![string is double $xbox] || ![string is double $ybox] || ![string is double $zbox] || 
        $xbox <= 0 || $ybox <= 0 || $zbox <= 0} {
        puts "ERROR: Invalid box dimensions - xbox:$xbox ybox:$ybox zbox:$zbox"
        set ::PETK::gui::poreValidityStatus "INVALID - Invalid box dimensions"
        return
    }
    
    # Perform comprehensive pore validation using the new validation system
    set validation_status [::PETK::gui::updatePoreValidityStatus $xbox $ybox $zbox]
    
    # Calculate box center and boundaries
    set center_x [expr $xbox / 2.0]
    set center_y [expr $ybox / 2.0]
    set center_z [expr $zbox / 2.0]
    
    set min_x [expr -$center_x]
    set max_x [expr $center_x]
    set min_y [expr -$center_y]
    set max_y [expr $center_y]
    set min_z [expr -$center_z]
    set max_z [expr $center_z]
    
    # Draw simulation box outline
    draw color cyan
    draw material Transparent
    
    # Bottom face
    draw line [list $min_x $min_y $min_z] [list $max_x $min_y $min_z]
    draw line [list $max_x $min_y $min_z] [list $max_x $max_y $min_z]
    draw line [list $max_x $max_y $min_z] [list $min_x $max_y $min_z]
    draw line [list $min_x $max_y $min_z] [list $min_x $min_y $min_z]
    
    # Top face
    draw line [list $min_x $min_y $max_z] [list $max_x $min_y $max_z]
    draw line [list $max_x $min_y $max_z] [list $max_x $max_y $max_z]
    draw line [list $max_x $max_y $max_z] [list $min_x $max_y $max_z]
    draw line [list $min_x $max_y $max_z] [list $min_x $min_y $max_z]
    
    # Vertical edges
    draw line [list $min_x $min_y $min_z] [list $min_x $min_y $max_z]
    draw line [list $max_x $min_y $min_z] [list $max_x $min_y $max_z]
    draw line [list $max_x $max_y $min_z] [list $max_x $max_y $max_z]
    draw line [list $min_x $max_y $min_z] [list $min_x $max_y $max_z]
    
    # Get nanopore parameters
    set thickness $::PETK::gui::nanoporeThickness
    set membrane_half_thickness [expr $thickness / 2.0]
    set pore_center_z 0.0  ; # Center the pore in Z direction
    set pore_start_z [expr $pore_center_z - $membrane_half_thickness]
    set pore_end_z [expr $pore_center_z + $membrane_half_thickness]
    
    # Draw membrane extending to box edges with pore as a hole
    draw color purple
    draw material Opaque
    
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        # Get pore parameters
        set pore_radius [expr $::PETK::gui::cylindricalDiameter / 2.0]
        
        # Check if corner radius is defined and > 0 for chamfered edges
        set has_corner_radius [expr {[info exists ::PETK::gui::cornerRadius] && $::PETK::gui::cornerRadius > 0}]
        if {$has_corner_radius} {
            set corner_radius $::PETK::gui::cornerRadius
            set edge_radius [expr $pore_radius + $corner_radius]
            set chamfer_depth $corner_radius  ; # Default chamfer depth
            puts "Using corner radius: $corner_radius, edge radius: $edge_radius"
        } else {
            puts "No corner radius specified, using regular cylindrical pore"
        }
        
        # Add visual warning if configuration is invalid
        if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus]} {
            draw color red
            draw material Transparent
            # Draw warning outline around the problematic area
            set warn_radius [expr $pore_radius + 10.0]
            if {$has_corner_radius} {
                set warn_radius [expr $edge_radius + 10.0]
            }
            
            for {set angle 0} {$angle < 360} {incr angle 30} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + 30) * 3.14159 / 180.0]
                set x1 [expr $warn_radius * cos($angle_rad)]
                set y1 [expr $warn_radius * sin($angle_rad)]
                set x2 [expr $warn_radius * cos($next_angle_rad)]
                set y2 [expr $warn_radius * sin($next_angle_rad)]
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z]
            }
            
            # Reset to purple for membrane drawing
            draw color purple
            draw material Opaque
        }
        
        # Draw membrane as rectangular sheets with hole
        set grid_size 5.0  ; # Grid spacing for membrane representation
        
        # Draw membrane in XY plane at multiple Z levels to show chamfer effect
        set z_levels 10  ; # Number of Z levels to draw for chamfer visualization
        
        for {set z_level 0} {$z_level < $z_levels} {incr z_level} {
            # Calculate Z position for this level
            set z_frac [expr double($z_level) / ($z_levels - 1)]  ; # 0 to 1
            set z_pos [expr $pore_start_z + $z_frac * $thickness]
            
            # Calculate local pore radius at this Z level
            if {$has_corner_radius} {
                # Distance from membrane edges
                set z_edge_dist [expr $membrane_half_thickness - abs($z_pos)]
                
                # Apply chamfer only within chamfer_depth from edges
                if {$z_edge_dist < $chamfer_depth} {
                    # Linear interpolation in chamfer zone
                    set chamfer_progress [expr $z_edge_dist / $chamfer_depth]
                    set chamfer_progress [expr max(0.0, min(1.0, $chamfer_progress))]
                    set local_pore_radius [expr $edge_radius + ($pore_radius - $edge_radius) * $chamfer_progress]
                } else {
                    set local_pore_radius $pore_radius
                }
            } else {
                set local_pore_radius $pore_radius
            }
            
            # Draw membrane grid at this Z level
            for {set x $min_x} {$x < $max_x} {set x [expr $x + $grid_size]} {
                for {set y $min_y} {$y < $max_y} {set y [expr $y + $grid_size]} {
                    
                    # Calculate grid cell boundaries, ensuring they don't exceed box limits
                    set x1 $x
                    set x2 [expr min($x + $grid_size, $max_x)]
                    set y1 $y
                    set y2 [expr min($y + $grid_size, $max_y)]
                    
                    # Check if this grid cell is outside the pore area
                    set cell_center_x [expr ($x1 + $x2) / 2.0]
                    set cell_center_y [expr ($y1 + $y2) / 2.0]
                    set distance [expr sqrt($cell_center_x*$cell_center_x + $cell_center_y*$cell_center_y)]
                    
                    # Only draw if cell center is outside local pore radius
                    if {$distance > $local_pore_radius} {
                        # Draw only at top and bottom for cleaner visualization
                        if {$z_level == 0} {
                            # Bottom face
                            draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                            draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] [list $x1 $y2 $pore_start_z]
                        } elseif {$z_level == [expr $z_levels - 1]} {
                            # Top face
                            draw triangle [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] [list $x2 $y1 $pore_end_z]
                            draw triangle [list $x1 $y1 $pore_end_z] [list $x1 $y2 $pore_end_z] [list $x2 $y2 $pore_end_z]
                        }
                    }
                }
            }
        }
        
        # Draw pore outline to show the hole and chamfer effect
        draw color yellow
        set angle_step 15
        
        # If we have corner radius, draw multiple circles at different Z levels
        if {$has_corner_radius} {
            set outline_z_levels 5  ; # Fewer levels for outline
            for {set z_level 0} {$z_level < $outline_z_levels} {incr z_level} {
                set z_frac [expr double($z_level) / ($outline_z_levels - 1)]
                set z_pos [expr $pore_start_z + $z_frac * $thickness]
                
                # Calculate local pore radius at this Z level
                set z_edge_dist [expr $membrane_half_thickness - abs($z_pos)]
                if {$z_edge_dist < $chamfer_depth} {
                    set chamfer_progress [expr $z_edge_dist / $chamfer_depth]
                    set chamfer_progress [expr max(0.0, min(1.0, $chamfer_progress))]
                    set local_pore_radius [expr $edge_radius + ($pore_radius - $edge_radius) * $chamfer_progress]
                } else {
                    set local_pore_radius $pore_radius
                }
                
                # Draw circle at this Z level
                for {set angle 0} {$angle < 360} {incr angle $angle_step} {
                    set angle_rad [expr $angle * 3.14159 / 180.0]
                    set next_angle_rad [expr ($angle + $angle_step) * 3.14159 / 180.0]
                    set x1 [expr $local_pore_radius * cos($angle_rad)]
                    set y1 [expr $local_pore_radius * sin($angle_rad)]
                    set x2 [expr $local_pore_radius * cos($next_angle_rad)]
                    set y2 [expr $local_pore_radius * sin($next_angle_rad)]
                    
                    draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $z_pos]
                }
                
                # Draw connecting lines between levels
                if {$z_level < [expr $outline_z_levels - 1]} {
                    set next_z_frac [expr double($z_level + 1) / ($outline_z_levels - 1)]
                    set next_z_pos [expr $pore_start_z + $next_z_frac * $thickness]
                    
                    for {set angle 0} {$angle < 360} {incr angle [expr $angle_step * 2]} {
                        set angle_rad [expr $angle * 3.14159 / 180.0]
                        set x1 [expr $local_pore_radius * cos($angle_rad)]
                        set y1 [expr $local_pore_radius * sin($angle_rad)]
                        
                        # Calculate radius at next level
                        set next_z_edge_dist [expr $membrane_half_thickness - abs($next_z_pos)]
                        if {$next_z_edge_dist < $chamfer_depth} {
                            set next_chamfer_progress [expr $next_z_edge_dist / $chamfer_depth]
                            set next_chamfer_progress [expr max(0.0, min(1.0, $next_chamfer_progress))]
                            set next_local_pore_radius [expr $edge_radius + ($pore_radius - $edge_radius) * $next_chamfer_progress]
                        } else {
                            set next_local_pore_radius $pore_radius
                        }
                        
                        set x2 [expr $next_local_pore_radius * cos($angle_rad)]
                        set y2 [expr $next_local_pore_radius * sin($angle_rad)]
                        
                        draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $next_z_pos]
                    }
                }
            }
        } else {
            # Regular cylindrical pore - draw simple circles
            for {set angle 0} {$angle < 360} {incr angle $angle_step} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + $angle_step) * 3.14159 / 180.0]
                set x1 [expr $pore_radius * cos($angle_rad)]
                set y1 [expr $pore_radius * sin($angle_rad)]
                set x2 [expr $pore_radius * cos($next_angle_rad)]
                set y2 [expr $pore_radius * sin($next_angle_rad)]
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z]
                draw line [list $x1 $y1 $pore_start_z] [list $x1 $y1 $pore_end_z]
            }
        }
        
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        # Draw membrane with double cone pore hole
        set inner_radius [expr $::PETK::gui::innerDiameter / 2.0]
        set outer_radius [expr $::PETK::gui::outerDiameter / 2.0]
        set middle_z $pore_center_z
        
        # Add visual warning if configuration is invalid
        if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus]} {
            draw color red
            draw material Transparent
            # Draw warning outline around the problematic area
            set warn_radius [expr $outer_radius + 10.0]
            
            for {set angle 0} {$angle < 360} {incr angle 30} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + 30) * 3.14159 / 180.0]
                set x1 [expr $warn_radius * cos($angle_rad)]
                set y1 [expr $warn_radius * sin($angle_rad)]
                set x2 [expr $warn_radius * cos($next_angle_rad)]
                set y2 [expr $warn_radius * sin($next_angle_rad)]
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z]
            }
            
            # Reset to purple for membrane drawing
            draw color purple
            draw material Opaque
        }
        
        # Draw membrane as rectangular sheets with cone-shaped hole
        set grid_size 5.0
        
        # Draw membrane exactly within box boundaries
        for {set x $min_x} {$x < $max_x} {set x [expr $x + $grid_size]} {
            for {set y $min_y} {$y < $max_y} {set y [expr $y + $grid_size]} {
                
                # Calculate grid cell boundaries, ensuring they don't exceed box limits
                set x1 $x
                set x2 [expr min($x + $grid_size, $max_x)]
                set y1 $y
                set y2 [expr min($y + $grid_size, $max_y)]
                
                # Check if this grid cell is outside the pore area
                set cell_center_x [expr ($x1 + $x2) / 2.0]
                set cell_center_y [expr ($y1 + $y2) / 2.0]
                set distance [expr sqrt($cell_center_x*$cell_center_x + $cell_center_y*$cell_center_y)]
                
                # Only draw if cell center is outside the outer pore radius
                if {$distance > $outer_radius} {
                    # Draw top and bottom membrane faces
                    draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                    draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] [list $x1 $y2 $pore_start_z]
                    draw triangle [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] [list $x2 $y1 $pore_end_z]
                    draw triangle [list $x1 $y1 $pore_end_z] [list $x1 $y2 $pore_end_z] [list $x2 $y2 $pore_end_z]
                }
            }
        }
        
        # Draw pore outline to show the double cone hole
        draw color yellow
        for {set angle 0} {$angle < 360} {incr angle 20} {
            set angle_rad [expr $angle * 3.14159 / 180.0]
            set next_angle_rad [expr ($angle + 20) * 3.14159 / 180.0]
            
            # Outer radius points
            set x1_out [expr $outer_radius * cos($angle_rad)]
            set y1_out [expr $outer_radius * sin($angle_rad)]
            set x2_out [expr $outer_radius * cos($next_angle_rad)]
            set y2_out [expr $outer_radius * sin($next_angle_rad)]
            
            # Inner radius points  
            set x1_in [expr $inner_radius * cos($angle_rad)]
            set y1_in [expr $inner_radius * sin($angle_rad)]
            set x2_in [expr $inner_radius * cos($next_angle_rad)]
            set y2_in [expr $inner_radius * sin($next_angle_rad)]
            
            # Draw cone outline
            draw line [list $x1_out $y1_out $pore_start_z] [list $x2_out $y2_out $pore_start_z]
            draw line [list $x1_in $y1_in $middle_z] [list $x2_in $y2_in $middle_z]
            draw line [list $x1_out $y1_out $pore_start_z] [list $x1_in $y1_in $middle_z]
            draw line [list $x1_in $y1_in $middle_z] [list $x1_out $y1_out $pore_end_z]
            draw line [list $x1_out $y1_out $pore_end_z] [list $x2_out $y2_out $pore_end_z]
        }
    }
    
    # Final status message
    if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus]} {
        puts "Pore visualization completed - WARNING: INVALID CONFIGURATION!"
        puts "See detailed validation report above."
    } else {
        if {$::PETK::gui::membraneType eq "cylindrical" && [info exists ::PETK::gui::cornerRadius] && $::PETK::gui::cornerRadius > 0} {
            puts "Pore visualization completed for chamfered cylindrical pore (corner radius: $::PETK::gui::cornerRadius) with purple membrane"
        } else {
            puts "Pore visualization completed for $::PETK::gui::membraneType pore with purple membrane"
        }
    }
    
    puts "Current validity status: $::PETK::gui::poreValidityStatus"
}

proc ::PETK::gui::updatePoreDiameter {} {
    
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        # For cylindrical pores, use the cylindrical diameter
        if {[info exists ::PETK::gui::cylindricalDiameter] && 
            [string is double $::PETK::gui::cylindricalDiameter] && 
            $::PETK::gui::cylindricalDiameter > 0} {
            set ::PETK::gui::poreDiameter $::PETK::gui::cylindricalDiameter
            puts "Pore diameter set to cylindrical diameter: $::PETK::gui::poreDiameter Å"
        } else {
            set ::PETK::gui::poreDiameter ""
            puts "Warning: Invalid or missing cylindrical diameter"
        }
        
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        # For double cone pores, use the smaller of inner and outer diameters
        if {[info exists ::PETK::gui::innerDiameter] && 
            [info exists ::PETK::gui::outerDiameter] &&
            [string is double $::PETK::gui::innerDiameter] && 
            [string is double $::PETK::gui::outerDiameter] &&
            $::PETK::gui::innerDiameter > 0 && 
            $::PETK::gui::outerDiameter > 0} {
            
            # Set to the minimum of inner and outer diameters
            if {$::PETK::gui::innerDiameter <= $::PETK::gui::outerDiameter} {
                set ::PETK::gui::poreDiameter $::PETK::gui::innerDiameter
                puts "Pore diameter set to inner diameter (smaller): $::PETK::gui::poreDiameter Å"
            } else {
                set ::PETK::gui::poreDiameter $::PETK::gui::outerDiameter
                puts "Pore diameter set to outer diameter (smaller): $::PETK::gui::poreDiameter Å"
            }
        } else {
            set ::PETK::gui::poreDiameter ""
            puts "Warning: Invalid or missing inner/outer diameters"
        }
        
    } else {
        # Unknown membrane type
        set ::PETK::gui::poreDiameter ""
        puts "Warning: Unknown membrane type: $::PETK::gui::membraneType"
    }
    
    # Trigger any dependent updates (like fit status)
    if {[info exists ::PETK::gui::poreDiameter] && $::PETK::gui::poreDiameter ne ""} {
        # Call fit status update if it exists
        if {[info procs ::PETK::gui::updateFitStatus] ne ""} {
            after idle ::PETK::gui::updateFitStatus
        }
    }
}

proc ::PETK::gui::validatePoreConfiguration {xbox ybox zbox} {
    
    set validation_messages {}
    set is_valid 1
    set critical_errors {}
    set warnings {}
    
    # Calculate box half-dimensions for radial checks
    set box_half_x [expr $xbox / 2.0]
    set box_half_y [expr $ybox / 2.0] 
    set box_half_z [expr $zbox / 2.0]
    set min_box_half_xy [expr min($box_half_x, $box_half_y)]
    
    # Get pore parameters
    set thickness $::PETK::gui::nanoporeThickness
    set membrane_half_thickness [expr $thickness / 2.0]
    
    # 1. THICKNESS VALIDATION - Check against Z dimension
    if {$thickness > $zbox} {
        set is_valid 0
        lappend critical_errors "Membrane thickness ($thickness Å) > box Z dimension ($zbox Å)"
        lappend critical_errors "  Membrane completely fills or exceeds the box height"
        lappend critical_errors "  Required: thickness < $zbox Å"
    } elseif {$thickness > [expr 0.8 * $zbox]} {
        lappend warnings "Membrane thickness ($thickness Å) > 80% of box Z dimension ($zbox Å)"
        lappend warnings "  Consider increasing box height for better simulation"
    }
    
    # 2. PORE TYPE SPECIFIC VALIDATIONS
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        set pore_radius [expr $::PETK::gui::cylindricalDiameter / 2.0]
        
        # Check basic pore radius
        if {$pore_radius > $min_box_half_xy} {
            set is_valid 0
            lappend critical_errors "Pore radius ($pore_radius Å) > minimum box half-dimension ($min_box_half_xy Å)"
            lappend critical_errors "  Box half-dimensions: X=$box_half_x Å, Y=$box_half_y Å"
            lappend critical_errors "  Required: pore radius < $min_box_half_xy Å"
        }
        
        # Check corner radius if present
        if {[info exists ::PETK::gui::cornerRadius] && $::PETK::gui::cornerRadius > 0} {
            set corner_radius $::PETK::gui::cornerRadius
            set edge_radius [expr $pore_radius + $corner_radius]
            
            if {$edge_radius > $min_box_half_xy} {
                set is_valid 0
                lappend critical_errors "Edge radius ($edge_radius Å) > minimum box half-dimension ($min_box_half_xy Å)"
                lappend critical_errors "  Edge radius = pore_radius ($pore_radius) + corner_radius ($corner_radius)"
                lappend critical_errors "  Box half-dimensions: X=$box_half_x Å, Y=$box_half_y Å"
                lappend critical_errors "  Required: edge_radius < $min_box_half_xy Å"
                lappend critical_errors "  Solutions:"
                lappend critical_errors "    1. Increase box dimensions (X,Y > [expr 2.0 * $edge_radius] Å)"
                lappend critical_errors "    2. Reduce corner radius (< [expr $min_box_half_xy - $pore_radius] Å)"
                lappend critical_errors "    3. Reduce pore radius"
            } elseif {$edge_radius > [expr 0.8 * $min_box_half_xy]} {
                lappend warnings "Edge radius ($edge_radius Å) > 80% of minimum box half-dimension"
                lappend warnings "  Consider increasing box size for better simulation quality"
            }
        } elseif {$pore_radius > [expr 0.8 * $min_box_half_xy]} {
            lappend warnings "Pore radius ($pore_radius Å) > 80% of minimum box half-dimension"
            lappend warnings "  Consider increasing box size for better simulation quality"
        }
        
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        set inner_radius [expr $::PETK::gui::innerDiameter / 2.0]
        set outer_radius [expr $::PETK::gui::outerDiameter / 2.0]
        
        # Check inner radius
        if {$inner_radius > $min_box_half_xy} {
            set is_valid 0
            lappend critical_errors "Inner radius ($inner_radius Å) > minimum box half-dimension ($min_box_half_xy Å)"
            lappend critical_errors "  Box half-dimensions: X=$box_half_x Å, Y=$box_half_y Å"
            lappend critical_errors "  Required: inner radius < $min_box_half_xy Å"
        }
        
        # Check outer radius (more critical)
        if {$outer_radius > $min_box_half_xy} {
            set is_valid 0
            lappend critical_errors "Outer radius ($outer_radius Å) > minimum box half-dimension ($min_box_half_xy Å)"
            lappend critical_errors "  Box half-dimensions: X=$box_half_x Å, Y=$box_half_y Å"
            lappend critical_errors "  Required: outer radius < $min_box_half_xy Å"
        } elseif {$outer_radius > [expr 0.8 * $min_box_half_xy]} {
            lappend warnings "Outer radius ($outer_radius Å) > 80% of minimum box half-dimension"
            lappend warnings "  Consider increasing box size for better simulation quality"
        }
        
    }
    
    # 3. ADDITIONAL GEOMETRIC VALIDATIONS
    
    # Check for minimum simulation space around pore
    if {$::PETK::gui::autoCalculateBoxDimensions == 0} {
        # Only check minimum spacing for manual box dimensions
        set min_spacing 10.0  ; # Minimum 10 Å spacing recommendation
        
        if {$::PETK::gui::membraneType eq "cylindrical"} {
            set effective_radius $pore_radius
            if {[info exists ::PETK::gui::cornerRadius] && $::PETK::gui::cornerRadius > 0} {
                set effective_radius [expr $pore_radius + $::PETK::gui::cornerRadius]
            }
            
            set available_space [expr $min_box_half_xy - $effective_radius]
            if {$available_space < $min_spacing} {
                lappend warnings "Limited space around pore: only $available_space Å"
                lappend warnings "  Recommended minimum: $min_spacing Å for proper solvation"
            }
        } elseif {$::PETK::gui::membraneType eq "doublecone"} {
            set available_space [expr $min_box_half_xy - $outer_radius]
            if {$available_space < $min_spacing} {
                lappend warnings "Limited space around pore: only $available_space Å"
                lappend warnings "  Recommended minimum: $min_spacing Å for proper solvation"
            }
        }
    }
    
    # 4. COMPILE VALIDATION RESULTS
    set detailed_message "=== PORE CONFIGURATION VALIDATION ===\n"
    append detailed_message "Box dimensions: $xbox × $ybox × $zbox Å\n"
    append detailed_message "Membrane thickness: $thickness Å\n"
    append detailed_message "Pore type: $::PETK::gui::membraneType\n"
    
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        append detailed_message "Pore diameter: $::PETK::gui::cylindricalDiameter Å\n"
        if {[info exists ::PETK::gui::cornerRadius] && $::PETK::gui::cornerRadius > 0} {
            append detailed_message "Corner radius: $::PETK::gui::cornerRadius Å\n"
        }
    } else {
        append detailed_message "Inner diameter: $::PETK::gui::innerDiameter Å\n"
        append detailed_message "Outer diameter: $::PETK::gui::outerDiameter Å\n"
    }
    
    if {$is_valid} {
        set status "Valid"
        append detailed_message "\n✓ CONFIGURATION IS VALID\n"
        ::PETK::gui::updatePoreDiameter
    } else {
        set status "INVALID"
        append detailed_message "\n✗ CONFIGURATION IS INVALID\n"
        append detailed_message "\nCRITICAL ERRORS:\n"
        foreach error $critical_errors {
            append detailed_message "  • $error\n"
        }
    }
    
    if {[llength $warnings] > 0} {
        append detailed_message "\nWARNINGS:\n"
        foreach warning $warnings {
            append detailed_message "  ⚠ $warning\n"
        }
    }
    
    append detailed_message "=====================================\n"
    
    return [list $status $detailed_message]
}


proc ::PETK::gui::updatePoreValidityStatus {xbox ybox zbox} {
    
    # Update the global pore validity status variable.
    
    #Args:
    #    xbox, ybox, zbox: Box dimensions in Angstroms
    
    
    set validation_result [::PETK::gui::validatePoreConfiguration $xbox $ybox $zbox]
    set status [lindex $validation_result 0]
    set detailed_message [lindex $validation_result 1]
    
    if {$status eq "Valid"} {
        set ::PETK::gui::poreValidityStatus "Valid"
    } else {
        set ::PETK::gui::poreValidityStatus "INVALID - Check console"
    }
    
    # Print detailed validation report
    puts $detailed_message
    
    return $status
}

####################################################
# Tab 2: Analyte Setup
####################################################


proc ::PETK::gui::buildTab2 {tab2} {

    # === CRITICAL: Configure tab to expand ===
    grid columnconfigure $tab2 0 -weight 1
    grid rowconfigure $tab2 0 -weight 1

    # Create scrollable canvas container
    canvas $tab2.canvas -highlightthickness 0
    ttk::scrollbar $tab2.vscroll -orient vertical -command [list $tab2.canvas yview]
    ttk::scrollbar $tab2.hscroll -orient horizontal -command [list $tab2.canvas xview]
    
    # Configure canvas scrolling
    $tab2.canvas configure -yscrollcommand [list $tab2.vscroll set]
    $tab2.canvas configure -xscrollcommand [list $tab2.hscroll set]
    
    # Create the actual content frame inside the canvas
    ttk::frame $tab2.canvas.content
    set canvas_window [$tab2.canvas create window 0 0 -anchor nw -window $tab2.canvas.content]
    
    # Grid the canvas and scrollbars with proper expansion
    grid $tab2.canvas -row 0 -column 0 -sticky nsew
    grid $tab2.vscroll -row 0 -column 1 -sticky ns
    grid $tab2.hscroll -row 1 -column 0 -sticky ew
    
    # Configure grid weights - CRITICAL for expansion
    grid rowconfigure $tab2 0 -weight 1
    grid columnconfigure $tab2 0 -weight 1

    # Now use the content frame as your container
    set container $tab2.canvas.content
    grid columnconfigure $container 0 -weight 1
    grid rowconfigure $container {3 4} -weight 1  ; # Make results and details expandable
    
    set row 0
    # === Analyte Input Frame ===
    ttk::labelframe $container.input -text "Analyte Input" -padding 10
    grid $container.input -column 0 -row $row -sticky ew -pady "0 10" -padx 10
    grid columnconfigure $container.input {1 2} -weight 1
    incr row

    ttk::label  $container.input.pdblbl -text "PDB file:" -width 12
    ttk::entry  $container.input.pdb -textvariable ::PETK::gui::analytePDB -width 35
    ttk::button $container.input.browse -text "Browse..." -command ::PETK::gui::browseAnalytePdb
    ttk::button $container.input.analyze -text "Analyze & Center" -command ::PETK::gui::analyzeAnalyte -style "Accent.TButton"

    ttk::label  $container.input.sellbl -text "Atom selection:" -width 12
    ttk::entry  $container.input.sel -textvariable ::PETK::gui::analyteSelection -width 35
    ttk::label  $container.input.selhelp -text "Examples: 'all', 'protein', 'not water'" -font {TkDefaultFont 9 italic} -foreground gray

    grid $container.input.pdblbl $container.input.pdb $container.input.browse $container.input.analyze -sticky ew -pady 3
    grid $container.input.sellbl $container.input.sel $container.input.selhelp - -sticky ew -pady 3

    # === Analysis Results Frame ===
    ttk::labelframe $container.results -text "Analysis Results" -padding 10
    grid $container.results -column 0 -row $row -sticky nsew -pady "0 10" -padx 10
    grid columnconfigure $container.results {0 1} -weight 1
    grid rowconfigure $container.results 0 -weight 1
    incr row

    # Results display in two columns
    # Left column - Basic measurements
    ttk::frame $container.results.left
    grid $container.results.left -column 0 -row 0 -sticky nsew -padx "0 10"
    grid rowconfigure $container.results.left 6 -weight 1

    ttk::label $container.results.left.title -text "Molecular Dimensions" -font {TkDefaultFont 11 bold}
    grid $container.results.left.title -sticky w -pady "0 10"

    # Bounding sphere radius
    ttk::label $container.results.left.radiuslbl -text "Bounding radius:" -width 15
    ttk::label $container.results.left.radius -textvariable ::PETK::gui::analyteDiameter -width 15 -anchor w -relief sunken -background white
    grid $container.results.left.radiuslbl $container.results.left.radius -sticky ew -pady 2

    # Approximate volume
    ttk::label $container.results.left.vollbl -text "Approx. volume:" -width 15
    ttk::label $container.results.left.vol -textvariable ::PETK::gui::analyteVolume -width 15 -anchor w -relief sunken -background white
    grid $container.results.left.vollbl $container.results.left.vol -sticky ew -pady 2

    # Extreme atom distance
    ttk::label $container.results.left.distlbl -text "Max distance:" -width 15
    ttk::label $container.results.left.dist -textvariable ::PETK::gui::analyteDistance -width 15 -anchor w -relief sunken -background white
    grid $container.results.left.distlbl $container.results.left.dist -sticky ew -pady 2

    # Fit status with nanopore
    ttk::label $container.results.left.statuslbl -text "Pore fit status:" -width 15
    ttk::label $container.results.left.statusval -textvariable ::PETK::gui::fitStatus -width 15 -anchor w -relief sunken -background white -wraplength 120
    grid $container.results.left.statuslbl $container.results.left.statusval -sticky ew -pady 2

    # Right column - Centering information
    ttk::frame $container.results.right
    grid $container.results.right -column 1 -row 0 -sticky nsew -padx "10 0"
    grid rowconfigure $container.results.right 6 -weight 1

    ttk::label $container.results.right.title -text "Centering & Alignment" -font {TkDefaultFont 11 bold}
    grid $container.results.right.title -sticky w -pady "0 10"

    # Verification score
    ttk::label $container.results.right.scorelbl -text "Quality score:" -width 15
    ttk::label $container.results.right.score -textvariable ::PETK::gui::verificationScore -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.scorelbl $container.results.right.score -sticky ew -pady 2

    # Centering status
    ttk::label $container.results.right.centerlbl -text "Centering:" -width 15
    ttk::label $container.results.right.center -textvariable ::PETK::gui::centeringStatus -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.centerlbl $container.results.right.center -sticky ew -pady 2

    # Surface alignment status
    ttk::label $container.results.right.alignlbl -text "Surface align:" -width 15
    ttk::label $container.results.right.align -textvariable ::PETK::gui::alignmentStatus -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.alignlbl $container.results.right.align -sticky ew -pady 2

    # Output file status
    ttk::label $container.results.right.outputlbl -text "Centered PDB:" -width 15
    ttk::label $container.results.right.output -textvariable ::PETK::gui::outputFileStatus -width 15 -anchor w -relief sunken -background white -wraplength 120
    grid $container.results.right.outputlbl $container.results.right.output -sticky ew -pady 2

    # === Detailed Information Frame (Expandable) ===
    ttk::labelframe $container.details -text "Detailed Information" -padding 10
    grid $container.details -column 0 -row $row -sticky nsew -pady "0 10" -padx 10
    grid columnconfigure $container.details 0 -weight 1
    grid rowconfigure $container.details 1 -weight 1
    incr row

    # Toggle button for detailed view
    ttk::frame $container.details.header
    grid $container.details.header -sticky ew -pady "0 5"
    grid columnconfigure $container.details.header 0 -weight 1

    ttk::button $container.details.header.toggle -text "▼ Show Details" -command ::PETK::gui::toggleDetailView
    ttk::button $container.details.header.export -text "Export Report" -command ::PETK::gui::exportAnalysisReport
    grid $container.details.header.toggle -sticky w
    grid $container.details.header.export -sticky e -row 0 -column 1

    # Scrollable text area for detailed information
    ttk::frame $container.details.content
    text $container.details.content.text -height 8 -width 80 -wrap word -state disabled \
        -yscrollcommand [list $container.details.content.scroll set] -font {TkFixedFont 9}
    ttk::scrollbar $container.details.content.scroll -orient vertical -command [list $container.details.content.text yview]

    grid $container.details.content.text $container.details.content.scroll -sticky nsew
    grid columnconfigure $container.details.content 0 -weight 1
    grid rowconfigure $container.details.content 0 -weight 1

    # Store reference to text widget and initially hide details
    set ::PETK::gui::detailsTextWidget $container.details.content.text
    set ::PETK::gui::detailsVisible 0

    # === Visualization Controls Frame ===
    ttk::labelframe $container.controls -text "Visualization & Actions" -padding 10
    grid $container.controls -column 0 -row $row -sticky ew -pady "0 10" -padx 10
    grid columnconfigure $container.controls {0 1 2 3} -weight 1
    incr row

    ttk::button $container.controls.show -text "Show Molecule" -command ::PETK::gui::showAnalyte
    ttk::button $container.controls.hide -text "Hide Molecule" -command ::PETK::gui::hideAnalyte
    ttk::button $container.controls.center -text "Center View" -command ::PETK::gui::centerAnalyteView
    ttk::button $container.controls.representations -text "Change Rep" -command ::PETK::gui::cycleAnalyteRepresentation

    grid $container.controls.show $container.controls.hide $container.controls.center $container.controls.representations -sticky ew -padx 5

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab2.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab2.canvas $canvas_window]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab2.canvas $canvas_window]
    
    # Store references for later use (with unique names for tab2)
    set ::PETK::gui::tab2MainCanvas $tab2.canvas
    set ::PETK::gui::tab2CanvasWindow $canvas_window
    set ::PETK::gui::tab2ContentContainer $container
    
    # Initialize result variables
    ::PETK::gui::initializeResultVariables
    
    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab2.canvas $canvas_window]
}

####################################################
# Tab 2 Function
####################################################

proc ::PETK::gui::initializeResultVariables {} {
    if {![info exists ::PETK::gui::analyteSelection]} {
        set ::PETK::gui::analyteSelection "all"
    }

    # Initialize all result display variables
    if {![info exists ::PETK::gui::analyteDiameter]} {
        set ::PETK::gui::analyteDiameter "Not analyzed"
    }
    if {![info exists ::PETK::gui::analyteVolume]} {
        set ::PETK::gui::analyteVolume "Not analyzed"
    }
    if {![info exists ::PETK::gui::analyteDistance]} {
        set ::PETK::gui::analyteDistance "Not analyzed"
    }
    if {![info exists ::PETK::gui::fitStatus]} {
        set ::PETK::gui::fitStatus "Set pore diameter to check"
    }
    if {![info exists ::PETK::gui::verificationScore]} {
        set ::PETK::gui::verificationScore "Not analyzed"
    }
    if {![info exists ::PETK::gui::centeringStatus]} {
        set ::PETK::gui::centeringStatus "Not analyzed"
    }
    if {![info exists ::PETK::gui::alignmentStatus]} {
        set ::PETK::gui::alignmentStatus "Not analyzed"
    }
    if {![info exists ::PETK::gui::outputFileStatus]} {
        set ::PETK::gui::outputFileStatus "Not created"
    }
    if {![info exists ::PETK::gui::detailsVisible]} {
        set ::PETK::gui::detailsVisible 0
    }
    if {![info exists ::PETK::gui::currentRepresentation]} {
        set ::PETK::gui::currentRepresentation 0
    }
}

####################################################
## Analyte Input
####################################################

proc ::PETK::gui::browseAnalytePdb {} {
    set pdbType { {{Protein Data Bank files} {.pdb}} {{All Files} *} }
    set tempfile [tk_getOpenFile -title "Select Analyte PDB File" -multiple 0 -filetypes $pdbType]
    if {![string eq $tempfile ""]} {
        set ::PETK::gui::analytePDB $tempfile
    }
}

####################################################
### Enhanced Analysis Function with Better Status Updates
####################################################

proc ::PETK::gui::analyzeAnalyte {} {
    puts "Starting complete analyte analysis with centering and alignment..."
    
    # Clear previous results
    ::PETK::gui::clearAnalysisResults
    
    # Update status
    set ::PETK::gui::verificationScore "Analyzing..."
    set ::PETK::gui::centeringStatus "Processing..."
    set ::PETK::gui::alignmentStatus "Processing..."
    set ::PETK::gui::outputFileStatus "Processing..."
    update
    
    # Step 1: Center and align the molecule
    if {[catch {::PETK::gui::analyzeAnalyteWithCentering} error]} {
        tk_messageBox -icon error -message "Analysis Error" -detail "Error during analysis:\n$error"
        ::PETK::gui::clearAnalysisResults
        return
    }
    
    # Step 2: Verify the quality of centering and alignment
    puts "\nVerifying centering and alignment quality..."
    set verification_score [::PETK::gui::verifyCentering]
    
    # Update detailed information
    ::PETK::gui::updateDetailedInformation
    
    puts "\nAnalysis complete. Verification score: $verification_score/12"
    
    return $verification_score
}

proc ::PETK::gui::analyzeAnalyteWithCentering {} {
    # Initialize variables if they don't exist
    if {![info exists ::PETK::gui::analyteMol]} {
        set ::PETK::gui::analyteMol ""
    }
    if {![info exists ::PETK::gui::analytePDB]} {
        tk_messageBox -icon error -message "Error!" -detail "No PDB file specified"
        return
    }
    if {![info exists ::PETK::gui::analyteSelection]} {
        set ::PETK::gui::analyteSelection "all"
    }
    
    # Clear previous results
    set ::PETK::gui::analyteDiameter ""
    set ::PETK::gui::analyteVolume ""
    set ::PETK::gui::fitStatus ""
    set ::PETK::gui::analyteDistance ""
    set ::PETK::gui::leftmostAtom ""
    set ::PETK::gui::rightmostAtom ""

    # Check if PDB file exists
    if {![file exists $::PETK::gui::analytePDB]} {
        tk_messageBox -icon error -message "Error!" -detail "PDB file not found: $::PETK::gui::analytePDB"
        return
    }

    # Load molecule if not already loaded
    if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
        catch {mol delete $::PETK::gui::analyteMol}
    }

    if {[catch {mol new $::PETK::gui::analytePDB waitfor all} molid]} {
        tk_messageBox -icon error -message "Error!" -detail "Could not load PDB file: $::PETK::gui::analytePDB"
        return
    }

    set ::PETK::gui::analyteMol $molid
    mol rename $molid "Analyte: [file tail $::PETK::gui::analytePDB]"

    # Create selection with error handling
    if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
        tk_messageBox -icon error -message "Error!" -detail "Could not create atom selection: $sel"
        return
    }
    
    # Check if selection is valid and has atoms
    if {![info exists sel] || $sel == "" || [catch {$sel num} num_atoms]} {
        tk_messageBox -icon error -message "Error!" -detail "Invalid atom selection created"
        if {[info exists sel] && $sel != ""} {
            catch {$sel delete}
        }
        return
    }
    
    if {$num_atoms == 0} {
        tk_messageBox -icon warning -message "Warning!" -detail "Selection contains 0 atoms!"
        $sel delete
        return
    }

    puts "=================================="
    puts "ANALYTE PROCESSING: CENTERING AND ALIGNMENT"
    puts "=================================="
    puts "Processing $num_atoms atoms with selection: $::PETK::gui::analyteSelection"

    # STEP 1: CENTER THE MOLECULE AT ORIGIN
    puts "\nSTEP 1: CENTERING MOLECULE"
    puts "-------------------------"
    set initial_center [measure center $sel]
    puts [format "Initial center of mass: (%.6f, %.6f, %.6f)" \
        [lindex $initial_center 0] [lindex $initial_center 1] [lindex $initial_center 2]]
    
    set move_vector [vecscale -1.0 $initial_center]
    $sel moveby $move_vector
    
    set new_center [measure center $sel]
    puts [format "New center of mass: (%.6f, %.6f, %.6f)" \
        [lindex $new_center 0] [lindex $new_center 1] [lindex $new_center 2]]
    puts "✓ Molecule centered at origin"

    # STEP 2: ALIGN Z-AXIS PERPENDICULAR TO MOLECULAR SURFACE
    puts "\nSTEP 2: ALIGNING Z-AXIS PERPENDICULAR TO SURFACE"
    puts "------------------------------------------------"
    
    # Method: Use coordinate variance to find surface normal
    # CORRECTED LOGIC: Largest variance = surface direction, smallest variance = normal direction
    
    set coords [$sel get {x y z}]
    set x_coords {}
    set y_coords {}
    set z_coords {}
    
    foreach coord $coords {
        lappend x_coords [lindex $coord 0]
        lappend y_coords [lindex $coord 1]
        lappend z_coords [lindex $coord 2]
    }
    
    set x_var [::PETK::gui::calculateVariance $x_coords]
    set y_var [::PETK::gui::calculateVariance $y_coords]
    set z_var [::PETK::gui::calculateVariance $z_coords]
    
    puts "\nSurface analysis using coordinate variance method:"
    puts [format "  X variance: %.3f (spread in X direction)" $x_var]
    puts [format "  Y variance: %.3f (spread in Y direction)" $y_var]
    puts [format "  Z variance: %.3f (spread in Z direction)" $z_var]
    
    # Find the direction with SMALLEST variance (this is perpendicular to the surface)
    set min_var $x_var
    set surface_normal {1 0 0}
    set normal_direction "X"
    set surface_directions "YZ"
    
    if {$y_var < $min_var} {
        set min_var $y_var
        set surface_normal {0 1 0}
        set normal_direction "Y"
        set surface_directions "XZ"
    }
    
    if {$z_var < $min_var} {
        set min_var $z_var
        set surface_normal {0 0 1}
        set normal_direction "Z"
        set surface_directions "XY"
    }
    
    puts [format "  Smallest variance: %.3f in %s direction (surface normal)" $min_var $normal_direction]
    puts [format "  Surface lies in %s plane" $surface_directions]
    puts [format "  Surface normal vector: (%8.5f, %8.5f, %8.5f)" \
          [lindex $surface_normal 0] [lindex $surface_normal 1] [lindex $surface_normal 2]]
    
    # Calculate rotation matrix to align surface normal with Z-axis
    set target_normal {0 0 1}  ; # We want Z to be perpendicular to surface
    set R [::PETK::gui::calculateRotationMatrix $surface_normal $target_normal]
    
    puts "\nApplying surface alignment rotation:"
    for {set i 0} {$i < 3} {incr i} {
        set row [lindex $R $i]
        puts [format "  \[%8.5f %8.5f %8.5f\]" [lindex $row 0] [lindex $row 1] [lindex $row 2]]
    }
    
    # Apply rotation to all atoms
    set new_coords {}
    foreach coord $coords {
        set rotated_coord [::PETK::gui::matrix_vector_multiply_3x3 $R $coord]
        lappend new_coords $rotated_coord
    }
    
    # Set new coordinates
    $sel set {x y z} $new_coords
    puts "✓ Z-axis aligned perpendicular to molecular surface"
    
    # Verify the surface alignment
    puts "\nVerifying surface alignment:"
    set coords [$sel get {x y z}]
    set x_coords {}
    set y_coords {}
    set z_coords {}
    
    foreach coord $coords {
        lappend x_coords [lindex $coord 0]
        lappend y_coords [lindex $coord 1]
        lappend z_coords [lindex $coord 2]
    }
    
    set x_var_new [::PETK::gui::calculateVariance $x_coords]
    set y_var_new [::PETK::gui::calculateVariance $y_coords]
    set z_var_new [::PETK::gui::calculateVariance $z_coords]
    
    puts [format "  New coordinate variances: X=%.3f, Y=%.3f, Z=%.3f" $x_var_new $y_var_new $z_var_new]
    
    if {$z_var_new < $x_var_new && $z_var_new < $y_var_new} {
        puts "✓ SUCCESS: Z-axis is now perpendicular to the flattest surface"
        puts "  Surface is now in the XY plane, Z points perpendicular to surface"
    } elseif {$z_var_new < ($x_var_new + $y_var_new) / 2.0} {
        puts "✓ GOOD: Z-axis alignment improved significantly"
    } else {
        puts "⚠ WARNING: Surface orientation may need manual adjustment"
        puts "  Consider checking molecule geometry or selection"
    }

    # STEP 3: ANALYZE GEOMETRIC PROPERTIES
    puts "\nSTEP 3: ANALYZING GEOMETRIC PROPERTIES"
    puts "--------------------------------------"
    
    # Refresh coordinates after alignment
    set coords [$sel get {x y z}]
    set atomNames [$sel get name]
    set resNames [$sel get resname]
    set resIds [$sel get resid]
    
    # Find extreme atoms along each axis
    set x_coords {}
    set y_coords {}
    set z_coords {}
    
    foreach coord $coords {
        lappend x_coords [lindex $coord 0]
        lappend y_coords [lindex $coord 1]
        lappend z_coords [lindex $coord 2]
    }
    
    # Find leftmost and rightmost atoms based on x-coordinate
    set leftmostX [lindex $x_coords 0]
    set rightmostX $leftmostX
    set leftmostIndex 0
    set rightmostIndex 0
    
    for {set i 0} {$i < [llength $coords]} {incr i} {
        set x [lindex $x_coords $i]
        if {$x < $leftmostX} {
            set leftmostX $x
            set leftmostIndex $i
        }
        if {$x > $rightmostX} {
            set rightmostX $x
            set rightmostIndex $i
        }
    }
    
    # Get extreme atom information
    set leftmostCoord [lindex $coords $leftmostIndex]
    set rightmostCoord [lindex $coords $rightmostIndex]
    set leftmostName [lindex $atomNames $leftmostIndex]
    set rightmostName [lindex $atomNames $rightmostIndex]
    set leftmostRes [lindex $resNames $leftmostIndex]
    set rightmostRes [lindex $resNames $rightmostIndex]
    set leftmostResId [lindex $resIds $leftmostIndex]
    set rightmostResId [lindex $resIds $rightmostIndex]
    
    # Calculate distance between extreme atoms
    set dx [expr {[lindex $rightmostCoord 0] - [lindex $leftmostCoord 0]}]
    set dy [expr {[lindex $rightmostCoord 1] - [lindex $leftmostCoord 1]}]
    set dz [expr {[lindex $rightmostCoord 2] - [lindex $leftmostCoord 2]}]
    set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
    
    # Calculate coordinate ranges
    set x_range [expr {[::tcl::mathfunc::max {*}$x_coords] - [::tcl::mathfunc::min {*}$x_coords]}]
    set y_range [expr {[::tcl::mathfunc::max {*}$y_coords] - [::tcl::mathfunc::min {*}$y_coords]}]
    set z_range [expr {[::tcl::mathfunc::max {*}$z_coords] - [::tcl::mathfunc::min {*}$z_coords]}]
    
    puts [format "Molecular dimensions after alignment:"]
    puts [format "  X dimension: %.3f Å (surface width)" $x_range]
    puts [format "  Y dimension: %.3f Å (surface length)" $y_range]
    puts [format "  Z dimension: %.3f Å (surface thickness)" $z_range]
    
    # Store extreme atom information
    set ::PETK::gui::leftmostAtom [format "%s %s%d (%.3f, %.3f, %.3f)" \
        $leftmostName $leftmostRes $leftmostResId \
        [lindex $leftmostCoord 0] [lindex $leftmostCoord 1] [lindex $leftmostCoord 2]]
    
    set ::PETK::gui::rightmostAtom [format "%s %s%d (%.3f, %.3f, %.3f)" \
        $rightmostName $rightmostRes $rightmostResId \
        [lindex $rightmostCoord 0] [lindex $rightmostCoord 1] [lindex $rightmostCoord 2]]
    
    set ::PETK::gui::analyteDistance [format "%.3f Å" $distance]

    # Calculate bounding sphere
    set final_center [measure center $sel]
    set maxDist 0.0
    foreach coord $coords {
        set dist [veclength [vecsub $coord $final_center]]
        if {$dist > $maxDist} {
            set maxDist $dist
        }
    }

    # Calculate volume (approximate as sphere)
    set volume [expr {4.0/3.0 * 3.14159 * pow($maxDist, 3)}]

    # Store results
    set ::PETK::gui::analyteDiameter [format "%.2f Å" $maxDist]
    set ::PETK::gui::analyteVolume [format "%.2f Å³" $volume]

    # Check fit with nanopore
    if {$::PETK::gui::poreDiameter != ""} {
        set poreRad [expr {double($::PETK::gui::poreDiameter)}]
        if {$distance < $poreRad} {
            set ::PETK::gui::fitStatus "FITS (clearance: [format "%.2f" [expr {$poreRad - $distance}]] Å)"
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval]} {
                $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval configure -foreground green
            }
        } else {
            set ::PETK::gui::fitStatus "TOO LARGE (excess: [format "%.2f" [expr {$distance - $poreRad}]] Å)"
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval]} {
                $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval configure -foreground red
            }
        }
    } else {
        set ::PETK::gui::fitStatus "Set pore Diameter to check fit"
        if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval]} {
            $::PETK::gui::window.hlf.nb.tab1.analyte.results.statusval configure -foreground black
        }
    }

    $sel delete

    # Print summary
    puts "\nPROCESSING COMPLETE:"
    puts "==================="
    puts "✓ Molecule centered at origin"
    puts "✓ Z-axis aligned perpendicular to molecular surface"
    puts "✓ Surface lies in XY plane"
    puts [format "Surface dimensions: %.1f × %.1f Å (XY), thickness: %.1f Å (Z)" $x_range $y_range $z_range]
    puts "Distance between extreme atoms: $::PETK::gui::analyteDistance"
    puts "Bounding sphere Diameter: $::PETK::gui::analyteDiameter"
    puts "Approximate volume: $::PETK::gui::analyteVolume"
    puts "=================================="

    # Auto-show the molecule
    ::PETK::gui::showAnalyte
    
    # Set a good viewing angle
    if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
        mol top $::PETK::gui::analyteMol
        display resetview
        scale by 1.2
    }
}

proc ::PETK::gui::showAnalyte {} {
    if {$::PETK::gui::analyteMol != ""} {
        mol on $::PETK::gui::analyteMol
        mol top $::PETK::gui::analyteMol
        display resetview
        scale by 0.8
    }
}

proc ::PETK::gui::hideAnalyte {} {
    if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
        catch {mol off $::PETK::gui::analyteMol}
    }
}

proc ::PETK::gui::centerAnalyteView {} {
    if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
        mol top $::PETK::gui::analyteMol
        display resetview
        scale by 1.0
    }
}

proc ::PETK::gui::cycleAnalyteRepresentation {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        return
    }
    
    if {![info exists ::PETK::gui::currentRepresentation]} {
        set ::PETK::gui::currentRepresentation 0
    }
    
    set representations {"Lines" "VDW" "CPK" "Licorice" "NewCartoon"}
    set rep_names {"Lines" "VDW (spheres)" "CPK" "Licorice (sticks)" "Cartoon"}
    
    set ::PETK::gui::currentRepresentation [expr {($::PETK::gui::currentRepresentation + 1) % [llength $representations]}]
    set new_rep [lindex $representations $::PETK::gui::currentRepresentation]
    set rep_name [lindex $rep_names $::PETK::gui::currentRepresentation]
    
    # Update the representation
    catch {
        mol modstyle 0 $::PETK::gui::analyteMol $new_rep
        puts "Changed representation to: $rep_name"
    }
    
    # Update button text to show current representation
    if {[winfo exists $::PETK::gui::window.hlf.nb.tab1.controls.representations]} {
        $::PETK::gui::window.hlf.nb.tab1.controls.representations configure -text "Rep: $rep_name"
    }
}


####################################################
### Enhanced Verification Function with Status Updates
####################################################

proc ::PETK::gui::verifyCentering {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        puts "No molecule loaded for verification"
        return 0
    }
    
    set molid $::PETK::gui::analyteMol
    
    # Check if molecule exists
    if {[lsearch [molinfo list] $molid] == -1} {
        puts "Error: Molecule $molid does not exist"
        return 0
    }
    
    # Initialize analyteSelection if it doesn't exist
    if {![info exists ::PETK::gui::analyteSelection]} {
        set ::PETK::gui::analyteSelection "all"
    }
    
    # Create selection with error handling
    if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
        puts "Error creating selection: $sel"
        return 0
    }
    
    # Check if selection is valid and has atoms
    if {![info exists sel] || $sel == "" || [catch {$sel num} num_atoms]} {
        puts "Error: Invalid atom selection created"
        if {[info exists sel] && $sel != ""} {
            catch {$sel delete}
        }
        return 0
    }
    
    if {$num_atoms == 0} {
        puts "No atoms in selection for verification"
        $sel delete
        return 0
    }
    
    puts "\n=========================================="
    puts "CENTERING AND ALIGNMENT VERIFICATION"
    puts "=========================================="
    puts "Molecule ID: $molid"
    
    # Get molecule name safely
    if {[catch {molinfo $molid get name} mol_name]} {
        set mol_name "Unknown"
    }
    puts "Molecule name: $mol_name"
    puts "Selection: $::PETK::gui::analyteSelection"
    puts "Number of atoms: $num_atoms"
    puts ""
    
    # 1. CENTER OF MASS VERIFICATION
    puts "1. CENTER OF MASS VERIFICATION:"
    puts "------------------------------"
    
    set com [measure center $sel]
    set com_x [lindex $com 0]
    set com_y [lindex $com 1]
    set com_z [lindex $com 2]
    
    puts [format "Center of Mass: (%.8f, %.8f, %.8f)" $com_x $com_y $com_z]
    
    set com_distance [veclength $com]
    puts [format "Distance from origin: %.8f Å" $com_distance]
    
    # Centering quality assessment
    if {$com_distance < 1e-6} {
        puts "✓ EXCELLENT: Perfect centering (< 1e-6 Å)"
        set centering_score 4
        set ::PETK::gui::centeringStatus "Excellent"
    } elseif {$com_distance < 1e-4} {
        puts "✓ EXCELLENT: Near-perfect centering (< 1e-4 Å)"
        set centering_score 3
        set ::PETK::gui::centeringStatus "Excellent"
    } elseif {$com_distance < 1e-2} {
        puts "✓ GOOD: Good centering (< 0.01 Å)"
        set centering_score 2
        set ::PETK::gui::centeringStatus "Good"
    } elseif {$com_distance < 0.1} {
        puts "⚠ FAIR: Reasonable centering (< 0.1 Å)"
        set centering_score 1
        set ::PETK::gui::centeringStatus "Fair"
    } else {
        puts "✗ POOR: Poor centering (> 0.1 Å)"
        set centering_score 0
        set ::PETK::gui::centeringStatus "Poor"
    }
    
    # 2. SURFACE ALIGNMENT VERIFICATION
    puts "\n2. SURFACE ALIGNMENT VERIFICATION:"
    puts "----------------------------------"
    
    # Check coordinate variances to verify surface alignment
    set coords [$sel get {x y z}]
    set x_coords {}
    set y_coords {}
    set z_coords {}
    
    foreach coord $coords {
        lappend x_coords [lindex $coord 0]
        lappend y_coords [lindex $coord 1]
        lappend z_coords [lindex $coord 2]
    }
    
    set x_var [::PETK::gui::calculateVariance $x_coords]
    set y_var [::PETK::gui::calculateVariance $y_coords]
    set z_var [::PETK::gui::calculateVariance $z_coords]
    
    puts [format "Coordinate variances:"]
    puts [format "  X variance: %.3f" $x_var]
    puts [format "  Y variance: %.3f" $y_var]
    puts [format "  Z variance: %.3f" $z_var]
    
    # Surface alignment scoring
    if {$z_var < $x_var && $z_var < $y_var} {
        set surface_ratio [expr {$z_var / (($x_var + $y_var) / 2.0)}]
        
        if {$surface_ratio < 0.1} {
            puts "✓ EXCELLENT: Z-axis perfectly perpendicular to surface"
            set surface_score 4
            set ::PETK::gui::alignmentStatus "Excellent"
        } elseif {$surface_ratio < 0.3} {
            puts "✓ GOOD: Z-axis well aligned perpendicular to surface"
            set surface_score 3
            set ::PETK::gui::alignmentStatus "Good"
        } elseif {$surface_ratio < 0.5} {
            puts "✓ FAIR: Z-axis reasonably aligned"
            set surface_score 2
            set ::PETK::gui::alignmentStatus "Fair"
        } else {
            puts "⚠ POOR: Z-axis alignment needs improvement"
            set surface_score 1
            set ::PETK::gui::alignmentStatus "Poor"
        }
        puts [format "  Surface flatness ratio: %.3f" $surface_ratio]
    } else {
        puts "✗ POOR: Z-axis is not the flattest dimension"
        set surface_score 0
        set ::PETK::gui::alignmentStatus "Poor"
    }
    
    # 3. PRINCIPAL AXES ALIGNMENT VERIFICATION
    puts "\n3. PRINCIPAL AXES ANALYSIS:"
    puts "---------------------------"

    # Get current inertia information
    set inertia_result [measure inertia $sel moments]
    set principal_axes [lindex $inertia_result 1]
    set moments_tensor [lindex $inertia_result 2]
    set moments [lindex $inertia_result 0]
    
    # Extract principal axes
    set axis1 [lindex $principal_axes 0]
    set axis2 [lindex $principal_axes 1]
    set axis3 [lindex $principal_axes 2]

    puts "Current Principal Axes:"
    puts [format "  Axis 1: (%8.6f, %8.6f, %8.6f) - Moment: %.1f" \
          [lindex $axis1 0] [lindex $axis1 1] [lindex $axis1 2] [lindex $moments 0]]
    puts [format "  Axis 2: (%8.6f, %8.6f, %8.6f) - Moment: %.1f" \
          [lindex $axis2 0] [lindex $axis2 1] [lindex $axis2 2] [lindex $moments 1]]
    puts [format "  Axis 3: (%8.6f, %8.6f, %8.6f) - Moment: %.1f" \
          [lindex $axis3 0] [lindex $axis3 1] [lindex $axis3 2] [lindex $moments 2]]

    # Check individual axis alignment with coordinate axes
    puts "\nAxis Alignment Check:"
    set coord_axes {{1 0 0} {0 1 0} {0 0 1}}
    set axis_names {"X" "Y" "Z"}
    
    set total_alignment 0.0
    foreach principal_axis [list $axis1 $axis2 $axis3] axis_num {1 2 3} {
        set best_alignment 0.0
        set best_coord_axis ""
        
        foreach coord_axis $coord_axes coord_name $axis_names {
            set dot_product [expr {[lindex $principal_axis 0] * [lindex $coord_axis 0] + \
                                   [lindex $principal_axis 1] * [lindex $coord_axis 1] + \
                                   [lindex $principal_axis 2] * [lindex $coord_axis 2]}]
            set alignment [expr {abs($dot_product)}]
            
            if {$alignment > $best_alignment} {
                set best_alignment $alignment
                set best_coord_axis $coord_name
            }
        }
        
        set total_alignment [expr {$total_alignment + $best_alignment}]
        puts [format "  Principal axis %d → %s-axis: %.6f" \
              $axis_num $best_coord_axis $best_alignment]
    }
    
    set avg_alignment [expr {$total_alignment / 3.0}]
    puts [format "  Average alignment: %.6f" $avg_alignment]
    
    # Alignment scoring
    if {$avg_alignment > 0.99} {
        puts "✓ EXCELLENT: Principal axes very well aligned"
        set alignment_score 4
    } elseif {$avg_alignment > 0.95} {
        puts "✓ GOOD: Principal axes well aligned"
        set alignment_score 3
    } elseif {$avg_alignment > 0.8} {
        puts "✓ FAIR: Principal axes reasonably aligned"
        set alignment_score 2
    } else {
        puts "⚠ POOR: Principal axes poorly aligned"
        set alignment_score 1
    }
    
    # 4. CALCULATE TOTAL SCORE AND UPDATE GUI
    set max_score 12  ; # 4 for centering + 4 for surface + 4 for alignment
    set total_score [expr {$centering_score + $surface_score + $alignment_score}]
    set percentage [expr {($total_score * 100.0) / $max_score}]
    
    # Update GUI with scores
    set ::PETK::gui::verificationScore [format "%d/12 (%.0f%%)" $total_score $percentage]
    
    puts "Quality Scores:"
    puts [format "  Centering: %d/4 (%s)" $centering_score $::PETK::gui::centeringStatus]
    puts [format "  Surface alignment: %d/4 (%s)" $surface_score $::PETK::gui::alignmentStatus]
    puts [format "  Principal axes: %d/4" $alignment_score]
    puts [format "  Total: %d/%d (%.1f%%)" $total_score $max_score $percentage]
    
    # 5. SAVE CENTERED PDB IF VERIFICATION PASSED
    set verification_passed 0
    if {$percentage >= 75.0} {
        puts "\n✓ OVERALL: EXCELLENT/VERY GOOD - Verification PASSED"
        set verification_passed 1
    } elseif {$percentage >= 50.0} {
        puts "\n✓ OVERALL: GOOD - Verification PASSED"
        set verification_passed 1
    } elseif {$percentage >= 33.0} {
        puts "\n⚠ OVERALL: FAIR - Verification MARGINALLY PASSED"
        set verification_passed 1
    } else {
        puts "\n✗ OVERALL: POOR - Verification FAILED"
        set verification_passed 0
    }
    
    # Save centered PDB if verification passed
    if {$verification_passed} {
        puts "\n6. SAVING CENTERED PDB:"
        puts "----------------------"
        
        # Determine output path and filename
        set original_filename [file tail $::PETK::gui::analytePDB]
        set name_without_ext [file rootname $original_filename]
        set output_filename "centered_${name_without_ext}.pdb"
        set ::PETK::gui::semAnalytePDB "centered_${name_without_ext}.pdb"
        # Use work directory if it exists, otherwise use same directory as input
        if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir != ""} {
            set output_path [file join $::PETK::gui::workdir $output_filename]
        } else {
            set input_dir [file dirname $::PETK::gui::analytePDB]
            set output_path [file join $input_dir $output_filename]
        }
        
        puts "Original file: $::PETK::gui::analytePDB"
        puts "Output file: $output_path"
        
        # Save the entire molecule (all atoms) with new coordinates
        set all_atoms [atomselect $molid "all"]
        
        if {[catch {$all_atoms writepdb $output_path} error]} {
            puts "✗ ERROR: Failed to save PDB file - $error"
            set ::PETK::gui::outputFileStatus "Save failed"
        } else {
            puts "✓ SUCCESS: Centered PDB saved successfully"
            puts "  File saved: $output_path"
            puts "  Contains: [$all_atoms num] atoms"
            
            set ::PETK::gui::outputFileStatus "Saved: [file tail $output_path]"
            
            # Verify the saved file exists and has content
            if {[file exists $output_path]} {
                set file_size [file size $output_path]
                puts [format "  File size: %.1f KB" [expr {$file_size / 1024.0}]]
            }
        }
        
        $all_atoms delete
        
    } else {
        puts "\n⚠ PDB NOT SAVED: Verification failed"
        puts "  Improve centering/alignment before saving"
        set ::PETK::gui::outputFileStatus "Not saved (low quality)"
    }
    
    puts "=========================================="
    
    $sel delete
    
    return $total_score
}

####################################################
### Information Display and Control Functions
####################################################

proc ::PETK::gui::clearAnalysisResults {} {
    set ::PETK::gui::analyteDiameter "Not analyzed"
    set ::PETK::gui::analyteVolume "Not analyzed"
    set ::PETK::gui::analyteDistance "Not analyzed"
    set ::PETK::gui::fitStatus "Set pore diameter to check"
    set ::PETK::gui::verificationScore "Not analyzed"
    set ::PETK::gui::centeringStatus "Not analyzed"
    set ::PETK::gui::alignmentStatus "Not analyzed"
    set ::PETK::gui::outputFileStatus "Not created"
    
    # Clear detailed information
    if {[info exists ::PETK::gui::detailsTextWidget] && [winfo exists $::PETK::gui::detailsTextWidget]} {
        $::PETK::gui::detailsTextWidget configure -state normal
        $::PETK::gui::detailsTextWidget delete 1.0 end
        $::PETK::gui::detailsTextWidget configure -state disabled
    }
}

proc ::PETK::gui::updateFitStatus {} {
    if {![info exists ::PETK::gui::poreDiameter] || $::PETK::gui::poreDiameter eq ""} {
        set ::PETK::gui::fitStatus "Set pore diameter to check"
        return
    }
    
    if {![info exists ::PETK::gui::analyteDiameter] || $::PETK::gui::analyteDiameter eq "Not analyzed"} {
        set ::PETK::gui::fitStatus "Analyze molecule first"
        return
    }
    
    # Extract numeric value from analyteDiameter (remove " Å")
    if {[regexp {^([\d.]+)} $::PETK::gui::analyteDiameter match radius]} {
        set poreRad [expr {double($::PETK::gui::poreDiameter) / 2.0}]
        set molRad [expr {double($radius)}]
        
        if {$molRad < $poreRad} {
            set clearance [expr {$poreRad - $molRad}]
            set ::PETK::gui::fitStatus "FITS (clearance: [format "%.1f" $clearance] Å)"
            # Update color if widget exists
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab1.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab1.results.left.statusval configure -foreground green
            }
        } else {
            set excess [expr {$molRad - $poreRad}]
            set ::PETK::gui::fitStatus "TOO LARGE (excess: [format "%.1f" $excess] Å)"
            # Update color if widget exists
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab1.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab1.results.left.statusval configure -foreground red
            }
        }
    }
}


proc ::PETK::gui::toggleDetailView {} {
    # Initialize detailsVisible if it doesn't exist
    if {![info exists ::PETK::gui::detailsVisible]} {
        set ::PETK::gui::detailsVisible 0
    }
    
    # Toggle the visibility state
    set ::PETK::gui::detailsVisible [expr {!$::PETK::gui::detailsVisible}]
    
    # Get references to the widgets with correct paths
    set button .petk_main_window.hlf.nb.tab2.canvas.content.details.header.toggle
    set content .petk_main_window.hlf.nb.tab2.canvas.content.details.content
    
    if {$::PETK::gui::detailsVisible} {
        # Show the details content
        grid $content -sticky nsew -pady "5 0"
        $button configure -text "▲ Hide Details"
        # Update the detailed information
        ::PETK::gui::updateDetailedInformation
    } else {
        # Hide the details content
        grid forget $content
        $button configure -text "▼ Show Details"
    }
}

proc ::PETK::gui::updateDetailedInformation {} {
    if {![info exists ::PETK::gui::detailsTextWidget] || ![winfo exists $::PETK::gui::detailsTextWidget]} {
        return
    }
    
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        return
    }
    
    set widget $::PETK::gui::detailsTextWidget
    $widget configure -state normal
    $widget delete 1.0 end
    
    # Add detailed analysis information
    $widget insert end "DETAILED ANALYSIS REPORT\n" {title}
    $widget insert end [string repeat "=" 50] 
    $widget insert end "\n\n"
    
    # File information
    $widget insert end "Input File Information:\n" {header}
    if {[info exists ::PETK::gui::analytePDB]} {
        $widget insert end "  File: [file tail $::PETK::gui::analytePDB]\n"
        $widget insert end "  Path: $::PETK::gui::analytePDB\n"
    }
    if {[info exists ::PETK::gui::analyteSelection]} {
        $widget insert end "  Selection: $::PETK::gui::analyteSelection\n"
    }
    
    # Get molecular information
    if {[catch {
        set molid $::PETK::gui::analyteMol
        set sel [atomselect $molid $::PETK::gui::analyteSelection]
        set num_atoms [$sel num]
        set coords [$sel get {x y z}]
        
        $widget insert end "  Number of atoms: $num_atoms\n\n"
        
        # Molecular dimensions
        $widget insert end "Molecular Dimensions:\n" {header}
        
        set x_coords {}
        set y_coords {}
        set z_coords {}
        
        foreach coord $coords {
            lappend x_coords [lindex $coord 0]
            lappend y_coords [lindex $coord 1]
            lappend z_coords [lindex $coord 2]
        }
        
        set x_min [::tcl::mathfunc::min {*}$x_coords]
        set x_max [::tcl::mathfunc::max {*}$x_coords]
        set y_min [::tcl::mathfunc::min {*}$y_coords]
        set y_max [::tcl::mathfunc::max {*}$y_coords]
        set z_min [::tcl::mathfunc::min {*}$z_coords]
        set z_max [::tcl::mathfunc::max {*}$z_coords]
        
        set x_range [expr {$x_max - $x_min}]
        set y_range [expr {$y_max - $y_min}]
        set z_range [expr {$z_max - $z_min}]
        
        $widget insert end [format "  X dimension: %.2f Å (%.2f to %.2f)\n" $x_range $x_min $x_max]
        $widget insert end [format "  Y dimension: %.2f Å (%.2f to %.2f)\n" $y_range $y_min $y_max]
        $widget insert end [format "  Z dimension: %.2f Å (%.2f to %.2f)\n" $z_range $z_min $z_max]
        
        # Center of mass
        set com [measure center $sel]
        $widget insert end "\nCenter of Mass:\n" {header}
        $widget insert end [format "  Position: (%.6f, %.6f, %.6f)\n" \
            [lindex $com 0] [lindex $com 1] [lindex $com 2]]
        $widget insert end [format "  Distance from origin: %.6f Å\n" [veclength $com]]
        
        # Coordinate variances (surface analysis)
        $widget insert end "\nSurface Analysis:\n" {header}
        set x_var [::PETK::gui::calculateVariance $x_coords]
        set y_var [::PETK::gui::calculateVariance $y_coords]
        set z_var [::PETK::gui::calculateVariance $z_coords]
        
        $widget insert end [format "  X variance: %.3f\n" $x_var]
        $widget insert end [format "  Y variance: %.3f\n" $y_var]
        $widget insert end [format "  Z variance: %.3f\n" $z_var]
        
        # Determine surface orientation
        set min_var $x_var
        set surface_normal "X"
        if {$y_var < $min_var} {
            set min_var $y_var
            set surface_normal "Y"
        }
        if {$z_var < $min_var} {
            set min_var $z_var
            set surface_normal "Z"
        }
        
        $widget insert end "  Surface normal direction: $surface_normal\n"
        
        # Principal axes information
        $widget insert end "\nPrincipal Axes Analysis:\n" {header}
        set inertia_result [measure inertia $sel moments]
        set principal_axes [lindex $inertia_result 1]
        set moments [lindex $inertia_result 0]
        
        for {set i 0} {$i < 3} {incr i} {
            set axis [lindex $principal_axes $i]
            set moment [lindex $moments $i]
            $widget insert end [format "  Axis %d: (%.6f, %.6f, %.6f) - Moment: %.1f\n" \
                [expr {$i+1}] [lindex $axis 0] [lindex $axis 1] [lindex $axis 2] $moment]
        }
        
        # Fit analysis
        if {[info exists ::PETK::gui::poreDiameter] && $::PETK::gui::poreDiameter ne ""} {
            $widget insert end "\nNanopore Fit Analysis:\n" {header}
            $widget insert end "  Target pore diameter: $::PETK::gui::poreDiameter Å\n"
            $widget insert end "  Target pore radius: [expr {$::PETK::gui::poreDiameter / 2.0}] Å\n"
            if {[info exists ::PETK::gui::analyteDiameter]} {
                $widget insert end "  Molecule bounding radius: $::PETK::gui::analyteDiameter\n"
            }
            $widget insert end "  Fit status: $::PETK::gui::fitStatus\n"
        }
        
        $sel delete
        
    } error]} {
        $widget insert end "Error retrieving detailed information: $error\n"
    }
    
    # Configure text tags for formatting
    $widget tag configure title -font {TkDefaultFont 12 bold}
    $widget tag configure header -font {TkDefaultFont 10 bold}
    
    $widget configure -state disabled
}

proc ::PETK::gui::exportAnalysisReport {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        tk_messageBox -icon warning -message "No Analysis" -detail "No molecule analysis to export."
        return
    }
    
    # Get output filename
    set filename [tk_getSaveFile -defaultextension ".txt" \
        -filetypes {{"Text files" ".txt"} {"All files" "*"}} \
        -title "Export Analysis Report"]
    
    if {$filename eq ""} {
        return
    }
    
    # Generate report content
    if {[catch {
        set report_file [open $filename w]
        
        puts $report_file "NANOPORE ANALYTE ANALYSIS REPORT"
        puts $report_file [string repeat "=" 60]
        puts $report_file "Generated: [clock format [clock seconds]]"
        puts $report_file ""
        
        # Basic information
        puts $report_file "INPUT INFORMATION:"
        puts $report_file [string repeat "-" 20]
        if {[info exists ::PETK::gui::analytePDB]} {
            puts $report_file "PDB File: $::PETK::gui::analytePDB"
        }
        if {[info exists ::PETK::gui::analyteSelection]} {
            puts $report_file "Selection: $::PETK::gui::analyteSelection"
        }
        puts $report_file ""
        
        # Analysis results
        puts $report_file "ANALYSIS RESULTS:"
        puts $report_file [string repeat "-" 20]
        puts $report_file "Bounding Radius: $::PETK::gui::analyteDiameter"
        puts $report_file "Approximate Volume: $::PETK::gui::analyteVolume"
        puts $report_file "Max Distance: $::PETK::gui::analyteDistance"
        puts $report_file "Verification Score: $::PETK::gui::verificationScore"
        puts $report_file "Centering Quality: $::PETK::gui::centeringStatus"
        puts $report_file "Alignment Quality: $::PETK::gui::alignmentStatus"
        puts $report_file "Output File Status: $::PETK::gui::outputFileStatus"
        puts $report_file ""
        
        # Pore fit analysis
        if {[info exists ::PETK::gui::poreDiameter] && $::PETK::gui::poreDiameter ne ""} {
            puts $report_file "NANOPORE FIT ANALYSIS:"
            puts $report_file [string repeat "-" 25]
            puts $report_file "Target Pore Diameter: $::PETK::gui::poreDiameter Å"
            puts $report_file "Fit Status: $::PETK::gui::fitStatus"
            puts $report_file ""
        }
        
        # Detailed molecular information
        if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
            set molid $::PETK::gui::analyteMol
            set sel [atomselect $molid $::PETK::gui::analyteSelection]
            
            puts $report_file "DETAILED MOLECULAR INFORMATION:"
            puts $report_file [string repeat "-" 35]
            puts $report_file "Number of atoms: [$sel num]"
            
            set coords [$sel get {x y z}]
            set x_coords {}
            set y_coords {}
            set z_coords {}
            
            foreach coord $coords {
                lappend x_coords [lindex $coord 0]
                lappend y_coords [lindex $coord 1]
                lappend z_coords [lindex $coord 2]
            }
            
            set x_min [::tcl::mathfunc::min {*}$x_coords]
            set x_max [::tcl::mathfunc::max {*}$x_coords]
            set y_min [::tcl::mathfunc::min {*}$y_coords]
            set y_max [::tcl::mathfunc::max {*}$y_coords]
            set z_min [::tcl::mathfunc::min {*}$z_coords]
            set z_max [::tcl::mathfunc::max {*}$z_coords]
            
            puts $report_file [format "X range: %.3f to %.3f Å (span: %.3f Å)" $x_min $x_max [expr {$x_max - $x_min}]]
            puts $report_file [format "Y range: %.3f to %.3f Å (span: %.3f Å)" $y_min $y_max [expr {$y_max - $y_min}]]
            puts $report_file [format "Z range: %.3f to %.3f Å (span: %.3f Å)" $z_min $z_max [expr {$z_max - $z_min}]]
            
            set com [measure center $sel]
            puts $report_file [format "Center of mass: (%.6f, %.6f, %.6f)" \
                [lindex $com 0] [lindex $com 1] [lindex $com 2]]
            puts $report_file [format "Distance from origin: %.6f Å" [veclength $com]]
            
            $sel delete
        }
        
        close $report_file
        
        tk_messageBox -icon info -message "Export Successful" \
            -detail "Analysis report exported to:\n$filename"
        
    } error]} {
        tk_messageBox -icon error -message "Export Failed" \
            -detail "Error exporting report:\n$error"
        if {[info exists report_file]} {
            catch {close $report_file}
        }
    }
}
####################################################
### Helper functions for matrix operations and surface alignment
####################################################
proc ::PETK::gui::calculateVariance {values} {
    set n [llength $values]
    if {$n == 0} {return 0}
    
    # Calculate mean
    set sum 0.0
    foreach val $values {
        set sum [expr {$sum + $val}]
    }
    set mean [expr {$sum / $n}]
    
    # Calculate variance
    set var_sum 0.0
    foreach val $values {
        set diff [expr {$val - $mean}]
        set var_sum [expr {$var_sum + $diff * $diff}]
    }
    
    return [expr {$var_sum / $n}]
}

proc ::PETK::gui::calculateRotationMatrix {from_vec to_vec} {
    # Normalize vectors
    set from_norm [::PETK::gui::normalizeVector $from_vec]
    set to_norm [::PETK::gui::normalizeVector $to_vec]
    
    # Calculate cross product (rotation axis)
    set cross [::PETK::gui::crossProduct $from_norm $to_norm]
    
    # Calculate dot product (cosine of angle)
    set dot [::PETK::gui::dotProduct $from_norm $to_norm]
    
    # Handle special cases
    if {$dot > 0.999999} {
        # Vectors are already aligned
        return {{1 0 0} {0 1 0} {0 0 1}}
    }
    
    if {$dot < -0.999999} {
        # Vectors are opposite - need 180 degree rotation
        # Find a perpendicular vector
        set perp [::PETK::gui::findPerpendicularVector $from_norm]
        return [::PETK::gui::rotationMatrix180 $perp]
    }
    
    # General case - use Rodrigues' rotation formula
    set angle [expr {acos($dot)}]
    set axis_norm [::PETK::gui::normalizeVector $cross]
    
    return [::PETK::gui::rodriguesRotation $axis_norm $angle]
}

proc ::PETK::gui::normalizeVector {vec} {
    set x [lindex $vec 0]
    set y [lindex $vec 1]
    set z [lindex $vec 2]
    
    set length [expr {sqrt($x*$x + $y*$y + $z*$z)}]
    
    if {$length < 1e-10} {
        return {1 0 0}  ; # Default to x-axis if zero vector
    }
    
    return [list [expr {$x/$length}] [expr {$y/$length}] [expr {$z/$length}]]
}

proc ::PETK::gui::crossProduct {a b} {
    set ax [lindex $a 0]; set ay [lindex $a 1]; set az [lindex $a 2]
    set bx [lindex $b 0]; set by [lindex $b 1]; set bz [lindex $b 2]
    
    set cx [expr {$ay*$bz - $az*$by}]
    set cy [expr {$az*$bx - $ax*$bz}]
    set cz [expr {$ax*$by - $ay*$bx}]
    
    return [list $cx $cy $cz]
}

proc ::PETK::gui::dotProduct {a b} {
    set ax [lindex $a 0]; set ay [lindex $a 1]; set az [lindex $a 2]
    set bx [lindex $b 0]; set by [lindex $b 1]; set bz [lindex $b 2]
    
    return [expr {$ax*$bx + $ay*$by + $az*$bz}]
}

proc ::PETK::gui::findPerpendicularVector {vec} {
    set x [lindex $vec 0]
    set y [lindex $vec 1]
    set z [lindex $vec 2]
    
    # Find the component with smallest absolute value
    if {abs($x) < abs($y) && abs($x) < abs($z)} {
        return [list 0 [expr {-$z}] $y]
    } elseif {abs($y) < abs($z)} {
        return [list [expr {-$z}] 0 $x]
    } else {
        return [list [expr {-$y}] $x 0]
    }
}

proc ::PETK::gui::rodriguesRotation {axis angle} {
    set ux [lindex $axis 0]; set uy [lindex $axis 1]; set uz [lindex $axis 2]
    set cos_a [expr {cos($angle)}]
    set sin_a [expr {sin($angle)}]
    set one_minus_cos [expr {1.0 - $cos_a}]
    
    # Rodrigues' rotation matrix
    set r00 [expr {$cos_a + $ux*$ux*$one_minus_cos}]
    set r01 [expr {$ux*$uy*$one_minus_cos - $uz*$sin_a}]
    set r02 [expr {$ux*$uz*$one_minus_cos + $uy*$sin_a}]
    
    set r10 [expr {$uy*$ux*$one_minus_cos + $uz*$sin_a}]
    set r11 [expr {$cos_a + $uy*$uy*$one_minus_cos}]
    set r12 [expr {$uy*$uz*$one_minus_cos - $ux*$sin_a}]
    
    set r20 [expr {$uz*$ux*$one_minus_cos - $uy*$sin_a}]
    set r21 [expr {$uz*$uy*$one_minus_cos + $ux*$sin_a}]
    set r22 [expr {$cos_a + $uz*$uz*$one_minus_cos}]
    
    return [list [list $r00 $r01 $r02] [list $r10 $r11 $r12] [list $r20 $r21 $r22]]
}

proc ::PETK::gui::rotationMatrix180 {axis} {
    set ux [lindex $axis 0]; set uy [lindex $axis 1]; set uz [lindex $axis 2]
    
    # 180-degree rotation matrix
    set r00 [expr {2*$ux*$ux - 1}]
    set r01 [expr {2*$ux*$uy}]
    set r02 [expr {2*$ux*$uz}]
    
    set r10 [expr {2*$uy*$ux}]
    set r11 [expr {2*$uy*$uy - 1}]
    set r12 [expr {2*$uy*$uz}]
    
    set r20 [expr {2*$uz*$ux}]
    set r21 [expr {2*$uz*$uy}]
    set r22 [expr {2*$uz*$uz - 1}]
    
    return [list [list $r00 $r01 $r02] [list $r10 $r11 $r12] [list $r20 $r21 $r22]]
}

# Helper functions for matrix operations
proc ::PETK::gui::matrix_transpose_3x3 {matrix} {
    set row0 [lindex $matrix 0]
    set row1 [lindex $matrix 1] 
    set row2 [lindex $matrix 2]
    
    return [list \
        [list [lindex $row0 0] [lindex $row1 0] [lindex $row2 0]] \
        [list [lindex $row0 1] [lindex $row1 1] [lindex $row2 1]] \
        [list [lindex $row0 2] [lindex $row1 2] [lindex $row2 2]]]
}

proc ::PETK::gui::matrix_multiply_3x3 {A B} {
    set result {}
    for {set i 0} {$i < 3} {incr i} {
        set row {}
        for {set j 0} {$j < 3} {incr j} {
            set sum 0.0
            for {set k 0} {$k < 3} {incr k} {
                set sum [expr {$sum + [lindex [lindex $A $i] $k] * [lindex [lindex $B $k] $j]}]
            }
            lappend row $sum
        }
        lappend result $row
    }
    return $result
}

proc ::PETK::gui::matrix_vector_multiply_3x3 {matrix vector} {
    set result {}
    for {set i 0} {$i < 3} {incr i} {
        set sum 0.0
        for {set j 0} {$j < 3} {incr j} {
            set sum [expr {$sum + [lindex [lindex $matrix $i] $j] * [lindex $vector $j]}]
        }
        lappend result $sum
    }
    return $result
}


####################################################
# Tab 3 SEM Setup
####################################################
proc ::PETK::gui::buildTab3 {tab3} {
    # Initialize any required variables
    ::PETK::gui::initializeSEMVariables

    # === CRITICAL: Configure tab to expand ===
    grid columnconfigure $tab3 0 -weight 1
    grid rowconfigure $tab3 0 -weight 1

    # Create scrollable canvas container
    canvas $tab3.canvas -highlightthickness 0
    ttk::scrollbar $tab3.vscroll -orient vertical -command [list $tab3.canvas yview]
    ttk::scrollbar $tab3.hscroll -orient horizontal -command [list $tab3.canvas xview]
    
    # Configure canvas scrolling
    $tab3.canvas configure -yscrollcommand [list $tab3.vscroll set]
    $tab3.canvas configure -xscrollcommand [list $tab3.hscroll set]
    
    # Create the actual content frame inside the canvas
    ttk::frame $tab3.canvas.content
    set canvas_window [$tab3.canvas create window 0 0 -anchor nw -window $tab3.canvas.content]
    
    # Grid the canvas and scrollbars with proper expansion
    grid $tab3.canvas -row 0 -column 0 -sticky nsew
    grid $tab3.vscroll -row 0 -column 1 -sticky ns
    grid $tab3.hscroll -row 1 -column 0 -sticky ew
    
    # Configure grid weights - CRITICAL for expansion
    grid rowconfigure $tab3 0 -weight 1
    grid columnconfigure $tab3 0 -weight 1

    # Now use the content frame as your container
    set container $tab3.canvas.content
    grid columnconfigure $container 0 -weight 1
    grid rowconfigure $container {4 5} -weight 1  ; # Make summary and parameter sections expandable
    
    set row 0

    # === SIMULATION PARAMETERS SECTION ===
    ttk::labelframe $container.simulation -text "Simulation Parameters" -padding 10
    grid $container.simulation -row $row -column 0 -sticky ew -padx 10 -pady "10 5"
    grid columnconfigure $container.simulation {1 3 5 7} -weight 1
    incr row

    # Applied voltage
    ttk::label $container.simulation.voltagelbl -text "Applied voltage (mV):" -width 18
    ttk::entry $container.simulation.voltage -textvariable ::PETK::gui::appliedVoltage -width 12 -justify center

    grid $container.simulation.voltagelbl $container.simulation.voltage - - - - -sticky ew -pady 3 -padx 5

    # Bulk and membrane conductivity
    ttk::label $container.simulation.bulklbl -text "Bulk conductivity (S/m):" -width 20
    ttk::entry $container.simulation.bulk -textvariable ::PETK::gui::bulkConductivity -width 12 -justify center
    ttk::label $container.simulation.memblbl -text "Membrane conductivity (S/m):" -width 25
    ttk::entry $container.simulation.memb -textvariable ::PETK::gui::membraneConductivity -width 12 -justify center

    grid $container.simulation.bulklbl $container.simulation.bulk $container.simulation.memblbl $container.simulation.memb -sticky ew -pady 3 -padx 5

    # === GRID SETTINGS SECTION ===
    ttk::labelframe $container.grid -text "Grid Settings" -padding 10
    grid $container.grid -row $row -column 0 -sticky ew -padx 10 -pady 5
    grid columnconfigure $container.grid {1 3 5} -weight 1
    incr row

    # Grid resolution
    ttk::label $container.grid.reslbl -text "Grid resolution (Å):" -width 18
    ttk::entry $container.grid.res -textvariable ::PETK::gui::gridResolution -width 12 -justify center
    ttk::label $container.grid.pointslbl -text "Estimated points:" -width 15
    ttk::label $container.grid.pointsval -textvariable ::PETK::gui::estimatedGridPoints -width 15 -anchor w -relief sunken -background white
    ttk::button $container.grid.update -text "Update Estimate" -command {::PETK::gui::estimateGridPoints}

    grid $container.grid.reslbl $container.grid.res $container.grid.pointslbl $container.grid.pointsval $container.grid.update - -sticky ew -pady 3

    # VdW radii settings
    ttk::label $container.grid.vdwlbl -text "Use VdW radii:" -width 15
    ttk::checkbutton $container.grid.vdw -variable ::PETK::gui::semUseVdWRadii
    ttk::label $container.grid.defaultlbl -text "Default radius (Å):" -width 18
    ttk::entry $container.grid.default -textvariable ::PETK::gui::semDefaultRadius -width 12 -justify center

    grid $container.grid.vdwlbl $container.grid.vdw $container.grid.defaultlbl $container.grid.default -sticky ew -pady 3

    # === MOVEMENT PARAMETERS SECTION ===
    ttk::labelframe $container.movement -text "Analyte Movement Parameters" -padding 10
    grid $container.movement -row $row -column 0 -sticky ew -padx 10 -pady "5 5"
    grid columnconfigure $container.movement {1 3 5} -weight 1
    incr row

    # Movement range
    ttk::label $container.movement.startlbl -text "Z start position (Å):" -width 18
    ttk::entry $container.movement.start -textvariable ::PETK::gui::zStartRange -width 12 -justify center
    ttk::label $container.movement.endlbl -text "Z end position (Å):" -width 18
    ttk::entry $container.movement.end -textvariable ::PETK::gui::zEndRange -width 12 -justify center
    ttk::label $container.movement.steplbl -text "Z step size (Å):" -width 15
    ttk::entry $container.movement.step -textvariable ::PETK::gui::zStep -width 12 -justify center

    grid $container.movement.startlbl $container.movement.start $container.movement.endlbl $container.movement.end $container.movement.steplbl $container.movement.step -sticky ew -pady 3

    # === PYTHON ENVIRONMENT SECTION ===
    ttk::labelframe $container.python -text "Python Environment" -padding 10
    grid $container.python -row $row -column 0 -sticky ew -padx 10 -pady "10 5"
    grid columnconfigure $container.python {1 2} -weight 1
    incr row

    # Conda environment selection
    ttk::label $container.python.envlbl -text "Conda Environment:" -width 20
    ttk::combobox $container.python.envcombo -textvariable ::PETK::gui::condaEnvironment -width 25
    ttk::button $container.python.refreshenv -text "Refresh" -command {::PETK::gui::refreshCondaEnvironments}
    ttk::button $container.python.createenv -text "Create SEM Env" -command {::PETK::gui::createSEMEnvironment}

    grid $container.python.envlbl $container.python.envcombo $container.python.refreshenv $container.python.createenv -sticky ew -pady 3

    # STORE WIDGET REFERENCE FOR LATER USE
    set ::PETK::gui::condaEnvComboWidget $container.python.envcombo

    # Environment status display - IMPROVED VERSION
    ttk::label $container.python.statuslbl -text "Status:" -width 20

    # Create a frame to hold the text widget and scrollbar
    ttk::frame $container.python.statusframe
    grid $container.python.statuslbl $container.python.statusframe - - -sticky ew -pady 3

    # Create text widget with scrollbar
    text $container.python.statusframe.text -height 4 -width 50 \
        -wrap word -font {TkDefaultFont 9} \
        -state disabled -cursor arrow \
        -relief sunken -bd 1

    ttk::scrollbar $container.python.statusframe.scroll -orient vertical \
        -command [list $container.python.statusframe.text yview]

    $container.python.statusframe.text configure -yscrollcommand [list $container.python.statusframe.scroll set]

    # Pack the text and scrollbar
    pack $container.python.statusframe.text -side left -fill both -expand 1
    pack $container.python.statusframe.scroll -side right -fill y

    # Configure the status frame to expand
    grid columnconfigure $container.python.statusframe 0 -weight 1

    # Store widget reference for status updates
    set ::PETK::gui::statusTextWidget $container.python.statusframe.text

    # Python executable path (auto-detected)
    ttk::label $container.python.pathlbl -text "Python Executable:" -width 20
    ttk::entry $container.python.path -textvariable ::PETK::gui::pythonExecutable -width 40 -state readonly
    ttk::button $container.python.testpython -text "Test Python" -command {::PETK::gui::testPythonEnvironment}

    grid $container.python.pathlbl $container.python.path $container.python.testpython - -sticky ew -pady 3
    
    # === PARAMETER SUMMARY SECTION ===
    ttk::labelframe $container.summary -text "Parameter Summary" -padding 10
    grid $container.summary -row $row -column 0 -sticky nsew -padx 10 -pady "5 5"
    grid columnconfigure $container.summary 0 -weight 1
    grid rowconfigure $container.summary 1 -weight 1
    incr row

    # Summary text widget
    ttk::frame $container.summary.content
    grid $container.summary.content -sticky nsew
    grid columnconfigure $container.summary.content 0 -weight 1
    grid rowconfigure $container.summary.content 0 -weight 1

    text $container.summary.content.text -height 8 -width 80 -wrap word -state disabled \
        -yscrollcommand [list $container.summary.content.scroll set] -font {TkFixedFont 9}
    ttk::scrollbar $container.summary.content.scroll -orient vertical -command [list $container.summary.content.text yview]

    grid $container.summary.content.text $container.summary.content.scroll -sticky nsew
    grid columnconfigure $container.summary.content 0 -weight 1
    grid rowconfigure $container.summary.content 0 -weight 1

    # Store reference to summary text widget
    set ::PETK::gui::semSummaryTextWidget $container.summary.content.text

    # Update button
    ttk::button $container.summary.update -text "Update Summary" -command {::PETK::gui::updateSEMParameterSummary}
    grid $container.summary.update -sticky ew -pady "5 0"

    # === ACTION BUTTONS SECTION ===
    ttk::labelframe $container.actions -text "Actions" -padding 10
    grid $container.actions -row $row -column 0 -sticky ew -padx 10 -pady "5 10"
    grid columnconfigure $container.actions {0 1 2} -weight 1
    incr row

    # Action buttons
    ttk::button $container.actions.validate -text "Validate Parameters" -command {::PETK::gui::validateSEMSetup}
    ttk::button $container.actions.generate -text "Preview Simulation" -command {::PETK::gui::runPreviewFromGUI}
    ttk::button $container.actions.visualize -text "Run Simulation" -command {::PETK::gui::runCalculationFromGUI}

    grid $container.actions.validate $container.actions.generate $container.actions.visualize -sticky ew -padx 3

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab3.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab3.canvas $canvas_window]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab3.canvas $canvas_window]
    
    # Store references for later use (with unique names for tab3)
    set ::PETK::gui::tab3MainCanvas $tab3.canvas
    set ::PETK::gui::tab3CanvasWindow $canvas_window
    set ::PETK::gui::tab3ContentContainer $container
    
    # Initialize SEM variables
    ::PETK::gui::initializeSEMVariables
    
    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab3.canvas $canvas_window]
}

####################################################
# Tab 3 Function
####################################################

####################################################
## SEM VARIABLE INITIALIZATION
####################################################

proc ::PETK::gui::initializeSEMVariables {} {
    
    # Simulation parameters
    if {![info exists ::PETK::gui::appliedVoltage]} {
        set ::PETK::gui::appliedVoltage "100"
    }
    if {![info exists ::PETK::gui::bulkConductivity]} {
        set ::PETK::gui::bulkConductivity "1.0"
    }
    if {![info exists ::PETK::gui::membraneConductivity]} {
        set ::PETK::gui::membraneConductivity "0.001"
    }
    
    # Grid settings
    if {![info exists ::PETK::gui::gridResolution]} {
        set ::PETK::gui::gridResolution "2.0"
    }
    if {![info exists ::PETK::gui::estimatedGridPoints]} {
        set ::PETK::gui::estimatedGridPoints "Not calculated"
    }
    if {![info exists ::PETK::gui::semUseVdWRadii]} {
        set ::PETK::gui::semUseVdWRadii 1
    }
    if {![info exists ::PETK::gui::semDefaultRadius]} {
        set ::PETK::gui::semDefaultRadius "1.5"
    }
    
    # Movement parameters
    if {![info exists ::PETK::gui::zStartRange]} {
        set ::PETK::gui::zStartRange "-50.0"
    }
    if {![info exists ::PETK::gui::zEndRange]} {
        set ::PETK::gui::zEndRange "50.0"
    }
    if {![info exists ::PETK::gui::zStep]} {
        set ::PETK::gui::zStep "2.0"
    }
    
    # Python environment
    if {![info exists ::PETK::gui::condaEnvironment]} {
        set ::PETK::gui::condaEnvironment "base"
    }
    if {![info exists ::PETK::gui::pythonExecutable]} {
        set ::PETK::gui::pythonExecutable ""
    }
    
}

proc ::PETK::gui::estimateGridPoints {} {
    set resolution 1.0
    if {[string is double $::PETK::gui::gridResolution]} {
        set resolution $::PETK::gui::gridResolution
    }
    
    if {$::PETK::gui::autoCalculateBoxDimensions} {
        # Use auto-calculated dimensions
        if {[info exists ::PETK::gui::calculatedBoxSizeX]} {
            set nx [expr {int($::PETK::gui::calculatedBoxSizeX / $resolution)}]
            set ny [expr {int($::PETK::gui::calculatedBoxSizeY / $resolution)}]
            set nz [expr {int($::PETK::gui::calculatedBoxSizeZ / $resolution)}]
        } else {
            set ::PETK::gui::estimatedGridPoints "Calculate dimensions first"
            return
        }
    } else {
        # Use manual dimensions
        if {[string is double $::PETK::gui::boxSizeX] && [string is double $::PETK::gui::boxSizeY] && [string is double $::PETK::gui::boxSizeZ]} {
            set nx [expr {int($::PETK::gui::boxSizeX / $resolution)}]
            set ny [expr {int($::PETK::gui::boxSizeY / $resolution)}]
            set nz [expr {int($::PETK::gui::boxSizeZ / $resolution)}]
        } else {
            set ::PETK::gui::estimatedGridPoints "Invalid dimensions"
            return
        }
    }
    
    set total_points [expr {$nx * $ny * $nz}]
    
    if {$total_points > 1000000000} {
        set ::PETK::gui::estimatedGridPoints [format "%.1fB points (!!)" [expr {$total_points / 1000000000.0}]]
    } elseif {$total_points > 1000000} {
        set ::PETK::gui::estimatedGridPoints [format "%.1fM points" [expr {$total_points / 1000000.0}]]
    } elseif {$total_points > 1000} {
        set ::PETK::gui::estimatedGridPoints [format "%.1fK points" [expr {$total_points / 1000.0}]]
    } else {
        set ::PETK::gui::estimatedGridPoints [format "%d points" $total_points]
    }
}


####################################################
# Parameter Summary Functions
####################################################

proc ::PETK::gui::updateSEMParameterSummary {} {
    if {![info exists ::PETK::gui::semSummaryTextWidget] || ![winfo exists $::PETK::gui::semSummaryTextWidget]} {
        return
    }
    
    set widget $::PETK::gui::semSummaryTextWidget
    $widget configure -state normal
    $widget delete 1.0 end
    
    # Add parameter summary
    $widget insert end "SEM CALCULATION PARAMETERS\n" {title}
    $widget insert end [string repeat "=" 50]
    $widget insert end "\n\n"
    
    # Input Configuration
    $widget insert end "INPUT CONFIGURATION:\n" {header}
    $widget insert end "  Analyte PDB: $::PETK::gui::semAnalytePDB\n"
    $widget insert end "  Python Environment: $::PETK::gui::condaEnvironment\n"
    $widget insert end "\n"
    
    # Pore Geometry
    $widget insert end "PORE GEOMETRY:\n" {header}
    $widget insert end "  Pore Type: $::PETK::gui::membraneType\n"
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        $widget insert end "  Pore Diameter: $::PETK::gui::cylindricalDiameter Å\n"
        $widget insert end "  Corner Radius: $::PETK::gui::cornerRadius Å\n"
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        $widget insert end "  Inner Diameter: $::PETK::gui::innerDiameter Å\n"
        $widget insert end "  Outer Diameter: $::PETK::gui::outerDiameter Å\n"
    }
    $widget insert end "  Membrane Thickness: $::PETK::gui::nanoporeThickness Å\n"
    $widget insert end "\n"
    
    # Box Dimensions
    $widget insert end "BOX DIMENSIONS:\n" {header}
    $widget insert end "  Distance Cutoff: $::PETK::gui::sysCutoff Å\n"

    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        $widget insert end "  Mode: Auto-calculated\n"
        
        # Display auto-calculated dimensions
        if {[info exists ::PETK::gui::autoBoxX] && $::PETK::gui::autoBoxX ne ""} {
            $widget insert end "  X dimension: $::PETK::gui::autoBoxX\n"
        } else {
            $widget insert end "  X dimension: Not calculated\n" {warning}
        }
        
        if {[info exists ::PETK::gui::autoBoxY] && $::PETK::gui::autoBoxY ne ""} {
            $widget insert end "  Y dimension: $::PETK::gui::autoBoxY\n"
        } else {
            $widget insert end "  Y dimension: Not calculated\n" {warning}
        }
        
        if {[info exists ::PETK::gui::autoBoxZ] && $::PETK::gui::autoBoxZ ne ""} {
            $widget insert end "  Z dimension: $::PETK::gui::autoBoxZ\n"
        } else {
            $widget insert end "  Z dimension: Not calculated\n" {warning}
        }
        
        # Extract and display total box volume if possible
        if {[info exists ::PETK::gui::autoBoxX] && [info exists ::PETK::gui::autoBoxY] && [info exists ::PETK::gui::autoBoxZ]} {
            set volume_calculated 0
            # Try to extract sizes from auto-calculated strings (format: "min to max (size)")
            if {[regexp {\(([0-9.]+)\)} $::PETK::gui::autoBoxX match x_size] && 
                [regexp {\(([0-9.]+)\)} $::PETK::gui::autoBoxY match y_size] && 
                [regexp {\(([0-9.]+)\)} $::PETK::gui::autoBoxZ match z_size]} {
                set volume [expr {$x_size * $y_size * $z_size}]
                $widget insert end "  Total Volume: [format "%.2e" $volume] Å³\n"
                set volume_calculated 1
            }
            if {!$volume_calculated} {
                $widget insert end "  Total Volume: Cannot calculate (dimensions not parsed)\n" {warning}
            }
        }
        
    } else {
        $widget insert end "  Mode: Manual\n"
        
        # Display manual dimensions
        if {[info exists ::PETK::gui::boxSizeX] && $::PETK::gui::boxSizeX ne ""} {
            $widget insert end "  X size: $::PETK::gui::boxSizeX Å\n"
        } else {
            $widget insert end "  X size: Not specified\n" {warning}
        }
        
        if {[info exists ::PETK::gui::boxSizeY] && $::PETK::gui::boxSizeY ne ""} {
            $widget insert end "  Y size: $::PETK::gui::boxSizeY Å\n"
        } else {
            $widget insert end "  Y size: Not specified\n" {warning}
        }
        
        if {[info exists ::PETK::gui::boxSizeZ] && $::PETK::gui::boxSizeZ ne ""} {
            $widget insert end "  Z size: $::PETK::gui::boxSizeZ Å\n"
        } else {
            $widget insert end "  Z size: Not specified\n" {warning}
        }
        
        # Calculate and display total volume and aspect ratio for manual dimensions
        if {[info exists ::PETK::gui::boxSizeX] && [info exists ::PETK::gui::boxSizeY] && [info exists ::PETK::gui::boxSizeZ]} {
            if {[string is double $::PETK::gui::boxSizeX] && [string is double $::PETK::gui::boxSizeY] && [string is double $::PETK::gui::boxSizeZ]} {
                set volume [expr {$::PETK::gui::boxSizeX * $::PETK::gui::boxSizeY * $::PETK::gui::boxSizeZ}]
                $widget insert end "  Total Volume: [format "%.2e" $volume] Å³\n"
                
                # Calculate aspect ratio
                set max_dim [expr {max($::PETK::gui::boxSizeX, $::PETK::gui::boxSizeY, $::PETK::gui::boxSizeZ)}]
                set min_dim [expr {min($::PETK::gui::boxSizeX, $::PETK::gui::boxSizeY, $::PETK::gui::boxSizeZ)}]
                set aspect_ratio [expr {$max_dim / $min_dim}]
                $widget insert end "  Aspect Ratio: [format "%.2f" $aspect_ratio]:1"
                
                if {$aspect_ratio > 10.0} {
                    $widget insert end " (WARNING: Extreme aspect ratio)" {warning}
                }
                $widget insert end "\n"
            }
        }
    }
    $widget insert end "\n"

    # Simulation Parameters
    $widget insert end "SIMULATION PARAMETERS:\n" {header}
    $widget insert end "  Applied Voltage: $::PETK::gui::appliedVoltage mV\n"
    $widget insert end "  Bulk Conductivity: $::PETK::gui::bulkConductivity S/m\n"
    $widget insert end "  Grid Resolution: $::PETK::gui::gridResolution Å\n"
    $widget insert end "  Use VdW Radii: [expr {$::PETK::gui::semUseVdWRadii ? "Yes" : "No"}]\n"
    $widget insert end "  Default Radius: $::PETK::gui::semDefaultRadius Å\n"
    $widget insert end "\n"
    
    # Movement Parameters
    $widget insert end "MOVEMENT PARAMETERS:\n" {header}
    $widget insert end "  Z Start: $::PETK::gui::zStartRange Å\n"
    $widget insert end "  Z End: $::PETK::gui::zEndRange Å\n"
    $widget insert end "  Z Step: $::PETK::gui::zStep Å\n"
    
    # Calculate number of steps
    if {[string is double $::PETK::gui::zStartRange] && [string is double $::PETK::gui::zEndRange] && [string is double $::PETK::gui::zStep]} {
        set range [expr {abs($::PETK::gui::zEndRange - $::PETK::gui::zStartRange)}]
        set num_steps [expr {int($range / $::PETK::gui::zStep) + 1}]
        $widget insert end "  Number of Steps: $num_steps\n"
        $widget insert end "  Total Range: $range Å\n"
    }
    $widget insert end "\n"
    
    # Output Settings
    $widget insert end "OUTPUT SETTINGS:\n" {header}
    $widget insert end "  Output Directory: $::PETK::gui::workdir\n"
    $widget insert end "  Output Prefix: $::PETK::gui::outputPrefix\n"
    $widget insert end "  Preview Frames: $::PETK::gui::semPreviewFrames\n"
    
    # Configure text tags for formatting
    $widget tag configure title -font {TkDefaultFont 12 bold}
    $widget tag configure header -font {TkDefaultFont 10 bold}
    
    $widget configure -state disabled
}

####################################################
# Conda Functions 
####################################################
####################################################
# Helper function to update status in text widget
####################################################

proc ::PETK::gui::updateStatusDisplay {message} {
    # Update both the variable (for compatibility) and the text widget
    set ::PETK::gui::pythonEnvStatus $message
    
    if {[info exists ::PETK::gui::statusTextWidget] && [winfo exists $::PETK::gui::statusTextWidget]} {
        $::PETK::gui::statusTextWidget configure -state normal
        $::PETK::gui::statusTextWidget delete 1.0 end
        $::PETK::gui::statusTextWidget insert 1.0 $message
        $::PETK::gui::statusTextWidget configure -state disabled
        
        # Auto-scroll to bottom if there's new content
        $::PETK::gui::statusTextWidget see end
    }
}

proc ::PETK::gui::appendStatusDisplay {message} {
    # Append message to existing status (useful for detailed logs)
    if {[info exists ::PETK::gui::statusTextWidget] && [winfo exists $::PETK::gui::statusTextWidget]} {
        $::PETK::gui::statusTextWidget configure -state normal
        $::PETK::gui::statusTextWidget insert end "\n$message"
        $::PETK::gui::statusTextWidget configure -state disabled
        $::PETK::gui::statusTextWidget see end
    }
}

####################################################
# Updated Conda Functions with improved status display
####################################################

proc ::PETK::gui::refreshCondaEnvironments {} {
    ::PETK::gui::updateStatusDisplay "Checking conda environments..."
    update
    
    # Check if conda is available
    if {[catch {exec conda env list} conda_output]} {
        ::PETK::gui::updateStatusDisplay "❌ Conda not found in PATH"
        return
    }
    
    # Parse conda environments
    set env_list {}
    set lines [split $conda_output "\n"]
    foreach line $lines {
        # Skip comment lines and empty lines
        if {[string match "#*" $line] || [string trim $line] eq ""} {
            continue
        }
        
        # Extract environment name (first word)
        set env_name [lindex [split $line] 0]
        if {$env_name ne ""} {
            lappend env_list $env_name
        }
    }
    
    # Update combobox using stored widget reference
    if {[info exists ::PETK::gui::condaEnvComboWidget] && [winfo exists $::PETK::gui::condaEnvComboWidget]} {
        $::PETK::gui::condaEnvComboWidget configure -values $env_list
        
        # Set default to first available SEM-related environment or base
        set default_env "base"
        foreach env $env_list {
            if {[string match "*sem*" [string tolower $env]] || 
                [string match "*fenics*" [string tolower $env]] ||
                [string match "*mdanalysis*" [string tolower $env]]} {
                set default_env $env
                break
            }
        }
        
        if {$::PETK::gui::condaEnvironment eq ""} {
            set ::PETK::gui::condaEnvironment $default_env
        }
        
        ::PETK::gui::updateStatusDisplay "✓ Found [llength $env_list] conda environments\nEnvironments: [join $env_list {, }]"
        
        # Debug output
        puts "Conda environments found: $env_list"
        puts "Selected environment: $::PETK::gui::condaEnvironment"
        
    } else {
        ::PETK::gui::updateStatusDisplay "❌ Could not find environment combobox widget"
        puts "ERROR: Widget reference not found"
    }
    
    # Update Python executable path
    ::PETK::gui::updatePythonExecutable
}

proc ::PETK::gui::updatePythonExecutable {} {
    if {$::PETK::gui::condaEnvironment eq ""} {
        set ::PETK::gui::pythonExecutable ""
        return
    }
    
    # Get conda info to find the path
    if {[catch {exec conda info --envs} conda_info]} {
        set ::PETK::gui::pythonExecutable ""
        return
    }
    
    # Find the path for the selected environment
    set lines [split $conda_info "\n"]
    foreach line $lines {
        if {[string match "*$::PETK::gui::condaEnvironment *" $line]} {
            set parts [split $line]
            if {[llength $parts] >= 2} {
                set env_path [lindex $parts end]
                if {[string eq $::tcl_platform(platform) "windows"]} {
                    set ::PETK::gui::pythonExecutable [file join $env_path "python.exe"]
                } else {
                    set ::PETK::gui::pythonExecutable [file join $env_path "bin" "python"]
                }
                return
            }
        }
    }
    
    # Fallback: try to construct path
    if {[catch {exec conda info --base} conda_base]} {
        set ::PETK::gui::pythonExecutable ""
        return
    }
    
    set conda_base [string trim $conda_base]
    if {$::PETK::gui::condaEnvironment eq "base"} {
        set env_path $conda_base
    } else {
        set env_path [file join $conda_base "envs" $::PETK::gui::condaEnvironment]
    }
    
    if {[string eq $::tcl_platform(platform) "windows"]} {
        set ::PETK::gui::pythonExecutable [file join $env_path "python.exe"]
    } else {
        set ::PETK::gui::pythonExecutable [file join $env_path "bin" "python"]
    }
}

proc ::PETK::gui::testPythonEnvironment {} {
    if {$::PETK::gui::pythonExecutable eq ""} {
        ::PETK::gui::updateStatusDisplay "❌ No Python executable specified"
        return
    }
    
    ::PETK::gui::updateStatusDisplay "Testing Python environment..."
    update
    
    # Test basic Python
    if {[catch {exec $::PETK::gui::pythonExecutable --version} python_version]} {
        ::PETK::gui::updateStatusDisplay "❌ Python test failed:\n$python_version"
        return
    }
    
    ::PETK::gui::appendStatusDisplay "✓ Python version: $python_version"
    
    # Test required packages
    set test_script "
import sys
packages = \['numpy', 'matplotlib', 'MDAnalysis', 'scipy'\]
missing = \[\]
available = \[\]

for package in packages:
    try:
        __import__(package)
        available.append(package)
        print(f'{package} OK')
    except ImportError:
        missing.append(package)
        print(f'{package} MISSING')

# Test FEniCS
try:
    from fenics import *
    available.append('fenics')
    print('fenics OK')
except Exception as e:
    missing.append('fenics')
    print(f'fenics FAILED: {e}')

print(f'AVAILABLE: {available}')
if missing:
    print(f'MISSING: {missing}')
else:
    print('ALL_OK')
"
    
    set exit_code [catch {exec conda run -n $::PETK::gui::condaEnvironment python -c $test_script 2>@1} result]

    # Parse the results regardless of warnings
    if {[string match "*ALL_OK*" $result]} {
        ::PETK::gui::appendStatusDisplay "✅ All required packages available!"
        ::PETK::gui::appendStatusDisplay "Ready for SEM simulation."
    } elseif {[regexp {MISSING: \[(.*)\]} $result match missing_list]} {
        ::PETK::gui::appendStatusDisplay "⚠️ Missing packages: $missing_list"
    } else {
        ::PETK::gui::appendStatusDisplay "❌ Environment test failed:\n$result"
        return
    }
    
    # Show available packages
    if {[regexp {AVAILABLE: \[(.*)\]} $result match available_list]} {
        ::PETK::gui::appendStatusDisplay "Available packages: $available_list"
    }
}

proc ::PETK::gui::createSEMEnvironment {} {
    set result [tk_messageBox -icon question -type yesno -title "Create SEM Environment" \
        -message "This will create a new conda environment called 'petk-sem' with all required packages.\n\nThis may take several minutes. Continue?"]
    
    if {$result eq "no"} {
        return
    }
    
    ::PETK::gui::updateStatusDisplay "Creating SEM environment 'petk-sem'..."
    ::PETK::gui::appendStatusDisplay "This may take several minutes, please wait..."
    update
    
    # Create environment with complete FEniCS ecosystem
    set create_cmd [list conda create -n petk-sem -y -c conda-forge \
        python=3.9 \
        fenics=2019.1.0 \
        pkg-config=0.29.2 \
        mpich \
        petsc \
        slepc \
        petsc4py \
        slepc4py \
        mpi4py \
        numpy=1.26.4 \
        scipy \
        matplotlib-base \
        mdanalysis \
        h5py]
    
    if {[catch {exec {*}$create_cmd} create_result]} {
        ::PETK::gui::updateStatusDisplay "❌ Environment creation failed"
        ::PETK::gui::appendStatusDisplay "Error details:\n$create_result"
        tk_messageBox -icon error -title "Environment Creation Error" \
            -message "Failed to create SEM environment:\n$create_result"
        return
    }
    
    ::PETK::gui::appendStatusDisplay "✅ Environment packages installed successfully!"
    ::PETK::gui::appendStatusDisplay "Setting up FEniCS environment scripts..."
    update
    
    # Set up FEniCS environment activation/deactivation scripts
    if {[catch {
        # Get the conda environment path
        set env_info [exec conda info --envs]
        set conda_prefix ""
        
        foreach line [split $env_info "\n"] {
            if {[string match "*petk-sem*" $line]} {
                set parts [regexp -all -inline {\S+} $line]
                set conda_prefix [lindex $parts end]
                break
            }
        }
        
        if {$conda_prefix eq ""} {
            error "Could not find petk-sem environment path"
        }
        
        ::PETK::gui::appendStatusDisplay "Environment path: $conda_prefix"
        
        # Create activation/deactivation directories
        file mkdir [file join $conda_prefix etc conda activate.d]
        file mkdir [file join $conda_prefix etc conda deactivate.d]
        
        # Create activation script
        set activate_script_path [file join $conda_prefix etc conda activate.d fenics_env.sh]
        set activate_script [open $activate_script_path w]
        puts $activate_script "#!/bin/bash"
        puts $activate_script "export PKG_CONFIG_PATH=\$CONDA_PREFIX/lib/pkgconfig:\$CONDA_PREFIX/share/pkgconfig:\$PKG_CONFIG_PATH"
        close $activate_script
        
        # Create deactivation script
        set deactivate_script_path [file join $conda_prefix etc conda deactivate.d fenics_env.sh]
        set deactivate_script [open $deactivate_script_path w]
        puts $deactivate_script "#!/bin/bash"
        puts $deactivate_script "export PKG_CONFIG_PATH=\$(echo \$PKG_CONFIG_PATH | sed \"s|\$CONDA_PREFIX/lib/pkgconfig:||g\" | sed \"s|\$CONDA_PREFIX/share/pkgconfig:||g\")"
        close $deactivate_script
        
        # Make scripts executable (only on Unix-like systems)
        if {$::tcl_platform(platform) ne "windows"} {
            exec chmod +x $activate_script_path
            exec chmod +x $deactivate_script_path
        }
        
        ::PETK::gui::appendStatusDisplay "✅ FEniCS environment scripts created successfully!"
        
    } script_error]} {
        ::PETK::gui::appendStatusDisplay "⚠️ Warning: Failed to create FEniCS environment scripts"
        ::PETK::gui::appendStatusDisplay "Error: $script_error"
        ::PETK::gui::appendStatusDisplay "Environment may still work but PKG_CONFIG_PATH may need manual setup"
    }
    
    # Final setup
    ::PETK::gui::updateStatusDisplay "✅ SEM environment 'petk-sem' created successfully!"
    ::PETK::gui::appendStatusDisplay "Environment is ready for use with FEniCS configuration."
    ::PETK::gui::appendStatusDisplay "Please restart your terminal or reactivate the environment to use the new scripts."
    
    set ::PETK::gui::condaEnvironment "petk-sem"
    ::PETK::gui::refreshCondaEnvironments
    ::PETK::gui::testPythonEnvironment
    
    tk_messageBox -icon info -title "Environment Created" \
        -message "SEM environment 'petk-sem' created successfully!\n\nFEniCS activation scripts have been configured.\nRestart your terminal or reactivate the environment to use the new configuration."
}
####################################################
# Action button function
####################################################
proc ::PETK::gui::validateSEMSetup {} {
    set ::PETK::gui::semCurrentStatus "Validating setup..."
    update
    
    set errors {}
    set warnings {}
    
    # Check input file
    if {$::PETK::gui::semAnalytePDB eq ""} {
        lappend errors "No analyte PDB file specified"
    } elseif {![file exists $::PETK::gui::semAnalytePDB]} {
        lappend errors "Analyte PDB file does not exist: $::PETK::gui::semAnalytePDB"
    }
    
    # Check numeric parameters that must be positive
    set positive_params {
        {cylindricalDiameter "Pore diameter"}
        {nanoporeThickness "Membrane thickness"}
        {zStep "Z step size"}
        {appliedVoltage "Applied voltage"}
        {bulkConductivity "Bulk conductivity"}
        {sysCutoff "Distance cutoff"}
        {gridResolution "Grid resolution"}
    }
    
    foreach param_info $positive_params {
        set var_name [lindex $param_info 0]
        set display_name [lindex $param_info 1]
        set var_value [set ::PETK::gui::$var_name]
        
        if {![string is double $var_value]} {
            lappend errors "$display_name is not a valid number: $var_value"
        } elseif {$var_value < 0} {
            lappend errors "$display_name must be positive: $var_value"
        }
    }
    
    # Check Z coordinate parameters (can be negative, just need to be valid numbers)
    set z_coord_params {
        {zStartRange "Z start position"}
        {zEndRange "Z end position"}
    }
    
    foreach param_info $z_coord_params {
        set var_name [lindex $param_info 0]
        set display_name [lindex $param_info 1]
        set var_value [set ::PETK::gui::$var_name]
        
        if {![string is double $var_value]} {
            lappend errors "$display_name is not a valid number: $var_value"
        }
    }
    
    # Check that Z range is larger than step size
    if {[string is double $::PETK::gui::zStartRange] && [string is double $::PETK::gui::zEndRange] && [string is double $::PETK::gui::zStep]} {
        set z_range [expr {abs($::PETK::gui::zStartRange - $::PETK::gui::zEndRange)}]
        if {$z_range <= $::PETK::gui::zStep} {
            lappend errors "Z range (|Zend - Zstart| = $z_range) must be greater than Z step size ($::PETK::gui::zStep)"
        }
    }

    # ========== BOX DIMENSION VALIDATION ==========
    # Check if auto-calculate is enabled or manual mode
    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        # Auto-calculate mode - validate auto-calculated dimensions exist and are reasonable
        if {[info exists ::PETK::gui::autoBoxX] && $::PETK::gui::autoBoxX ne ""} {
            # Parse auto-calculated values (assuming format like "min to max (size)")
            if {[regexp {([0-9.-]+)\s+to\s+([0-9.-]+)\s+\(([0-9.]+)\)} $::PETK::gui::autoBoxX match min max size]} {
                if {$size < 10.0} {
                    lappend warnings "Auto-calculated X dimension is very small ($size Å)"
                } elseif {$size > 1000.0} {
                    lappend warnings "Auto-calculated X dimension is very large ($size Å)"
                }
            } else {
                lappend warnings "Auto-calculated X dimension format not recognized: $::PETK::gui::autoBoxX"
            }
        } else {
            lappend errors "Auto-calculate mode enabled but X dimension not calculated"
        }
        
        if {[info exists ::PETK::gui::autoBoxY] && $::PETK::gui::autoBoxY ne ""} {
            if {[regexp {([0-9.-]+)\s+to\s+([0-9.-]+)\s+\(([0-9.]+)\)} $::PETK::gui::autoBoxY match min max size]} {
                if {$size < 10.0} {
                    lappend warnings "Auto-calculated Y dimension is very small ($size Å)"
                } elseif {$size > 1000.0} {
                    lappend warnings "Auto-calculated Y dimension is very large ($size Å)"
                }
            } else {
                lappend warnings "Auto-calculated Y dimension format not recognized: $::PETK::gui::autoBoxY"
            }
        } else {
            lappend errors "Auto-calculate mode enabled but Y dimension not calculated"
        }
        
        if {[info exists ::PETK::gui::autoBoxZ] && $::PETK::gui::autoBoxZ ne ""} {
            if {[regexp {([0-9.-]+)\s+to\s+([0-9.-]+)\s+\(([0-9.]+)\)} $::PETK::gui::autoBoxZ match min max size]} {
                if {$size < 10.0} {
                    lappend warnings "Auto-calculated Z dimension is very small ($size Å)"
                } elseif {$size > 1000.0} {
                    lappend warnings "Auto-calculated Z dimension is very large ($size Å)"
                }
            } else {
                lappend warnings "Auto-calculated Z dimension format not recognized: $::PETK::gui::autoBoxZ"
            }
        } else {
            lappend errors "Auto-calculate mode enabled but Z dimension not calculated"
        }
        
    } else {
        # Manual mode - validate manually entered box dimensions
        set box_params {
            {boxSizeX "Box X dimension"}
            {boxSizeY "Box Y dimension"}
            {boxSizeZ "Box Z dimension"}
        }
        
        foreach param_info $box_params {
            set var_name [lindex $param_info 0]
            set display_name [lindex $param_info 1]
            
            if {[info exists ::PETK::gui::$var_name]} {
                set var_value [set ::PETK::gui::$var_name]
                
                if {$var_value eq ""} {
                    lappend errors "$display_name is empty"
                } elseif {![string is double $var_value]} {
                    lappend errors "$display_name is not a valid number: $var_value"
                } elseif {$var_value <= 0} {
                    lappend errors "$display_name must be positive: $var_value"
                } elseif {$var_value < 10.0} {
                    lappend warnings "$display_name is very small ($var_value Å)"
                } elseif {$var_value > 1000.0} {
                    lappend warnings "$display_name is very large ($var_value Å)"
                }
            } else {
                lappend errors "$display_name variable not found"
            }
        }
        
        # Check box dimension ratios for reasonableness
        if {[info exists ::PETK::gui::boxSizeX] && [info exists ::PETK::gui::boxSizeY] && [info exists ::PETK::gui::boxSizeZ]} {
            if {[string is double $::PETK::gui::boxSizeX] && [string is double $::PETK::gui::boxSizeY] && [string is double $::PETK::gui::boxSizeZ]} {
                set max_dim [expr {max($::PETK::gui::boxSizeX, $::PETK::gui::boxSizeY, $::PETK::gui::boxSizeZ)}]
                set min_dim [expr {min($::PETK::gui::boxSizeX, $::PETK::gui::boxSizeY, $::PETK::gui::boxSizeZ)}]
                
                if {$max_dim / $min_dim > 10.0} {
                    lappend warnings "Box dimensions have extreme aspect ratio (max/min = [format "%.1f" [expr {$max_dim / $min_dim}]])"
                }
            }
        }
    }
    
    # Check output directory
    if {![file exists $::PETK::gui::workdir]} {
        lappend warnings "Output directory does not exist: $::PETK::gui::workdir"
    } elseif {![file writable $::PETK::gui::workdir]} {
        lappend errors "Output directory is not writable: $::PETK::gui::workdir"
    }
    
    # Check for Python availability
    if {[catch {exec python3 --version} python_version]} {
        if {[catch {exec python --version} python_version]} {
            lappend errors "Python is not available in the system PATH"
        } else {
            lappend warnings "Using 'python' instead of 'python3'"
        }
    }
    
    # Calculate estimated calculation time
    if {[string is double $::PETK::gui::zStartRange] && [string is double $::PETK::gui::zEndRange] && [string is double $::PETK::gui::zStep]} {
        set z_range [expr {abs($::PETK::gui::zStartRange - $::PETK::gui::zEndRange)}]
        set num_steps [expr {int($z_range / $::PETK::gui::zStep) + 1}]
        set estimated_minutes [expr {$num_steps * 0.5}]  ; # Rough estimate: 30 seconds per step
        set ::PETK::gui::semEstimatedTime [format "~%.1f minutes (%d steps)" $estimated_minutes $num_steps]
    }
    
    # Display results
    set message "Validation Results:\n\n"
    
    if {[llength $errors] > 0} {
        append message "ERRORS (must be fixed):\n"
        foreach error $errors {
            append message "• $error\n"
        }
        append message "\n"
    }
    
    if {[llength $warnings] > 0} {
        append message "WARNINGS:\n"
        foreach warning $warnings {
            append message "• $warning\n"
        }
        append message "\n"
    }
    
    if {[llength $errors] == 0} {
        append message "✓ Setup validation PASSED\n"
        append message "Ready to run SEM calculations\n"
        set ::PETK::gui::semValidationPassed 1
        set ::PETK::gui::semCurrentStatus "Validation passed - ready to run"
        set icon "info"
        ::PETK::gui::copyModuleToOutputDir
        ::PETK::gui::outputParametersToConfig
    } else {
        append message "✗ Setup validation FAILED\n"
        append message "Please fix errors before proceeding\n"
        set ::PETK::gui::semCurrentStatus "Validation failed - fix errors"
        set icon "error"
    }
    
    tk_messageBox -icon $icon -title "SEM Setup Validation" -message $message
    
    return [expr {[llength $errors] == 0}]
}

proc ::PETK::gui::copyModuleToOutputDir {} {
    set tcl_script_dir [file dirname [info script]]
    set module_file [file join $tcl_script_dir "vertical_movement_sem.py"]
    set output_module [file join $::PETK::gui::workdir "vertical_movement_sem.py"]
    
    if {[file exists $module_file]} {
        if {[catch {file copy -force $module_file $output_module} error]} {
            tk_messageBox -icon error -title "Module Copy Error" \
                -message "Failed to copy module file:\n\n$error"
            return 0
        }
        puts "Copied module to output directory: $output_module"
        
        # Optional: Show success dialog (uncomment if desired)
        # tk_messageBox -icon info -title "Module Copied" \
        #     -message "Module copied successfully to:\n$output_module"
        
        return 1
    } else {
        tk_messageBox -icon error -title "Module Not Found" \
            -message "Module file not found:\n$module_file\n\nPlease ensure vertical_movement_sem.py is in the same directory as the Tcl script."
        return 0
    }
}

proc ::PETK::gui::outputParametersToConfig {{output_file ""}} {
    
    # Determine output file location
    if {$output_file eq ""} {
        # Check if workdir variable exists and has a value, similar to copyModuleToOutputDir
        if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
            set output_file [file join $::PETK::gui::workdir "config.json"]
        } else {
            set output_file [file join [pwd] "config.json"]
        }
    }
    
    # Convert boolean values to JSON format
    set use_vdw_json [expr {$::PETK::gui::semUseVdWRadii ? "true" : "false"}]
    # set preview_only_json [expr {$::PETK::gui::semPreviewOnly ? "true" : "false"}]
    
    # Get relative path for PDB file (relative to output directory)
    set output_dir [file dirname $output_file]
    if {$::PETK::gui::semAnalytePDB ne ""} {
        set pdb_relative [file tail $::PETK::gui::semAnalytePDB]
        set pdb_path "../$pdb_relative"
    }
    
    # Helper function to parse box dimension ranges
    proc parseBoxRange {range_string} {
        # Handle formats like "(-150.0, 150.0)" or "-150.0 to 150.0"
        set range_string [string trim $range_string]
        
        # Try to extract two numbers from the string
        if {[regexp {\(?\s*(-?[\d\.]+)\s*,?\s*to?\s*(-?[\d\.]+)\s*\)?} $range_string -> min_val max_val]} {
            return "\[$min_val, $max_val\]"
        }
        
        # If parsing fails, return a default range
        return "\[-150.0, 150.0\]"
    }
    
    # Build JSON content
    set json_content "{\n"
    
    # Metadata section
    append json_content "  \"metadata\": {\n"
    append json_content "    \"generated_by\": \"PETK GUI\",\n"
    append json_content "    \"timestamp\": \"[clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]\",\n"
    append json_content "    \"version\": \"1.0\"\n"
    append json_content "  },\n"
    
    # Input section
    append json_content "  \"input\": {\n"
    append json_content "    \"moving_pdb\": \"$pdb_path\"\n"
    append json_content "  },\n"
    
    # Pore geometry section
    append json_content "  \"pore_geometry\": {\n"
    append json_content "    \"pore_type\": \"$::PETK::gui::currentPoreType\",\n"
    
    # Add pore-specific parameters
    if {$::PETK::gui::membraneType eq "cylindrical"} {
        set pore_radius [expr {$::PETK::gui::cylindricalDiameter / 2.0}]
        append json_content "    \"pore_radius\": $pore_radius,\n"
        append json_content "    \"corner_radius\": $::PETK::gui::cornerRadius,\n"
    } elseif {$::PETK::gui::membraneType eq "double_cone"} {
        set inner_radius [expr {$::PETK::gui::innerDiameter / 2.0}]
        set outer_radius [expr {$::PETK::gui::outerDiameter / 2.0}]
        append json_content "    \"inner_radius\": $inner_radius,\n"
        append json_content "    \"outer_radius\": $outer_radius,\n"
    }
    
    append json_content "    \"membrane_thickness\": $::PETK::gui::nanoporeThickness\n"
    append json_content "  },\n"
    
    # Simulation section
    append json_content "  \"simulation\": {\n"
    append json_content "    \"voltage\": $::PETK::gui::appliedVoltage,\n"
    append json_content "    \"bulk_conductivity\": $::PETK::gui::bulkConductivity,\n"
    append json_content "    \"grid_resolution\": $::PETK::gui::gridResolution,\n"
    append json_content "    \"use_vdw_radii\": $use_vdw_json,\n"
    append json_content "    \"default_radius\": $::PETK::gui::semDefaultRadius\n"
    append json_content "  },\n"
    
    # Movement section
    append json_content "  \"movement\": {\n"
    append json_content "    \"z_start\": $::PETK::gui::zStartRange,\n"
    append json_content "    \"z_end\": $::PETK::gui::zEndRange,\n"
    append json_content "    \"z_step\": $::PETK::gui::zStep\n"
    append json_content "  },\n"
    
    # Output section
    append json_content "  \"output\": {\n"
    append json_content "    \"output_prefix\": \"$::PETK::gui::outputPrefix\",\n"
    # append json_content "    \"preview_only\": $preview_only_json,\n"
    append json_content "    \"preview_frames\": $::PETK::gui::semPreviewFrames\n"
    append json_content "  },\n"
    
    # Box dimensions section
    append json_content "  \"box_dimensions\": {\n"
    
    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        # Auto-calculated mode - parse range strings
        set x_range [parseBoxRange $::PETK::gui::autoBoxX]
        set y_range [parseBoxRange $::PETK::gui::autoBoxY] 
        set z_range [parseBoxRange $::PETK::gui::autoBoxZ]
        
        append json_content "    \"x\": $x_range,\n"
        append json_content "    \"y\": $y_range,\n"
        append json_content "    \"z\": $z_range\n"
    } else {
        # Manual mode - convert sizes to ranges (assume centered at origin)
        set x_half [expr {$::PETK::gui::boxSizeX / 2.0}]
        set y_half [expr {$::PETK::gui::boxSizeY / 2.0}]
        set z_half [expr {$::PETK::gui::boxSizeZ / 2.0}]
        
        append json_content "    \"x\": \[-$x_half, $x_half\],\n"
        append json_content "    \"y\": \[-$y_half, $y_half\],\n"
        append json_content "    \"z\": \[-$z_half, $z_half\]\n"
    }
    
    append json_content "  }\n"
    append json_content "}"
    
    # Write the JSON file
    if {[catch {
        set fp [open $output_file w]
        puts $fp $json_content
        close $fp
    } error]} {
        tk_messageBox -icon error -title "Config Export Error" \
            -message "Failed to write config.json:\n\n$error"
        return 0
    }
    
    puts "Created config.json: $output_file"
    return $output_file
}

proc ::PETK::gui::runPreviewFromGUI {} {
    # Run a preview simulation by generating a DCD trajectory showing
    # the analyte moving from zStartRange to zEndRange.
    # Uses corrected VMD animate syntax.
    
    puts "=== SEM PREVIEW SIMULATION ==="
    
    # Use the stored analyte molecule ID
    set analyte_molid $::PETK::gui::analyteMol
    puts "Using stored analyte molecule ID: $analyte_molid"
    puts "Analyte name: [molinfo $analyte_molid get name]"
    
    # Get movement parameters
    set z_start [expr double($::PETK::gui::zStartRange)]
    set z_end [expr double($::PETK::gui::zEndRange)]
    set z_step [expr double($::PETK::gui::zStep)]
    
    # Calculate trajectory parameters
    set num_frames [expr int(abs($z_end - $z_start) / $z_step) + 1]
    set direction [expr {$z_end > $z_start ? 1 : -1}]
    
    puts "Generating trajectory: $num_frames frames from $z_start to $z_end Å"
    puts "Direction: [expr {$direction > 0 ? "forward (+Z)" : "reverse (-Z)"}]"
    puts "Step size: $z_step Å"
    
    # Generate the trajectory using corrected syntax
    if {[::PETK::gui::generateSEMTrajectory $analyte_molid $z_start $z_end $z_step]} {
        # Load and play the trajectory
        ::PETK::gui::loadAndPlaySEMTrajectory $analyte_molid
        
        # Update the molecule name to indicate it now has a trajectory
        mol rename $analyte_molid "Analyte: [file tail $::PETK::gui::analytePDB] (SEM Preview)"
        
        tk_messageBox -type ok -icon info -title "Preview Ready" \
            -message "SEM preview trajectory generated successfully!\n\nFrames: $num_frames\nRange: $z_start to $z_end Å\nStep: $z_step Å\n\nUse the playback controls to view the simulation."
    }
}

proc ::PETK::gui::getAnalyteMoleculeID {} {
    # Find the molecule ID of the loaded analyte using the stored variable.
    # This is much more reliable than searching by filename.
    #
    # Returns:
    #   Molecule ID or -1 if not found/invalid
    
    # First check if we have the stored molecule ID
    if {[info exists ::PETK::gui::analyteMol]} {
        set molid $::PETK::gui::analyteMol
        
        # Verify the molecule still exists in VMD
        if {[lsearch [molinfo list] $molid] != -1} {
            # Double-check that it's actually loaded with frames
            if {[molinfo $molid get numframes] > 0} {
                puts "Found analyte molecule: ID $molid ([molinfo $molid get name])"
                return $molid
            } else {
                puts "Warning: Analyte molecule $molid exists but has no frames"
            }
        } else {
            puts "Warning: Stored analyte molecule ID $molid no longer exists in VMD"
        }
    } else {
        puts "Warning: No analyte molecule ID stored in ::PETK::gui::analyteMol"
    }
    
    # Fallback: Look for a molecule with "Analyte:" in the name
    set mollist [molinfo list]
    foreach molid $mollist {
        set molname [molinfo $molid get name]
        if {[string match "Analyte:*" $molname]} {
            puts "Found analyte by name match: ID $molid ($molname)"
            return $molid
        }
    }
    
    # Final fallback: Look for any molecule that matches the PDB filename
    if {[info exists ::PETK::gui::analytePDB] && $::PETK::gui::analytePDB ne ""} {
        foreach molid $mollist {
            set filename [molinfo $molid get filename]
            if {[string match "*[file tail $::PETK::gui::analytePDB]*" $filename]} {
                puts "Found analyte by filename match: ID $molid"
                return $molid
            }
        }
    }
    
    puts "Error: Could not find analyte molecule in VMD"
    return -1
}

proc ::PETK::gui::generateSEMTrajectory {molid z_start z_end z_step} {
    # Generate a DCD trajectory file with the analyte moving vertically.
    # Uses correct VMD animate command syntax based on official documentation.
    #
    # Args:
    #   molid: VMD molecule ID
    #   z_start: Starting Z position
    #   z_end: Ending Z position  
    #   z_step: Step size
    #
    # Returns:
    #   1 if successful, 0 if failed
    
    # Use success flag to avoid catch/return interaction
    set success 0
    
    if {[catch {
        # Create output directory if it doesn't exist
        set output_dir [file join $::PETK::gui::workdir "sem_preview"]
        file mkdir $output_dir
        
        # Output files
        set dcd_file [file join $output_dir "analyte_trajectory.dcd"]
        set pdb_file [file join $output_dir "analyte_trajectory.pdb"]
        
        puts "Generating trajectory files:"
        puts "  PDB: $pdb_file"
        puts "  DCD: $dcd_file"
        
        # Get original coordinates
        set all_atoms [atomselect $molid "all"]
        set original_coords [$all_atoms get {x y z}]
        set num_atoms [$all_atoms num]
        
        # Calculate center of mass in Z
        set com_z [lindex [measure center $all_atoms] 2]
        puts "Original center of mass Z: $com_z Å"
        
        # Calculate trajectory parameters
        set direction [expr {$z_end > $z_start ? 1 : -1}]
        set num_frames [expr int(abs($z_end - $z_start) / $z_step) + 1]
        
        puts "Trajectory parameters:"
        puts "  Frames: $num_frames"
        puts "  Direction: [expr {$direction > 0 ? "forward" : "reverse"}]"
        puts "  Step size: $z_step Å"
        
        # Delete existing frames except the first one
        set current_frames [molinfo $molid get numframes]
        if {$current_frames > 1} {
            animate delete beg 1 end -1 $molid
        }
        
        # Generate frames
        for {set frame 0} {$frame < $num_frames} {incr frame} {
            # Calculate Z position for this frame
            set current_z [expr $z_start + ($frame * $z_step * $direction)]
            set z_offset [expr $current_z - $com_z]
            
            if {$frame > 0} {
                # Add a new frame by duplicating the current frame
                animate dup $molid
            }
            
            # Move to the new frame
            animate goto $frame
            
            # Get atom selection for current frame
            set frame_atoms [atomselect $molid "all" frame $frame]
            
            # Reset to original coordinates first
            $frame_atoms set {x y z} $original_coords
            
            # Apply Z translation
            $frame_atoms moveby [list 0 0 $z_offset]
            
            $frame_atoms delete
            
            puts "  Frame [expr $frame + 1]/$num_frames: Z = [format "%.2f" $current_z] Å"
        }
        
        # Write trajectory files
        puts "Writing trajectory files..."
        
        # Write PDB (frame 0 as reference)
        animate goto 0
        $all_atoms writepdb $pdb_file
        
        # Write DCD trajectory
        animate write dcd $dcd_file beg 0 end -1 $molid
        
        $all_atoms delete
        
        # Store trajectory info for later use
        set ::PETK::gui::currentTrajectoryDCD $dcd_file
        set ::PETK::gui::currentTrajectoryPDB $pdb_file
        set ::PETK::gui::currentTrajectoryFrames $num_frames
        
        puts "Trajectory generation completed successfully!"
        
        # Set success flag instead of returning
        set success 1
        
    } error]} {
        puts "ERROR generating trajectory: $error"
        tk_messageBox -type ok -icon error -title "Trajectory Generation Failed" \
            -message "Failed to generate trajectory:\n\n$error"
        set success 0
    }
    
    # Return success status outside the catch block
    return $success
}

proc ::PETK::gui::loadAndPlaySEMTrajectory {molid} {
    # Load the generated trajectory and set up for playback.
    # Uses correct VMD syntax.
    #
    # Args:
    #   molid: VMD molecule ID
    
    if {![info exists ::PETK::gui::currentTrajectoryDCD]} {
        puts "No trajectory file to load"
        return
    }
    
    if {[catch {
        # First, delete any existing trajectory frames from the molecule (keep frame 0)
        set current_frames [molinfo $molid get numframes]
        if {$current_frames > 1} {
            animate delete beg 1 end -1 $molid
        }
        
        # Load the DCD file
        animate read dcd $::PETK::gui::currentTrajectoryDCD $molid
        
        puts "Trajectory loaded: [molinfo $molid get numframes] frames"
        
        # Set up display
        animate goto 0
        display update
        
        # Center the view
        display resetview
        
        # Set up nice representations if needed
        set num_reps [molinfo $molid get numreps]
        if {$num_reps == 0} {
            mol representation NewCartoon
            mol color ColorID 0
            mol addrep $molid
            
            mol representation VDW 0.3
            mol color Element
            mol selection "not protein"
            mol addrep $molid
        }
        
        puts "Trajectory display configured"
        
    } error]} {
        puts "ERROR loading trajectory: $error"
        tk_messageBox -type ok -icon error -title "Trajectory Load Failed" \
            -message "Failed to load trajectory:\n\n$error"
    }
}
proc ::PETK::gui::runCalculationFromGUI {} {
    # Generate config
    ::PETK::gui::outputParametersToConfig
    
    # Run calculations
    set config_file [file join $::PETK::gui::workdir "config.json"]
    ::PETK::gui::run_sem_calculation $config_file
}

# Main procedure for running SEM calculations
proc ::PETK::gui::run_sem_calculation {config_file} {
    # Check if validation has passed
    if {![info exists ::PETK::gui::semValidationPassed] || $::PETK::gui::semValidationPassed != 1} {
        set error_msg "Cannot run calculations: SEM setup validation has not passed.\nPlease run validation first and fix any errors."
        puts $error_msg
        tk_messageBox -icon error -title "Calculation Error" -message $error_msg
        return -code error $error_msg
    }
    
    # Create results folder in workdir
    set results_dir [file join $::PETK::gui::workdir "results"]
    
    # Create the results directory if it doesn't exist
    if {![file exists $results_dir]} {
        file mkdir $results_dir
        puts "Created results directory: $results_dir"
    } else {
        puts "Results directory already exists: $results_dir"
    }
    
    # Change to results directory
    set original_dir [pwd]
    cd $results_dir
    puts "Changed to results directory: [pwd]"
    
    # Copy config file to results directory if it's not already there
    set config_basename [file tail $config_file]
    set local_config [file join $results_dir $config_basename]
    
    if {![file exists $local_config] || [file mtime $config_file] > [file mtime $local_config]} {
        file copy -force $config_file $local_config
        puts "Copied config file to results directory"
    }
    
    # Construct the command to run full SEM calculation
    set sem_script_path [file join $::PETK::gui::workdir "vertical_movement_sem.py"]
    # Note: Removed "preview_only" flag for full calculation
    set cmd [list conda run -n $::PETK::gui::condaEnvironment python $sem_script_path $config_basename run]
    
    puts "Running SEM calculation with command: $cmd"
    puts "========================================="
    puts "Note: This may take some time to complete..."
    
    # Update GUI to show calculation is running (optional)
    # You might want to disable buttons or show progress here
    
    # Execute the command and capture output
    if {[catch {
        set result [exec {*}$cmd 2>&1]
        puts $result
        puts "========================================="
        puts "SEM calculation completed successfully!"
        
        # List generated files
        set result_files [glob -nocomplain "*.png" "*.dat" "*.txt" "*.csv" "*.out"]
        if {[llength $result_files] > 0} {
            puts "Generated files:"
            foreach file $result_files {
                puts "  - $file"
            }
        } else {
            puts "No result files found"
        }
        
        # Optional: Show completion message to user
        tk_messageBox -icon info -title "Calculation Complete" -message "SEM calculation completed successfully!\nResults saved in: $results_dir"
        
    } error]} {
        puts "Error running SEM calculation:"
        puts $error
        cd $original_dir
        tk_messageBox -icon error -title "Calculation Error" -message "Error running SEM calculation:\n$error"
        return -code error $error
    }
    
    # Return to original directory
    cd $original_dir
    puts "Returned to original directory: [pwd]"
    
    return $results_dir
}



####################################################
# Scroll Region Update Function (from Tab 2)
####################################################

proc ::PETK::gui::onCanvasConfigured {canvas canvas_window} {
    set canvas_width [winfo width $canvas]
    set canvas_height [winfo height $canvas]
    
    # Skip if canvas not ready
    if {$canvas_width <= 1 || $canvas_height <= 1} {
        return
    }
    
    # Update content frame sizing
    ::PETK::gui::resizeCanvasContent $canvas $canvas_window
}

proc ::PETK::gui::onContentConfigured {canvas canvas_window} {
    after idle [list ::PETK::gui::resizeCanvasContent $canvas $canvas_window]
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

proc ::PETK::gui::forceInitialCanvasResize {canvas canvas_window} {
    # Force initial update
    update idletasks
    
    # Call resize function
    ::PETK::gui::resizeCanvasContent $canvas $canvas_window
    
    # Schedule additional updates to ensure proper sizing
    after 100 [list ::PETK::gui::resizeCanvasContent $canvas $canvas_window]
    after 500 [list ::PETK::gui::resizeCanvasContent $canvas $canvas_window]
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