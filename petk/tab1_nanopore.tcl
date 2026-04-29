#
# PETK GUI - Tab 1 (Nanopore Setup) module
#
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
    ttk::label $container.boxdim.cutofflbl -text "Distance padding (Å):" -width 18
    ttk::entry $container.boxdim.cutoff -textvariable ::PETK::gui::sysPadding -width 12 -justify center

    grid $container.boxdim.auto $container.boxdim.cutofflbl $container.boxdim.cutoff -columnspan 5 -sticky w -pady "0 5"

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

    ttk::label $container.boxdim.auto_display.aspectlbl -text "Bulk height vs. lateral:" -width 22
    ttk::label $container.boxdim.auto_display.aspectval -textvariable ::PETK::gui::autoBoxAspectSummary \
        -width 40 -anchor w -relief sunken -background white
    grid $container.boxdim.auto_display.aspectlbl $container.boxdim.auto_display.aspectval - - - - -sticky ew -pady 2

    ttk::checkbutton $container.boxdim.auto_display.golden \
        -text "Enforce golden aspect ratio (H/L = 1.2)" \
        -variable ::PETK::gui::useGoldenAspectRatio \
        -command {::PETK::gui::calculateBoxDimensions}
    grid $container.boxdim.auto_display.golden -columnspan 6 -sticky w -pady "6 0"

    ttk::button $container.boxdim.auto_display.calc -text "Recalculate Box" -command {::PETK::gui::calculateBoxDimensions}
    grid $container.boxdim.auto_display.calc -row 3 -column 10 -columnspan 2 -sticky w -pady "5 0"
    
    # Grid the appropriate frame based on mode
    if {$::PETK::gui::autoCalculateBoxDimensions} {
        grid $container.boxdim.auto_display -row 1 -column 0 -columnspan 6 -sticky ew
    } else {
        grid $container.boxdim.manual -row 1 -column 0 -columnspan 6 -sticky ew
    }

    # === PORE TYPE SELECTION (SOLID-STATE vs BIOLOGICAL) ===
    ttk::labelframe $container.poreoption -text "Pore Type Selection" -padding 10
    grid $container.poreoption -row $row -column 0 -sticky ew -padx 10 -pady 5
    grid columnconfigure $container.poreoption {0 1} -weight 1
    incr row

    # Initialize pore option variable if not exists
    if {![info exists ::PETK::gui::poreOption]} {
        set ::PETK::gui::poreOption "solid-state"
    }

    ttk::radiobutton $container.poreoption.solidstate -text "Solid-State Pore" -value "solid-state" \
        -variable ::PETK::gui::poreOption -command {::PETK::gui::transitionToPoreOption "solid-state"}
    ttk::radiobutton $container.poreoption.biological -text "Biological Pore" -value "biological" \
        -variable ::PETK::gui::poreOption -command {::PETK::gui::transitionToPoreOption "biological"}

    grid $container.poreoption.solidstate $container.poreoption.biological -sticky w -pady 3 -padx 20

    # === BIOLOGICAL PORE SELECTION SECTION ===
    ttk::labelframe $container.biopore -text "Biological Pore Selection" -padding 10
    # Don't grid it initially - will be shown/hidden by updatePoreOptionDisplay
    grid columnconfigure $container.biopore {1 2} -weight 1

    # Pre-prepared biological pores from local folder
    ttk::label $container.biopore.preselectlbl -text "Select pre-prepared pore:" -width 20
    ttk::combobox $container.biopore.preselect -textvariable ::PETK::gui::selectedBioPore -width 30 -state readonly
    ttk::button $container.biopore.refresh -text "Refresh List" -command {::PETK::gui::refreshBioPoreList}

    grid $container.biopore.preselectlbl $container.biopore.preselect $container.biopore.refresh -sticky ew -pady 3

    # Custom biological pore upload
    ttk::separator $container.biopore.sep1 -orient horizontal
    grid $container.biopore.sep1 -columnspan 3 -sticky ew -pady 10

    ttk::label $container.biopore.customlbl -text "Or upload custom pore:" -width 20
    ttk::entry $container.biopore.custompath -textvariable ::PETK::gui::customBioPorePath -width 35
    ttk::button $container.biopore.browse -text "Browse..." -command {::PETK::gui::selectCustomBioPore}

    grid $container.biopore.customlbl $container.biopore.custompath $container.biopore.browse -sticky ew -pady 3
    
    # Pore information display
    ttk::separator $container.biopore.sep2 -orient horizontal
    grid $container.biopore.sep2 -columnspan 3 -sticky ew -pady 10

    ttk::label $container.biopore.infolbl -text "Pore Configuration:" -font {TkDefaultFont 10 bold}
    grid $container.biopore.infolbl -columnspan 3 -sticky w -pady "0 5"

    text $container.biopore.info -height 6 -width 60 -wrap word -state disabled \
        -font {TkDefaultFont 9}
    ttk::scrollbar $container.biopore.infoscroll -orient vertical -command [list $container.biopore.info yview]
    $container.biopore.info configure -yscrollcommand [list $container.biopore.infoscroll set]

    grid $container.biopore.info -row 6 -column 0 -columnspan 2 -sticky ew -pady 2
    grid $container.biopore.infoscroll -row 6 -column 2 -sticky ns -pady 2

    # === BIOLOGICAL MEMBRANE PARAMETERS SECTION ===
    ttk::labelframe $container.biomembrane -text "Solid-State Membrane Parameters" -padding 10
    # This should be shown only when biological pore is selected
    grid columnconfigure $container.biomembrane {0 1} -weight 1

    # Information label about fixed parameters
    ttk::label $container.biomembrane.fixedlbl -text "Box size and membrane thickness are automatically set based on pore type" \
        -font {TkDefaultFont 8} -foreground blue -wraplength 400
    grid $container.biomembrane.fixedlbl -columnspan 2 -sticky ew -pady "0 10" -row 0

    # Membrane thickness control (read-only for biological pores)
    ttk::label $container.biomembrane.thicklbl -text "Membrane thickness (Å):" -width 20
    ttk::spinbox $container.biomembrane.thick -textvariable ::PETK::gui::nanoporeThickness \
        -from 5 -to 100 -increment 1.0 -width 15 -state readonly
    ttk::label $container.biomembrane.thickinfo -text "Fixed - optimized for this biological pore" \
        -font {TkDefaultFont 8} -foreground gray

    grid $container.biomembrane.thicklbl $container.biomembrane.thick -sticky ew -pady 2 -row 1
    grid $container.biomembrane.thickinfo -columnspan 2 -sticky w -pady "0 5" -row 2

    # Membrane Z-offset control (user-modifiable)
    ttk::label $container.biomembrane.zoffsetlbl -text "Membrane Z-offset (Å):" -width 20
    ttk::spinbox $container.biomembrane.zoffset -textvariable ::PETK::gui::membraneZOffset \
        -from -100 -to 100 -increment 1.0 -width 15 \
        -command {::PETK::gui::updatePoreVmdVisualization}
    ttk::label $container.biomembrane.zoffsetinfo -text "Fixed - optimized for this biological pore" \
        -font {TkDefaultFont 8} -foreground gray

    grid $container.biomembrane.zoffsetlbl $container.biomembrane.zoffset -sticky ew -pady 2 -row 3
    grid $container.biomembrane.zoffsetinfo -columnspan 2 -sticky w -pady "0 5" -row 4

    # Membrane visualization options
    ttk::separator $container.biomembrane.sep -orient horizontal
    grid $container.biomembrane.sep -columnspan 2 -sticky ew -pady 5 -row 5

    # ttk::checkbutton $container.biomembrane.showmem -text "Show membrane outline" \
    #     -variable ::PETK::gui::showMembraneOutline \
    #     -command {::PETK::gui::updatePoreVmdVisualization}
    # ttk::checkbutton $container.biomembrane.showrad -text "Show local radius calculation" \
    #     -variable ::PETK::gui::showRadiusCalculation \
    #     -command {::PETK::gui::updatePoreVmdVisualization}

    # grid $container.biomembrane.showmem -columnspan 2 -sticky w -pady 2 -row 6
    #grid $container.biomembrane.showrad -columnspan 2 -sticky w -pady 2 -row 7

    # Membrane analysis button
    ttk::button $container.biomembrane.analyze -text "Update Prview" \
        -command {::PETK::gui::updatePoreVisualization}
    grid $container.biomembrane.analyze -columnspan 2 -pady 5 -row 8

    # === PORE TYPE SELECTION (NANOPORE GEOMETRY - only for solid-state) ===
    ttk::labelframe $container.poretype -text "Solid-State Nanopore Type Selection" -padding 10
    # Don't grid initially - will be shown/hidden by updatePoreOptionDisplay
    grid columnconfigure $container.poretype {0 1} -weight 1

    ttk::radiobutton $container.poretype.cylindrical -text "Cylindrical Pore" -value "cylindrical" \
        -variable ::PETK::gui::membraneType -command {::PETK::gui::updateMembraneTypeDisplay}
    ttk::radiobutton $container.poretype.doublecone -text "Double Cone Pore" -value "doublecone" \
        -variable ::PETK::gui::membraneType -command {::PETK::gui::updateMembraneTypeDisplay}
    ttk::radiobutton $container.poretype.conical -text "Conical Pore" -value "conical" \
        -variable ::PETK::gui::membraneType -command {::PETK::gui::updateMembraneTypeDisplay}

    grid $container.poretype.cylindrical $container.poretype.doublecone $container.poretype.conical -sticky w -pady 3 -padx 20

    # === PORE PARAMETERS SECTION (only for solid-state) ===
    ttk::labelframe $container.params -text "Pore Parameters" -padding 10
    # Don't grid initially - will be shown/hidden by updatePoreOptionDisplay
    grid columnconfigure $container.params {1 3 5} -weight 1

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

    # Conical (single truncated cone / frustum) pore parameters
    ttk::frame $container.params.conical
    ttk::label $container.params.conical.toplbl -text "Top diameter (Å):" -width 20
    ttk::entry $container.params.conical.top -textvariable ::PETK::gui::topDiameter -width 12 -justify center
    ttk::label $container.params.conical.botlbl -text "Bottom diameter (Å):" -width 18
    ttk::entry $container.params.conical.bot -textvariable ::PETK::gui::bottomDiameter -width 12 -justify center

    grid $container.params.conical.toplbl $container.params.conical.top $container.params.conical.botlbl $container.params.conical.bot -sticky ew -pady 3 -padx 5

    # === VISUALIZATION SECTION ===
    ttk::labelframe $container.visualization -text "Pore Geometry Visualization" -padding 10
    grid $container.visualization -row 7 -column 0 -sticky nsew -padx 10 -pady 5
    grid columnconfigure $container.visualization {0 1} -weight 1
    grid rowconfigure $container.visualization 1 -weight 1

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
    bind $container.visualization.image.canvas <Configure> \
        [list ::PETK::gui::onPoreImageCanvasConfigure $container.visualization.image.canvas]
    
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
    ttk::button $container.visualization.summary.update -text "Update Preview" -command {::PETK::gui::updatePreviewAndMovement}
    grid $container.visualization.summary.update -sticky ew -pady "10 0"

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab1.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab1.canvas $canvas_window $container]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab1.canvas $canvas_window $container]
    
    # Store references for later use (Tab 1 specific)
    set ::PETK::gui::tab1MainCanvas $tab1.canvas
    set ::PETK::gui::tab1CanvasWindow $canvas_window
    set ::PETK::gui::tab1ContentContainer $container
    # Maintain legacy variable names used elsewhere
    set ::PETK::gui::mainCanvas $tab1.canvas
    set ::PETK::gui::canvasWindow $canvas_window
    set ::PETK::gui::contentContainer $container
    
    # === CRITICAL: Initialize biological pore configuration system ===
    ::PETK::gui::initializeBiologicalPoreConfigurations
    
    # Add binding for biological pore selection - CRITICAL for auto-configuration
    bind $container.biopore.preselect <<ComboboxSelected>> {::PETK::gui::onBiologicalPoreSelected}
    
    # Initialize variables and display
    ::PETK::gui::initializeBiologicalPoreVariables
    ::PETK::gui::updateMembraneTypeDisplay
    ::PETK::gui::loadPoreImages
    ::PETK::gui::updatePoreOptionDisplay
    ::PETK::gui::refreshBioPoreList
    
    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab1.canvas $canvas_window $container]
    after idle [list ::PETK::gui::updateTab1ScrollRegion]
}

####################################################
# Variable Initialization
####################################################
proc ::PETK::gui::initializeBioPoreVariables {} {
    if {![info exists ::PETK::gui::poreOption]} {
        set ::PETK::gui::poreOption "solid-state"
    }
    if {![info exists ::PETK::gui::selectedBioPore]} {
        set ::PETK::gui::selectedBioPore ""
    }
    if {![info exists ::PETK::gui::customBioPorePath]} {
        set ::PETK::gui::customBioPorePath ""
    }
}

proc ::PETK::gui::initializeMembraneVariables {} {
    if {[info exists ::PETK::gui::membraneVarsInitialized] && $::PETK::gui::membraneVarsInitialized} {
        return
    }
    # Initialize pore type
    if {![info exists ::PETK::gui::membraneType]} {
        set ::PETK::gui::membraneType "cylindrical"
    }
    if {![info exists ::PETK::gui::useGoldenAspectRatio]} {
        set ::PETK::gui::useGoldenAspectRatio 0
    }
    if {![info exists ::PETK::gui::autoBoxAspectSummary]} {
        set ::PETK::gui::autoBoxAspectSummary "L=—, H=— (H/L=—)"
    }

    if {![info exists ::PETK::gui::poreImageRedrawDelay]} {
        set ::PETK::gui::poreImageRedrawDelay 90
    }
    if {![info exists ::PETK::gui::poreImageScaleGranularity]} {
        set ::PETK::gui::poreImageScaleGranularity 8
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

    # Initialize conical (single frustum) parameters
    if {![info exists ::PETK::gui::topDiameter]} {
        set ::PETK::gui::topDiameter "240.0"
    }
    if {![info exists ::PETK::gui::bottomDiameter]} {
        set ::PETK::gui::bottomDiameter "120.0"
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

    # Define solid-state presets and storage for user overrides
    set ::PETK::gui::solidStateDefaults(cylindrical) [dict create \
        autoCalculateBox 0 \
        boxSizeX "300.0" \
        boxSizeY "300.0" \
        boxSizeZ "300.0" \
        nanoporeThickness "200.0" \
        cylindricalDiameter "200.0" \
        cornerRadius "50.0"]

    set ::PETK::gui::solidStateDefaults(doublecone) [dict create \
        autoCalculateBox 0 \
        boxSizeX "300.0" \
        boxSizeY "300.0" \
        boxSizeZ "300.0" \
        nanoporeThickness "200.0" \
        innerDiameter "100.0" \
        outerDiameter "300.0"]

    set ::PETK::gui::solidStateDefaults(conical) [dict create \
        autoCalculateBox 0 \
        boxSizeX "300.0" \
        boxSizeY "300.0" \
        boxSizeZ "300.0" \
        nanoporeThickness "200.0" \
        topDiameter "240.0" \
        bottomDiameter "120.0"]

    if {![info exists ::PETK::gui::solidStateState(cylindrical)]} {
        set ::PETK::gui::solidStateState(cylindrical) ""
    }
    if {![info exists ::PETK::gui::solidStateState(doublecone)]} {
        set ::PETK::gui::solidStateState(doublecone) ""
    }
    if {![info exists ::PETK::gui::solidStateState(conical)]} {
        set ::PETK::gui::solidStateState(conical) ""
    }
    if {![info exists ::PETK::gui::previousMembraneType]} {
        set ::PETK::gui::previousMembraneType $::PETK::gui::membraneType
    }
    if {![info exists ::PETK::gui::forceSolidStateDefaults]} {
        set ::PETK::gui::forceSolidStateDefaults 0
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
    
    set ::PETK::gui::membraneVarsInitialized 1
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

proc ::PETK::gui::transitionToPoreOption {new_option} {
    set old_option ""
    if {[info exists ::PETK::gui::poreOption]} {
        set old_option $::PETK::gui::poreOption
    }

    # Persist current solid-state configuration when leaving the mode
    if {$old_option eq "solid-state" && $new_option ne "solid-state"} {
        ::PETK::gui::persistSolidStateParameters
        set ::PETK::gui::previousMembraneType ""
    }

    # Decide whether we need to reset to the solid-state defaults
    set apply_preset 0
    set preset_type ""
    if {$old_option eq "biological" && $new_option eq "solid-state"} {
        set apply_preset 1
        if {![info exists ::PETK::gui::membraneType] || $::PETK::gui::membraneType eq "" || $::PETK::gui::membraneType ni {cylindrical doublecone conical}} {
            set ::PETK::gui::membraneType "cylindrical"
        }
        set preset_type $::PETK::gui::membraneType
    }
    
    # Update the option
    set ::PETK::gui::poreOption $new_option
    
    puts "Transitioning from '$old_option' to '$new_option'"
    
    if {$apply_preset} {
        if {$preset_type eq ""} {
            set preset_type "cylindrical"
        }
        set ::PETK::gui::membraneType $preset_type
        set ::PETK::gui::forceSolidStateDefaults 1
    }
    
    # Perform transition-specific cleanup
    if {$old_option eq "biological" && $new_option eq "solid-state"} {
        # Clean up biological pore resources
        ::PETK::gui::cleanupBiologicalPoreResources
        # Restore solid-state controls
        ::PETK::gui::restoreSolidStateGUIControls
    } elseif {$old_option eq "solid-state" && $new_option eq "biological"} {
        # Clean up solid-state resources if needed
        # No specific cleanup needed for solid-state currently
    }
    
    # Update the display
    ::PETK::gui::updatePoreOptionDisplay
    if {[info procs ::PETK::gui::updateBulkElectrolyteButtonState] ne ""} {
        ::PETK::gui::updateBulkElectrolyteButtonState
    }
    
    # Remember the new selection for future transitions
    set ::PETK::gui::previousPoreOption $new_option
    
    # Force complete layout refresh (coalesced)
    ::PETK::gui::requestCanvasLayoutRefresh
}

proc ::PETK::gui::updatePoreOptionDisplay {} {
    if {![info exists ::PETK::gui::contentContainer]} {
        return
    }
    
    set container $::PETK::gui::contentContainer
    set __petkPaused 0
    if {[info procs ::PETK::gui::pausePoreVisualization] ne ""} {
        ::PETK::gui::pausePoreVisualization
        set __petkPaused 1
    }
    
    # Manage all sections and their row positions based on pore type
    if {$::PETK::gui::poreOption eq "biological"} {
        # === BIOLOGICAL PORE LAYOUT ===
        # Show biological pore section at row 3
        grid $container.biopore -row 3 -column 0 -sticky ew -padx 10 -pady 5
        
        # Show biological membrane parameters at row 4
        grid $container.biomembrane -row 4 -column 0 -sticky ew -padx 10 -pady 5
        
        # Hide solid-state specific sections
        grid remove $container.poretype  
        grid remove $container.params
        
        # HIDE visualization section for biological pores
        grid remove $container.visualization
        
        # Configure row weights for biological layout (no visualization)
        grid rowconfigure $container 3 -weight 0  ;# biopore section
        grid rowconfigure $container 4 -weight 1  ;# biomembrane section (expandable)
        grid rowconfigure $container 5 -weight 0  ;# clear unused rows
        grid rowconfigure $container 6 -weight 0
        grid rowconfigure $container 7 -weight 0
        
        # Update section titles for biological mode
        if {[winfo exists $container.biopore]} {
            $container.biopore configure -text "Biological Pore Selection"
        }
        if {[winfo exists $container.biomembrane]} {
            $container.biomembrane configure -text "Biological Pore Membrane Parameters"
        }
        
        puts "Biological pore mode: Visualization section hidden"
        
    } else {
        # === SOLID-STATE PORE LAYOUT ===
        # Hide biological pore sections
        grid remove $container.biopore
        grid remove $container.biomembrane
        
        # Show solid-state sections
        grid $container.poretype -row 3 -column 0 -sticky ew -padx 10 -pady 5
        grid $container.params -row 4 -column 0 -sticky ew -padx 10 -pady 5
        
        # SHOW visualization section for solid-state pores
        grid $container.visualization -row 5 -column 0 -sticky nsew -padx 10 -pady 5
        
        # Configure row weights for solid-state layout (with visualization)
        grid rowconfigure $container 3 -weight 0  ;# poretype section
        grid rowconfigure $container 4 -weight 0  ;# params section
        grid rowconfigure $container 5 -weight 1  ;# visualization (expandable)
        grid rowconfigure $container 6 -weight 0  ;# clear unused rows
        grid rowconfigure $container 7 -weight 0
        
        # Update section titles for solid-state mode
        if {[winfo exists $container.poretype]} {
            $container.poretype configure -text "Solid-State Nanopore Type Selection"
        }
        if {[winfo exists $container.params]} {
            $container.params configure -text "Pore Parameters"
        }
        if {[winfo exists $container.visualization]} {
            $container.visualization configure -text "Pore Geometry Visualization"
        }
        
        # Ensure solid-state GUI controls are restored
        ::PETK::gui::restoreSolidStateGUIControls
        
        puts "Solid-state pore mode: Visualization section visible"
    }
    
    # Force layout update after idle to avoid repeated synchronous refresh
    ::PETK::gui::requestCanvasLayoutRefresh
    
    # Update the membrane type display if we're in solid-state mode
    if {$::PETK::gui::poreOption eq "solid-state"} {
        ::PETK::gui::updateMembraneTypeDisplay
    }
    
    # Update visualization (only matters for solid-state now)
    if {[info procs "::PETK::gui::updatePoreVisualization"] ne ""} {
        ::PETK::gui::updatePoreVisualization
    }
    
    if {$__petkPaused} {
        ::PETK::gui::resumePoreVisualization
    }
    
    puts "Pore option display updated: $::PETK::gui::poreOption mode"
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
    set movementUpdate 0
    if {[info exists ::PETK::gui::forceMovementUpdate] && $::PETK::gui::forceMovementUpdate} {
        set movementUpdate 1
    }
    
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "biological"} {
        set box_x 0.0
        set box_y 0.0
        set box_z 0.0
        
        if {[info exists ::PETK::gui::customBioBoundingInfo] && $::PETK::gui::customBioBoundingInfo ne ""} {
            if {[catch {set box_x [expr {double($::PETK::gui::boxSizeX)}]}]} { set box_x 0.0 }
            if {[catch {set box_y [expr {double($::PETK::gui::boxSizeY)}]}]} { set box_y 0.0 }
            if {[catch {set box_z [expr {double($::PETK::gui::boxSizeZ)}]}]} { set box_z 0.0 }
            
            ::PETK::gui::updateAutoBoxDisplayFromDimensions $box_x $box_y $box_z "" "" $movementUpdate
            puts "Auto-calculated box dimensions (custom biological pore):"
            puts [format "  Using stored custom bounding box sizes: %.1f × %.1f × %.1f Å" $box_x $box_y $box_z]
            puts [format "  Movement presets: %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
            return
        } elseif {[info exists ::PETK::gui::boxSizeX] && [info exists ::PETK::gui::boxSizeY] && [info exists ::PETK::gui::boxSizeZ]} {
            if {[catch {set box_x [expr {double($::PETK::gui::boxSizeX)}]}]} { set box_x 0.0 }
            if {[catch {set box_y [expr {double($::PETK::gui::boxSizeY)}]}]} { set box_y 0.0 }
            if {[catch {set box_z [expr {double($::PETK::gui::boxSizeZ)}]}]} { set box_z 0.0 }
            
            ::PETK::gui::updateAutoBoxDisplayFromDimensions $box_x $box_y $box_z "" "" $movementUpdate
            puts "Auto-calculated box dimensions (biological pore):"
            puts [format "  Using current biological pore box: %.1f × %.1f × %.1f Å" $box_x $box_y $box_z]
            puts [format "  Movement presets: %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
            return
        }
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
    } elseif {$::PETK::gui::membraneType eq "conical"} {
        # Use the larger of top/bottom diameter as the bounding radius
        set top_d 0.0
        set bot_d 0.0
        if {[string is double $::PETK::gui::topDiameter]} {
            set top_d [expr {$::PETK::gui::topDiameter}]
        }
        if {[string is double $::PETK::gui::bottomDiameter]} {
            set bot_d [expr {$::PETK::gui::bottomDiameter}]
        }
        set pore_radius [expr {max($top_d, $bot_d) / 2.0}]
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

    set z_padding 0.0
    if {[string is double $::PETK::gui::sysPadding]} {
        set z_padding [expr {max(0.0, double($::PETK::gui::sysPadding))}]
    }
    set z_min [expr {$z_min - $z_padding}]
    set z_max [expr {$z_max + $z_padding}]
    
    set box_x [expr {2.0 * $xy_size}]
    set box_y [expr {2.0 * $xy_size}]
    set box_z [expr {$z_max - $z_min}]
    if {$box_z < 0} {
        set box_z [expr {abs($box_z)}]
    }
    
    lassign [::PETK::gui::maybeApplyGoldenAspectRatio $box_x $box_y $box_z $z_min $z_max] box_x box_y box_z z_min z_max
    ::PETK::gui::updateAutoBoxDisplayFromDimensions $box_x $box_y $box_z $z_min $z_max $movementUpdate
        
    set half_x [expr {$box_x / 2.0}]
    set half_y [expr {$box_y / 2.0}]

    puts "Auto-calculated box dimensions:"
    puts [format "  X: %.1f to %.1f Å" [expr {-$half_x}] $half_x]
    puts [format "  Y: %.1f to %.1f Å" [expr {-$half_y}] $half_y]
    puts [format "  Z: %.1f to %.1f Å" $z_min $z_max]
    puts "  Pore radius: $pore_radius Å"
    puts "  Movement range: $z_end to $z_start Å"
    if {![info exists ::PETK::gui::semCalculationMode] || $::PETK::gui::semCalculationMode eq "run"} {
        puts [format "  Movement presets: %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
    }
}

proc ::PETK::gui::maybeApplyGoldenAspectRatio {box_x box_y box_z z_min z_max} {
    if {![info exists ::PETK::gui::useGoldenAspectRatio] || !$::PETK::gui::useGoldenAspectRatio} {
        return [list $box_x $box_y $box_z $z_min $z_max]
    }
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "biological"} {
        return [list $box_x $box_y $box_z $z_min $z_max]
    }

    set lateral [expr {max($box_x, $box_y)}]
    if {$lateral <= 0} {
        return [list $box_x $box_y $box_z $z_min $z_max]
    }

    set thickness 0.0
    if {[string is double $::PETK::gui::nanoporeThickness]} {
        set thickness [expr {double($::PETK::gui::nanoporeThickness)}]
    }

    set bulk_height [expr {$box_z - $thickness}]
    if {$bulk_height <= 0} {
        return [list $box_x $box_y $box_z $z_min $z_max]
    }

    set alpha 1.2
    set tolerance 0.01
    set ratio [expr {$bulk_height / double($lateral)}]

    set new_box_x $box_x
    set new_box_y $box_y
    set new_box_z $box_z

    if {$ratio > $alpha + $tolerance} {
        set needed_lateral [expr {$bulk_height / $alpha}]
        if {$needed_lateral > $lateral} {
            set new_box_x $needed_lateral
            set new_box_y $needed_lateral
        }
    } elseif {$ratio < $alpha - $tolerance} {
        set needed_height [expr {$alpha * $lateral}]
        if {$needed_height > $bulk_height} {
            set delta [expr {$needed_height - $bulk_height}]
            set new_box_z [expr {$box_z + $delta}]
            set bulk_height $needed_height
        }
    }

    set lateral_final [expr {max($new_box_x, $new_box_y)}]
    if {$lateral_final > 0} {
        set new_box_x $lateral_final
        set new_box_y $lateral_final
    }

    if {![string is double -strict $z_min] || ![string is double -strict $z_max]} {
        set center_z 0.0
    } else {
        set center_z [expr {($z_max + $z_min) / 2.0}]
    }
    set half_z [expr {$new_box_z / 2.0}]
    set z_min_new [expr {$center_z - $half_z}]
    set z_max_new [expr {$center_z + $half_z}]

    return [list $new_box_x $new_box_y $new_box_z $z_min_new $z_max_new]
}

####################################################
# Biological Pore
####################################################
proc ::PETK::gui::refreshBioPoreList {} {
    set container $::PETK::gui::contentContainer
    
    set search_dirs [::PETK::gui::getBioPoreSearchDirs]
    
    $container.biopore.preselect configure -values {}
    set ::PETK::gui::selectedBioPore ""
    set ::PETK::gui::bioPorePathMap [dict create]
    
    set file_types {*.pdb *.gro *.xyz *.mol2 *.sdf}
    set pore_names {}
    set path_map [dict create]
    
    foreach dir $search_dirs {
        foreach pattern $file_types {
            foreach file [glob -nocomplain -directory $dir $pattern] {
                set name [file tail $file]
                if {![dict exists $path_map $name]} {
                    dict set path_map $name $file
                    lappend pore_names $name
                }
            }
        }
    }
    
    set ::PETK::gui::bioPorePathMap $path_map
    set pore_names [lsort $pore_names]
    $container.biopore.preselect configure -values $pore_names
    
    if {[llength $pore_names] > 0} {
        set ::PETK::gui::selectedBioPore [lindex $pore_names 0]
        ::PETK::gui::updateBioPoreInfo
    } else {
        set search_str [join $search_dirs ", "]
        if {$search_str eq ""} {
            set search_str "bio_pore/"
        }
        ::PETK::gui::updateBioPoreInfo "No biological pore files found in: $search_str"
    }
}

proc ::PETK::gui::getBioPoreSearchDirs {} {
    set dirs {}
    set user_dir [file join [pwd] "bio_pore"]
    if {![file exists $user_dir]} {
        if {[catch {file mkdir $user_dir} err]} {
            puts "Warning: Cannot create user bio_pore directory at $user_dir: $err"
        }
    }
    if {[file isdirectory $user_dir]} {
        lappend dirs $user_dir
    }
    
    set asset_dir [::PETK::gui::resourcePath bio_pore]
    if {[file isdirectory $asset_dir]} {
        if {[lsearch -exact $dirs $asset_dir] < 0} {
            lappend dirs $asset_dir
        }
    }
    return $dirs
}

proc ::PETK::gui::resolveBioPoreFile {pore_id} {
    if {$pore_id eq ""} {
        return ""
    }
    if {[file exists $pore_id]} {
        return $pore_id
    }
    if {[info exists ::PETK::gui::bioPorePathMap] && [dict exists $::PETK::gui::bioPorePathMap $pore_id]} {
        return [dict get $::PETK::gui::bioPorePathMap $pore_id]
    }
    foreach dir [::PETK::gui::getBioPoreSearchDirs] {
        set candidate [file join $dir $pore_id]
        if {[file exists $candidate]} {
            return $candidate
        }
    }
    return ""
}

proc ::PETK::gui::updateBioPoreInfo {{custom_message ""}} {
    set container $::PETK::gui::contentContainer
    set info_widget $container.biopore.info
    
    # Enable text widget for editing
    $info_widget configure -state normal
    $info_widget delete 1.0 end
    
    if {$custom_message ne ""} {
        $info_widget insert end $custom_message
    } else {
        # Determine which file to analyze
        set filepath ""
        if {$::PETK::gui::selectedBioPore ne ""} {
            set filepath [::PETK::gui::resolveBioPoreFile $::PETK::gui::selectedBioPore]
        } elseif {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
            set filepath $::PETK::gui::customBioPorePath
        }
        
        if {$filepath ne "" && [file exists $filepath]} {
            # Get basic file information
            set filesize [file size $filepath]
            set fileext [file extension $filepath]
            
            $info_widget insert end "File: [file tail $filepath]\n"
            $info_widget insert end "Type: [string toupper [string trimleft $fileext .]]\n"
            $info_widget insert end "Size: $filesize bytes\n"
            
            # Try to get more specific information based on file type
            if {[catch {
                set info [::PETK::gui::analyzeBioPoreFile $filepath]
                $info_widget insert end $info
            } err]} {
                $info_widget insert end "Additional analysis failed: $err\n"
            }
            
        } else {
            $info_widget insert end "No valid pore file selected.\n"
        }
    }
    
    # Disable text widget to make it read-only
    $info_widget configure -state disabled
}

# Function to analyze biological pore files (basic implementation)
proc ::PETK::gui::analyzeBioPoreFile {filepath} {
    set info ""
    set fileext [string tolower [file extension $filepath]]
    
    switch $fileext {
        ".pdb" {
            if {[catch {open $filepath r} fh]} {
                return "Error: Cannot open file\n"
            }
            
            set atom_count 0
            set has_hetatm 0
            set chain_ids {}
            
            while {[gets $fh line] >= 0} {
                if {[string match "ATOM*" $line]} {
                    incr atom_count
                    set chain_id [string range $line 21 21]
                    if {$chain_id ni $chain_ids} {
                        lappend chain_ids $chain_id
                    }
                } elseif {[string match "HETATM*" $line]} {
                    set has_hetatm 1
                }
            }
            close $fh
            
            append info "Atoms: $atom_count\n"
            append info "Chains: [join $chain_ids {, }]\n"
            if {$has_hetatm} {
                append info "Contains hetero atoms: Yes\n"
            }
        }
        
        ".gro" {
            if {[catch {open $filepath r} fh]} {
                return "Error: Cannot open file\n"
            }
            
            # Skip title line
            gets $fh
            # Get atom count
            gets $fh atom_line
            set atom_count [string trim $atom_line]
            close $fh
            
            append info "Atoms: $atom_count\n"
            append info "Format: GROMACS structure\n"
        }
        
        default {
            append info "File format: [string toupper [string trimleft $fileext .]]\n"
            append info "Basic file analysis only\n"
        }
    }
    
    return $info
}


####################################################
# Pore Type Management Functions
####################################################

proc ::PETK::gui::updateMembraneTypeDisplay {} {
    # FIX: Use consistent container reference
    if {![info exists ::PETK::gui::contentContainer]} {
        return
    }
    set container $::PETK::gui::contentContainer

    if {![info exists ::PETK::gui::previousMembraneType]} {
        set ::PETK::gui::previousMembraneType ""
    }
    
    if {$::PETK::gui::poreOption eq "solid-state"} {
        set new_type $::PETK::gui::membraneType
        if {$new_type eq "" || $new_type ni {cylindrical doublecone conical}} {
            set new_type "cylindrical"
            set ::PETK::gui::membraneType $new_type
        }
        set prior_type $::PETK::gui::previousMembraneType
        if {$prior_type ne "" && $prior_type ne $new_type} {
            ::PETK::gui::persistSolidStateParameters $prior_type
        }
        set use_defaults 0
        if {[info exists ::PETK::gui::forceSolidStateDefaults] && $::PETK::gui::forceSolidStateDefaults} {
            set use_defaults 1
            set ::PETK::gui::forceSolidStateDefaults 0
        }
        ::PETK::gui::applySolidStatePreset $new_type $use_defaults
        set ::PETK::gui::previousMembraneType $new_type
    }

    # Hide all parameter frames first
    catch {grid forget $container.params.cyl}
    catch {grid forget $container.params.cone}
    catch {grid forget $container.params.conical}

    if {$::PETK::gui::membraneType eq "cylindrical"} {
        grid $container.params.cyl -row 1 -column 0 -columnspan 6 -sticky ew -pady 5
        set ::PETK::gui::currentPoreType "Cylindrical"
        ::PETK::gui::updateMembraneStatus "Switched to cylindrical pore mode"
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        grid $container.params.cone -row 1 -column 0 -columnspan 6 -sticky ew -pady 5
        set ::PETK::gui::currentPoreType "Double Cone"
        ::PETK::gui::updateMembraneStatus "Switched to double cone pore mode"
    } elseif {$::PETK::gui::membraneType eq "conical"} {
        grid $container.params.conical -row 1 -column 0 -columnspan 6 -sticky ew -pady 5
        set ::PETK::gui::currentPoreType "Conical"
        ::PETK::gui::updateMembraneStatus "Switched to conical (frustum) pore mode"
    }
    
    # FIX: Force canvas and content refresh
    ::PETK::gui::requestCanvasLayoutRefresh
    ::PETK::gui::updateParameterDisplay
    ::PETK::gui::calculateBoxDimensions
}

proc ::PETK::gui::updateTab1ScrollRegion {} {
    if {![info exists ::PETK::gui::tab1MainCanvas] || ![winfo exists $::PETK::gui::tab1MainCanvas]} {
        return
    }
    if {![info exists ::PETK::gui::tab1ContentContainer] || ![winfo exists $::PETK::gui::tab1ContentContainer]} {
        return
    }

    set canvas $::PETK::gui::tab1MainCanvas
    set container $::PETK::gui::tab1ContentContainer

    update idletasks

    set content_width [winfo reqwidth $container]
    set content_height [winfo reqheight $container]
    set canvas_width [winfo width $canvas]

    if {$canvas_width > 1} {
        set content_width $canvas_width
    }

    if {[info exists ::PETK::gui::tab1CanvasWindow]} {
        $canvas itemconfig $::PETK::gui::tab1CanvasWindow -width $content_width
    }

    set scroll_height [expr {$content_height + 20}]
    $canvas configure -scrollregion [list 0 0 $content_width $scroll_height]

    puts "Updated Tab 1 scroll region: ${content_width}x${scroll_height}"
}

proc ::PETK::gui::refreshCanvasLayout {} {
    ::PETK::gui::updateTab1ScrollRegion
}

proc ::PETK::gui::requestCanvasLayoutRefresh {} {
    if {[info exists ::PETK::gui::pendingLayoutRefresh] && $::PETK::gui::pendingLayoutRefresh ne ""} {
        return
    }

    set ::PETK::gui::pendingLayoutRefresh [after idle {
        if {[info exists ::PETK::gui::pendingLayoutRefresh] && $::PETK::gui::pendingLayoutRefresh ne ""} {
            unset ::PETK::gui::pendingLayoutRefresh
        }
        catch {::PETK::gui::refreshCanvasLayout}
    }]
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
    } elseif {$::PETK::gui::membraneType eq "conical"} {
        set ::PETK::gui::displayParam1 "Top: $::PETK::gui::topDiameter Å"
        set ::PETK::gui::displayParam2 "Bottom: $::PETK::gui::bottomDiameter Å"

        # Volume of a single truncated cone (frustum):
        # V = π/3 * h * (r_top² + r_top*r_bot + r_bot²)
        if {[string is double $::PETK::gui::topDiameter] && [string is double $::PETK::gui::bottomDiameter] && [string is double $::PETK::gui::nanoporeThickness]} {
            set rt [expr {$::PETK::gui::topDiameter / 2.0}]
            set rb [expr {$::PETK::gui::bottomDiameter / 2.0}]
            set h $::PETK::gui::nanoporeThickness
            set volume [expr {3.14159/3.0 * $h * ($rt*$rt + $rt*$rb + $rb*$rb)}]
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
    # Load pore shape images from module assets
    ::PETK::gui::clearScaledImageCache
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

proc ::PETK::gui::updatePoreVisualization {} {
    if {[info exists ::PETK::gui::pendingPoreCanvasRedraw]} {
        unset ::PETK::gui::pendingPoreCanvasRedraw
    }

    if {[info exists ::PETK::gui::poreVisualizationSuspendCount] && $::PETK::gui::poreVisualizationSuspendCount > 0} {
        set ::PETK::gui::poreVisualizationPending 1
        return
    }

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
        set canvas_width [winfo width $canvas]
        set canvas_height [winfo height $canvas]
        if {$canvas_width <= 1} {set canvas_width [winfo reqwidth $canvas]}
        if {$canvas_height <= 1} {set canvas_height [winfo reqheight $canvas]}
        if {$canvas_width <= 1} {set canvas_width 300}
        if {$canvas_height <= 1} {set canvas_height 250}
        
        # Paint a white backdrop so transparent pixels in the pore image
        # don't reveal the default canvas background.
        $canvas create rectangle 0 0 $canvas_width $canvas_height \
            -fill white -outline ""
        
        set img_width [image width $image]
        set img_height [image height $image]
        
        # Calculate scaling to fit canvas while maintaining aspect ratio
        set scale_x [expr {double($canvas_width) / $img_width}]
        set scale_y [expr {double($canvas_height) / $img_height}]
        set scale [expr {min($scale_x, $scale_y) * 0.8}]  ; # 80% of canvas size
        
        # Create scaled image if needed
        if {$scale < 1.0} {
            set new_width [expr {max(1, int(round($img_width * $scale)))}]
            set new_height [expr {max(1, int(round($img_height * $scale)))}]
            set new_width [::PETK::gui::quantizePoreImageLength $new_width]
            set new_height [::PETK::gui::quantizePoreImageLength $new_height]
            set draw_image [::PETK::gui::getScaledPoreImage $image $new_width $new_height]
        } else {
            set draw_image $image
        }

        set draw_width [image width $draw_image]
        set draw_height [image height $draw_image]

        # Center whichever image we’re drawing
        set x [expr {($canvas_width - $draw_width) / 2}]
        set y [expr {($canvas_height - $draw_height) / 2}]
        $canvas create image $x $y -anchor nw -image $draw_image
        
    } else {
        # Draw a simple schematic if no image available
        ::PETK::gui::drawPoreSchematic $canvas
    }
    
    # Add parameter labels on the image
    # ::PETK::gui::addParameterLabels $canvas
    ::PETK::gui::updateParameterDisplay
    ::PETK::gui::updatePoreVmdVisualization 
    if {[info procs ::PETK::gui::syncMovementRangeToBox] ne ""} {
        catch {::PETK::gui::syncMovementRangeToBox}
    }
}

proc ::PETK::gui::updatePreviewAndMovement {} {
    set ::PETK::gui::forceMovementUpdate 1
    catch {::PETK::gui::calculateBoxDimensions}
    set ::PETK::gui::forceMovementUpdate 0
    ::PETK::gui::updatePoreVisualization
}

proc ::PETK::gui::pausePoreVisualization {} {
    if {![info exists ::PETK::gui::poreVisualizationSuspendCount]} {
        set ::PETK::gui::poreVisualizationSuspendCount 0
    }
    incr ::PETK::gui::poreVisualizationSuspendCount
}

proc ::PETK::gui::resumePoreVisualization {} {
    if {![info exists ::PETK::gui::poreVisualizationSuspendCount]} {
        return
    }

    incr ::PETK::gui::poreVisualizationSuspendCount -1
    if {$::PETK::gui::poreVisualizationSuspendCount < 0} {
        set ::PETK::gui::poreVisualizationSuspendCount 0
    }

    if {$::PETK::gui::poreVisualizationSuspendCount == 0 && [info exists ::PETK::gui::poreVisualizationPending] && $::PETK::gui::poreVisualizationPending} {
        unset ::PETK::gui::poreVisualizationPending
        after idle {catch {::PETK::gui::updatePoreVisualization}}
    }
}

proc ::PETK::gui::onPoreImageCanvasConfigure {canvas} {
    if {![winfo exists $canvas]} {
        return
    }

    set width [winfo width $canvas]
    set height [winfo height $canvas]

    if {$width <= 1 || $height <= 1} {
        return
    }

    set new_size [list $width $height]
    if {[info exists ::PETK::gui::lastPoreCanvasSize] && $::PETK::gui::lastPoreCanvasSize eq $new_size} {
        return
    }

    set ::PETK::gui::lastPoreCanvasSize $new_size

    if {[info exists ::PETK::gui::pendingPoreCanvasRedraw] && $::PETK::gui::pendingPoreCanvasRedraw ne ""} {
        catch {after cancel $::PETK::gui::pendingPoreCanvasRedraw}
    }
    set delay 90
    if {[info exists ::PETK::gui::poreImageRedrawDelay]} {
        set delay $::PETK::gui::poreImageRedrawDelay
    }
    set ::PETK::gui::pendingPoreCanvasRedraw \
        [after $delay {catch {::PETK::gui::updatePoreVisualization}}]
}

proc ::PETK::gui::quantizePoreImageLength {value} {
    set val [expr {int($value)}]
    if {$val <= 1} {
        return 1
    }

    set gran 1
    if {[info exists ::PETK::gui::poreImageScaleGranularity] && $::PETK::gui::poreImageScaleGranularity > 1} {
        set gran $::PETK::gui::poreImageScaleGranularity
    }

    if {$gran <= 1 || $val <= $gran} {
        return $val
    }

    set quantized [expr {int(floor(double($val) / double($gran))) * $gran}]
    if {$quantized < 1} {
        set quantized 1
    }
    return $quantized
}

proc ::PETK::gui::clearScaledImageCache {} {
    if {[info exists ::PETK::gui::scaledImageCache]} {
        foreach key [array names ::PETK::gui::scaledImageCache] {
            set cached $::PETK::gui::scaledImageCache($key)
            if {$cached ne ""} {
                catch {image delete $cached}
            }
            unset ::PETK::gui::scaledImageCache($key)
        }
        array unset ::PETK::gui::scaledImageCache
    }

    if {[info exists ::PETK::gui::photoColorDepth]} {
        array unset ::PETK::gui::photoColorDepth
    }
    if {[info exists ::PETK::gui::photoTransparencySupport]} {
        array unset ::PETK::gui::photoTransparencySupport
    }
}

proc ::PETK::gui::getScaledPoreImage {image target_width target_height} {
    set target_width [expr {max(1, int($target_width))}]
    set target_height [expr {max(1, int($target_height))}]

    set src_w [image width $image]
    set src_h [image height $image]

    if {$target_width >= $src_w && $target_height >= $src_h} {
        return $image
    }

    set cache_key [format "%s:%dx%d" $image $target_width $target_height]
    if {[info exists ::PETK::gui::scaledImageCache($cache_key)]} {
        set cached $::PETK::gui::scaledImageCache($cache_key)
        if {$cached ne "" && ![catch {image width $cached}]} {
            return $cached
        } else {
            unset ::PETK::gui::scaledImageCache($cache_key)
        }
    }

    set scaled_image [::PETK::gui::resamplePhoto $image $target_width $target_height]
    set ::PETK::gui::scaledImageCache($cache_key) $scaled_image
    return $scaled_image
}

proc ::PETK::gui::resamplePhoto {image target_width target_height} {
    set target_width [expr {max(1, int($target_width))}]
    set target_height [expr {max(1, int($target_height))}]

    set src_w [image width $image]
    set src_h [image height $image]

    if {$target_width >= $src_w && $target_height >= $src_h} {
        return $image
    }

    set scaled [image create photo -width $target_width -height $target_height]
    set max_channel [::PETK::gui::getPhotoColorDepth $image]
    set white_value $max_channel
    if {$white_value <= 0} {
        set white_value 255
    }
    set color_divisor [expr {$max_channel > 255 ? 257.0 : 1.0}]
    set transparency_lookup [::PETK::gui::photoSupportsTransparencyQuery $image]

    set x_scale [expr {double($src_w) / $target_width}]
    set y_scale [expr {double($src_h) / $target_height}]

    # Pre-compute source ranges for each destination pixel
    for {set dx 0} {$dx < $target_width} {incr dx} {
        set x_start($dx) [expr {int(floor($dx * $x_scale))}]
        set x_end($dx) [expr {int(ceil(($dx + 1) * $x_scale))}]
        if {$x_end($dx) > $src_w} {set x_end($dx) $src_w}
        if {$x_end($dx) <= $x_start($dx)} {set x_end($dx) [expr {$x_start($dx) + 1}]}
    }
    for {set dy 0} {$dy < $target_height} {incr dy} {
        set y_start($dy) [expr {int(floor($dy * $y_scale))}]
        set y_end($dy) [expr {int(ceil(($dy + 1) * $y_scale))}]
        if {$y_end($dy) > $src_h} {set y_end($dy) $src_h}
        if {$y_end($dy) <= $y_start($dy)} {set y_end($dy) [expr {$y_start($dy) + 1}]}
    }

    for {set dy 0} {$dy < $target_height} {incr dy} {
        set row_colors {}
        for {set dx 0} {$dx < $target_width} {incr dx} {
            set r 0
            set g 0
            set b 0
            set a 0
            set samples 0

            for {set sy $y_start($dy)} {$sy < $y_end($dy)} {incr sy} {
                for {set sx $x_start($dx)} {$sx < $x_end($dx)} {incr sx} {
                    set pixel [$image get $sx $sy]
                    lassign $pixel pr pg pb pa
                    set pixel_len [llength $pixel]
                    if {$pixel_len < 4} {
                        if {$transparency_lookup} {
                            set transparent [$image transparency get $sx $sy]
                            if {$transparent} {
                                set pa 0
                            } else {
                                set pa $white_value
                            }
                        } else {
                            set pa $white_value
                        }
                    }
                    incr r $pr
                    incr g $pg
                    incr b $pb
                    incr a $pa
                    incr samples
                }
            }

            if {$samples == 0} {
                lappend row_colors "#000000"
                continue
            }

            set avg_r [expr {$r / double($samples)}]
            set avg_g [expr {$g / double($samples)}]
            set avg_b [expr {$b / double($samples)}]
            set avg_a [expr {$a / double($samples)}]

            foreach component {avg_r avg_g avg_b avg_a} {
                if {[set $component] < 0} {set $component 0}
                if {[set $component] > 65535} {set $component 65535}
            }

            # Composite semi-transparent pixels over white so Tk color strings remain valid
            set alpha_frac [expr {$avg_a / double($white_value)}]
            if {$alpha_frac < 0} {set alpha_frac 0}
            if {$alpha_frac > 1} {set alpha_frac 1}

            set blended_r [expr {int(round($alpha_frac * $avg_r + (1.0 - $alpha_frac) * $white_value))}]
            set blended_g [expr {int(round($alpha_frac * $avg_g + (1.0 - $alpha_frac) * $white_value))}]
            set blended_b [expr {int(round($alpha_frac * $avg_b + (1.0 - $alpha_frac) * $white_value))}]

            foreach component {blended_r blended_g blended_b} {
                if {[set $component] < 0} {set $component 0}
                if {[set $component] > $white_value} {set $component $white_value}
            }

            # Convert to 8-bit hex color for photo put
            set color [format "#%02x%02x%02x" \
                [expr {int(round($blended_r / $color_divisor))}] \
                [expr {int(round($blended_g / $color_divisor))}] \
                [expr {int(round($blended_b / $color_divisor))}]]
            lappend row_colors $color
        }
        $scaled put [list $row_colors] -to 0 $dy
    }

    return $scaled
}

proc ::PETK::gui::getPhotoColorDepth {image} {
    if {[info exists ::PETK::gui::photoColorDepth($image)]} {
        return $::PETK::gui::photoColorDepth($image)
    }

    set width [image width $image]
    set height [image height $image]

    if {$width <= 0 || $height <= 0} {
        set ::PETK::gui::photoColorDepth($image) 65535
        return 65535
    }

    set mid_x [expr {int($width / 2)}]
    set mid_y [expr {int($height / 2)}]
    set coords [list \
        [list 0 0] \
        [list [expr {$width - 1}] 0] \
        [list 0 [expr {$height - 1}]] \
        [list [expr {$width - 1}] [expr {$height - 1}]] \
        [list $mid_x $mid_y]]

    set uses_high_depth 0
    foreach coord $coords {
        lassign $coord sx sy
        if {$sx < 0 || $sy < 0} {
            continue
        }
        set pixel [$image get $sx $sy]
        foreach value $pixel {
            if {$value > 255} {
                set uses_high_depth 1
                break
            }
        }
        if {$uses_high_depth} {
            break
        }
    }

    set depth [expr {$uses_high_depth ? 65535 : 255}]
    set ::PETK::gui::photoColorDepth($image) $depth
    return $depth
}

proc ::PETK::gui::photoSupportsTransparencyQuery {image} {
    if {[info exists ::PETK::gui::photoTransparencySupport($image)]} {
        return $::PETK::gui::photoTransparencySupport($image)
    }

    set width [image width $image]
    set height [image height $image]
    if {$width <= 0 || $height <= 0} {
        set ::PETK::gui::photoTransparencySupport($image) 0
        return 0
    }

    set supports 0
    if {[catch {$image transparency get 0 0} result] == 0} {
        set supports 1
    }

    set ::PETK::gui::photoTransparencySupport($image) $supports
    return $supports
}

proc ::PETK::gui::registerPreviewMolid {molid} {
    if {$molid eq ""} {
        return
    }

    if {![string is integer -strict $molid] || $molid < 0} {
        return
    }

    if {![info exists ::PETK::gui::petk_draw_molids]} {
        set ::PETK::gui::petk_draw_molids {}
    }

    if {[lsearch -exact $::PETK::gui::petk_draw_molids $molid] == -1} {
        lappend ::PETK::gui::petk_draw_molids $molid
    }
}

proc ::PETK::gui::clearPreviewGraphics {} {
    if {[info commands draw] ne ""} {
        catch {draw delete all}
    }

    if {![info exists ::PETK::gui::petk_draw_molids]} {
        return
    }

    if {[info commands graphics] eq ""} {
        return
    }

    set remaining {}
    foreach molid $::PETK::gui::petk_draw_molids {
        if {$molid eq ""} {
            continue
        }

        if {![string is integer -strict $molid] || $molid < 0} {
            continue
        }

        if {[catch {graphics $molid delete all}]} {
            # Molecule might no longer exist; skip keeping it in the list
            continue
        }

        lappend remaining $molid
    }

    set ::PETK::gui::petk_draw_molids $remaining
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

    } elseif {$::PETK::gui::membraneType eq "conical"} {
        # Draw conical (single frustum) pore schematic.
        # Top edge has half-width derived from topDiameter, bottom from bottomDiameter.
        set max_half 40
        set top_d 0.0
        set bot_d 0.0
        if {[string is double $::PETK::gui::topDiameter]} {
            set top_d [expr {$::PETK::gui::topDiameter}]
        }
        if {[string is double $::PETK::gui::bottomDiameter]} {
            set bot_d [expr {$::PETK::gui::bottomDiameter}]
        }
        set max_d [expr {max($top_d, $bot_d, 1.0)}]
        set top_half [expr {int($max_half * $top_d / $max_d)}]
        set bot_half [expr {int($max_half * $bot_d / $max_d)}]
        if {$top_half < 6} { set top_half 6 }
        if {$bot_half < 6} { set bot_half 6 }

        # Membrane (gray rectangles flanking the pore)
        set max_pore_half [expr {max($top_half, $bot_half)}]
        $canvas create rectangle 20 [expr {$cy-40}] [expr {$cx-$max_pore_half}] [expr {$cy+40}] -fill gray70 -outline black
        $canvas create rectangle [expr {$cx+$max_pore_half}] [expr {$cy-40}] [expr {$width-20}] [expr {$cy+40}] -fill gray70 -outline black

        # Frustum (trapezoid): top-edge at cy-40, bottom-edge at cy+40
        set coords [list [expr {$cx-$top_half}] [expr {$cy-40}] \
                         [expr {$cx+$top_half}] [expr {$cy-40}] \
                         [expr {$cx+$bot_half}] [expr {$cy+40}] \
                         [expr {$cx-$bot_half}] [expr {$cy+40}]]
        $canvas create polygon $coords -fill white -outline black

        # Labels
        $canvas create text $cx [expr {$cy-60}] -text "Conical Pore" -font {TkDefaultFont 10 bold}
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
    } elseif {$::PETK::gui::membraneType eq "conical"} {
        $canvas create text [expr {$width-10}] 20 -text "Top = $::PETK::gui::topDiameter Å" -anchor ne -font {TkDefaultFont 9} -fill blue
        $canvas create text [expr {$width-10}] 35 -text "Bottom = $::PETK::gui::bottomDiameter Å" -anchor ne -font {TkDefaultFont 9} -fill blue
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
    # Clear any existing graphics we previously drew
    ::PETK::gui::clearPreviewGraphics

    # Track the molecule that will receive new draw commands
    if {[info commands molinfo] ne ""} {
        if {![catch {molinfo top} current_molid]} {
            ::PETK::gui::registerPreviewMolid $current_molid
        }
    }

    if {[info exists ::PETK::gui::biological_pore_molid]} {
        catch {mol delete $::PETK::gui::biological_pore_molid}
        unset ::PETK::gui::biological_pore_molid
    }    
    
    # Check if we're dealing with biological or solid-state pore
    if {![info exists ::PETK::gui::poreOption]} {
        set ::PETK::gui::poreOption "solid-state"
    }
    
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
    draw line [list $min_x $min_y $min_z] [list $max_x $min_y $min_z] width 10
    draw line [list $max_x $min_y $min_z] [list $max_x $max_y $min_z] width 10
    draw line [list $max_x $max_y $min_z] [list $min_x $max_y $min_z] width 10
    draw line [list $min_x $max_y $min_z] [list $min_x $min_y $min_z] width 10
    
    # Top face
    draw line [list $min_x $min_y $max_z] [list $max_x $min_y $max_z] width 10
    draw line [list $max_x $min_y $max_z] [list $max_x $max_y $max_z] width 10
    draw line [list $max_x $max_y $max_z] [list $min_x $max_y $max_z] width 10
    draw line [list $min_x $max_y $max_z] [list $min_x $min_y $max_z] width 10
    
    # Vertical edges
    draw line [list $min_x $min_y $min_z] [list $min_x $min_y $max_z] width 10
    draw line [list $max_x $min_y $min_z] [list $max_x $min_y $max_z] width 10
    draw line [list $max_x $max_y $min_z] [list $max_x $max_y $max_z] width 10
    draw line [list $min_x $max_y $min_z] [list $min_x $max_y $max_z] width 10
    
    # Handle visualization based on pore type
    if {$::PETK::gui::poreOption eq "biological"} {
        ::PETK::gui::drawBiologicalPore $min_x $max_x $min_y $max_y $min_z $max_z
    } else {
        ::PETK::gui::drawSolidStatePore $min_x $max_x $min_y $max_y $min_z $max_z $xbox $ybox $zbox
    }
}

# Biological Pore Configuration System
# This creates default configurations for specific biological pores

# Define default configurations for biological pores
proc ::PETK::gui::initializeBiologicalPoreConfigurations {} {
    if {[info exists ::PETK::gui::bioPoreConfigsInitialized] && $::PETK::gui::bioPoreConfigsInitialized} {
        return
    }
    # Create global array to store pore configurations
    global ::PETK::gui::bioPoreConfigs
    
    # Configuration for centered_6MRT.pdb
    set ::PETK::gui::bioPoreConfigs(centered_6MRT.pdb) {
        boxSizeX 202.0
        boxSizeY 202.0  
        boxSizeZ 215.0
        membraneThickness 40.0
        membraneZOffset -53.0
        description "ClyA nanopore"
        autoCalculateBox 0
    }
    
    # Configuration for other pores (you can add more)
    set ::PETK::gui::bioPoreConfigs(centered_1UUN.pdb) {
        boxSizeX 140.0
        boxSizeY 140.0
        boxSizeZ 172.0
        membraneThickness 40.0
        membraneZOffset -34.0
        description "Mspa nanopore"
        autoCalculateBox 0
    }
    
    set ::PETK::gui::bioPoreConfigs(centered_7AHL.pdb) {
        boxSizeX 150.0
        boxSizeY 150.0
        boxSizeZ 198.0
        membraneThickness 40.0
        membraneZOffset -43.0
        description "alpha-hemolysin nanopore"
        autoCalculateBox 0
    }

    set ::PETK::gui::bioPoreConfigs(centered_3X2R.pdb) {
        boxSizeX 130.0
        boxSizeY 128.0
        boxSizeZ 120.0
        membraneThickness 40.0
        membraneZOffset -32.0
        description "CsgG nanopore"
        autoCalculateBox 0
    }

    # Add more pore configurations as needed
    # set ::PETK::gui::bioPoreConfigs(pore_name.pdb) { ... }
    
    puts "Initialized biological pore configurations for [array size ::PETK::gui::bioPoreConfigs] pore types"
    set ::PETK::gui::bioPoreConfigsInitialized 1
}

# Get configuration for a specific biological pore
proc ::PETK::gui::getBiologicalPoreConfig {pore_filename} {
    global ::PETK::gui::bioPoreConfigs
    
    # Extract just the filename from full path
    set filename [file tail $pore_filename]
    
    if {[info exists ::PETK::gui::bioPoreConfigs($filename)]} {
        return $::PETK::gui::bioPoreConfigs($filename)
    } else {
        # Return default configuration for unknown pores
        return {
            boxSizeX 150.0
            boxSizeY 150.0
            boxSizeZ 150.0
            membraneThickness 20.0
            membraneZOffset 0.0
            description "Unknown pore - using default settings"
            autoCalculateBox 0
        }
    }
}

# Apply configuration when a biological pore is loaded
proc ::PETK::gui::applyBiologicalPoreConfiguration {pore_filename} {
    set config [::PETK::gui::getBiologicalPoreConfig $pore_filename]
    set filename [file tail $pore_filename]
    
    puts "=== APPLYING BIOLOGICAL PORE CONFIGURATION ==="
    puts "Pore file: $filename"
    
    
    # Extract configuration values
    array set cfg $config
    
    puts "Configuration:"
    puts "  Box size: $cfg(boxSizeX) x $cfg(boxSizeY) x $cfg(boxSizeZ) Å"
    puts "  Membrane thickness: $cfg(membraneThickness) Å"
    puts "  Membrane Z offset: $cfg(membraneZOffset) Å"
    puts "  Description: $cfg(description)"
    
    # Apply box dimensions
    set ::PETK::gui::boxSizeX $cfg(boxSizeX)
    set ::PETK::gui::boxSizeY $cfg(boxSizeY)
    set ::PETK::gui::boxSizeZ $cfg(boxSizeZ)
    
    # Apply membrane parameters
    set ::PETK::gui::nanoporeThickness $cfg(membraneThickness)
    set ::PETK::gui::membraneZOffset $cfg(membraneZOffset)

    # Set auto-calculate box to false (use fixed dimensions)
    set ::PETK::gui::autoCalculateBoxDimensions $cfg(autoCalculateBox)
    
    # Force the GUI to show manual mode (not auto-calculate mode)
    ::PETK::gui::ensureManualBoxMode
    
    # Update GUI elements to reflect the new values
    ::PETK::gui::updateBiologicalPoreGUIElements $config
    
    # Update pore information display
    ::PETK::gui::updateBiologicalPoreInfo $config
    
    puts "Configuration applied successfully"
    puts "============================================="
}

# Ensure manual box dimension mode is displayed
proc ::PETK::gui::ensureManualBoxMode {} {
    # Get the container path from the stored global variable
    if {[info exists ::PETK::gui::contentContainer]} {
        set container $::PETK::gui::contentContainer
    } else {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    }
    
    # Make sure we're in manual mode for biological pores
    set ::PETK::gui::autoCalculateBoxDimensions 0
    
    # Call the existing toggle function to update the display
    if {[info procs ::PETK::gui::toggleBoxDimensionMode] ne ""} {
        ::PETK::gui::toggleBoxDimensionMode
    } else {
        # Fallback: manually grid the correct frame
        if {[winfo exists $container.boxdim.manual] && [winfo exists $container.boxdim.auto_display]} {
            grid remove $container.boxdim.auto_display
            grid $container.boxdim.manual -row 1 -column 0 -columnspan 6 -sticky ew
        }
    }
}

# Update GUI elements when biological pore configuration is applied
proc ::PETK::gui::updateBiologicalPoreGUIElements {config} {
    array set cfg $config
    
    # Get the container path from the stored global variable
    if {[info exists ::PETK::gui::contentContainer]} {
        set container $::PETK::gui::contentContainer
    } else {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    }
    
    # Disable auto-calculate checkbox for biological pores (use fixed dimensions)
    if {[winfo exists $container.boxdim.auto]} {
        $container.boxdim.auto configure -state disabled
        set ::PETK::gui::autoCalculateBoxDimensions 0
    }
    
    # Disable box dimension manual entry controls (make them read-only for biological pores)
    if {[winfo exists $container.boxdim.manual.x]} {
        $container.boxdim.manual.x configure -state readonly
    }
    if {[winfo exists $container.boxdim.manual.y]} {
        $container.boxdim.manual.y configure -state readonly  
    }
    if {[winfo exists $container.boxdim.manual.z]} {
        $container.boxdim.manual.z configure -state readonly
    }
    
    # Disable distance padding control
    if {[winfo exists $container.boxdim.cutoff]} {
        $container.boxdim.cutoff configure -state readonly
    }
    
    # Disable recalculate button
    if {[winfo exists $container.boxdim.auto_display.calc]} {
        $container.boxdim.auto_display.calc configure -state disabled
    }
    
    # Ensure manual mode is displayed (not auto-calculate mode)
    ::PETK::gui::ensureManualBoxMode
    
    # Disable membrane thickness control (read-only for biological pores)
    if {[winfo exists $container.biomembrane.thick]} {
        $container.biomembrane.thick configure -state readonly
    }
    if {[winfo exists $container.biomembrane.thickinfo]} {
        $container.biomembrane.thickinfo configure \
            -text "Fixed - optimized for this biological pore" \
            -foreground gray
    }
    
    # Keep Z-offset control enabled (user can modify this)
    if {[winfo exists $container.biomembrane.zoffset]} {
        $container.biomembrane.zoffset configure -state readonly
    }
    if {[winfo exists $container.biomembrane.zoffsetinfo]} {
        $container.biomembrane.zoffsetinfo configure \
            -text "Fixed - optimized for this biological pore" \
            -foreground gray
    }
    
    puts "GUI elements updated for biological pore:"
    puts "  - Box dimensions are now FIXED (read-only)"
    puts "  - Membrane thickness is FIXED (read-only)" 
    puts "  - Auto-calculate is DISABLED"
    puts "  - Z-offset is (read-only)"
}

# Restore GUI controls when switching away from biological pores
proc ::PETK::gui::restoreSolidStateGUIControls {} {
    # Get the container path from the stored global variable
    if {[info exists ::PETK::gui::contentContainer]} {
        set container $::PETK::gui::contentContainer
    } else {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    }
    
    # Re-enable auto-calculate checkbox
    if {[winfo exists $container.boxdim.auto]} {
        $container.boxdim.auto configure -state normal
    }
    
    # Re-enable box dimension manual entry controls
    if {[winfo exists $container.boxdim.manual.x]} {
        $container.boxdim.manual.x configure -state normal
    }
    if {[winfo exists $container.boxdim.manual.y]} {
        $container.boxdim.manual.y configure -state normal
    }
    if {[winfo exists $container.boxdim.manual.z]} {
        $container.boxdim.manual.z configure -state normal
    }
    
    # Re-enable distance padding control
    if {[winfo exists $container.boxdim.cutoff]} {
        $container.boxdim.cutoff configure -state normal
    }
    
    # Re-enable recalculate button
    if {[winfo exists $container.boxdim.auto_display.calc]} {
        $container.boxdim.auto_display.calc configure -state normal
    }
    
    # Re-enable membrane thickness control
    if {[winfo exists $container.biomembrane.thick]} {
        $container.biomembrane.thick configure -state normal
    }
    
    puts "Restored solid-state GUI controls:"
    puts "  - Box dimensions are now USER-ADJUSTABLE"
    puts "  - Membrane thickness is USER-ADJUSTABLE"
    puts "  - Auto-calculate is ENABLED"
    puts "  - All parameters are user-modifiable"
}

proc ::PETK::gui::persistSolidStateParameters {{targetType ""}} {
    if {$targetType eq ""} {
        if {![info exists ::PETK::gui::membraneType]} {
            return
        }
        set targetType $::PETK::gui::membraneType
    }
    
    if {$targetType ni {cylindrical doublecone conical}} {
        return
    }
    
    if {![info exists ::PETK::gui::autoCalculateBoxDimensions]} {
        set ::PETK::gui::autoCalculateBoxDimensions 0
    }
    
    set state [dict create \
        boxSizeX $::PETK::gui::boxSizeX \
        boxSizeY $::PETK::gui::boxSizeY \
        boxSizeZ $::PETK::gui::boxSizeZ \
        nanoporeThickness $::PETK::gui::nanoporeThickness \
        autoCalculateBox $::PETK::gui::autoCalculateBoxDimensions]
    
    if {$targetType eq "cylindrical"} {
        dict set state cylindricalDiameter $::PETK::gui::cylindricalDiameter
        dict set state cornerRadius $::PETK::gui::cornerRadius
    } elseif {$targetType eq "doublecone"} {
        dict set state innerDiameter $::PETK::gui::innerDiameter
        dict set state outerDiameter $::PETK::gui::outerDiameter
    } elseif {$targetType eq "conical"} {
        dict set state topDiameter $::PETK::gui::topDiameter
        dict set state bottomDiameter $::PETK::gui::bottomDiameter
    }
    
    set ::PETK::gui::solidStateState($targetType) $state
    puts "Persisted solid-state parameters for $targetType mode"
}

proc ::PETK::gui::applySolidStatePreset {{profile_type ""} {useDefaults 0}} {
    # Resolve target pore type
    if {$profile_type eq ""} {
        if {[info exists ::PETK::gui::membraneType] && $::PETK::gui::membraneType ne ""} {
            set profile_type $::PETK::gui::membraneType
        } else {
            set profile_type "cylindrical"
        }
    }
    
    if {$profile_type ni {cylindrical doublecone conical}} {
        puts "Warning: Unknown solid-state pore type '$profile_type' - falling back to cylindrical defaults"
        set profile_type "cylindrical"
    }
    
    if {![info exists ::PETK::gui::solidStateDefaults($profile_type)]} {
        puts "Warning: No solid-state defaults defined for '$profile_type'"
        return
    }
    
    set defaults $::PETK::gui::solidStateDefaults($profile_type)
    set source ""
    
    if {!$useDefaults && [info exists ::PETK::gui::solidStateState($profile_type)]} {
        set stored $::PETK::gui::solidStateState($profile_type)
        if {$stored ne ""} {
            set source $stored
        }
    }
    
    if {$source eq ""} {
        set source $defaults
        set useDefaults 1
    }
    
    # Apply auto/manual configuration
    if {[dict exists $source autoCalculateBox]} {
        set auto_val [dict get $source autoCalculateBox]
    } else {
        set auto_val [dict get $defaults autoCalculateBox]
    }
    set ::PETK::gui::autoCalculateBoxDimensions [expr {$auto_val ? 1 : 0}]
    
    # Common dimension parameters
    foreach key {boxSizeX boxSizeY boxSizeZ nanoporeThickness} {
        if {[dict exists $source $key]} {
            set value [dict get $source $key]
        } else {
            set value [dict get $defaults $key]
        }
        set ::PETK::gui::$key $value
    }
    
    # Type-specific parameters
    if {$profile_type eq "cylindrical"} {
        if {[dict exists $source cylindricalDiameter]} {
            set ::PETK::gui::cylindricalDiameter [dict get $source cylindricalDiameter]
        } else {
            set ::PETK::gui::cylindricalDiameter [dict get $defaults cylindricalDiameter]
        }
        if {[dict exists $source cornerRadius]} {
            set ::PETK::gui::cornerRadius [dict get $source cornerRadius]
        } else {
            set ::PETK::gui::cornerRadius [dict get $defaults cornerRadius]
        }
    } elseif {$profile_type eq "conical"} {
        if {[dict exists $source topDiameter]} {
            set ::PETK::gui::topDiameter [dict get $source topDiameter]
        } else {
            set ::PETK::gui::topDiameter [dict get $defaults topDiameter]
        }
        if {[dict exists $source bottomDiameter]} {
            set ::PETK::gui::bottomDiameter [dict get $source bottomDiameter]
        } else {
            set ::PETK::gui::bottomDiameter [dict get $defaults bottomDiameter]
        }
    } else {
        if {[dict exists $source innerDiameter]} {
            set ::PETK::gui::innerDiameter [dict get $source innerDiameter]
        } else {
            set ::PETK::gui::innerDiameter [dict get $defaults innerDiameter]
        }
        if {[dict exists $source outerDiameter]} {
            set ::PETK::gui::outerDiameter [dict get $source outerDiameter]
        } else {
            set ::PETK::gui::outerDiameter [dict get $defaults outerDiameter]
        }
    }
    
    # Update the GUI to reflect manual vs auto mode
    if {[info procs ::PETK::gui::toggleBoxDimensionMode] ne "" && [info exists ::PETK::gui::window]} {
        ::PETK::gui::toggleBoxDimensionMode
    }
    
    # Refresh displayed parameters to match the applied values
    ::PETK::gui::updateParameterDisplay
    
    # Update the VMD preview to reflect the reset values
    if {[info procs ::PETK::gui::updatePoreVisualization] ne ""} {
        ::PETK::gui::updatePoreVisualization
    }
    
    if {$useDefaults} {
        puts "Applied solid-state defaults for $profile_type pore"
    } else {
        puts "Restored stored solid-state parameters for $profile_type pore"
    }
    
    set ::PETK::gui::previousMembraneType $profile_type
}

# Update biological pore information display
proc ::PETK::gui::updateBiologicalPoreInfo {config} {
    array set cfg $config
    
    # Get the container path from the stored global variable
    if {[info exists ::PETK::gui::contentContainer]} {
        set container $::PETK::gui::contentContainer
    } else {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    }
    
    # Create info text for display
    set info_text "=== BIOLOGICAL PORE CONFIGURATION ===\n"
    append info_text "Pore: [file tail $::PETK::gui::selectedBioPore]\n"
    append info_text "Description: $cfg(description)\n\n"
    append info_text "FIXED PARAMETERS (optimized for this pore):\n"
    append info_text "• Box dimensions: $cfg(boxSizeX) × $cfg(boxSizeY) × $cfg(boxSizeZ) Å\n"
    append info_text "• Membrane thickness: $cfg(membraneThickness) Å\n"
    append info_text "• Auto-calculate box: DISABLED\n\n"
    append info_text "NOTES:\n"
    append info_text "• Box size and membrane thickness are read-only\n"
    append info_text "• These values are optimized for this specific pore\n"
    append info_text "• Switch to solid-state mode to regain full control\n"
    
    # Update info display in GUI
    if {[winfo exists $container.biopore.info]} {
        $container.biopore.info configure -state normal
        $container.biopore.info delete 1.0 end
        $container.biopore.info insert 1.0 $info_text
        $container.biopore.info configure -state disabled
    }
}

# Update information display for custom biological pores
proc ::PETK::gui::updateCustomBiologicalPoreInfo {filename} {
    # Get the container path from the stored global variable
    set container ""
    if {[info exists ::PETK::gui::contentContainer] && $::PETK::gui::contentContainer ne ""} {
        set container $::PETK::gui::contentContainer
    } elseif {[winfo exists ".petk_main_window.hlf.nb.tab1.canvas.content"]} {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    } else {
        puts "Warning: Unable to locate Tab 1 container for custom pore info update"
        return
    }
    
    set short_name [file tail $filename]
    
    # Create info text for display
    set info_text "=== CUSTOM BIOLOGICAL PORE ===\n"
    append info_text "File: $short_name\n"
    append info_text "Path: $filename\n\n"
    append info_text "ALL PARAMETERS ARE USER-ADJUSTABLE:\n"
    append info_text "• Box size: $::PETK::gui::boxSizeX × $::PETK::gui::boxSizeY × $::PETK::gui::boxSizeZ Å\n"
    append info_text "• Membrane thickness: $::PETK::gui::nanoporeThickness Å\n"
    append info_text "• Membrane Z-offset: $::PETK::gui::membraneZOffset Å\n\n"
    append info_text "NOTES:\n"
    append info_text "• Custom pores use default settings initially\n"
    append info_text "• Adjust parameters as needed for your system\n"
    append info_text "• Use analysis tool to optimize settings\n"
    
    # Update info display in GUI
    if {[winfo exists $container.biopore.info]} {
        $container.biopore.info configure -state normal
        $container.biopore.info delete 1.0 end
        $container.biopore.info insert 1.0 $info_text
        $container.biopore.info configure -state disabled
    }
}

# Enhanced pore selection handler
proc ::PETK::gui::onBiologicalPoreSelected {} {
    # Check if a pore is selected
    if {![info exists ::PETK::gui::selectedBioPore] || $::PETK::gui::selectedBioPore eq ""} {
        return
    }
    
    set pore_filename [::PETK::gui::resolveBioPoreFile $::PETK::gui::selectedBioPore]
    if {$pore_filename eq ""} {
        puts "Warning: Unable to locate biological pore file for '$::PETK::gui::selectedBioPore'"
        return
    }
    puts "Biological pore selected: [file tail $pore_filename]"
    
    set ::PETK::gui::customBioBoundingInfo ""
    
    # Apply configuration for this pore
    ::PETK::gui::applyBiologicalPoreConfiguration $pore_filename
    
    # Clear any custom pore path when pre-prepared pore is selected
    if {[info exists ::PETK::gui::customBioPorePath]} {
        set ::PETK::gui::customBioPorePath ""
    }
    
    # Update visualization
    ::PETK::gui::updatePoreVisualization
}

# Custom biological pore selection handler
proc ::PETK::gui::selectCustomBioPore {} {
    set filename [tk_getOpenFile \
        -title "Select Biological Pore Structure File" \
        -filetypes {
            {"PDB files" {.pdb}}
            {"GRO files" {.gro}}
            {"XYZ files" {.xyz}}
            {"MOL2 files" {.mol2}}
            {"SDF files" {.sdf}}
            {"All files" {*}}
        }]
    
    if {$filename ne ""} {
        set ::PETK::gui::customBioPorePath $filename
        
        # Clear pre-prepared selection
        set ::PETK::gui::selectedBioPore ""
        
        # Apply default configuration for custom pore
        set applied_ok [::PETK::gui::applyCustomBiologicalPoreConfiguration $filename]
        
        puts "Custom biological pore selected: [file tail $filename]"
        
        if {!$applied_ok} {
            puts "Warning: Failed to fully configure custom biological pore - using fallback parameters"
        }
        
        ::PETK::gui::updateParameterDisplay
        # Update visualization
        ::PETK::gui::updatePoreVisualization
    }
}

# Apply default configuration for custom biological pores
proc ::PETK::gui::applyCustomBiologicalPoreConfiguration {pore_filename} {
    set filename [file tail $pore_filename]
    
    puts "=== APPLYING CUSTOM BIOLOGICAL PORE CONFIGURATION ==="
    puts "Custom pore file: $filename"
    
    set ::PETK::gui::customBioBoundingInfo ""
    
    # Default to auto-calculated box display for custom pores
    set ::PETK::gui::autoCalculateBoxDimensions 1
    set ::PETK::gui::nanoporeThickness "40.0"
    set ::PETK::gui::membraneZOffset "0.0"
    set thickness_val 40.0
    if {[string is double -strict $::PETK::gui::nanoporeThickness]} {
        set thickness_val [expr {double($::PETK::gui::nanoporeThickness)}]
    }
    
    # Determine padding to apply around the biological structure
    set padding 20.0
    if {[info exists ::PETK::gui::sysPadding] && [string is double -strict $::PETK::gui::sysPadding]} {
        set padding [expr {double($::PETK::gui::sysPadding)}]
    }
    set padding [expr {min(25.0, max(10.0, $padding))}]
    
    set load_success [::PETK::gui::loadBiologicalPoreStructure $pore_filename]
    if {!$load_success} {
        puts "Warning: Failed to load custom biological pore structure. Using fallback dimensions."
        set fallback_box 150.0
        set ::PETK::gui::boxSizeX [format "%.1f" $fallback_box]
        set ::PETK::gui::boxSizeY [format "%.1f" $fallback_box]
        set ::PETK::gui::boxSizeZ [format "%.1f" $fallback_box]
        set fallback_half [expr {$fallback_box / 2.0}]
        ::PETK::gui::updateAutoBoxDisplayFromDimensions $fallback_box $fallback_box $fallback_box [expr {-$fallback_half}] [expr {$fallback_half}]
        puts [format "Auto movement presets (fallback): %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]

        ::PETK::gui::enableCustomPoreGUIControls
        if {[info procs ::PETK::gui::toggleBoxDimensionMode] ne "" && [info exists ::PETK::gui::window]} {
            ::PETK::gui::toggleBoxDimensionMode
        }
        ::PETK::gui::updateCustomBiologicalPoreInfo $pore_filename
        ::PETK::gui::updateParameterDisplay
        puts "============================================="
        return 0
    }
    
    set molid $::PETK::gui::biological_pore_molid
    set bbox [::PETK::gui::calculateStructureBoundingBox $molid]
    
    if {$bbox eq ""} {
        puts "Warning: Unable to determine bounding box for custom pore. Using fallback dimensions."
        set fallback_box 150.0
        set ::PETK::gui::boxSizeX [format "%.1f" $fallback_box]
        set ::PETK::gui::boxSizeY [format "%.1f" $fallback_box]
        set ::PETK::gui::boxSizeZ [format "%.1f" $fallback_box]
        set fallback_half [expr {$fallback_box / 2.0}]
        ::PETK::gui::updateAutoBoxDisplayFromDimensions $fallback_box $fallback_box $fallback_box [expr {-$fallback_half}] [expr {$fallback_half}]
        puts [format "Auto movement presets (fallback): %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
    } else {
        set ::PETK::gui::customBioBoundingInfo $bbox
        
        set x_range [dict get $bbox x_range]
        set y_range [dict get $bbox y_range]
        set z_range [dict get $bbox z_range]
        
        set box_x [expr {$x_range + 2.0 * $padding}]
        set box_y [expr {$y_range + 2.0 * $padding}]
        set box_z [expr {$z_range + 2.0 * $padding}]
        
        set min_z [expr {$thickness_val + 20.0}]
        set computed_box_x [expr {max(40.0, $box_x)}]
        set computed_box_y [expr {max(40.0, $box_y)}]
        set computed_box_z [expr {max($min_z, $box_z)}]
        
        set ::PETK::gui::boxSizeX [format "%.1f" $computed_box_x]
        set ::PETK::gui::boxSizeY [format "%.1f" $computed_box_y]
        set ::PETK::gui::boxSizeZ [format "%.1f" $computed_box_z]
        
        set computed_half_z [expr {$computed_box_z / 2.0}]
        ::PETK::gui::updateAutoBoxDisplayFromDimensions $computed_box_x $computed_box_y $computed_box_z [expr {-$computed_half_z}] [expr {$computed_half_z}]
        puts [format "Auto movement presets: %s to %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
        
        puts "Calculated custom pore extents:"
        puts [format "  X: %.2f to %.2f (%.2f Å range)" [dict get $bbox x_min] [dict get $bbox x_max] $x_range]
        puts [format "  Y: %.2f to %.2f (%.2f Å range)" [dict get $bbox y_min] [dict get $bbox y_max] $y_range]
        puts [format "  Z: %.2f to %.2f (%.2f Å range)" [dict get $bbox z_min] [dict get $bbox z_max] $z_range]
        puts [format "  Atoms considered: %d" [dict get $bbox num_atoms]]
        puts [format "  Applied padding: %.1f Å" $padding]
    }
    
    puts "Applied simulation parameters:"
    puts "  Box size: $::PETK::gui::boxSizeX × $::PETK::gui::boxSizeY × $::PETK::gui::boxSizeZ Å"
    puts "  Membrane thickness: $::PETK::gui::nanoporeThickness Å"
    puts "  Membrane Z-offset: $::PETK::gui::membraneZOffset Å"
    
    ::PETK::gui::enableCustomPoreGUIControls
    if {[info procs ::PETK::gui::toggleBoxDimensionMode] ne "" && [info exists ::PETK::gui::window]} {
        ::PETK::gui::toggleBoxDimensionMode
    }
    ::PETK::gui::updateCustomBiologicalPoreInfo $pore_filename
    ::PETK::gui::updateParameterDisplay
    
    puts "Custom pore configuration applied"
    puts "============================================="
    return 1
}

# Enable all controls for custom biological pores
proc ::PETK::gui::enableCustomPoreGUIControls {} {
    # Get the container path from the stored global variable
    if {[info exists ::PETK::gui::contentContainer]} {
        set container $::PETK::gui::contentContainer
    } else {
        # Fallback to the known path structure
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    }
    
    # Enable auto-calculate checkbox for custom pores
    if {[winfo exists $container.boxdim.auto]} {
        $container.boxdim.auto configure -state normal
    }
    
    # Enable box dimension manual entry controls for custom pores
    if {[winfo exists $container.boxdim.manual.x]} {
        $container.boxdim.manual.x configure -state normal
    }
    if {[winfo exists $container.boxdim.manual.y]} {
        $container.boxdim.manual.y configure -state normal  
    }
    if {[winfo exists $container.boxdim.manual.z]} {
        $container.boxdim.manual.z configure -state normal
    }
    
    # Enable distance padding control
    if {[winfo exists $container.boxdim.cutoff]} {
        $container.boxdim.cutoff configure -state normal
    }
    
    # Enable recalculate button
    if {[winfo exists $container.boxdim.auto_display.calc]} {
        $container.boxdim.auto_display.calc configure -state normal
    }
    
    # Enable membrane thickness control for custom pores
    if {[winfo exists $container.biomembrane.thick]} {
        $container.biomembrane.thick configure -state normal
    }
    if {[winfo exists $container.biomembrane.thickinfo]} {
        $container.biomembrane.thickinfo configure \
            -text "Adjustable for custom biological pore" \
            -foreground green4
    }
    
    # Z-offset remains enabled
    if {[winfo exists $container.biomembrane.zoffset]} {
        $container.biomembrane.zoffset configure -state normal
    }
    if {[winfo exists $container.biomembrane.zoffsetinfo]} {
        $container.biomembrane.zoffsetinfo configure \
            -text "Adjustable for custom biological pore" \
            -foreground green4
    }
    
    puts "Custom pore GUI controls enabled:"
    puts "  - Box dimensions are USER-ADJUSTABLE"
    puts "  - Membrane thickness is USER-ADJUSTABLE"
    puts "  - Auto-calculate is AVAILABLE"
    puts "  - All parameters are user-modifiable"
}

# Update information display for custom biological pores
proc ::PETK::gui::updateCustomBiologicalPoreInfo {filename} {
    # Resolve container widget (fallback to known path if global ref missing)
    set container ""
    if {[info exists ::PETK::gui::contentContainer] && $::PETK::gui::contentContainer ne ""} {
        set container $::PETK::gui::contentContainer
    } elseif {[winfo exists ".petk_main_window.hlf.nb.tab1.canvas.content"]} {
        set container ".petk_main_window.hlf.nb.tab1.canvas.content"
    } else {
        puts "Warning: Unable to locate Tab 1 container for custom pore info update"
        return
    }
    
    set short_name [file tail $filename]
    
    # Create info text for display
    set info_text "=== CUSTOM BIOLOGICAL PORE ===\n"
    append info_text "File: $short_name\n"
    append info_text "Path: $filename\n\n"
    
    append info_text "Computed Structure Extents (Å):\n"
    if {[info exists ::PETK::gui::customBioBoundingInfo] && $::PETK::gui::customBioBoundingInfo ne ""} {
        set bbox $::PETK::gui::customBioBoundingInfo
        append info_text [format "• X: %.2f to %.2f (%.2f)\n" \
            [dict get $bbox x_min] [dict get $bbox x_max] [dict get $bbox x_range]]
        append info_text [format "• Y: %.2f to %.2f (%.2f)\n" \
            [dict get $bbox y_min] [dict get $bbox y_max] [dict get $bbox y_range]]
        append info_text [format "• Z: %.2f to %.2f (%.2f)\n" \
            [dict get $bbox z_min] [dict get $bbox z_max] [dict get $bbox z_range]]
        if {[dict exists $bbox num_atoms]} {
            append info_text [format "• Atoms considered: %d\n" [dict get $bbox num_atoms]]
        }
        if {[dict exists $bbox center]} {
            set center [dict get $bbox center]
            append info_text [format "• Approx. center: (%.2f, %.2f, %.2f)\n" \
                [lindex $center 0] [lindex $center 1] [lindex $center 2]]
        }
        append info_text "\n"
    } else {
        append info_text "• Extents unavailable (using fallback dimensions)\n\n"
    }
    
    append info_text "Applied Simulation Box (Å):\n"
    append info_text "• Box size: $::PETK::gui::boxSizeX × $::PETK::gui::boxSizeY × $::PETK::gui::boxSizeZ Å\n"
    append info_text "• Membrane thickness: $::PETK::gui::nanoporeThickness Å\n"
    append info_text "• Membrane Z-offset: $::PETK::gui::membraneZOffset Å\n\n"
    
    append info_text "Notes:\n"
    append info_text "• These values were derived from the loaded structure with padding.\n"
    append info_text "• Adjust parameters as needed for your system.\n"
    append info_text "• Auto-calculate is disabled by default for custom pores.\n"
    
    # Update info display in GUI
    if {[winfo exists $container.biopore.info]} {
        $container.biopore.info configure -state normal
        $container.biopore.info delete 1.0 end
        $container.biopore.info insert 1.0 $info_text
        $container.biopore.info configure -state disabled
    }
}

# Procedure to load biological pore structure in VMD
proc ::PETK::gui::loadBiologicalPoreStructure {pore_file} {
    # Remove any existing biological pore molecules
    if {[info exists ::PETK::gui::biological_pore_molid]} {
        catch {mol delete $::PETK::gui::biological_pore_molid}
        unset ::PETK::gui::biological_pore_molid
    }
    
    # Determine file type and appropriate loader
    set file_ext [string tolower [file extension $pore_file]]
    
    # Try to load the structure
    set molid -1
    if {[catch {
        switch $file_ext {
            ".pdb" {
                set molid [mol new $pore_file type pdb]
            }
            ".gro" {
                set molid [mol new $pore_file type gro]
            }
            ".xyz" {
                set molid [mol new $pore_file type xyz]
            }
            ".mol2" {
                set molid [mol new $pore_file type mol2]
            }
            ".sdf" {
                set molid [mol new $pore_file type sdf]
            }
            default {
                # Try to load as PDB by default
                set molid [mol new $pore_file type pdb]
            }
        }
    } load_error]} {
        puts "Error loading file: $load_error"
        return 0
    }
    
    # Check if molecule was loaded successfully
    if {$molid == -1} {
        puts "Failed to create molecule - mol new returned -1"
        return 0
    }
    
    # Verify the molecule has atoms
    set num_atoms [molinfo $molid get numatoms]
    if {$num_atoms == 0} {
        puts "Loaded molecule has no atoms"
        catch {mol delete $molid}
        return 0
    }
    
    # Store the molecule ID for later reference
    set ::PETK::gui::biological_pore_molid $molid
    ::PETK::gui::registerPreviewMolid $molid
    puts "Loaded biological pore structure with $num_atoms atoms"
    
    # Set up visualization style for biological pore with better error handling
    # Check if there are any representations first
    set num_reps [molinfo $molid get numreps]
    if {$num_reps > 0} {
        # Remove default representation - ignore errors
        catch {mol delrep 0 $molid}
    }
    
    # Add protein representation
    catch {
        mol representation NewCartoon
        mol color ColorID 3
        mol selection "protein"
        mol material Opaque
        mol addrep $molid
        puts "Added protein representation"
    }
    
    # Add nucleic acid representation if present
    catch {
        set nucleic_sel [atomselect $molid "nucleic"]
        set nucleic_count [$nucleic_sel num]
        $nucleic_sel delete
        
        if {$nucleic_count > 0} {
            mol representation NewCartoon
            mol color ColorID 7
            mol selection "nucleic"
            mol material Opaque
            mol addrep $molid
            puts "Added nucleic acid representation ($nucleic_count atoms)"
        }
    }
    
    # Add hetero atom representation if present
    catch {
        set hetero_sel [atomselect $molid "hetero and not water"]
        set hetero_count [$hetero_sel num]
        $hetero_sel delete
        
        if {$hetero_count > 0} {
            mol representation CPK
            mol color Element
            mol selection "hetero and not water"
            mol material Opaque
            mol addrep $molid
            puts "Added hetero atom representation ($hetero_count atoms)"
        }
    }
    
    puts "Biological pore structure loaded successfully"
    return 1
}

proc ::PETK::gui::calculateStructureBoundingBox {molid {selection "all"}} {
    if {$molid eq ""} {
        puts "Warning: No molecule id available for bounding box calculation"
        return ""
    }
    
    if {[catch {atomselect $molid $selection} sel]} {
        puts "Warning: Unable to create atom selection '$selection' for molid $molid"
        return ""
    }
    
    if {[catch {$sel num} num_atoms] || $num_atoms <= 0} {
        puts "Warning: Selection '$selection' for molid $molid has no atoms"
        catch {$sel delete}
        return ""
    }
    
    set bbox ""
    if {[catch {
        set minmax [measure minmax $sel]
        set min_coord [lindex $minmax 0]
        set max_coord [lindex $minmax 1]
        
        set x_min [lindex $min_coord 0]
        set y_min [lindex $min_coord 1]
        set z_min [lindex $min_coord 2]
        
        set x_max [lindex $max_coord 0]
        set y_max [lindex $max_coord 1]
        set z_max [lindex $max_coord 2]
        
        set center [measure center $sel]
        
        set bbox [dict create \
            num_atoms $num_atoms \
            x_min $x_min \
            x_max $x_max \
            x_range [expr {$x_max - $x_min}] \
            y_min $y_min \
            y_max $y_max \
            y_range [expr {$y_max - $y_min}] \
            z_min $z_min \
            z_max $z_max \
            z_range [expr {$z_max - $z_min}] \
            center $center]
    } err]} {
        puts "Warning: Failed to measure bounding box for selection '$selection': $err"
    }
    
    catch {$sel delete}
    return $bbox
}

proc ::PETK::gui::updateAutoBoxDisplayFromDimensions {box_x box_y box_z {min_z ""} {max_z ""} {updateMovement 1}} {
    if {![string is double $box_x]} { set box_x 0.0 }
    if {![string is double $box_y]} { set box_y 0.0 }
    if {![string is double $box_z]} { set box_z 0.0 }
    
    # Ensure positive dimensions
    if {$box_x < 0} { set box_x [expr {abs($box_x)}] }
    if {$box_y < 0} { set box_y [expr {abs($box_y)}] }
    if {$box_z < 0} { set box_z [expr {abs($box_z)}] }
    
    set half_x [expr {$box_x / 2.0}]
    set half_y [expr {$box_y / 2.0}]
    set half_z [expr {$box_z / 2.0}]
    
    set min_x [expr {floor(-$half_x)}]
    set max_x [expr {ceil($half_x)}]
    set min_y [expr {floor(-$half_y)}]
    set max_y [expr {ceil($half_y)}]
    
    if {$min_z eq "" || $max_z eq ""} {
        set min_z [expr {floor(-$half_z)}]
        set max_z [expr {ceil($half_z)}]
    } else {
        set min_z [expr {floor($min_z)}]
        set max_z [expr {ceil($max_z)}]
    }
    
    set rounded_box_x [expr {$max_x - $min_x}]
    set rounded_box_y [expr {$max_y - $min_y}]
    set rounded_box_z [expr {$max_z - $min_z}]
    
    set ::PETK::gui::autoBoxX [format "%.0f to %.0f (%.1f)" $min_x $max_x $rounded_box_x]
    set ::PETK::gui::autoBoxY [format "%.0f to %.0f (%.1f)" $min_y $max_y $rounded_box_y]
    set ::PETK::gui::autoBoxZ [format "%.0f to %.0f (%.1f)" $min_z $max_z $rounded_box_z]
    
    set ::PETK::gui::calculatedBoxSizeX $rounded_box_x
    set ::PETK::gui::calculatedBoxSizeY $rounded_box_y
    set ::PETK::gui::calculatedBoxSizeZ $rounded_box_z

    if {![info exists ::PETK::gui::poreOption]} {
        set ::PETK::gui::poreOption "solid-state"
    }
    if {![info exists ::PETK::gui::membraneType]} {
        set ::PETK::gui::membraneType "cylindrical"
    }

    if {$::PETK::gui::poreOption eq "solid-state" && $::PETK::gui::membraneType in {"cylindrical" "doublecone" "conical"}} {
        set lateral_size [expr {max($box_x, $box_y)}]
        set bulk_height [expr {$box_z}]
        set thickness 0.0
        if {[string is double $::PETK::gui::nanoporeThickness]} {
            set thickness [expr {double($::PETK::gui::nanoporeThickness)}]
        }
        set bulk_height [expr {$bulk_height - $thickness}]
        if {$bulk_height < 0} {
            set bulk_height 0
        }
        if {$lateral_size > 0} {
            set aspect_ratio [expr {$bulk_height / double($lateral_size)}]
        } else {
            set aspect_ratio 0
        }
        set ::PETK::gui::autoBoxAspectSummary [format "L=%.1f Å, H=%.1f Å  (H/L=%.2f)" $lateral_size $bulk_height $aspect_ratio]
    } else {
        set ::PETK::gui::autoBoxAspectSummary "Bulk height ratio only applies to cylindrical or double-cone pores"
    }
    
    if {$updateMovement && (![info exists ::PETK::gui::semCalculationMode] || $::PETK::gui::semCalculationMode eq "run")} {
        ::PETK::gui::updateMovementRangeFromBounds $min_z $max_z
    }
}

proc ::PETK::gui::updateMovementRangeFromBounds {min_z max_z} {
    if {![string is double -strict $min_z] || ![string is double -strict $max_z]} {
        return
    }
    
    if {$min_z > $max_z} {
        set temp $min_z
        set min_z $max_z
        set max_z $temp
    }
    
    set lower_bound [expr {ceil($min_z)}]
    set upper_bound [expr {floor($max_z)}]
    
    if {$lower_bound >= $upper_bound} {
        set lower_bound $min_z
        set upper_bound $max_z
    }
    
    set ::PETK::gui::zStartRange [::PETK::gui::formatMovementValue $upper_bound]
    set ::PETK::gui::zEndRange [::PETK::gui::formatMovementValue $lower_bound]
}

proc ::PETK::gui::formatMovementValue {value} {
    if {![string is double $value]} {
        return $value
    }
    set formatted [format "%.1f" $value]
    if {$formatted eq "-0.0"} {
        set formatted "0.0"
    }
    if {[string match "*\.0" $formatted]} {
        set formatted [string range $formatted 0 end-2]
    }
    if {$formatted eq "-0"} {
        set formatted "0"
    }
    return $formatted
}

# Procedure to handle biological pore visualization with membrane
proc ::PETK::gui::drawBiologicalPore {min_x max_x min_y max_y min_z max_z} {
    puts "Drawing biological pore visualization..."
    
    # Determine which biological pore file to load
    set pore_file ""
    set pore_source ""
    
    if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
        set pore_file [::PETK::gui::resolveBioPoreFile $::PETK::gui::selectedBioPore]
        set pore_source "pre-prepared"
    } elseif {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
        set pore_file $::PETK::gui::customBioPorePath
        set pore_source "custom"
    }
    
    if {$pore_file eq "" || ![file exists $pore_file]} {
        # No valid biological pore file - draw placeholder
        ::PETK::gui::drawBiologicalPorePlaceholder $min_x $max_x $min_y $max_y $min_z $max_z
        set ::PETK::gui::poreValidityStatus "WARNING - No biological pore file selected"
        return
    }
    
    # Try to load the biological pore structure
    set load_success [::PETK::gui::loadBiologicalPoreStructure $pore_file]
    
    if {$load_success} {
        puts "Successfully loaded biological pore from $pore_source: [file tail $pore_file]"
        
        # Apply configuration for this pore (important!)
        if {$pore_source eq "pre-prepared"} {
            ::PETK::gui::applyBiologicalPoreConfiguration $pore_file
        } else {
            puts "Custom pore loaded - using current parameters"
        }
        
        # Perform comprehensive validation (like solid-state)
        set validation_status [::PETK::gui::validateBiologicalPoreConfiguration]
        
        # Draw the solid state membrane around the biological pore
        ::PETK::gui::drawSolidStateMembrane $min_x $max_x $min_y $max_y $min_z $max_z        
    } else {
        # Failed to load - draw placeholder with error indication
        ::PETK::gui::drawBiologicalPorePlaceholder $min_x $max_x $min_y $max_y $min_z $max_z "error"
        set ::PETK::gui::poreValidityStatus "ERROR - Failed to load biological pore file"
        puts "Failed to load biological pore structure from: $pore_file"
    }
}

# Enhanced placeholder procedure with fixed positioning
proc ::PETK::gui::drawBiologicalPorePlaceholder {min_x max_x min_y max_y min_z max_z {mode "normal"}} {
    if {$mode eq "error"} {
        draw color red
        set status_text "ERROR: LOAD FAILED"
    } else {
        draw color orange
        set status_text "NO PORE SELECTED"
    }
    
    draw material Transparent
    
    # Draw a generic pore-like structure as placeholder
    set center_x 0.0
    set center_y 0.0
    set center_z 0.0
    set placeholder_radius 15.0
    # Placeholder height
    set placeholder_height 30.0
    
    # Draw cylindrical placeholder
    for {set angle 0} {$angle < 360} {incr angle 30} {
        set angle_rad [expr $angle * 3.14159 / 180.0]
        set next_angle_rad [expr ($angle + 30) * 3.14159 / 180.0]
        set x1 [expr $placeholder_radius * cos($angle_rad)]
        set y1 [expr $placeholder_radius * sin($angle_rad)]
        set x2 [expr $placeholder_radius * cos($next_angle_rad)]
        set y2 [expr $placeholder_radius * sin($next_angle_rad)]
        
        # Draw outline circles
        draw line [list $x1 $y1 [expr -$placeholder_height/2]] [list $x2 $y2 [expr -$placeholder_height/2]] width 4
        draw line [list $x1 $y1 [expr $placeholder_height/2]] [list $x2 $y2 [expr $placeholder_height/2]] width 4
        
        # Draw vertical lines
        draw line [list $x1 $y1 [expr -$placeholder_height/2]] [list $x1 $y1 [expr $placeholder_height/2]] width 4
    }
    
    # Fixed text positioning - centered within the box, away from edges
    draw color white
    set text_x [expr ($min_x + $max_x) / 2.0]
    set text_y [expr ($min_y + $max_y) / 2.0] 
    # Position text 20 Å from top, well within box
    set text_z [expr $max_z - 20.0]
    draw text [list $text_x $text_y $text_z] $status_text size 1.2
    
    puts "Drew biological pore placeholder ($mode mode)"
}

# New procedure to draw solid state membrane around biological pore
proc ::PETK::gui::drawSolidStateMembrane {min_x max_x min_y max_y min_z max_z} {
    if {![info exists ::PETK::gui::biological_pore_molid]} {
        puts "No biological pore loaded - cannot create membrane"
        return
    }
    
    # Get membrane parameters (use nanopore thickness and any membrane offset)
    set membrane_thickness $::PETK::gui::nanoporeThickness
    set membrane_half_thickness [expr $membrane_thickness / 2.0]
    
    # Get membrane Z offset (default to 0 if not set)
    if {[info exists ::PETK::gui::membraneZOffset]} {
        set membrane_z_offset $::PETK::gui::membraneZOffset
    } else {
        set membrane_z_offset 0.0
    }
    
    # Calculate membrane boundaries using EXACTLY same naming as solid-state
    # Center the pore at offset
    set pore_center_z $membrane_z_offset
    set pore_start_z [expr $pore_center_z - $membrane_half_thickness]
    set pore_end_z [expr $pore_center_z + $membrane_half_thickness]
    
    puts "Creating solid state membrane around biological pore:"
    puts "  Membrane thickness: $membrane_thickness Å"
    puts "  Membrane Z-offset: $membrane_z_offset Å"
    puts "  Membrane boundaries: [format %.1f $pore_start_z] to [format %.1f $pore_end_z] Å"
    
    # Add validation warnings first (EXACTLY same as solid-state)
    ::PETK::gui::addBiologicalPoreValidationWarning $min_x $max_x $min_y $max_y $min_z $max_z
    
    # Calculate local pore radius at different Z positions
    # 2 Å steps for Z sampling
    set z_step 2.0
    set z_positions [list]
    set local_radii [list]
    
    # Sample Z positions through the membrane region
    for {set z $pore_start_z} {$z <= $pore_end_z} {set z [expr $z + $z_step]} {
        lappend z_positions $z
        
        # Get local radius at this Z position
        set local_radius [::PETK::gui::calculateLocalPoreRadius $z $membrane_z_offset $membrane_half_thickness]
        lappend local_radii $local_radius
        
        puts "  Z = [format %.1f $z] Å: local radius = [format %.1f $local_radius] Å"
    }
    
    # Draw the membrane structure using grid-based approach
    ::PETK::gui::visualizeMembrane $z_positions $local_radii $membrane_z_offset $membrane_half_thickness $min_x $max_x $min_y $max_y
    
    # Final status message (EXACTLY like solid-state pattern)
    if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus] || 
        [string match "*ERROR*" $::PETK::gui::poreValidityStatus]} {
        puts "Pore visualization completed - WARNING: INVALID CONFIGURATION!"
        puts "See detailed validation report above."
    } else {
        puts "Pore visualization completed for biological pore with purple membrane"
        puts "  Variable pore radius adapts to biological structure"
        puts "  Membrane positioned at Z-offset: $membrane_z_offset Å"
    }
    
    puts "Current validity status: $::PETK::gui::poreValidityStatus"
}

# Calculate local pore radius at a given Z position
proc ::PETK::gui::calculateLocalPoreRadius {z_pos membrane_z_offset membrane_half_thickness} {
    # Check if biological pore molecule exists
    if {![info exists ::PETK::gui::biological_pore_molid] || $::PETK::gui::biological_pore_molid eq ""} {
        return 0.0
    }
    
    set molid $::PETK::gui::biological_pore_molid
    
    # Check if the molecule still exists in VMD
    if {[catch {molinfo $molid get numatoms} num_atoms]} {
        return 0.0
    }
    
    if {$num_atoms == 0} {
        return 0.0
    }
    
    # Default large radius for regions FAR outside membrane (with buffer zone)
    set default_radius 50.0
    
    # FIXED: Use strict inequality (>) instead of (>=) to include membrane boundaries
    # Add small buffer zone beyond membrane boundaries before using default radius
    set boundary_buffer 5.0  ; # Allow 5 Å beyond membrane boundaries for biological calculation
    
    # Check if Z position is WELL outside membrane region (not just at boundaries)
    if {[expr abs($z_pos - $membrane_z_offset)] > ($membrane_half_thickness + $boundary_buffer)} {
        puts "  Z=$z_pos is far outside membrane (>[expr $membrane_half_thickness + $boundary_buffer] Å from center), using default radius $default_radius"
        return $default_radius
    }
    
    # For Z positions within or near membrane boundaries, ALWAYS calculate biological pore radius
    puts "  Z=$z_pos is within/near membrane region, calculating biological pore radius..."
    
    # Find atoms near this Z position (within ±2 Å slice)
    set slice_tolerance 2.0
    
    # Create atomselection with error checking
    if {[catch {
        set sel [atomselect $molid "abs(z - $z_pos) <= $slice_tolerance"]
    } sel_error]} {
        puts "Warning: Failed to create atomselection in calculateLocalPoreRadius: $sel_error"
        return 0.0
    }
    
    # Check if selection is valid and has atoms
    if {[catch {set num_atoms_in_slice [$sel num]} num_error]} {
        puts "Warning: Failed to get atom count: $num_error"
        catch {$sel delete}
        return 0.0
    }
    
    # MINIMAL MODIFICATION: Show dialog and stop process when no atoms found
    if {$num_atoms_in_slice == 0} {
        $sel delete
        
        # Display information dialog and stop process
        set message "No atoms found in Z-slice during biological pore calculation!\n\n"
        append message "Z-position: $z_pos Å\n"
        append message "Slice range: [expr $z_pos - $slice_tolerance] to [expr $z_pos + $slice_tolerance] Å\n"
        append message "Membrane center: $membrane_z_offset Å\n"
        append message "Membrane thickness: [expr $membrane_half_thickness * 2] Å\n\n"
        append message "This might indicate:\n"
        append message "• Open pore region (normal for some biological pores)\n"
        append message "• Incorrect membrane positioning\n"
        append message "• Structure gaps in the biological pore\n\n"
        append message "Process will be stopped for investigation."
        
        tk_messageBox -icon warning -title "No Atoms Found in Slice" -message $message
        
        # Stop the process
        error "Process stopped: No atoms found at Z=$z_pos Å"
    }
    
    # Get atom positions and elements with error checking
    if {[catch {
        set positions [$sel get {x y z}]
        set elements [$sel get element]
    } get_error]} {
        puts "Warning: Failed to get atom properties: $get_error"
        $sel delete
        return 0.0
    }
    
    $sel delete
    
    # Calculate minimum local radius based on atom positions
    # Start with large value
    set min_extent 1000.0
    
    foreach pos $positions element $elements {
        set x [lindex $pos 0]
        set y [lindex $pos 1]
        set z [lindex $pos 2]
        
        # Calculate radial distance from Z-axis
        set radial_dist [expr sqrt($x*$x + $y*$y)]
        
        # Get approximate atom radius (simplified)
        set atom_radius [::PETK::gui::getAtomRadius $element]
        
        # Calculate extent (radial distance + atom radius)
        set extent [expr $radial_dist + $atom_radius]
        
        if {$extent < $min_extent} {
            set min_extent $extent
        }
    }
    
    # Ensure minimum radius for membrane (at least 2 Å from atoms)
    set min_membrane_radius [expr $min_extent + 2.0]
    
    puts "    Calculated biological pore radius: [format %.2f $min_membrane_radius] Å at Z=$z_pos"
    
    return $min_membrane_radius
}

# Get approximate atomic radius for common elements
proc ::PETK::gui::getAtomRadius {element} {
    switch -nocase $element {
        "H" {return 1.2}
        "C" {return 1.7}
        "N" {return 1.55}
        "O" {return 1.52}
        "P" {return 1.8}
        "S" {return 1.8}
        default {return 1.5}
    }
}

# Updated visualizeMembrane function for biological pores - Simple cylindrical approach
proc ::PETK::gui::visualizeMembrane {z_positions local_radii membrane_z_offset membrane_half_thickness min_x max_x min_y max_y} {
    
    # Set membrane visualization properties (same colors as solid-state)
    draw color purple
    draw material Opaque
    
    # Calculate membrane boundaries
    set pore_start_z [expr $membrane_z_offset - $membrane_half_thickness]
    set pore_end_z [expr $membrane_z_offset + $membrane_half_thickness]
    set thickness [expr $pore_end_z - $pore_start_z]
    
    # Grid-based membrane drawing (same grid approach but simpler logic)
    set grid_size 5.0
    
    # For biological pores, we use simple cylindrical approach (NO chamfer effects)
    # Number of Z levels for membrane visualization
    set z_levels 10
    
    puts "Drawing biological pore membrane with simple cylindrical approach (no chamfer)..."
    puts "  Membrane Z range: [format %.1f $pore_start_z] to [format %.1f $pore_end_z] Å"
    puts "  Grid size: $grid_size Å"
    puts "  Membrane type: Simple cylindrical (adapts to biological structure)"
    
    # Draw membrane at multiple Z levels using simple cylindrical logic
    for {set z_level 0} {$z_level < $z_levels} {incr z_level} {
        # Calculate Z position for this level
        set z_frac [expr double($z_level) / ($z_levels - 1)]
        set z_pos [expr $pore_start_z + $z_frac * $thickness]
        
        # Get local pore radius at this Z position (varies based on biological structure)
        set local_pore_radius [::PETK::gui::calculateLocalPoreRadius $z_pos $membrane_z_offset $membrane_half_thickness]
        
        # Skip levels with no valid radius
        if {$local_pore_radius <= 0} continue
        
        # Draw membrane grid at this Z level using simple cylindrical exclusion
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
                
                # Only draw if cell center is outside local pore radius (simple cylindrical test)
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
    
    # Draw pore outline to show the variable radius hole (simple cylindrical style)
    draw color yellow
    set angle_step 15
    # Fewer levels for outline
    set outline_z_levels 5
    
    puts "Drawing cylindrical pore outline (adapts to biological structure)..."
    
    for {set z_level 0} {$z_level < $outline_z_levels} {incr z_level} {
        set z_frac [expr double($z_level) / ($outline_z_levels - 1)]
        set z_pos [expr $pore_start_z + $z_frac * $thickness]
        
        # Get local pore radius at this Z level
        set local_pore_radius [::PETK::gui::calculateLocalPoreRadius $z_pos $membrane_z_offset $membrane_half_thickness]
        
        if {$local_pore_radius > 0} {
            # Draw simple circle at this Z level (NO chamfer effects - just regular cylindrical)
            for {set angle 0} {$angle < 360} {incr angle $angle_step} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + $angle_step) * 3.14159 / 180.0]
                set x1 [expr $local_pore_radius * cos($angle_rad)]
                set y1 [expr $local_pore_radius * sin($angle_rad)]
                set x2 [expr $local_pore_radius * cos($next_angle_rad)]
                set y2 [expr $local_pore_radius * sin($next_angle_rad)]
                
                draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $z_pos] width 4
            }
            
            # Draw connecting lines between levels (if not the last level)
            if {$z_level < [expr $outline_z_levels - 1]} {
                set next_z_frac [expr double($z_level + 1) / ($outline_z_levels - 1)]
                set next_z_pos [expr $pore_start_z + $next_z_frac * $thickness]
                set next_local_pore_radius [::PETK::gui::calculateLocalPoreRadius $next_z_pos $membrane_z_offset $membrane_half_thickness]
                
                if {$next_local_pore_radius > 0} {
                    # Draw simple vertical connecting lines (no chamfer effects)
                    for {set angle 0} {$angle < 360} {incr angle [expr $angle_step * 2]} {
                        set angle_rad [expr $angle * 3.14159 / 180.0]
                        set x1 [expr $local_pore_radius * cos($angle_rad)]
                        set y1 [expr $local_pore_radius * sin($angle_rad)]
                        set x2 [expr $next_local_pore_radius * cos($angle_rad)]
                        set y2 [expr $next_local_pore_radius * sin($angle_rad)]
                        
                        draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $next_z_pos] width 4
                    }
                }
            }
        }
    }
    
    puts "Biological pore membrane visualization complete:"
    puts "  - Simple cylindrical membrane (no chamfer effects)"
    puts "  - Variable radius adapts to biological pore structure"
    puts "  - Membrane embeds the biological pore"
}

# Add validation warnings for biological pore membrane (EXACTLY like solid-state)
proc ::PETK::gui::addBiologicalPoreValidationWarning {min_x max_x min_y max_y min_z max_z} {
    if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus] || 
        [string match "*ERROR*" $::PETK::gui::poreValidityStatus]} {
        
        draw color red
        draw material Transparent
        
        # Draw warning outline around the membrane region (EXACTLY like solid-state)
        if {[info exists ::PETK::gui::nanoporeThickness]} {
            set membrane_half_thickness [expr $::PETK::gui::nanoporeThickness / 2.0]
            set membrane_z_offset 0.0
            if {[info exists ::PETK::gui::membraneZOffset]} {
                set membrane_z_offset $::PETK::gui::membraneZOffset
            }
            
            set pore_start_z [expr $membrane_z_offset - $membrane_half_thickness]
            set pore_end_z [expr $membrane_z_offset + $membrane_half_thickness]
            
            # Get a reasonable warning radius based on typical biological pore size
            set warn_radius 30.0
            
            # Draw warning outline EXACTLY like solid-state (same pattern)
            for {set angle 0} {$angle < 360} {incr angle 30} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + 30) * 3.14159 / 180.0]
                set x1 [expr $warn_radius * cos($angle_rad)]
                set y1 [expr $warn_radius * sin($angle_rad)]
                set x2 [expr $warn_radius * cos($next_angle_rad)]
                set y2 [expr $warn_radius * sin($next_angle_rad)]
                
                # Use EXACTLY same thick line style as solid-state
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] width 10
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] width 10
            }
            
            # Reset to purple for membrane drawing (EXACTLY like solid-state)
            draw color purple
            draw material Opaque
        }
        
        puts "WARNING: Biological pore configuration issues detected!"
    }
}

# Add this to your GUI variables initialization
proc ::PETK::gui::initializeBiologicalPoreVariables {} {
    if {[info exists ::PETK::gui::biologicalVarsInitialized] && $::PETK::gui::biologicalVarsInitialized} {
        return
    }
    # Initialize biological pore specific variables
    if {![info exists ::PETK::gui::membraneZOffset]} {
        set ::PETK::gui::membraneZOffset 0.0
    }
    
    if {![info exists ::PETK::gui::biological_pore_molid]} {
        set ::PETK::gui::biological_pore_molid ""
    }
    
    #if {![info exists ::PETK::gui::showMembraneOutline]} {
    #    set ::PETK::gui::showMembraneOutline 1
    #}
    
    if {![info exists ::PETK::gui::showRadiusCalculation]} {
        set ::PETK::gui::showRadiusCalculation 0
    }
    
    # Initialize pore option tracking for cleanup
    if {![info exists ::PETK::gui::previousPoreOption]} {
        set ::PETK::gui::previousPoreOption ""
    }
    
    # Initialize pore selection variables
    if {![info exists ::PETK::gui::selectedBioPore]} {
        set ::PETK::gui::selectedBioPore ""
    }
    
    if {![info exists ::PETK::gui::customBioPorePath]} {
        set ::PETK::gui::customBioPorePath ""
    }
    
    if {![info exists ::PETK::gui::customBioBoundingInfo]} {
        set ::PETK::gui::customBioBoundingInfo ""
    }
    
    # Ensure nanopore thickness is set for biological pores
    if {![info exists ::PETK::gui::nanoporeThickness]} {
        set ::PETK::gui::nanoporeThickness 20.0
    }
    
    # Ensure box dimensions are initialized
    if {![info exists ::PETK::gui::boxSizeX]} {
        set ::PETK::gui::boxSizeX 150.0
    }
    if {![info exists ::PETK::gui::boxSizeY]} {
        set ::PETK::gui::boxSizeY 150.0
    }
    if {![info exists ::PETK::gui::boxSizeZ]} {
        set ::PETK::gui::boxSizeZ 150.0
    }
    
    # Initialize biological pore configurations
    ::PETK::gui::initializeBiologicalPoreConfigurations
    set ::PETK::gui::biologicalVarsInitialized 1
}

# Clean up biological pore resources when switching to solid-state mode
proc ::PETK::gui::cleanupBiologicalPoreResources {} {
    set prior_molid ""
    # Remove biological pore molecule if it exists
    if {[info exists ::PETK::gui::biological_pore_molid] && $::PETK::gui::biological_pore_molid ne ""} {
        set prior_molid $::PETK::gui::biological_pore_molid
        if {[catch {
            mol delete $::PETK::gui::biological_pore_molid
            puts "Cleaned up biological pore molecule (ID: $::PETK::gui::biological_pore_molid)"
        } cleanup_error]} {
            puts "Warning: Could not delete biological pore molecule: $cleanup_error"
        }
        set ::PETK::gui::biological_pore_molid ""
    }
    

    # Clear any PETK preview drawings associated with prior molecules
    ::PETK::gui::clearPreviewGraphics
    
    if {$prior_molid ne "" && [info exists ::PETK::gui::petk_draw_molids]} {
        set idx [lsearch -exact $::PETK::gui::petk_draw_molids $prior_molid]
        if {$idx != -1} {
            set ::PETK::gui::petk_draw_molids [lreplace $::PETK::gui::petk_draw_molids $idx $idx]
        }
    }

    set ::PETK::gui::customBioBoundingInfo ""
    
    puts "Biological pore resources cleaned up"
}

# Comprehensive validation for biological pore configuration (like solid-state)
proc ::PETK::gui::validateBiologicalPoreConfiguration {} {
    if {![info exists ::PETK::gui::biological_pore_molid] || $::PETK::gui::biological_pore_molid eq ""} {
        set ::PETK::gui::poreValidityStatus "INVALID - No biological pore loaded"
        return "INVALID"
    }
    
    set molid $::PETK::gui::biological_pore_molid
    set issues_found 0
    set warnings [list]
    
    puts "=== BIOLOGICAL PORE VALIDATION ==="
    
    # Get pore structure information
    set all_atoms [atomselect $molid "all"]
    set num_atoms [$all_atoms num]
    set minmax [measure minmax $all_atoms]
    $all_atoms delete
    
    if {$num_atoms == 0} {
        set ::PETK::gui::poreValidityStatus "INVALID - No biological pore loaded"
        puts "✗ No atoms found in biological pore structure"
        return "INVALID"
    }
    
    puts "Biological pore contains $num_atoms atoms"
    
    # Get membrane parameters
    if {![info exists ::PETK::gui::nanoporeThickness] || ![info exists ::PETK::gui::membraneZOffset]} {
        set ::PETK::gui::poreValidityStatus "VALID - Biological pore configuration OK"
        puts "✓ Validation passed - membrane parameters not set yet"
        return "VALID"
    }
    
    set membrane_thickness $::PETK::gui::nanoporeThickness
    set membrane_z_offset $::PETK::gui::membraneZOffset
    
    puts "Membrane thickness: $membrane_thickness Å"
    puts "Membrane Z-offset: $membrane_z_offset Å"
    
    # Get pore Z range
    set min_coord [lindex $minmax 0]
    set max_coord [lindex $minmax 1]
    set pore_z_min [lindex $min_coord 2]
    set pore_z_max [lindex $max_coord 2]
    set pore_z_span [expr $pore_z_max - $pore_z_min]
    
    # Get membrane Z range
    set membrane_half_thickness [expr $membrane_thickness / 2.0]
    set membrane_bottom [expr $membrane_z_offset - $membrane_half_thickness]
    set membrane_top [expr $membrane_z_offset + $membrane_half_thickness]
    
    puts "Pore Z span: [format %.1f $pore_z_min] to [format %.1f $pore_z_max] Å ([format %.1f $pore_z_span] Å)"
    puts "Membrane region: [format %.1f $membrane_bottom] to [format %.1f $membrane_top] Å"
    
    # ONLY CHECK: Does membrane extend outside pore Z range?
    set membrane_extends_below 0
    set membrane_extends_above 0
    
    if {$membrane_bottom < $pore_z_min} {
        set membrane_extends_below 1
        set extension_below [expr $pore_z_min - $membrane_bottom]
        lappend warnings "Membrane extends [format %.1f $extension_below] Å below pore structure"
        incr issues_found
    }
    
    if {$membrane_top > $pore_z_max} {
        set membrane_extends_above 1
        set extension_above [expr $membrane_top - $pore_z_max]
        lappend warnings "Membrane extends [format %.1f $extension_above] Å above pore structure"
        incr issues_found
    }
    
    # Show results
    if {$issues_found == 0} {
        set ::PETK::gui::poreValidityStatus "VALID - Biological pore configuration OK"
        puts "✓ Validation passed - membrane is within pore Z range"
        return "VALID"
    } else {
        set ::PETK::gui::poreValidityStatus "WARNING - Membrane extends outside pore range"
        puts "⚠ Validation warning - membrane positioning issue:"
        
        foreach warning $warnings {
            puts "  - $warning"
        }
        
        # Show warning dialog for membrane extension
        ::PETK::gui::showMembraneExtensionWarning $warnings $pore_z_min $pore_z_max $membrane_bottom $membrane_top $membrane_extends_below $membrane_extends_above
        
        return "WARNING"
    }
}

# SIMPLIFIED: Show warning only for membrane extension outside pore range
proc ::PETK::gui::showMembraneExtensionWarning {warnings pore_z_min pore_z_max membrane_bottom membrane_top extends_below extends_above} {
    
    # Build focused message
    set message "Membrane extends outside biological pore Z-range!\n\n"
    
    # Add positioning info
    append message "POSITIONING:\n"
    append message "• Pore Z-range: [format %.1f $pore_z_min] to [format %.1f $pore_z_max] Å\n"
    append message "• Membrane region: [format %.1f $membrane_bottom] to [format %.1f $membrane_top] Å\n\n"
    
    # Add specific extension details
    append message "MEMBRANE EXTENSIONS:\n"
    foreach warning $warnings {
        append message "• $warning\n"
    }
    
    # Add recommendations
    append message "\nRECOMMENDATIONS:\n"
    
    if {$extends_below && $extends_above} {
        # Membrane extends both ways - suggest centering
        set suggested_offset [expr ($pore_z_min + $pore_z_max) / 2.0]
        append message "• Center membrane on pore: set Z-offset to ~[format %.0f $suggested_offset] Å\n"
        
        # Calculate what thickness would fit
        set pore_span [expr $pore_z_max - $pore_z_min]
        set suggested_thickness [expr $pore_span * 0.8]  ;# 80% of pore span
        append message "• Reduce membrane thickness to ~[format %.0f $suggested_thickness] Å to fit within pore\n"
        
    } elseif {$extends_below} {
        # Only extends below - suggest moving up
        set needed_shift [expr $pore_z_min - $membrane_bottom + 5]  ;# 5 Å buffer
        set suggested_offset [expr $::PETK::gui::membraneZOffset + $needed_shift]
        append message "• Move membrane up: set Z-offset to ~[format %.0f $suggested_offset] Å\n"
        
    } elseif {$extends_above} {
        # Only extends above - suggest moving down
        set needed_shift [expr $membrane_top - $pore_z_max + 5]  ;# 5 Å buffer
        set suggested_offset [expr $::PETK::gui::membraneZOffset - $needed_shift]
        append message "• Move membrane down: set Z-offset to ~[format %.0f $suggested_offset] Å\n"
    }
    
    append message "\nNote: This is just a positioning warning. The simulation may still work fine.\n\n"
    append message "Continue with current configuration?"
    
    # Show warning dialog
    set user_choice [tk_messageBox -icon warning -title "Membrane Extension Warning" -message $message -type yesno]
    
    if {$user_choice eq "no"} {
        error "Process stopped by user due to membrane extension outside pore range"
    }
    
    puts "User chose to continue with membrane extending outside pore range."
}

# Original solid-state pore drawing procedure (extracted from main function)
proc ::PETK::gui::drawSolidStatePore {min_x max_x min_y max_y min_z max_z xbox ybox zbox} {
    # Perform comprehensive pore validation using the new validation system
    set validation_status [::PETK::gui::updatePoreValidityStatus $xbox $ybox $zbox]
    
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
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] width 10
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] width 10
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
                    
                    draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $z_pos] width 4
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
                        
                        draw line [list $x1 $y1 $z_pos] [list $x2 $y2 $next_z_pos] width 4
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
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] width 4
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] width 4
                draw line [list $x1 $y1 $pore_start_z] [list $x1 $y1 $pore_end_z] width 4
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
                
                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] width 4
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] width 4
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
            draw line [list $x1_out $y1_out $pore_start_z] [list $x2_out $y2_out $pore_start_z] width 4
            draw line [list $x1_in $y1_in $middle_z] [list $x2_in $y2_in $middle_z] width 4
            draw line [list $x1_out $y1_out $pore_start_z] [list $x1_in $y1_in $middle_z] width 4
            draw line [list $x1_in $y1_in $middle_z] [list $x1_out $y1_out $pore_end_z] width 4
            draw line [list $x1_out $y1_out $pore_end_z] [list $x2_out $y2_out $pore_end_z] width 4
        }

    } elseif {$::PETK::gui::membraneType eq "conical"} {
        # Draw membrane with single truncated-cone (frustum) hole.
        # top_radius is at z = pore_end_z (+half_thickness),
        # bottom_radius is at z = pore_start_z (-half_thickness).
        set top_radius [expr $::PETK::gui::topDiameter / 2.0]
        set bottom_radius [expr $::PETK::gui::bottomDiameter / 2.0]
        set max_radius [expr {max($top_radius, $bottom_radius)}]

        # Add visual warning if configuration is invalid
        if {[string match "*INVALID*" $::PETK::gui::poreValidityStatus]} {
            draw color red
            draw material Transparent
            set warn_radius [expr $max_radius + 10.0]

            for {set angle 0} {$angle < 360} {incr angle 30} {
                set angle_rad [expr $angle * 3.14159 / 180.0]
                set next_angle_rad [expr ($angle + 30) * 3.14159 / 180.0]
                set x1 [expr $warn_radius * cos($angle_rad)]
                set y1 [expr $warn_radius * sin($angle_rad)]
                set x2 [expr $warn_radius * cos($next_angle_rad)]
                set y2 [expr $warn_radius * sin($next_angle_rad)]

                draw line [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] width 4
                draw line [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] width 4
            }

            draw color purple
            draw material Opaque
        }

        # Draw membrane as rectangular sheets with a frustum hole.
        # On the bottom face (z = pore_start_z) the hole has radius bottom_radius.
        # On the top face (z = pore_end_z) the hole has radius top_radius.
        set grid_size 5.0
        for {set x $min_x} {$x < $max_x} {set x [expr $x + $grid_size]} {
            for {set y $min_y} {$y < $max_y} {set y [expr $y + $grid_size]} {
                set x1 $x
                set x2 [expr min($x + $grid_size, $max_x)]
                set y1 $y
                set y2 [expr min($y + $grid_size, $max_y)]

                set cx [expr ($x1 + $x2) / 2.0]
                set cy [expr ($y1 + $y2) / 2.0]
                set distance [expr sqrt($cx*$cx + $cy*$cy)]

                # Bottom face — hole has radius bottom_radius
                if {$distance > $bottom_radius} {
                    draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z]
                    draw triangle [list $x1 $y1 $pore_start_z] [list $x2 $y2 $pore_start_z] [list $x1 $y2 $pore_start_z]
                }

                # Top face — hole has radius top_radius
                if {$distance > $top_radius} {
                    draw triangle [list $x1 $y1 $pore_end_z] [list $x2 $y2 $pore_end_z] [list $x2 $y1 $pore_end_z]
                    draw triangle [list $x1 $y1 $pore_end_z] [list $x1 $y2 $pore_end_z] [list $x2 $y2 $pore_end_z]
                }
            }
        }

        # Draw the frustum outline in yellow.
        draw color yellow
        set angle_step 20
        for {set angle 0} {$angle < 360} {incr angle $angle_step} {
            set angle_rad [expr $angle * 3.14159 / 180.0]
            set next_angle_rad [expr ($angle + $angle_step) * 3.14159 / 180.0]

            set xt1 [expr $top_radius * cos($angle_rad)]
            set yt1 [expr $top_radius * sin($angle_rad)]
            set xt2 [expr $top_radius * cos($next_angle_rad)]
            set yt2 [expr $top_radius * sin($next_angle_rad)]

            set xb1 [expr $bottom_radius * cos($angle_rad)]
            set yb1 [expr $bottom_radius * sin($angle_rad)]
            set xb2 [expr $bottom_radius * cos($next_angle_rad)]
            set yb2 [expr $bottom_radius * sin($next_angle_rad)]

            # Top circle (at +half), bottom circle (at -half), and the slanted wall.
            draw line [list $xt1 $yt1 $pore_end_z]   [list $xt2 $yt2 $pore_end_z]   width 4
            draw line [list $xb1 $yb1 $pore_start_z] [list $xb2 $yb2 $pore_start_z] width 4
            draw line [list $xb1 $yb1 $pore_start_z] [list $xt1 $yt1 $pore_end_z]   width 4
        }
    }
    
    # Final status message for solid-state pores
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

    } elseif {$::PETK::gui::membraneType eq "conical"} {
        # For conical (single frustum) pores use the smaller of top/bottom
        # — this is the constriction the analyte must clear.
        if {[info exists ::PETK::gui::topDiameter] &&
            [info exists ::PETK::gui::bottomDiameter] &&
            [string is double $::PETK::gui::topDiameter] &&
            [string is double $::PETK::gui::bottomDiameter] &&
            $::PETK::gui::topDiameter > 0 &&
            $::PETK::gui::bottomDiameter > 0} {

            if {$::PETK::gui::topDiameter <= $::PETK::gui::bottomDiameter} {
                set ::PETK::gui::poreDiameter $::PETK::gui::topDiameter
                puts "Pore diameter set to top diameter (smaller): $::PETK::gui::poreDiameter Å"
            } else {
                set ::PETK::gui::poreDiameter $::PETK::gui::bottomDiameter
                puts "Pore diameter set to bottom diameter (smaller): $::PETK::gui::poreDiameter Å"
            }
        } else {
            set ::PETK::gui::poreDiameter ""
            puts "Warning: Invalid or missing top/bottom diameters"
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

    } elseif {$::PETK::gui::membraneType eq "conical"} {
        set top_radius [expr $::PETK::gui::topDiameter / 2.0]
        set bottom_radius [expr $::PETK::gui::bottomDiameter / 2.0]
        set max_radius [expr {max($top_radius, $bottom_radius)}]

        # The widest mouth is what matters for box-fit
        if {$max_radius > $min_box_half_xy} {
            set is_valid 0
            lappend critical_errors "Max pore radius ($max_radius Å) > minimum box half-dimension ($min_box_half_xy Å)"
            lappend critical_errors "  Box half-dimensions: X=$box_half_x Å, Y=$box_half_y Å"
            lappend critical_errors "  Required: max(top, bottom) radius < $min_box_half_xy Å"
        } elseif {$max_radius > [expr 0.8 * $min_box_half_xy]} {
            lappend warnings "Max pore radius ($max_radius Å) > 80% of minimum box half-dimension"
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
        } elseif {$::PETK::gui::membraneType eq "conical"} {
            set available_space [expr $min_box_half_xy - $max_radius]
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
    } elseif {$::PETK::gui::membraneType eq "doublecone"} {
        append detailed_message "Inner diameter: $::PETK::gui::innerDiameter Å\n"
        append detailed_message "Outer diameter: $::PETK::gui::outerDiameter Å\n"
    } elseif {$::PETK::gui::membraneType eq "conical"} {
        append detailed_message "Top diameter: $::PETK::gui::topDiameter Å\n"
        append detailed_message "Bottom diameter: $::PETK::gui::bottomDiameter Å\n"
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
    
    # Print detailed validation report only if it changed to avoid duplicates
    if {![info exists ::PETK::gui::lastValidationReport] || $::PETK::gui::lastValidationReport ne $detailed_message} {
        set ::PETK::gui::lastValidationReport $detailed_message
        puts $detailed_message
    }
    
    return $status
}
