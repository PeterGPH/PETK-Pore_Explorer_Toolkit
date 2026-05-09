#
# PETK GUI - Tab 2 (Analyte Setup) module
#

if {![info exists ::PETK::gui::aminoAcidFriendlyNames]} {
    set ::PETK::gui::aminoAcidFriendlyNames [dict create \
        ALA "Alanine" \
        ARG "Arginine" \
        ASN "Asparagine" \
        ASP "Aspartic acid" \
        CYS "Cysteine" \
        GLN "Glutamine" \
        GLU "Glutamic acid" \
        GLY "Glycine" \
        HIS "Histidine" \
        ILE "Isoleucine" \
        LEU "Leucine" \
        LYS "Lysine" \
        MET "Methionine" \
        PHE "Phenylalanine" \
        PRO "Proline" \
        SER "Serine" \
        THR "Threonine" \
        TRP "Tryptophan" \
        TYR "Tyrosine" \
        VAL "Valine" \
    ]
}

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
    grid rowconfigure $container {4 5} -weight 1  ; # Make results and details expandable
    
    set row 0

    if {![info exists ::PETK::gui::selectedCenteredAmino]} {
        set ::PETK::gui::selectedCenteredAmino ""
    }
    if {![info exists ::PETK::gui::centeredAminoMap]} {
        set ::PETK::gui::centeredAminoMap [dict create]
    }
    
    # === Analyte Source Frame ===
    ttk::labelframe $container.source -text "Analyte Source" -padding 10
    grid $container.source -column 0 -row $row -sticky ew -pady "0 10" -padx 10
    grid columnconfigure $container.source {1 3} -weight 1
    incr row

    ttk::label  $container.source.pdblbl -text "PDB file:" -width 18
    ttk::entry  $container.source.pdb -textvariable ::PETK::gui::analytePDB -width 35
    ttk::button $container.source.browse -text "Browse..." -command ::PETK::gui::browseAnalytePdb

    ttk::label  $container.source.sellbl -text "Atom selection:" -width 18
    ttk::entry  $container.source.sel -textvariable ::PETK::gui::analyteSelection -width 35
    ttk::label  $container.source.selhelp -text "Examples: 'all', 'protein', 'not water'" -font {TkDefaultFont 9 italic} -foreground gray

    ttk::label  $container.source.aflbl -text "AlphaFold UniProt:" -width 18
    ttk::entry  $container.source.afid -textvariable ::PETK::gui::alphafoldUniProt -width 35
    ttk::button $container.source.afbtn -text "Download from AlphaFold" -command ::PETK::gui::downloadAlphaFoldAnalyte
    ttk::label  $container.source.afhelp -text "Enter UniProt ID (AlphaFold) or 4-letter PDB code (RCSB) and click download" -font {TkDefaultFont 9 italic} -foreground gray

    ttk::label  $container.source.centeredlbl -text "Centered amino acids:" -width 18
    ttk::combobox $container.source.centeredcombo -textvariable ::PETK::gui::selectedCenteredAmino \
        -width 35 -state readonly
    ttk::button $container.source.centeredload -text "Load Sample" -command ::PETK::gui::loadSelectedCenteredAmino
    ttk::button $container.source.centeredrefresh -text "Refresh" -command ::PETK::gui::refreshCenteredAminoOptions
    ttk::label  $container.source.centeredhelp -text "Load any of the 20 canonical amino acids centered to the origin." \
        -font {TkDefaultFont 9 italic} -foreground gray

    set ::PETK::gui::centeredAminoCombo $container.source.centeredcombo

    ttk::label $container.source.note -text "Browse loads analyte with its original coordinates for review." \
        -font {TkDefaultFont 9 italic} -foreground gray

    grid $container.source.pdblbl $container.source.pdb $container.source.browse -sticky ew -pady 3
    grid $container.source.sellbl $container.source.sel $container.source.selhelp -sticky ew -pady 3
    grid $container.source.aflbl $container.source.afid $container.source.afbtn $container.source.afhelp -sticky ew -pady 3
    grid $container.source.centeredlbl $container.source.centeredcombo $container.source.centeredload $container.source.centeredrefresh -sticky ew -pady 3
    grid $container.source.centeredhelp -column 0 -columnspan 4 -row 4 -sticky w -pady "0 5"
    grid $container.source.note -column 0 -columnspan 4 -row 5 -sticky w -pady "5 0"

    ::PETK::gui::refreshCenteredAminoOptions

    # === Phase 1: Centering Frame ===
    ttk::labelframe $container.centering -text "Phase 1: Centering" -padding 10
    grid $container.centering -column 0 -row $row -sticky ew -pady "0 10" -padx 10
    grid columnconfigure $container.centering 0 -weight 1
    incr row

    ttk::label $container.centering.desc -text "Review measured center of mass (Å) and optionally override the target values before centering." \
        -wraplength 520
    grid $container.centering.desc -column 0 -row 0 -sticky w

    ttk::frame $container.centering.comtable
    grid $container.centering.comtable -column 0 -row 1 -sticky ew -pady "8 5"
    grid columnconfigure $container.centering.comtable {1 2} -weight 1

    ttk::label $container.centering.comtable.haxis -text "Axis" -font {TkDefaultFont 9 bold}
    ttk::label $container.centering.comtable.hcurrent -text "Current COM (Å)" -font {TkDefaultFont 9 bold}
    ttk::label $container.centering.comtable.htarget -text "Target COM (Å)" -font {TkDefaultFont 9 bold}
    grid $container.centering.comtable.haxis $container.centering.comtable.hcurrent $container.centering.comtable.htarget -sticky ew -pady 2

    ttk::label $container.centering.comtable.xaxis -text "X"
    ttk::label $container.centering.comtable.xcurrent -textvariable ::PETK::gui::currentComX -width 18 -relief sunken -background white
    ttk::entry $container.centering.comtable.xtarget -textvariable ::PETK::gui::targetComX -width 12 -justify center -state disabled
    grid $container.centering.comtable.xaxis $container.centering.comtable.xcurrent $container.centering.comtable.xtarget -sticky ew -pady 2

    ttk::label $container.centering.comtable.yaxis -text "Y"
    ttk::label $container.centering.comtable.ycurrent -textvariable ::PETK::gui::currentComY -width 18 -relief sunken -background white
    ttk::entry $container.centering.comtable.ytarget -textvariable ::PETK::gui::targetComY -width 12 -justify center -state disabled
    grid $container.centering.comtable.yaxis $container.centering.comtable.ycurrent $container.centering.comtable.ytarget -sticky ew -pady 2

    ttk::label $container.centering.comtable.zaxis -text "Z"
    ttk::label $container.centering.comtable.zcurrent -textvariable ::PETK::gui::currentComZ -width 18 -relief sunken -background white
    ttk::entry $container.centering.comtable.ztarget -textvariable ::PETK::gui::targetComZ -width 12 -justify center -state disabled
    grid $container.centering.comtable.zaxis $container.centering.comtable.zcurrent $container.centering.comtable.ztarget -sticky ew -pady 2

    set ::PETK::gui::manualCenterEntryWidgets [list \
        $container.centering.comtable.xtarget \
        $container.centering.comtable.ytarget \
        $container.centering.comtable.ztarget]

    ttk::checkbutton $container.centering.manualchk -text "Enable manual centering targets" \
        -variable ::PETK::gui::manualCenteringEnabled -command ::PETK::gui::updateManualCenteringState
    grid $container.centering.manualchk -column 0 -row 2 -sticky w -pady "5 0"

    ttk::label $container.centering.manualhint -text "When enabled, the target COM values are used instead of auto-centering to (0,0,0)." \
        -font {TkDefaultFont 9 italic} -foreground gray
    grid $container.centering.manualhint -column 0 -row 3 -sticky w

    ttk::frame $container.centering.actions
    grid $container.centering.actions -column 0 -row 4 -sticky ew -pady "8 0"
    grid columnconfigure $container.centering.actions 1 -weight 1

    ttk::button $container.centering.actions.analyze -text "1. Center Molecule" \
        -command ::PETK::gui::analyzeAnalyte -style "Accent.TButton"
    ttk::label $container.centering.actions.status -textvariable ::PETK::gui::centeringStatus -foreground darkgreen
    grid $container.centering.actions.analyze -column 0 -row 0 -sticky w
    grid $container.centering.actions.status -column 1 -row 0 -sticky w -padx "10 0"

    if {![info exists ::PETK::gui::alphafoldUniProt]} {
        set ::PETK::gui::alphafoldUniProt ""
    }

    # === Phase 2: Axis Selection Frame ===
    ttk::labelframe $container.alignment -text "Phase 2: Principal Axis Alignment" -padding 10
    grid $container.alignment -column 0 -row $row -sticky ew -pady "0 10" -padx 10
    grid columnconfigure $container.alignment {0 1 2 3} -weight 1
    incr row

    # Principal axes information display
    ttk::frame $container.alignment.info
    grid $container.alignment.info -column 0 -row 0 -columnspan 4 -sticky ew -pady "0 10"
    grid columnconfigure $container.alignment.info {1 2 3} -weight 1

    ttk::label $container.alignment.info.title -text "Principal Axes Information:" -font {TkDefaultFont 10 bold}
    grid $container.alignment.info.title -column 0 -row 0 -columnspan 4 -sticky w -pady "0 5"

    # Headers
    ttk::label $container.alignment.info.hdr1 -text "Axis" -font {TkDefaultFont 9 bold} -width 8
    ttk::label $container.alignment.info.hdr2 -text "Direction Vector" -font {TkDefaultFont 9 bold} -width 25
    ttk::label $container.alignment.info.hdr3 -text "Moment" -font {TkDefaultFont 9 bold} -width 12
    grid $container.alignment.info.hdr1 $container.alignment.info.hdr2 $container.alignment.info.hdr3 -sticky ew -pady 2

    # Axis 1 info
    ttk::label $container.alignment.info.ax1 -text "Axis 1" -width 8
    ttk::label $container.alignment.info.vec1 -textvariable ::PETK::gui::axis1Vector -width 25 -relief sunken -background white -font {TkFixedFont 9}
    ttk::label $container.alignment.info.mom1 -textvariable ::PETK::gui::axis1Moment -width 12 -relief sunken -background white
    grid $container.alignment.info.ax1 $container.alignment.info.vec1 $container.alignment.info.mom1 -sticky ew -pady 1

    # Axis 2 info
    ttk::label $container.alignment.info.ax2 -text "Axis 2" -width 8
    ttk::label $container.alignment.info.vec2 -textvariable ::PETK::gui::axis2Vector -width 25 -relief sunken -background white -font {TkFixedFont 9}
    ttk::label $container.alignment.info.mom2 -textvariable ::PETK::gui::axis2Moment -width 12 -relief sunken -background white
    grid $container.alignment.info.ax2 $container.alignment.info.vec2 $container.alignment.info.mom2 -sticky ew -pady 1

    # Axis 3 info
    ttk::label $container.alignment.info.ax3 -text "Axis 3" -width 8
    ttk::label $container.alignment.info.vec3 -textvariable ::PETK::gui::axis3Vector -width 25 -relief sunken -background white -font {TkFixedFont 9}
    ttk::label $container.alignment.info.mom3 -textvariable ::PETK::gui::axis3Moment -width 12 -relief sunken -background white
    grid $container.alignment.info.ax3 $container.alignment.info.vec3 $container.alignment.info.mom3 -sticky ew -pady 1

    # Axis selection controls with reset capability
    ttk::frame $container.alignment.controls
    grid $container.alignment.controls -column 0 -row 1 -columnspan 4 -sticky ew -pady "10 0"
    grid columnconfigure $container.alignment.controls 0 -weight 1

    # Selection label and radio buttons on first row
    ttk::label $container.alignment.controls.sellbl -text "Select axis to align to Z (axis with smallest moment is automated anligned with z axis):" -font {TkDefaultFont 10 bold}
    grid $container.alignment.controls.sellbl -column 0 -row 0 -sticky w -pady "0 5"

    # Radio buttons in a sub-frame for better control
    ttk::frame $container.alignment.controls.radioframe
    grid $container.alignment.controls.radioframe -column 0 -row 1 -sticky ew -pady "0 5"
    grid columnconfigure $container.alignment.controls.radioframe {0 1 2} -weight 1

    # Radio buttons for axis selection
    ttk::radiobutton $container.alignment.controls.radioframe.axis1 -text "Axis 1" -variable ::PETK::gui::selectedAxis -value 1
    ttk::radiobutton $container.alignment.controls.radioframe.axis2 -text "Axis 2" -variable ::PETK::gui::selectedAxis -value 2
    ttk::radiobutton $container.alignment.controls.radioframe.axis3 -text "Axis 3" -variable ::PETK::gui::selectedAxis -value 3

    grid $container.alignment.controls.radioframe.axis1 -column 0 -row 0 -sticky w -padx "0 10"
    grid $container.alignment.controls.radioframe.axis2 -column 1 -row 0 -sticky w -padx "0 10"
    grid $container.alignment.controls.radioframe.axis3 -column 2 -row 0 -sticky w -padx "0 10"

    # Button frame for align and reset buttons
    ttk::frame $container.alignment.controls.buttonframe
    grid $container.alignment.controls.buttonframe -column 0 -row 2 -sticky ew -pady "5 0"
    grid columnconfigure $container.alignment.controls.buttonframe {0 1} -weight 1

    # Align button
    ttk::button $container.alignment.controls.buttonframe.align -text "2. Align to Z-axis" -command ::PETK::gui::alignToZAxis -style "Accent.TButton" -state disabled
    grid $container.alignment.controls.buttonframe.align -column 0 -row 0 -sticky w -padx "0 5"

    # Reset button (initially disabled)
    ttk::button $container.alignment.controls.buttonframe.reset -text "Reset to Centered" -command ::PETK::gui::resetToCentered -state disabled
    grid $container.alignment.controls.buttonframe.reset -column 1 -row 0 -sticky w -padx "5 0"

    # Store references to alignment controls for enabling/disabling
    set ::PETK::gui::alignControls [list \
        $container.alignment.controls.radioframe.axis1 \
        $container.alignment.controls.radioframe.axis2 \
        $container.alignment.controls.radioframe.axis3 \
        $container.alignment.controls.buttonframe.align]

    # Store reset button reference
    set ::PETK::gui::resetButton $container.alignment.controls.buttonframe.reset

    # Add status indicator
    ttk::label $container.alignment.controls.status -text "Status: Complete Phase 1 to enable alignment" \
        -font {TkDefaultFont 9 italic} -foreground gray
    grid $container.alignment.controls.status -column 0 -row 3 -sticky w -pady "5 0"
    set ::PETK::gui::alignmentStatusLabel $container.alignment.controls.status

    # Add recommendation text on bottom row  
    ttk::label $container.alignment.controls.recommend -text "💡 Tip: You can try different alignments - use Reset to return to centered state" \
        -font {TkDefaultFont 9 italic} -foreground darkblue
    grid $container.alignment.controls.recommend -column 0 -row 4 -sticky w -pady "2 0"

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

    # Right column - Processing status
    ttk::frame $container.results.right
    grid $container.results.right -column 1 -row 0 -sticky nsew -padx "10 0"
    grid rowconfigure $container.results.right 6 -weight 1

    ttk::label $container.results.right.title -text "Processing Status" -font {TkDefaultFont 11 bold}
    grid $container.results.right.title -sticky w -pady "0 10"

    # Phase status
    ttk::label $container.results.right.phaselbl -text "Current phase:" -width 15
    ttk::label $container.results.right.phase -textvariable ::PETK::gui::currentPhase -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.phaselbl $container.results.right.phase -sticky ew -pady 2

    # Verification score
    ttk::label $container.results.right.scorelbl -text "Quality score:" -width 15
    ttk::label $container.results.right.score -textvariable ::PETK::gui::verificationScore -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.scorelbl $container.results.right.score -sticky ew -pady 2

    # Centering status
    ttk::label $container.results.right.centerlbl -text "Centering:" -width 15
    ttk::label $container.results.right.center -textvariable ::PETK::gui::centeringStatus -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.centerlbl $container.results.right.center -sticky ew -pady 2

    # Alignment status
    ttk::label $container.results.right.alignlbl -text "Alignment:" -width 15
    ttk::label $container.results.right.align -textvariable ::PETK::gui::alignmentStatus -width 15 -anchor w -relief sunken -background white
    grid $container.results.right.alignlbl $container.results.right.align -sticky ew -pady 2

    # Output file status
    ttk::label $container.results.right.outputlbl -text "Output file:" -width 15
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
    text $container.details.content.text -height 10 -width 80 -wrap word -state disabled \
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
    grid columnconfigure $container.controls {0 1 2 3 4} -weight 1
    incr row

    ttk::button $container.controls.show -text "Show Molecule" -command ::PETK::gui::showAnalyte
    ttk::button $container.controls.hide -text "Hide Molecule" -command ::PETK::gui::hideAnalyte
    ttk::button $container.controls.center -text "Center View" -command ::PETK::gui::centerAnalyteView
    ttk::button $container.controls.representations -text "Change Rep" -command ::PETK::gui::cycleAnalyteRepresentation
    ttk::button $container.controls.refresh -text "Refresh Analyte" -command ::PETK::gui::refreshAnalyte
    #ttk::button $container.controls.save -text "Save Centered PDB" -command ::PETK::gui::saveCenteredPDB -state disabled

    grid $container.controls.show $container.controls.hide $container.controls.center $container.controls.representations $container.controls.refresh -sticky ew -padx 3

    # Store reference to save button
    #set ::PETK::gui::saveButton $container.controls.save

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab2.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab2.canvas $canvas_window $container]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab2.canvas $canvas_window $container]
    
    # Store references for later use (with unique names for tab2)
    set ::PETK::gui::tab2MainCanvas $tab2.canvas
    set ::PETK::gui::tab2CanvasWindow $canvas_window
    set ::PETK::gui::tab2ContentContainer $container
    
    # Initialize result variables
    ::PETK::gui::initializeResultVariables
    ::PETK::gui::updateManualCenteringState

    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab2.canvas $canvas_window $container]
}

####################################################
# Tab 2 Function
####################################################

proc ::PETK::gui::initializeResultVariables {} {
    # Basic results
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
    if {![info exists ::PETK::gui::leftmostAtom]} {
        set ::PETK::gui::leftmostAtom ""
    }
    if {![info exists ::PETK::gui::rightmostAtom]} {
        set ::PETK::gui::rightmostAtom ""
    }
    
    # Status variables
    if {![info exists ::PETK::gui::currentPhase]} {
        set ::PETK::gui::currentPhase "Ready"
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

    # Principal axes display variables
    if {![info exists ::PETK::gui::axis1Vector]} {
        set ::PETK::gui::axis1Vector "Not calculated"
    }
    if {![info exists ::PETK::gui::axis2Vector]} {
        set ::PETK::gui::axis2Vector "Not calculated"
    }
    if {![info exists ::PETK::gui::axis3Vector]} {
        set ::PETK::gui::axis3Vector "Not calculated"
    }
    if {![info exists ::PETK::gui::axis1Moment]} {
        set ::PETK::gui::axis1Moment ""
    }
    if {![info exists ::PETK::gui::axis2Moment]} {
        set ::PETK::gui::axis2Moment ""
    }
    if {![info exists ::PETK::gui::axis3Moment]} {
        set ::PETK::gui::axis3Moment ""
    }
    

    # Selection variables 
    if {![info exists ::PETK::gui::selectedAxis]} {
        set ::PETK::gui::selectedAxis 3
    }  

    if {![info exists ::PETK::gui::analyteSelection]} {
        set ::PETK::gui::analyteSelection "all"
    }

    # Center-of-mass display variables
    if {![info exists ::PETK::gui::currentComX]} {
        set ::PETK::gui::currentComX "N/A"
    }
    if {![info exists ::PETK::gui::currentComY]} {
        set ::PETK::gui::currentComY "N/A"
    }
    if {![info exists ::PETK::gui::currentComZ]} {
        set ::PETK::gui::currentComZ "N/A"
    }

    if {![info exists ::PETK::gui::targetComX]} {
        set ::PETK::gui::targetComX 0.0
    }
    if {![info exists ::PETK::gui::targetComY]} {
        set ::PETK::gui::targetComY 0.0
    }
    if {![info exists ::PETK::gui::targetComZ]} {
        set ::PETK::gui::targetComZ 0.0
    }

    if {![info exists ::PETK::gui::manualCenteringEnabled]} {
        set ::PETK::gui::manualCenteringEnabled 0
    }
}

proc ::PETK::gui::updateManualCenteringState {} {
    if {![info exists ::PETK::gui::manualCenteringEnabled]} {
        set ::PETK::gui::manualCenteringEnabled 0
    }

    if {$::PETK::gui::manualCenteringEnabled} {
        set state normal
    } else {
        set state disabled
    }

    if {[info exists ::PETK::gui::manualCenterEntryWidgets]} {
        foreach widget $::PETK::gui::manualCenterEntryWidgets {
            if {[winfo exists $widget]} {
                $widget configure -state $state
            }
        }
    }
}

proc ::PETK::gui::resetManualCenterTargets {} {
    foreach axis {X Y Z} {
        set varname "::PETK::gui::targetCom$axis"
        set $varname 0.0
    }
}

proc ::PETK::gui::updateCurrentComDisplay {center_values} {
    if {[llength $center_values] != 3} {
        return
    }

    foreach axis {X Y Z} idx {0 1 2} {
        set varname "::PETK::gui::currentCom$axis"
        set value [lindex $center_values $idx]
        set $varname [format "%.3f" $value]
    }
}

proc ::PETK::gui::updateMovementRangeFromSelection {sel} {
    if {$sel eq ""} {
        return
    }
    if {[catch {$sel num} num_atoms] || $num_atoms <= 0} {
        return
    }
    if {[catch {measure minmax $sel} minmax]} {
        return
    }
    set min_coord [lindex $minmax 0]
    set max_coord [lindex $minmax 1]
    set min_z [lindex $min_coord 2]
    set max_z [lindex $max_coord 2]

    if {[info procs ::PETK::gui::updateMovementRangeFromBounds] ne ""} {
        ::PETK::gui::updateMovementRangeFromBounds $min_z $max_z
    } else {
        set start [format "%.1f" [expr {max($min_z,$max_z)}]]
        set end [format "%.1f" [expr {min($min_z,$max_z)}]]
        if {[info procs ::PETK::gui::formatMovementValue] ne ""} {
            set ::PETK::gui::zStartRange [::PETK::gui::formatMovementValue $start]
            set ::PETK::gui::zEndRange [::PETK::gui::formatMovementValue $end]
        } else {
            set ::PETK::gui::zStartRange $start
            set ::PETK::gui::zEndRange $end
        }
    }

    puts [format "Analyte movement range updated from geometry: %s → %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
}

proc ::PETK::gui::captureViewState {} {
    if {[info commands display] eq ""} {
        return ""
    }
    if {[catch {display get view} state]} {
        return ""
    }
    return $state
}

proc ::PETK::gui::restoreViewState {state} {
    if {$state eq ""} {
        return
    }
    if {[info commands display] eq ""} {
        return
    }
    after 50 [list ::PETK::gui::applyViewState $state]
}

proc ::PETK::gui::applyViewState {state} {
    if {[info commands display] eq ""} {
        return
    }
    catch {display set view $state}
}

proc ::PETK::gui::enableAlignmentControls {} {
    # Enable the alignment controls after Phase 1 is complete
    if {[info exists ::PETK::gui::alignControls]} {
        foreach widget $::PETK::gui::alignControls {
            if {[winfo exists $widget]} {
                $widget configure -state normal
            }
        }
    }
    
    # Update status
    set ::PETK::gui::currentPhase "Phase 1 complete"
    
    # Update status label
    if {[info exists ::PETK::gui::alignmentStatusLabel] && [winfo exists $::PETK::gui::alignmentStatusLabel]} {
        $::PETK::gui::alignmentStatusLabel configure -text "Status: Ready for alignment - select an axis and click Align" -foreground darkgreen
    }
    
    # Update principal axes display (previous code...)
    if {[info exists ::PETK::gui::principalAxes] && [info exists ::PETK::gui::inertiaMoments]} {
        set axes $::PETK::gui::principalAxes
        set moments $::PETK::gui::inertiaMoments
        
        # Format axis vectors and moments (previous formatting code)
        for {set i 0} {$i < 3} {incr i} {
            set axis [lindex $axes $i]
            set moment [lindex $moments $i]
            
            set var_name_vec "::PETK::gui::axis[expr {$i+1}]Vector"
            set var_name_mom "::PETK::gui::axis[expr {$i+1}]Moment"
            
            set $var_name_vec [format "(%.3f, %.3f, %.3f)" \
                [lindex $axis 0] [lindex $axis 1] [lindex $axis 2]]
            
            if {[string is double $moment]} {
                if {$moment < 0.001} {
                    set $var_name_mom [format "%.1e" $moment]
                } elseif {$moment < 1.0} {
                    set $var_name_mom [format "%.4f" $moment]
                } else {
                    set $var_name_mom [format "%.2f" $moment]
                }
            } else {
                set $var_name_mom "Invalid"
            }
        }
    }
}

# Enhanced alignment completion function
proc ::PETK::gui::onAlignmentComplete {} {
    # Enable reset button after alignment
    if {[info exists ::PETK::gui::resetButton] && [winfo exists $::PETK::gui::resetButton]} {
        $::PETK::gui::resetButton configure -state normal
    }
    
    # Update status
    if {[info exists ::PETK::gui::alignmentStatusLabel] && [winfo exists $::PETK::gui::alignmentStatusLabel]} {
        $::PETK::gui::alignmentStatusLabel configure -text "Status: Aligned - use Reset to try different alignment" -foreground blue
    }
    
    # Enable save button
    ::PETK::gui::enableSaveButton
    ::PETK::gui::refreshPorePreviewIfNeeded
}

# Enhanced reset completion function
proc ::PETK::gui::onResetComplete {} {
    # Update status
    if {[info exists ::PETK::gui::alignmentStatusLabel] && [winfo exists $::PETK::gui::alignmentStatusLabel]} {
        $::PETK::gui::alignmentStatusLabel configure -text "Status: Reset to centered - ready for new alignment" -foreground darkgreen
    }
    
    # Disable save button
    ::PETK::gui::disableSaveButton
    ::PETK::gui::refreshPorePreviewIfNeeded
}

proc ::PETK::gui::enableSaveButton {} {
    if {[info exists ::PETK::gui::saveButton] && [winfo exists $::PETK::gui::saveButton]} {
        $::PETK::gui::saveButton configure -state normal
    }
}

proc ::PETK::gui::disableSaveButton {} {
    if {[info exists ::PETK::gui::saveButton] && [winfo exists $::PETK::gui::saveButton]} {
        $::PETK::gui::saveButton configure -state disabled
    }
}

proc ::PETK::gui::testInertiaCalculation {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        puts "No molecule loaded for testing"
        return
    }
    
    set molid $::PETK::gui::analyteMol
    
    # Test with different selections
    set test_selections {"all" "name CA" "protein" "backbone"}
    
    foreach sel_text $test_selections {
        puts "\n=== Testing selection: $sel_text ==="
        
        if {[catch {atomselect $molid $sel_text} sel]} {
            puts "ERROR: Could not create selection '$sel_text': $sel"
            continue
        }
        
        set num_atoms [$sel num]
        puts "Number of atoms: $num_atoms"
        
        if {$num_atoms == 0} {
            puts "WARNING: Selection contains 0 atoms"
            $sel delete
            continue
        }
        
        # Test center of mass
        if {[catch {measure center $sel} com]} {
            puts "ERROR: Could not calculate center of mass: $com"
            $sel delete
            continue
        }
        puts "Center of mass: $com"
        
        # Test inertia calculation
        if {[catch {measure inertia $sel} inertia_basic]} {
            puts "ERROR: Basic inertia calculation failed: $inertia_basic"
            $sel delete
            continue
        }
        puts "Basic inertia result length: [llength $inertia_basic]"
        
        # Test inertia with moments
        if {[catch {measure inertia $sel moments} inertia_moments]} {
            puts "ERROR: Inertia moments calculation failed: $inertia_moments"
            $sel delete
            continue
        }
        
        puts "Inertia moments result:"
        puts "  Full result: $inertia_moments"
        puts "  Result length: [llength $inertia_moments]"
        
        if {[llength $inertia_moments] == 3} {
            set moments [lindex $inertia_moments 0]
            set axes [lindex $inertia_moments 1]
            set tensor [lindex $inertia_moments 2]
            
            puts "  Moments: $moments (length: [llength $moments])"
            puts "  Axes: $axes (length: [llength $axes])"
            
            # Check individual moment values
            for {set i 0} {$i < [llength $moments]} {incr i} {
                set moment [lindex $moments $i]
                puts "    Moment $i: $moment (valid number: [string is double $moment])"
            }
        }
        
        $sel delete
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
        ::PETK::gui::loadAnalyteOriginal
    }
}

proc ::PETK::gui::friendlyAminoName {code} {
    set key [string toupper $code]
    if {[info exists ::PETK::gui::aminoAcidFriendlyNames] && [dict exists $::PETK::gui::aminoAcidFriendlyNames $key]} {
        return [dict get $::PETK::gui::aminoAcidFriendlyNames $key]
    }
    return $key
}

proc ::PETK::gui::refreshCenteredAminoOptions {} {
    set dir [::PETK::gui::resourcePath analytes centered_amino_acids]
    set values {}
    set lookup [dict create]

    if {[file isdirectory $dir]} {
        set files [glob -nocomplain -directory $dir "centered_*.pdb"]
        foreach file [lsort -dictionary $files] {
            set base [file tail $file]
            if {![regexp {centered_([A-Z]{3})\.pdb} $base -> rescode]} {
                continue
            }
            set friendly [::PETK::gui::friendlyAminoName $rescode]
            set label [format "%s - %s" $rescode $friendly]
            dict set lookup $label $file
            lappend values $label
        }
    } else {
        puts "Centered amino acid library not found: $dir"
    }

    set ::PETK::gui::centeredAminoMap $lookup

    if {![info exists ::PETK::gui::selectedCenteredAmino]} {
        set ::PETK::gui::selectedCenteredAmino ""
    }

    if {[llength $values] == 0} {
        set ::PETK::gui::selectedCenteredAmino ""
    } elseif {$::PETK::gui::selectedCenteredAmino eq "" || ![dict exists $lookup $::PETK::gui::selectedCenteredAmino]} {
        set ::PETK::gui::selectedCenteredAmino [lindex $values 0]
    }

    if {[info exists ::PETK::gui::centeredAminoCombo] && [winfo exists $::PETK::gui::centeredAminoCombo]} {
        $::PETK::gui::centeredAminoCombo configure -values $values
    }

    return $values
}

proc ::PETK::gui::loadSelectedCenteredAmino {} {
    if {![info exists ::PETK::gui::selectedCenteredAmino] || $::PETK::gui::selectedCenteredAmino eq ""} {
        tk_messageBox -icon info -title "Centered Amino Acids" -message "Select a centered amino acid first."
        return
    }

    if {![info exists ::PETK::gui::centeredAminoMap] || ![dict exists $::PETK::gui::centeredAminoMap $::PETK::gui::selectedCenteredAmino]} {
        tk_messageBox -icon warning -title "Centered Amino Acids" \
            -message "The selected amino acid entry is unavailable. Refresh the list and try again."
        return
    }

    set path [dict get $::PETK::gui::centeredAminoMap $::PETK::gui::selectedCenteredAmino]

    if {![file exists $path]} {
        tk_messageBox -icon error -title "Centered Amino Acids" \
            -message "File not found for the selected amino acid." -detail $path
        ::PETK::gui::refreshCenteredAminoOptions
        return
    }

    set ::PETK::gui::analytePDB $path
    # Library amino acids are pre-centered, so the SEM-ready PDB is the same
    # file. Sync semAnalytePDB so tab 3 (parameter summary, validation, run
    # command) sees a populated path without requiring the centering workflow.
    set ::PETK::gui::semAnalytePDB $path
    if {![info exists ::PETK::gui::analyteSelection] || $::PETK::gui::analyteSelection eq ""} {
        set ::PETK::gui::analyteSelection "all"
    }

    puts "Loading centered amino acid sample: [file tail $path]"
    ::PETK::gui::loadAnalyteOriginal
}

proc ::PETK::gui::refreshAnalyte {} {
    if {![info exists ::PETK::gui::analytePDB] || $::PETK::gui::analytePDB eq ""} {
        tk_messageBox -icon info -title "Refresh Analyte" -message "Select an analyte PDB first."
        return
    }
    ::PETK::gui::loadAnalyteOriginal
}

proc ::PETK::gui::loadAnalyteOriginal {} {
    set viewState [::PETK::gui::captureViewState]
    set code [catch {
        set redrawPorePreview 0
        if {[info exists ::PETK::gui::petk_draw_molids] && [llength $::PETK::gui::petk_draw_molids] > 0} {
            set redrawPorePreview 1
        }

        if {![info exists ::PETK::gui::analytePDB] || $::PETK::gui::analytePDB eq ""} {
            tk_messageBox -icon warning -title "Load Analyte" -message "No PDB file selected."
            return
        }

        if {![file exists $::PETK::gui::analytePDB]} {
            tk_messageBox -icon error -title "Load Analyte" -message "PDB file not found:" -detail $::PETK::gui::analytePDB
            return
        }

        if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol ne ""} {
            if {[info exists ::PETK::gui::petk_draw_molids]} {
                set idx [lsearch -exact $::PETK::gui::petk_draw_molids $::PETK::gui::analyteMol]
                if {$idx != -1} {
                    set ::PETK::gui::petk_draw_molids [lreplace $::PETK::gui::petk_draw_molids $idx $idx]
                }
            }
            catch {mol delete $::PETK::gui::analyteMol}
        }

        if {[catch {mol new $::PETK::gui::analytePDB waitfor all} molid]} {
            tk_messageBox -icon error -title "Load Analyte" -message "Could not load PDB file" -detail $molid
            return
        }

        set ::PETK::gui::analyteMol $molid
        mol rename $molid "Analyte: [file tail $::PETK::gui::analytePDB]"

        if {![info exists ::PETK::gui::analyteSelection]} {
            set ::PETK::gui::analyteSelection "all"
        }

        if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
            puts "Analyte load: selection '$::PETK::gui::analyteSelection' failed ($sel), falling back to 'all'"
            if {[catch {atomselect $molid "all"} sel]} {
                tk_messageBox -icon error -title "Load Analyte" -message "Could not create atom selection" -detail $sel
                return
            }
        }

        if {[$sel num] == 0} {
            $sel delete
            set ::PETK::gui::currentComX "N/A"
            set ::PETK::gui::currentComY "N/A"
            set ::PETK::gui::currentComZ "N/A"
            tk_messageBox -icon warning -title "Load Analyte" -message "Selection contains 0 atoms." -detail "Adjust the atom selection text."
        } else {
        set com [measure center $sel]
        ::PETK::gui::updateCurrentComDisplay $com
        ::PETK::gui::updateMovementRangeFromSelection $sel
        $sel delete
        }

        ::PETK::gui::clearAnalysisResults
        set ::PETK::gui::centeringStatus "Original coordinates"
        set ::PETK::gui::alignmentStatus "Not aligned"
        set ::PETK::gui::currentPhase "Loaded (original coords)"

        catch {unset ::PETK::gui::centeredCoords}
        ::PETK::gui::resetManualCenterTargets
        set ::PETK::gui::manualCenteringEnabled 0
        ::PETK::gui::updateManualCenteringState

        if {[info exists ::PETK::gui::alignControls]} {
            foreach widget $::PETK::gui::alignControls {
                if {[winfo exists $widget]} {
                    $widget configure -state disabled
                }
            }
        }

        if {[info exists ::PETK::gui::resetButton] && [winfo exists $::PETK::gui::resetButton]} {
            $::PETK::gui::resetButton configure -state disabled
        }

        if {[info exists ::PETK::gui::alignmentStatusLabel] && [winfo exists $::PETK::gui::alignmentStatusLabel]} {
            $::PETK::gui::alignmentStatusLabel configure -text "Status: Complete Phase 1 to enable alignment" -foreground gray
        }

        ::PETK::gui::showAnalyte
        if {$redrawPorePreview} {
            ::PETK::gui::refreshPorePreviewIfNeeded
        }
        puts "Loaded analyte with original coordinates: $::PETK::gui::analytePDB"
    } result options]
    ::PETK::gui::restoreViewState $viewState
    if {$code} {
        return -code $code -options $options $result
    }
    return $result
}

proc ::PETK::gui::refreshPorePreviewIfNeeded {} {
    if {![info exists ::PETK::gui::petk_draw_molids] || [llength $::PETK::gui::petk_draw_molids] == 0} {
        return
    }
    if {[info procs ::PETK::gui::updatePoreVisualization] eq ""} {
        return
    }
    after idle {catch {::PETK::gui::updatePoreVisualization}}
}

proc ::PETK::gui::downloadAlphaFoldAnalyte {} {
    if {![info exists ::PETK::gui::alphafoldUniProt]} {
        set ::PETK::gui::alphafoldUniProt ""
    }
    set raw_id [string trim $::PETK::gui::alphafoldUniProt]
    if {$raw_id eq ""} {
        tk_messageBox -icon warning -title "AlphaFold Download"             -message "Please enter a UniProt accession (e.g., P69905)."
        return
    }

    set id_upper [string toupper $raw_id]

    if {[regexp {^[A-Z0-9]{4}$} $id_upper]} {
        set ::PETK::gui::alphafoldUniProt $id_upper
        return [::PETK::gui::downloadRCSBAnalyte $id_upper]
    }

    if {![regexp {^[A-Z0-9]+(?:-[A-Z0-9]+)?$} $id_upper]} {
        tk_messageBox -icon error -title "Analyte Download" \
            -message "Unrecognized ID: $raw_id\nEnter a UniProt accession (e.g., P69905) or a four-character PDB code (e.g., 1XYZ)."
        return
    }

    set uniprot $id_upper

    if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
        set dest_dir $::PETK::gui::workdir
    } else {
        set dest_dir [pwd]
    }
    if {![file isdirectory $dest_dir]} {
        catch {file mkdir $dest_dir}
    }

    # Prepare TLS search path for Homebrew installs
    set tls_paths {}
    if {[catch {exec brew --prefix tcl-tls} brew_prefix] == 0} {
        set brew_dir [string trim $brew_prefix]
        if {$brew_dir ne ""} {
            lappend tls_paths [file join $brew_dir "lib"]
            foreach sub [glob -nocomplain [file join $brew_dir "../Cellar/tcl-tls" * lib]] {
                lappend tls_paths $sub
            }
        }
    }
    foreach path {/opt/homebrew/opt/tcl-tls/lib /usr/local/opt/tcl-tls/lib} {
        if {[file isdirectory $path]} {
            lappend tls_paths $path
        }
    }
    foreach path $tls_paths {
        if {[file isdirectory $path] && [lsearch $auto_path $path] < 0} {
            lappend auto_path $path
        }
    }

    set packages_available 0
    set tls_err ""
    set http_err ""
    if {[catch {package require tls} tls_err] == 0 && [catch {package require http} http_err] == 0} {
        set packages_available 1
    }

    set versions {v10 v9 v8 v7 v6 v5 v4 v3 v2 v1 ""}
    set formats {{".pdb.gz" 1} {".pdb" 0}}
    set attempts {}
    set download_ok 0
    set outfile ""
    set final_url ""
    set needs_gunzip 0
    set final_error ""

    foreach ver $versions {
        if {$ver eq ""} {
            set ver_suffix ""
        } else {
            set ver_suffix "_$ver"
        }
        foreach fmt $formats {
            set ext [lindex $fmt 0]
            set gz_flag [lindex $fmt 1]
            set filename_base [format "AF-%s-F1-model%s" $uniprot $ver_suffix]
            set filename "${filename_base}${ext}"
            set url [format "https://alphafold.ebi.ac.uk/files/%s" $filename]
            lappend attempts $url

            set outfile_raw [file join $dest_dir $filename]
            catch {file delete -force $outfile_raw}

            set success 0
            set status_msg ""

            if {$packages_available} {
                ::http::register https 443 ::tls::socket
                if {[catch {set token [::http::geturl $url -binary 1 -timeout 60000]} err_http]} {
                    set status_msg $err_http
                } else {
                    set code [::http::ncode $token]
                    if {[::http::status $token] eq "ok" && $code == 200} {
                        set data [::http::data $token]
                        if {[string length $data] > 0} {
                            if {[catch {
                                set fh [open $outfile_raw w]
                                fconfigure $fh -translation binary -encoding binary
                                puts -nonewline $fh $data
                                close $fh
                            } err_write]} {
                                set status_msg $err_write
                            } else {
                                set success 1
                            }
                        } else {
                            set status_msg "Empty response"
                        }
                    } else {
                        set status_msg "HTTP $code"
                    }
                    ::http::cleanup $token
                }
                ::http::unregister https
            }

            if {!$success} {
                set curl_cmd [list curl -L -f -s -o $outfile_raw $url]
                puts "AlphaFold: attempting curl download: [join $curl_cmd { }]"
                if {[catch {eval exec $curl_cmd} curl_err]} {
                    set status_msg $curl_err
                } else {
                    set success 1
                }
            }

            if {$success} {
                set download_ok 1
                set outfile $outfile_raw
                set final_url $url
                set needs_gunzip $gz_flag
                set final_error ""
                break
            } else {
                set final_error $status_msg
            }
        }
        if {$download_ok} {
            break
        }
    }

    if {!$download_ok} {
        set msg "Failed to download AlphaFold model for $uniprot.\nTried URLs:\n[join $attempts "\n"]\n\nLast error: $final_error"
        if {!$packages_available} {
            append msg "\n\nTLS/http packages unavailable (tls: $tls_err, http: $http_err)."
        }
        tk_messageBox -icon error -title "AlphaFold Download" -message $msg
        return
    }

    if {$needs_gunzip} {
        set decompressed_file [file rootname $outfile]
        if {[catch {
            set in_fh [open $outfile r]
            fconfigure $in_fh -translation binary -encoding binary
            set gzdata [read $in_fh]
            close $in_fh
            set decompressed [zlib gunzip $gzdata]
            set out_fh [open $decompressed_file w]
            fconfigure $out_fh -translation binary -encoding binary
            puts -nonewline $out_fh $decompressed
            close $out_fh
            file delete -force $outfile
            set outfile $decompressed_file
        } zlib_err]} {
            if {[catch {exec gunzip -c $outfile} gunzip_data]} {
                set msg "Downloaded AlphaFold model from:\n$final_url\nbut failed to decompress (.pdb.gz).\n\nzlib error: $zlib_err\ngunzip error: $gunzip_data"
                tk_messageBox -icon error -title "AlphaFold Download" -message $msg
                return
            } else {
                set out_fh [open $decompressed_file w]
                fconfigure $out_fh -translation binary -encoding binary
                puts -nonewline $out_fh $gunzip_data
                close $out_fh
                file delete -force $outfile
                set outfile $decompressed_file
            }
        }
    }

    set ::PETK::gui::analytePDB $outfile
    set ::PETK::gui::alphafoldUniProt $uniprot
    tk_messageBox -icon info -title "AlphaFold Download" \
        -message "Downloaded AlphaFold model for $uniprot to:\n$outfile"
    puts "AlphaFold download completed: $outfile"
    ::PETK::gui::loadAnalyteOriginal
}

proc ::PETK::gui::downloadRCSBAnalyte {pdb_id} {
    set pdb_code [string toupper $pdb_id]
    set url [format "https://files.rcsb.org/download/%s.pdb" $pdb_code]

    if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
        set dest_dir $::PETK::gui::workdir
    } else {
        set dest_dir [pwd]
    }
    if {![file isdirectory $dest_dir]} {
        catch {file mkdir $dest_dir}
    }
    set outfile [file join $dest_dir [format "%s.pdb" $pdb_code]]
    catch {file delete -force $outfile}

    set download_ok 0
    set error_msg ""

    if {[catch {package require tls}]==0 && [catch {package require http}]==0} {
        ::http::register https 443 ::tls::socket
        if {[catch {set token [::http::geturl $url -binary 1 -timeout 30000]} err]} {
            set error_msg $err
        } else {
            if {[::http::status $token] eq "ok" && [::http::ncode $token] == 200} {
                set data [::http::data $token]
                if {[string length $data] > 0} {
                    if {[catch {
                        set fh [open $outfile w]
                        fconfigure $fh -translation binary -encoding binary
                        puts -nonewline $fh $data
                        close $fh
                    } write_err]} {
                        set error_msg $write_err
                    } else {
                        set download_ok 1
                    }
                } else {
                    set error_msg "Empty response"
                }
            } else {
                set error_msg [format "HTTP %s" [::http::ncode $token]]
            }
            ::http::cleanup $token
        }
        ::http::unregister https
    }

    if {!$download_ok} {
        set curl_cmd [list curl -L -f -s -o $outfile $url]
        puts "RCSB: attempting curl download: [join $curl_cmd { }]"
        if {[catch {eval exec $curl_cmd} curl_err]} {
            set error_msg $curl_err
        } else {
            set download_ok 1
        }
    }

    if {!$download_ok} {
        set msg "Failed to download PDB entry $pdb_code from RCSB.\nURL: $url\nError: $error_msg"
        tk_messageBox -icon error -title "RCSB Download" -message $msg
        return
    }

    set ::PETK::gui::analytePDB $outfile
    tk_messageBox -icon info -title "RCSB Download" \
        -message "Downloaded PDB entry $pdb_code to:\n$outfile"
    puts "RCSB download completed: $outfile"
    ::PETK::gui::loadAnalyteOriginal
}

proc ::PETK::gui::analyzeAnalyte {} {
    puts "Starting Phase 1: Centering and principal axis analysis..."
    
    # Clear previous results
    ::PETK::gui::clearAnalysisResults
    
    # Update status
    set ::PETK::gui::verificationScore "Analyzing..."
    set ::PETK::gui::centeringStatus "Processing..."
    set ::PETK::gui::alignmentStatus "Not aligned"
    set ::PETK::gui::outputFileStatus "Phase 1 only"
    update
    
    # Phase 1: Center molecule and calculate principal axes
    if {[catch {::PETK::gui::centerAndAnalyzePrincipalAxes} error]} {
        tk_messageBox -icon error -message "Analysis Error" -detail "Error during Phase 1 analysis:\n$error"
        ::PETK::gui::clearAnalysisResults
        return
    }
    
    puts "\nPhase 1 complete. Please select alignment axis and click 'Align to Z-axis'"
    return "phase1_complete"
}

proc ::PETK::gui::centerAndAnalyzePrincipalAxes {} {
    set viewState [::PETK::gui::captureViewState]
    set code [catch {
        # Initialize variables
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

        # Load molecule
        if {[info exists ::PETK::gui::analyteMol] && $::PETK::gui::analyteMol != ""} {
            catch {mol delete $::PETK::gui::analyteMol}
        }

        if {[catch {mol new $::PETK::gui::analytePDB waitfor all} molid]} {
            tk_messageBox -icon error -message "Error!" -detail "Could not load PDB file: $::PETK::gui::analytePDB"
            return
        }

        set ::PETK::gui::analyteMol $molid
        mol rename $molid "Analyte: [file tail $::PETK::gui::analytePDB]"

        # Create selection
        if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
            tk_messageBox -icon error -message "Error!" -detail "Could not create atom selection: $sel"
            return
        }
        
        if {[$sel num] == 0} {
            tk_messageBox -icon warning -message "Warning!" -detail "Selection contains 0 atoms!"
            $sel delete
            return
        }

        puts "=================================="
        puts "PHASE 1: CENTERING AND PRINCIPAL AXIS ANALYSIS"
        puts "=================================="
        puts "Processing [$sel num] atoms with selection: $::PETK::gui::analyteSelection"

        # STEP 1: CENTER THE MOLECULE AT ORIGIN
        puts "\nSTEP 1: CENTERING MOLECULE"
        puts "-------------------------"
        set initial_center [measure center $sel]
        puts [format "Initial center of mass: (%.6f, %.6f, %.6f)" \
            [lindex $initial_center 0] [lindex $initial_center 1] [lindex $initial_center 2]]
        ::PETK::gui::updateCurrentComDisplay $initial_center

        set target_center {0.0 0.0 0.0}
        if {[info exists ::PETK::gui::manualCenteringEnabled] && $::PETK::gui::manualCenteringEnabled} {
            set idx 0
            foreach axis {X Y Z} {
                set varname "::PETK::gui::targetCom$axis"
                set raw_value [set $varname]
                if {[catch {expr {double($raw_value)}} value]} {
                    $sel delete
                    error "Invalid target COM value for axis $axis: $raw_value"
                }
                lset target_center $idx $value
                incr idx
            }
            puts [format "Manual target COM: (%.6f, %.6f, %.6f)" \
                [lindex $target_center 0] [lindex $target_center 1] [lindex $target_center 2]]
        } else {
            puts "Manual centering disabled; targeting (0, 0, 0)"
        }

        set move_vector [vecsub $target_center $initial_center]
        $sel moveby $move_vector

    set new_center [measure center $sel]
    puts [format "New center of mass: (%.6f, %.6f, %.6f)" \
        [lindex $new_center 0] [lindex $new_center 1] [lindex $new_center 2]]
    ::PETK::gui::updateCurrentComDisplay $new_center
    ::PETK::gui::updateMovementRangeFromSelection $sel
        puts "✓ Molecule centered at origin"

        set ::PETK::gui::centeredCoords [$sel get {x y z}]
        puts "✓ Centered coordinates stored for reset capability"

        # STEP 2: CALCULATE PRINCIPAL AXES
        puts "\nSTEP 2: PRINCIPAL AXES ANALYSIS"
        puts "-------------------------------"
        
        # Get inertia information
        set inertia_result [measure inertia $sel moments]
        set moments [lindex $inertia_result 0]
        set principal_axes [lindex $inertia_result 1]
        set moments_tensor [lindex $inertia_result 2]
        
        # Store principal axes for later use
        set ::PETK::gui::principalAxes $principal_axes
        set ::PETK::gui::inertiaMoments $moments
        
        # Extract individual axes
        set axis1 [lindex $principal_axes 0]
        set axis2 [lindex $principal_axes 1] 
        set axis3 [lindex $principal_axes 2]

        puts "Principal Axes (sorted by moment of inertia):"
        puts [format "  Axis 1: (%8.5f, %8.5f, %8.5f) - Moment: %10.1f (largest)" \
              [lindex $axis1 0] [lindex $axis1 1] [lindex $axis1 2] [lindex $moments 0]]
        puts [format "  Axis 2: (%8.5f, %8.5f, %8.5f) - Moment: %10.1f (medium)" \
              [lindex $axis2 0] [lindex $axis2 1] [lindex $axis2 2] [lindex $moments 1]]
        puts [format "  Axis 3: (%8.5f, %8.5f, %8.5f) - Moment: %10.1f (smallest)" \
              [lindex $axis3 0] [lindex $axis3 1] [lindex $axis3 2] [lindex $moments 2]]

        # Calculate coordinate variances for comparison
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
        
        puts "\nCurrent coordinate spreads:"
        puts [format "  X variance: %.3f" $x_var]
        puts [format "  Y variance: %.3f" $y_var]
        puts [format "  Z variance: %.3f" $z_var]

        # Suggest which axis to align for flat surface perpendicular to Z
        puts "\nAlignment recommendations:"
        puts "  For flat surface ⊥ Z-axis: Align SMALLEST moment axis to Z"
        puts "  For elongated molecule: Align LARGEST moment axis to Z"
        puts "  For custom orientation: Choose based on your preference"

        # Calculate basic geometric properties
        puts "\nSTEP 3: BASIC GEOMETRIC ANALYSIS"
        puts "--------------------------------"
        
        # Find extreme atoms
        set x_coords_only [lmap coord $coords {lindex $coord 0}]
        set leftmostX [::tcl::mathfunc::min {*}$x_coords_only]
        set rightmostX [::tcl::mathfunc::max {*}$x_coords_only]
        set leftmostIndex [lsearch $x_coords_only $leftmostX]
        set rightmostIndex [lsearch $x_coords_only $rightmostX]
        
        # Get atom information
        set atomNames [$sel get name]
        set resNames [$sel get resname]
        set resIds [$sel get resid]
        
        set leftmostCoord [lindex $coords $leftmostIndex]
        set rightmostCoord [lindex $coords $rightmostIndex]
        set distance [veclength [vecsub $rightmostCoord $leftmostCoord]]
        
        # Store extreme atom information
        set ::PETK::gui::leftmostAtom [format "%s %s%d (%.3f, %.3f, %.3f)" \
            [lindex $atomNames $leftmostIndex] [lindex $resNames $leftmostIndex] [lindex $resIds $leftmostIndex] \
            [lindex $leftmostCoord 0] [lindex $leftmostCoord 1] [lindex $leftmostCoord 2]]
        
        set ::PETK::gui::rightmostAtom [format "%s %s%d (%.3f, %.3f, %.3f)" \
            [lindex $atomNames $rightmostIndex] [lindex $resNames $rightmostIndex] [lindex $resIds $rightmostIndex] \
            [lindex $rightmostCoord 0] [lindex $rightmostCoord 1] [lindex $rightmostCoord 2]]
        
        set ::PETK::gui::analyteDistance [format "%.3f Å" $distance]

        # Calculate bounding sphere
        set center [measure center $sel]
        set maxDist 0.0
        foreach coord $coords {
            set dist [veclength [vecsub $coord $center]]
            if {$dist > $maxDist} {
                set maxDist $dist
            }
        }

        # Calculate volume
        set volume [expr {4.0/3.0 * 3.14159 * pow($maxDist, 3)}]
        set ::PETK::gui::analyteDiameter [format "%.2f Å" $maxDist]
        set ::PETK::gui::analyteVolume [format "%.2f Å³" $volume]

        # Update fit status
        ::PETK::gui::updateFitStatus

        # Set centering status
        set ::PETK::gui::centeringStatus "Centered"
        set ::PETK::gui::alignmentStatus "Ready to align"
         puts "\nVerifying alignment quality..."
        set verification_score [::PETK::gui::verifyCentering]

        $sel delete

        puts "\nPHASE 1 COMPLETE:"
        puts "================="
        puts "✓ Molecule centered at origin"
        puts "✓ Principal axes calculated"
        puts "✓ Ready for user selection and alignment"
        puts "=================================="

        # Auto-show the molecule
        ::PETK::gui::showAnalyte
        
        # Update the detailed information display
        ::PETK::gui::updateDetailedInformation
        
        # Enable alignment controls in GUI
        ::PETK::gui::enableAlignmentControls
        ::PETK::gui::refreshPorePreviewIfNeeded
    } result options]
    ::PETK::gui::restoreViewState $viewState
    if {$code} {
        return -code $code -options $options $result
    }
    return $result
}

proc ::PETK::gui::alignToZAxis {} {
    puts "Starting Phase 2: Aligning selected axis to Z..."
    
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        tk_messageBox -icon error -message "Error!" -detail "No molecule loaded. Run Phase 1 analysis first."
        return
    }
    
    if {![info exists ::PETK::gui::selectedAxis]} {
        tk_messageBox -icon error -message "Error!" -detail "No axis selected. Please select an axis to align to Z."
        return
    }
    
    # Update status
    set ::PETK::gui::alignmentStatus "Aligning..."
    update
    
    if {[catch {::PETK::gui::performAxisAlignment} error]} {
        tk_messageBox -icon error -message "Alignment Error" -detail "Error during alignment:\n$error"
        set ::PETK::gui::alignmentStatus "Alignment failed"
        return
    }
    
    # Verify the alignment
    puts "\nVerifying alignment quality..."
    set verification_score [::PETK::gui::verifyCentering]
    
    # Update detailed information
    ::PETK::gui::updateDetailedInformation
    
    puts "\nPhase 2 complete. Verification score: $verification_score/12"
    return $verification_score
}

proc ::PETK::gui::performAxisAlignment {} {
    set molid $::PETK::gui::analyteMol
    
    # Create selection
    if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
        error "Could not create atom selection: $sel"
    }
    
    if {[$sel num] == 0} {
        $sel delete
        error "Selection contains 0 atoms!"
    }

    puts "=================================="
    puts "PHASE 2: AXIS ALIGNMENT TO Z"
    puts "=================================="

    # Get the selected principal axis
    set selected_axis_vector [lindex $::PETK::gui::principalAxes [expr {$::PETK::gui::selectedAxis - 1}]]
    set axis_names {1 2 3}
    set axis_name [lindex $axis_names [expr {$::PETK::gui::selectedAxis - 1}]]
    
    puts "Selected axis: Principal Axis $axis_name"
    puts [format "Axis vector: (%8.5f, %8.5f, %8.5f)" \
          [lindex $selected_axis_vector 0] [lindex $selected_axis_vector 1] [lindex $selected_axis_vector 2]]

    # Calculate rotation matrix to align selected axis with Z-axis
    set target_axis {0 0 1}  ; # Z-axis
    set R [::PETK::gui::calculateRotationMatrix $selected_axis_vector $target_axis]
    
    puts "\nRotation matrix to align axis to Z:"
    for {set i 0} {$i < 3} {incr i} {
        set row [lindex $R $i]
        puts [format "  \[%8.5f %8.5f %8.5f\]" [lindex $row 0] [lindex $row 1] [lindex $row 2]]
    }
    
    # Apply rotation to all atoms
    set coords [$sel get {x y z}]
    set new_coords {}
    foreach coord $coords {
        set rotated_coord [::PETK::gui::matrix_vector_multiply_3x3 $R $coord]
        lappend new_coords $rotated_coord
    }
    
    # Set new coordinates
    $sel set {x y z} $new_coords
    puts "✓ Selected principal axis aligned to Z-axis"
    
    # Verify the alignment
    puts "\nVerifying axis alignment:"
    
    # Recalculate principal axes after rotation
    set new_inertia_result [measure inertia $sel moments]
    set new_principal_axes [lindex $new_inertia_result 1]
    set new_axis_aligned_to_z [lindex $new_principal_axes [expr {$::PETK::gui::selectedAxis - 1}]]
    
    # Check alignment with Z-axis
    set z_axis {0 0 1}
    set dot_product [expr {[lindex $new_axis_aligned_to_z 0] * 0 + \
                           [lindex $new_axis_aligned_to_z 1] * 0 + \
                           [lindex $new_axis_aligned_to_z 2] * 1}]
    set alignment_quality [expr {abs($dot_product)}]
    
    puts [format "Aligned axis: (%8.5f, %8.5f, %8.5f)" \
          [lindex $new_axis_aligned_to_z 0] [lindex $new_axis_aligned_to_z 1] [lindex $new_axis_aligned_to_z 2]]
    puts [format "Alignment with Z-axis: %.6f (1.0 = perfect)" $alignment_quality]
    
    if {$alignment_quality > 0.999} {
        puts "✓ EXCELLENT: Perfect alignment achieved"
        set ::PETK::gui::alignmentStatus "Excellent"
    } elseif {$alignment_quality > 0.99} {
        puts "✓ VERY GOOD: Near-perfect alignment"
        set ::PETK::gui::alignmentStatus "Very Good"
    } elseif {$alignment_quality > 0.95} {
        puts "✓ GOOD: Good alignment"
        set ::PETK::gui::alignmentStatus "Good"
    } else {
        puts "⚠ WARNING: Alignment may need improvement"
        set ::PETK::gui::alignmentStatus "Fair"
    }
    
    # Update coordinate ranges after alignment
    set coords [$sel get {x y z}]
    set x_coords {}
    set y_coords {}
    set z_coords {}
    
    foreach coord $coords {
        lappend x_coords [lindex $coord 0]
        lappend y_coords [lindex $coord 1]
        lappend z_coords [lindex $coord 2]
    }
    
    set x_range [expr {[::tcl::mathfunc::max {*}$x_coords] - [::tcl::mathfunc::min {*}$x_coords]}]
    set y_range [expr {[::tcl::mathfunc::max {*}$y_coords] - [::tcl::mathfunc::min {*}$y_coords]}]
    set z_range [expr {[::tcl::mathfunc::max {*}$z_coords] - [::tcl::mathfunc::min {*}$z_coords]}]
    
    puts [format "\nMolecular dimensions after alignment:"]
    puts [format "  X dimension: %.3f Å" $x_range]
    puts [format "  Y dimension: %.3f Å" $y_range]
    puts [format "  Z dimension: %.3f Å" $z_range]

    ::PETK::gui::updateMovementRangeFromSelection $sel

    $sel delete

    puts "\nPHASE 2 COMPLETE:"
    puts "================="
    puts "✓ Selected principal axis aligned to Z"
    puts "✓ Molecule properly oriented"
    puts "✓ Ready for verification and saving"
    puts "=================================="
    ::PETK::gui::onAlignmentComplete
}

proc ::PETK::gui::showAnalyte {} {
    if {$::PETK::gui::analyteMol != ""} {
        mol on $::PETK::gui::analyteMol
        mol top $::PETK::gui::analyteMol
        scale by 1.0
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

proc ::PETK::gui::resetToCentered {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol == ""} {
        tk_messageBox -icon error -message "Error!" -detail "No molecule loaded."
        return
    }
    
    if {![info exists ::PETK::gui::centeredCoords]} {
        tk_messageBox -icon error -message "Error!" -detail "No centered coordinates stored. Please run Phase 1 first."
        return
    }
    
    set molid $::PETK::gui::analyteMol
    
    if {[catch {atomselect $molid $::PETK::gui::analyteSelection} sel]} {
        tk_messageBox -icon error -message "Error!" -detail "Could not create atom selection."
        return
    }
    
    # Reset coordinates
    $sel set {x y z} $::PETK::gui::centeredCoords
    
    # Verify center
    set center [measure center $sel]
    set center_dist [veclength $center]
    
    puts [format "Reset complete. Center: (%.6f, %.6f, %.6f), distance: %.6f" \
          [lindex $center 0] [lindex $center 1] [lindex $center 2] $center_dist]
    
    $sel delete
    
    # Update status
    set ::PETK::gui::alignmentStatus "Reset to centered"
    set ::PETK::gui::moleculeAligned 0
    ::PETK::gui::disableSaveButton
    ::PETK::gui::onResetComplete
    tk_messageBox -icon info -message "Reset Complete" -detail "Molecule reset to centered state.\nYou can now select a different axis alignment."
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
    # First check what type of pore we're dealing with
    if {![info exists ::PETK::gui::poreOption]} {
        set ::PETK::gui::fitStatus "Pore type not specified"
        return
    }
    
    # Check if molecule has been analyzed
    if {![info exists ::PETK::gui::analyteDiameter] || $::PETK::gui::analyteDiameter eq "Not analyzed"} {
        set ::PETK::gui::fitStatus "Analyze molecule first"
        return
    }
    
    # Extract numeric value from analyteDiameter
    if {![regexp {^([\d.]+)} $::PETK::gui::analyteDiameter match molRadius]} {
        set ::PETK::gui::fitStatus "Invalid molecule radius"
        return
    }
    set molRad [expr {double($molRadius)}]
    
    # Handle different pore types
    if {$::PETK::gui::poreOption eq "biological"} {
        # For biological pores, calculate minimum local radius
        if {![info exists ::PETK::gui::biological_pore_molid] || $::PETK::gui::biological_pore_molid eq ""} {
            set ::PETK::gui::fitStatus "Load biological pore first"
            return
        }
        
        # Check if biological pore molecule still exists
        if {[catch {molinfo $::PETK::gui::biological_pore_molid get numatoms} num_atoms]} {
            set ::PETK::gui::fitStatus "Biological pore not loaded"
            return
        }
        
        if {$num_atoms == 0} {
            set ::PETK::gui::fitStatus "Empty biological pore structure"
            return
        }
        
        # Get membrane parameters
        set membrane_thickness $::PETK::gui::nanoporeThickness
        set membrane_half_thickness [expr $membrane_thickness / 2.0]
        
        # Get membrane Z offset (default to 0 if not set)
        if {[info exists ::PETK::gui::membraneZOffset]} {
            set membrane_z_offset $::PETK::gui::membraneZOffset
        } else {
            set membrane_z_offset 0.0
        }
        
        # Calculate membrane boundaries
        set pore_start_z [expr $membrane_z_offset - $membrane_half_thickness]
        set pore_end_z [expr $membrane_z_offset + $membrane_half_thickness]
        
        # Sample pore radius at multiple Z positions to find minimum
        set z_step 2.0  ;# Sample every 2 Å
        set min_pore_radius 1000.0  ;# Start with large value
        set valid_samples 0
        
        for {set z $pore_start_z} {$z <= $pore_end_z} {set z [expr $z + $z_step]} {
            if {[catch {
                set local_radius [::PETK::gui::calculateLocalPoreRadius $z $membrane_z_offset $membrane_half_thickness]
                if {$local_radius > 0 && $local_radius < $min_pore_radius} {
                    set min_pore_radius $local_radius
                }
                incr valid_samples
            }]} {
                # Skip this sample if calculation fails
                continue
            }
        }
        
        if {$valid_samples == 0} {
            set ::PETK::gui::fitStatus "Cannot calculate pore radius"
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval configure -foreground orange
            }
            return
        }
        
        # Compare minimum pore radius with molecule radius
        if {$molRad < $min_pore_radius} {
            set clearance [expr {$min_pore_radius - $molRad}]
            set ::PETK::gui::fitStatus "FITS (clearance: [format "%.1f" $clearance] Å, min pore: [format "%.1f" $min_pore_radius] Å)"
            # Update color to green
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval configure -foreground green
            }
        } else {
            set excess [expr {$molRad - $min_pore_radius}]
            set ::PETK::gui::fitStatus "TOO LARGE (excess: [format "%.1f" $excess] Å, min pore: [format "%.1f" $min_pore_radius] Å)"
            # Update color to red
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval configure -foreground red
            }
        }
        
    } elseif {$::PETK::gui::poreOption eq "solid-state"} {
        # For solid-state pores, use the original diameter-based logic
        if {![info exists ::PETK::gui::poreDiameter] || $::PETK::gui::poreDiameter eq ""} {
            set ::PETK::gui::fitStatus "Set pore diameter to check"
            return
        }
        
        set poreRad [expr {double($::PETK::gui::poreDiameter) / 2.0}]
        
        if {$molRad < $poreRad} {
            set clearance [expr {$poreRad - $molRad}]
            set ::PETK::gui::fitStatus "FITS (clearance: [format "%.1f" $clearance] Å)"
            # Update color if widget exists
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval configure -foreground green
            }
        } else {
            set excess [expr {$molRad - $poreRad}]
            set ::PETK::gui::fitStatus "TOO LARGE (excess: [format "%.1f" $excess] Å)"
            # Update color if widget exists
            if {[info exists ::PETK::gui::window] && [winfo exists $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval]} {
                $::PETK::gui::window.hlf.nb.tab2.canvas.content.results.left.statusval configure -foreground red
            }
        }
        
    } else {
        set ::PETK::gui::fitStatus "Unknown pore type: $::PETK::gui::poreOption"
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
