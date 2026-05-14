#
# PETK GUI - Tab 3 (SEM Setup) module
#
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

    # Structure preparation (pdb2pqr only)
    ttk::label $container.grid.preplbl -text "Structure preparation:" -width 18
    ttk::label $container.grid.prepval -text "Run pdb2pqr (CHARMM)" -anchor w
    ttk::label $container.grid.defaultlbl -text "Default radius (Å):" -width 18
    ttk::entry $container.grid.default -textvariable ::PETK::gui::semDefaultRadius -width 12 -justify center

    grid $container.grid.preplbl $container.grid.prepval $container.grid.defaultlbl $container.grid.default -sticky ew -pady 3

    # === MESH MODE: uniform (single grid resolution) vs multigrid (gmsh fine + coarse) ===
    # A separator frame inside Grid Settings, vertically below the rows above.
    ttk::separator $container.grid.sep -orient horizontal
    grid $container.grid.sep -columnspan 6 -sticky ew -pady "8 4"

    ttk::label $container.grid.meshmodelbl -text "Mesh mode:" -width 18
    ttk::radiobutton $container.grid.meshmode_uniform -text "Uniform" \
        -variable ::PETK::gui::meshMode -value "uniform" \
        -command {::PETK::gui::onMeshModeChanged}
    ttk::radiobutton $container.grid.meshmode_graded -text "Multigrid (gmsh)" \
        -variable ::PETK::gui::meshMode -value "multigrid" \
        -command {::PETK::gui::onMeshModeChanged}
    ttk::label $container.grid.meshmodehelp \
        -text "Multigrid uses a fine cell size in a box around the pore and coarser cells in the bulk." \
        -foreground gray35
    grid $container.grid.meshmodelbl $container.grid.meshmode_uniform \
        $container.grid.meshmode_graded -sticky w -pady 3
    grid $container.grid.meshmodehelp -columnspan 6 -sticky w -pady {0 3}

    # Fine and coarse cell sizes (multigrid mode only).
    ttk::label $container.grid.finelbl -text "Fine cell size (Å):" -width 18
    ttk::entry $container.grid.fine -textvariable ::PETK::gui::gmshFineSize \
        -width 12 -justify center
    ttk::label $container.grid.coarselbl -text "Coarse cell size (Å):" -width 18
    ttk::entry $container.grid.coarse -textvariable ::PETK::gui::gmshCoarseSize \
        -width 12 -justify center
    grid $container.grid.finelbl $container.grid.fine \
        $container.grid.coarselbl $container.grid.coarse -sticky ew -pady 3

    # Fine-box sizing: explicit X/Y/Z or auto-fit-to-pore with margin.
    ttk::label $container.grid.fineboxmodelbl -text "Fine box:" -width 18
    ttk::radiobutton $container.grid.fineboxmode_manual -text "Explicit X/Y/Z" \
        -variable ::PETK::gui::gmshFineBoxMode -value "manual" \
        -command {::PETK::gui::onMeshModeChanged}
    ttk::radiobutton $container.grid.fineboxmode_auto \
        -text "Auto-fit to pore (margin Å)" \
        -variable ::PETK::gui::gmshFineBoxMode -value "auto" \
        -command {::PETK::gui::onMeshModeChanged}
    grid $container.grid.fineboxmodelbl $container.grid.fineboxmode_manual \
        $container.grid.fineboxmode_auto -sticky w -pady 3

    # Explicit X/Y/Z entries (manual fine-box mode).
    ttk::label $container.grid.fbxlbl -text "Fine box X (Å):" -width 18
    ttk::entry $container.grid.fbx -textvariable ::PETK::gui::gmshFineBoxX \
        -width 10 -justify center
    ttk::label $container.grid.fbylbl -text "Y (Å):" -width 6
    ttk::entry $container.grid.fby -textvariable ::PETK::gui::gmshFineBoxY \
        -width 10 -justify center
    ttk::label $container.grid.fbzlbl -text "Z (Å):" -width 6
    ttk::entry $container.grid.fbz -textvariable ::PETK::gui::gmshFineBoxZ \
        -width 10 -justify center
    grid $container.grid.fbxlbl $container.grid.fbx \
        $container.grid.fbylbl $container.grid.fby \
        $container.grid.fbzlbl $container.grid.fbz -sticky ew -pady 3

    # Auto-margin entry (auto fine-box mode).
    ttk::label $container.grid.autolbl -text "Auto margin (Å):" -width 18
    ttk::entry $container.grid.auto -textvariable ::PETK::gui::gmshAutoMargin \
        -width 12 -justify center
    ttk::label $container.grid.autohelp \
        -text "Fine box = (pore widest + 2·margin) × (pore widest + 2·margin) × (membrane + 2·margin)." \
        -foreground gray35 -wraplength 480
    grid $container.grid.autolbl $container.grid.auto -sticky ew -pady 3
    grid $container.grid.autohelp -columnspan 6 -sticky w -pady {0 3}

    # Stash references so the mode-change handler can toggle states.
    set ::PETK::gui::tab3MeshGradedWidgets [list \
        $container.grid.fine \
        $container.grid.coarse \
        $container.grid.fineboxmode_manual \
        $container.grid.fineboxmode_auto \
    ]
    set ::PETK::gui::tab3MeshFineBoxManualWidgets [list \
        $container.grid.fbx $container.grid.fby $container.grid.fbz \
    ]
    set ::PETK::gui::tab3MeshFineBoxAutoWidgets [list \
        $container.grid.auto \
    ]

    # === CALCULATION MODE SECTION ===
    # Each mode is a radio button on its own row, followed by a one-line
    # description indented underneath. Keeps the three options visually
    # distinct and gives the long "Hybrid" label room to breathe.
    ttk::labelframe $container.calcmode -text "Calculation Mode" -padding 10
    grid $container.calcmode -row $row -column 0 -sticky ew -padx 10 -pady "5 5"
    grid columnconfigure $container.calcmode 0 -weight 1
    incr row

    # Description-style label: one shared style so the three captions render
    # consistently even if Tk picks different defaults on macOS / Linux / Windows.
    ttk::style configure PETK.CalcModeDesc.TLabel \
        -font {TkDefaultFont 9} -foreground gray35

    ttk::radiobutton $container.calcmode.run -text "Vertical Passage" \
        -variable ::PETK::gui::semCalculationMode -value "run" \
        -command {::PETK::gui::onSemCalculationModeChanged}
    ttk::label $container.calcmode.run_desc \
        -style PETK.CalcModeDesc.TLabel \
        -text "Translocate the analyte along Z at a fixed orientation." \
        -wraplength 520 -justify left

    ttk::radiobutton $container.calcmode.rotation -text "Rotation Scan" \
        -variable ::PETK::gui::semCalculationMode -value "rotation_scan" \
        -command {::PETK::gui::onSemCalculationModeChanged}
    ttk::label $container.calcmode.rotation_desc \
        -style PETK.CalcModeDesc.TLabel \
        -text "Sample random orientations at a fixed Z position (no Z sweep)." \
        -wraplength 520 -justify left

    ttk::radiobutton $container.calcmode.hybrid -text "Hybrid (Translocation + Rotation)" \
        -variable ::PETK::gui::semCalculationMode -value "hybrid" \
        -command {::PETK::gui::onSemCalculationModeChanged}
    ttk::label $container.calcmode.hybrid_desc \
        -style PETK.CalcModeDesc.TLabel \
        -text "Run the full Z sweep at each sampled orientation." \
        -wraplength 520 -justify left

    # Stack vertically: radio, indented description, radio, indented description, …
    grid $container.calcmode.run          -row 0 -column 0 -sticky w -pady "2 0"
    grid $container.calcmode.run_desc     -row 1 -column 0 -sticky w -padx {22 0} -pady "0 6"
    grid $container.calcmode.rotation     -row 2 -column 0 -sticky w -pady "2 0"
    grid $container.calcmode.rotation_desc -row 3 -column 0 -sticky w -padx {22 0} -pady "0 6"
    grid $container.calcmode.hybrid       -row 4 -column 0 -sticky w -pady "2 0"
    grid $container.calcmode.hybrid_desc  -row 5 -column 0 -sticky w -padx {22 0} -pady "0 2"

    # === MOVEMENT PARAMETERS (RUN MODE) ===
    ttk::labelframe $container.movement -text "Analyte Movement Parameters" -padding 10
    grid columnconfigure $container.movement {1 3 5} -weight 1

    ttk::label $container.movement.startlbl -text "Z start position (Å):" -width 18
    ttk::entry $container.movement.start -textvariable ::PETK::gui::zStartRange -width 12 -justify center
    ttk::label $container.movement.endlbl -text "Z end position (Å):" -width 18
    ttk::entry $container.movement.end -textvariable ::PETK::gui::zEndRange -width 12 -justify center
    ttk::label $container.movement.steplbl -text "Z step size (Å):" -width 15
    ttk::entry $container.movement.step -textvariable ::PETK::gui::zStep -width 12 -justify center

    grid $container.movement.startlbl $container.movement.start $container.movement.endlbl $container.movement.end \
        $container.movement.steplbl $container.movement.step -sticky ew -pady 3

    # === ROTATION SCAN PARAMETERS ===
    ttk::labelframe $container.rotation -text "Rotation Scan Parameters" -padding 10
    grid columnconfigure $container.rotation 1 -weight 1

    ttk::label $container.rotation.sampleslbl -text "Number of samples:" -width 20
    ttk::entry $container.rotation.samples -textvariable ::PETK::gui::rotationSamples -width 12 -justify center
    ttk::label $container.rotation.seedlbl -text "Random seed:" -width 20
    ttk::entry $container.rotation.seed -textvariable ::PETK::gui::rotationSeed -width 12 -justify center
    ttk::frame $container.rotation.options
    ttk::checkbutton $container.rotation.options.bulkchk -text "In Bulk Electrolyte" \
        -variable ::PETK::gui::useBulkElectrolyte -command ::PETK::gui::applyBulkElectrolyteSetting
    ttk::checkbutton $container.rotation.options.reusemesh -text "Reuse mesh" \
        -variable ::PETK::gui::reuseSemMesh
    pack $container.rotation.options.bulkchk $container.rotation.options.reusemesh -side left -padx 5

    grid $container.rotation.sampleslbl $container.rotation.samples -sticky ew -pady 3
    grid $container.rotation.seedlbl $container.rotation.seed -sticky ew -pady 3
    grid $container.rotation.options -column 0 -columnspan 2 -sticky w -pady "0 3"
    set ::PETK::gui::bulkElectrolyteButton $container.rotation.options.bulkchk
    if {[info procs ::PETK::gui::updateBulkElectrolyteButtonState] ne ""} {
        ::PETK::gui::updateBulkElectrolyteButtonState
    }

    set ::PETK::gui::tab3MovementFrame $container.movement
    set ::PETK::gui::tab3RotationFrame $container.rotation
    set ::PETK::gui::tab3MovementRow $row
    incr row
    set ::PETK::gui::tab3RotationRow $row
    incr row

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

    incr row
    # === ARBD EXPORT SECTION ===
    ttk::labelframe $container.arbd -text "ARBD Grid Export (optional)" -padding 10
    grid $container.arbd -row $row -column 0 -sticky ew -padx 10 -pady "5 5"
    grid columnconfigure $container.arbd {1 3} -weight 1
    incr row

    ttk::checkbutton $container.arbd.enable \
        -text "Export phi / per-ion / steric DX grids per step" \
        -variable ::PETK::gui::semArbdEnabled
    grid $container.arbd.enable -row 0 -column 0 -columnspan 4 -sticky w -pady "0 4"

    ttk::label $container.arbd.ionslbl -text "Ions (name:valence, comma-separated):" -width 32 -anchor w
    ttk::entry $container.arbd.ions -textvariable ::PETK::gui::semArbdIons -width 24 -justify left
    ttk::label $container.arbd.stridelbl -text "Stride (0 = every step):" -width 22 -anchor w
    ttk::entry $container.arbd.stride -textvariable ::PETK::gui::semArbdStride -width 8 -justify center
    grid $container.arbd.ionslbl $container.arbd.ions $container.arbd.stridelbl $container.arbd.stride -sticky ew -pady 3 -padx 5

    ttk::label $container.arbd.walllbl -text "Steric wall height (kcal/mol):" -width 32 -anchor w
    ttk::entry $container.arbd.wall -textvariable ::PETK::gui::semArbdWallHeight -width 12 -justify center
    ttk::label $container.arbd.templbl -text "Temperature (K):" -width 22 -anchor w
    ttk::entry $container.arbd.temp -textvariable ::PETK::gui::semArbdTemperature -width 8 -justify center
    grid $container.arbd.walllbl $container.arbd.wall $container.arbd.templbl $container.arbd.temp -sticky ew -pady 3 -padx 5

    ttk::label $container.arbd.help -anchor w -wraplength 700 \
        -text "Produces {prefix}_{pore_type}_z{...}_open_pore_phi.dx (volts), one DX per listed ion (kcal/mol), and *_steric.dx, for each solved step. Output goes to the workdir."
    grid $container.arbd.help -columnspan 4 -sticky w -pady "4 0"

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
    set ::PETK::gui::tab3ActionsValidateButton $container.actions.validate
    set ::PETK::gui::tab3ActionsGenerateButton $container.actions.generate
    set ::PETK::gui::tab3ActionsRunButton $container.actions.visualize

    incr row
    # === PLOTTING SECTION (STATIC) ===
    ttk::labelframe $container.plotting -text "Result Plotting" -padding 10
    grid $container.plotting -row $row -column 0 -sticky ew -padx 10 -pady "5 10"
    grid columnconfigure $container.plotting 0 -weight 1
    incr row

    ttk::label $container.plotting.desc -text "Plot stored SEM results using current workdir/output prefix." -wraplength 520
    grid $container.plotting.desc -sticky w -pady "0 5"

    ttk::label $container.plotting.serieslbl -text "Select series to plot:" -font {TkDefaultFont 10 bold}
    grid $container.plotting.serieslbl -sticky w

    ttk::frame $container.plotting.series
    grid $container.plotting.series -sticky w -pady "2 8"
    grid columnconfigure $container.plotting.series {0 1} -weight 1

    ttk::radiobutton $container.plotting.series.current -text "Current (nA)" \
        -variable ::PETK::gui::plot::selectedSeries -value "current" \
        -command {::PETK::gui::plot::setSeries "current"}
    ttk::radiobutton $container.plotting.series.blockage -text "Δ Conductance (nS)" \
        -variable ::PETK::gui::plot::selectedSeries -value "conductance" \
        -command {::PETK::gui::plot::setSeries "conductance"}
    grid $container.plotting.series.current $container.plotting.series.blockage -sticky w -padx 3

    ttk::frame $container.plotting.buttons
    grid $container.plotting.buttons -sticky ew -pady "5 0"
    grid columnconfigure $container.plotting.buttons {0 1} -weight 1

    ttk::button $container.plotting.buttons.current -text "Plot Current Data" -command {
        set data_dir $::PETK::gui::workdir
        set results_dir [file join $::PETK::gui::workdir "results"]
        if {[file isdirectory $results_dir]} {
            set data_dir $results_dir
        }
        set alldata [::PETK::gui::plot::parseDirectory $data_dir $::PETK::gui::outputPrefix]
        set title "PETK Results - [file tail $data_dir]"
        ::PETK::gui::plot::createPlot $alldata $title $data_dir
    }
    ttk::button $container.plotting.buttons.rotation -text "Plot Rotation Scan" -command ::PETK::gui::plot::plotRotationHistogram

    grid $container.plotting.buttons.current $container.plotting.buttons.rotation -sticky ew -padx 3

    ttk::label $container.plotting.tip -text "Tip: Run a simulation first so the latest data appears in results/." \
        -foreground darkblue -font {TkDefaultFont 9 italic}
    grid $container.plotting.tip -sticky w -pady "5 0"

    # Add some bottom padding to ensure scrolling works properly
    ttk::frame $container.bottom_padding
    grid $container.bottom_padding -sticky ew -pady 10

    # === ENHANCED RESIZE HANDLING ===
    
    # Bind resize events for proper canvas expansion
    bind $tab3.canvas <Configure> [list ::PETK::gui::onCanvasConfigured $tab3.canvas $canvas_window $container]
    bind $container <Configure> [list ::PETK::gui::onContentConfigured $tab3.canvas $canvas_window $container]
    
    # Store references for later use (with unique names for tab3)
    set ::PETK::gui::tab3MainCanvas $tab3.canvas
    set ::PETK::gui::tab3CanvasWindow $canvas_window
    set ::PETK::gui::tab3ContentContainer $container
    
    # Initialize SEM variables
    ::PETK::gui::initializeSEMVariables
    ::PETK::gui::onSemCalculationModeChanged
    ::PETK::gui::onMeshModeChanged

    # Force proper sizing after everything is created
    after idle [list ::PETK::gui::updateTab3ScrollRegion]
    after idle [list ::PETK::gui::forceInitialCanvasResize $tab3.canvas $canvas_window $container]
}

####################################################
# Tab 3 Function
####################################################
proc ::PETK::gui::updateTab3ScrollRegion {} {
    # Update scroll region specifically for Tab 3
    if {[info exists ::PETK::gui::tab3MainCanvas] && [winfo exists $::PETK::gui::tab3MainCanvas]} {
        set canvas $::PETK::gui::tab3MainCanvas
        set container $::PETK::gui::tab3ContentContainer
        
        if {[winfo exists $container]} {
            # Force geometry update
            update idletasks
            
            # Get actual content dimensions
            set content_width [winfo reqwidth $container]
            set content_height [winfo reqheight $container]
            set canvas_width [winfo width $canvas]
            
            # Ensure content width matches canvas width
            if {$canvas_width > 1} {
                set content_width $canvas_width
            }
            
            # Update the canvas window item and scroll region
            if {[info exists ::PETK::gui::tab3CanvasWindow]} {
                $canvas itemconfig $::PETK::gui::tab3CanvasWindow -width $content_width
            }
            
            # Set scroll region with some padding
            set scroll_height [expr {$content_height + 20}]  ;# Add 20px padding at bottom
            $canvas configure -scrollregion [list 0 0 $content_width $scroll_height]
            
            puts "Updated Tab 3 scroll region: ${content_width}x${scroll_height}"
        }
    }
}
####################################################
## SEM VARIABLE INITIALIZATION
####################################################

proc ::PETK::gui::initializeSEMVariables {} {
    if {![info exists ::PETK::gui::semCalculationMode]} {
        set ::PETK::gui::semCalculationMode "run"
    }
    if {![info exists ::PETK::gui::rotationSamples]} {
        set ::PETK::gui::rotationSamples "10"
    }
    if {![info exists ::PETK::gui::rotationSeed]} {
        set ::PETK::gui::rotationSeed "42"
    }
    # Always exists so the parameter summary, validation, and run-command
    # builder can read it without an "info exists" guard. Tab 2 syncs this when
    # an analyte is loaded; before that it stays empty.
    if {![info exists ::PETK::gui::semAnalytePDB]} {
        set ::PETK::gui::semAnalytePDB ""
    }

    # Mesh mode (uniform / graded) and graded-mesh parameters.
    if {![info exists ::PETK::gui::meshMode]}        { set ::PETK::gui::meshMode "uniform" }
    if {![info exists ::PETK::gui::gmshFineSize]}    { set ::PETK::gui::gmshFineSize "1.0" }
    if {![info exists ::PETK::gui::gmshCoarseSize]}  { set ::PETK::gui::gmshCoarseSize "5.0" }
    if {![info exists ::PETK::gui::gmshFineBoxMode]} { set ::PETK::gui::gmshFineBoxMode "auto" }
    if {![info exists ::PETK::gui::gmshFineBoxX]}    { set ::PETK::gui::gmshFineBoxX "100.0" }
    if {![info exists ::PETK::gui::gmshFineBoxY]}    { set ::PETK::gui::gmshFineBoxY "100.0" }
    if {![info exists ::PETK::gui::gmshFineBoxZ]}    { set ::PETK::gui::gmshFineBoxZ "100.0" }
    if {![info exists ::PETK::gui::gmshAutoMargin]}  { set ::PETK::gui::gmshAutoMargin "50.0" }
    if {![info exists ::PETK::gui::useBulkElectrolyte]} {
        set ::PETK::gui::useBulkElectrolyte 0
    }
    if {![info exists ::PETK::gui::reuseSemMesh]} {
        set ::PETK::gui::reuseSemMesh 0
    }
    # Simulation parameters
    if {![info exists ::PETK::gui::appliedVoltage]} {
        set ::PETK::gui::appliedVoltage "100"
    }
    if {![info exists ::PETK::gui::bulkConductivity]} {
        set ::PETK::gui::bulkConductivity "11.2"
    }
    if {![info exists ::PETK::gui::membraneConductivity]} {
        set ::PETK::gui::membraneConductivity "0.0000001"
    }
    if {![info exists ::PETK::gui::originalMembraneConductivity]} {
        set ::PETK::gui::originalMembraneConductivity $::PETK::gui::membraneConductivity
    }
    
    # Grid settings
    if {![info exists ::PETK::gui::gridResolution]} {
        set ::PETK::gui::gridResolution "1.0"
    }
    if {![info exists ::PETK::gui::estimatedGridPoints]} {
        set ::PETK::gui::estimatedGridPoints "Not calculated"
    }
    if {![info exists ::PETK::gui::semUseVdWRadii]} {
        set ::PETK::gui::semUseVdWRadii 0
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

    # ARBD-compatible DX export (phi + per-ion + steric grids).
    # When semArbdEnabled is 1, outputParametersToConfig emits an "arbd_export"
    # block under output{} that sem.cli reads and passes to VerticalMovementSEM.
    # See sem/arbd_export.py for output details.
    if {![info exists ::PETK::gui::semArbdEnabled]}    { set ::PETK::gui::semArbdEnabled 0 }
    if {![info exists ::PETK::gui::semArbdIons]}       { set ::PETK::gui::semArbdIons "POT:1, CLA:-1" }
    if {![info exists ::PETK::gui::semArbdStride]}     { set ::PETK::gui::semArbdStride "0" }
    if {![info exists ::PETK::gui::semArbdWallHeight]} { set ::PETK::gui::semArbdWallHeight "100.0" }
    if {![info exists ::PETK::gui::semArbdTemperature]} { set ::PETK::gui::semArbdTemperature "295.0" }

    # Python environment
    if {![info exists ::PETK::gui::condaEnvironment]} {
        set ::PETK::gui::condaEnvironment "base"
    }
    if {![info exists ::PETK::gui::pythonExecutable]} {
        set ::PETK::gui::pythonExecutable ""
    }

    # Plotting variables
    # Open pore current reference line
    if {![info exists ::PETK::gui::plot::showOpenPoreLine]} {
        set ::PETK::gui::plot::showOpenPoreLine 1
    }
    if {![info exists ::PETK::gui::plot::openPoreCurrent]} {
        set ::PETK::gui::plot::openPoreCurrent ""
    }
    if {![info exists ::PETK::gui::plot::openPoreCurrentMode]} {
        set ::PETK::gui::plot::openPoreCurrentMode "auto"
    }
}

proc ::PETK::gui::applyBulkElectrolyteSetting {} {
    if {![info exists ::PETK::gui::useBulkElectrolyte]} {
        set ::PETK::gui::useBulkElectrolyte 0
    }
    if {![info exists ::PETK::gui::originalMembraneConductivity]} {
        set ::PETK::gui::originalMembraneConductivity $::PETK::gui::membraneConductivity
    }

    if {![info exists ::PETK::gui::poreOption] || $::PETK::gui::poreOption ne "solid-state"} {
        set ::PETK::gui::useBulkElectrolyte 0
        tk_messageBox -icon info -title "Bulk Electrolyte" \
            -message "The bulk electrolyte option is only available for solid-state nanopores."
        if {[info procs ::PETK::gui::updateBulkElectrolyteButtonState] ne ""} {
            ::PETK::gui::updateBulkElectrolyteButtonState
        }
        return
    }

    if {!$::PETK::gui::useBulkElectrolyte} {
        if {[info exists ::PETK::gui::originalMembraneConductivity]} {
            set ::PETK::gui::membraneConductivity $::PETK::gui::originalMembraneConductivity
        }
        if {[info procs ::PETK::gui::updateSEMParameterSummary] ne ""} {
            ::PETK::gui::updateSEMParameterSummary
        }
        return
    }

    if {![info exists ::PETK::gui::bulkConductivity] || $::PETK::gui::bulkConductivity eq ""} {
        tk_messageBox -icon warning -title "Bulk Electrolyte" \
            -message "Enter a bulk conductivity value before enabling the bulk electrolyte option."
        set ::PETK::gui::useBulkElectrolyte 0
        if {[info procs ::PETK::gui::updateBulkElectrolyteButtonState] ne ""} {
            ::PETK::gui::updateBulkElectrolyteButtonState
        }
        return
    }

    set ::PETK::gui::membraneConductivity $::PETK::gui::bulkConductivity
    puts "Membrane conductivity set to bulk conductivity ($::PETK::gui::bulkConductivity S/m)."

    if {[info procs ::PETK::gui::updateSEMParameterSummary] ne ""} {
        ::PETK::gui::updateSEMParameterSummary
    }
}

proc ::PETK::gui::updateBulkElectrolyteButtonState {} {
    if {![info exists ::PETK::gui::bulkElectrolyteButton]} {
        return
    }
    set widget $::PETK::gui::bulkElectrolyteButton
    if {![winfo exists $widget]} {
        return
    }

    set state disabled
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "solid-state"} {
        set state normal
    } else {
        if {[info exists ::PETK::gui::useBulkElectrolyte]} {
            set ::PETK::gui::useBulkElectrolyte 0
        }
    }
    $widget configure -state $state
}

proc ::PETK::gui::onSemCalculationModeChanged {} {
    puts "SEM calculation mode switched to $::PETK::gui::semCalculationMode"
    ::PETK::gui::updateSemCalculationModeUI
    if {[info procs ::PETK::gui::updateSEMParameterSummary] ne ""} {
        ::PETK::gui::updateSEMParameterSummary
    }
}

proc ::PETK::gui::onMeshModeChanged {} {
    set mode "uniform"
    if {[info exists ::PETK::gui::meshMode]} {
        set mode $::PETK::gui::meshMode
    }
    set fbmode "auto"
    if {[info exists ::PETK::gui::gmshFineBoxMode]} {
        set fbmode $::PETK::gui::gmshFineBoxMode
    }
    set graded_state   [expr {$mode eq "multigrid"                          ? "normal" : "disabled"}]
    set fb_manual_state [expr {$mode eq "multigrid" && $fbmode eq "manual" ? "normal" : "disabled"}]
    set fb_auto_state   [expr {$mode eq "multigrid" && $fbmode eq "auto"   ? "normal" : "disabled"}]

    if {[info exists ::PETK::gui::tab3MeshGradedWidgets]} {
        foreach w $::PETK::gui::tab3MeshGradedWidgets {
            if {[winfo exists $w]} { $w configure -state $graded_state }
        }
    }
    if {[info exists ::PETK::gui::tab3MeshFineBoxManualWidgets]} {
        foreach w $::PETK::gui::tab3MeshFineBoxManualWidgets {
            if {[winfo exists $w]} { $w configure -state $fb_manual_state }
        }
    }
    if {[info exists ::PETK::gui::tab3MeshFineBoxAutoWidgets]} {
        foreach w $::PETK::gui::tab3MeshFineBoxAutoWidgets {
            if {[winfo exists $w]} { $w configure -state $fb_auto_state }
        }
    }
    if {[info procs ::PETK::gui::updateSEMParameterSummary] ne ""} {
        ::PETK::gui::updateSEMParameterSummary
    }
}

proc ::PETK::gui::updateSemCalculationModeUI {} {
    if {![info exists ::PETK::gui::tab3ContentContainer]} {
        return
    }
    set container $::PETK::gui::tab3ContentContainer
    set mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set mode $::PETK::gui::semCalculationMode
    }
    set movement_row 0
    if {[info exists ::PETK::gui::tab3MovementRow]} {
        set movement_row $::PETK::gui::tab3MovementRow
    }
    set rotation_row [expr {$movement_row + 1}]
    if {[info exists ::PETK::gui::tab3RotationRow]} {
        set rotation_row $::PETK::gui::tab3RotationRow
    }

    set movement_frame {}
    if {[info exists ::PETK::gui::tab3MovementFrame]} {
        set movement_frame $::PETK::gui::tab3MovementFrame
    }
    set rotation_frame {}
    if {[info exists ::PETK::gui::tab3RotationFrame]} {
        set rotation_frame $::PETK::gui::tab3RotationFrame
    }

    set show_movement [expr {$mode eq "run"           || $mode eq "hybrid"}]
    set show_rotation [expr {$mode eq "rotation_scan" || $mode eq "hybrid"}]

    if {$movement_frame ne "" && [winfo exists $movement_frame]} {
        set is_gridded [expr {[winfo manager $movement_frame] eq "grid"}]
        if {$show_movement && !$is_gridded} {
            grid $movement_frame -row $movement_row -column 0 -sticky ew -padx 10 -pady "5 5"
        } elseif {!$show_movement && $is_gridded} {
            grid remove $movement_frame
        }
    }
    if {$rotation_frame ne "" && [winfo exists $rotation_frame]} {
        set is_gridded [expr {[winfo manager $rotation_frame] eq "grid"}]
        if {$show_rotation && !$is_gridded} {
            grid $rotation_frame -row $rotation_row -column 0 -sticky ew -padx 10 -pady "5 5"
        } elseif {!$show_rotation && $is_gridded} {
            grid remove $rotation_frame
        }
    }

    # Toggle widget states based on which frame is active
    set movement_state [expr {$show_movement ? "normal" : "disabled"}]
    set rotation_state [expr {$show_rotation ? "normal" : "disabled"}]

    if {$movement_frame ne ""} {
        foreach w {start end step} {
            set path [format "%s.%s" $movement_frame $w]
            if {[winfo exists $path]} {
                $path configure -state $movement_state
            }
        }
    }

    if {$rotation_frame ne ""} {
        foreach w {samples seed boxx boxy boxz} {
            set path [format "%s.%s" $rotation_frame $w]
            if {[winfo exists $path]} {
                $path configure -state $rotation_state
            }
        }
    }

    # In hybrid mode the bulk-electrolyte option doesn't apply (full FEM is forced).
    if {[info exists ::PETK::gui::bulkElectrolyteButton] && \
            [winfo exists $::PETK::gui::bulkElectrolyteButton]} {
        if {$mode eq "hybrid"} {
            $::PETK::gui::bulkElectrolyteButton configure -state disabled
        } else {
            $::PETK::gui::bulkElectrolyteButton configure -state normal
        }
    }

    # Generate (preview) button: enabled when a translocation z-sweep is part of
    # the run, i.e. for "run" and "hybrid" modes.
    if {[info exists ::PETK::gui::tab3ActionsGenerateButton] && [winfo exists $::PETK::gui::tab3ActionsGenerateButton]} {
        $::PETK::gui::tab3ActionsGenerateButton configure -state [expr {$show_movement ? "normal" : "disabled"}]
    }

    if {[info exists ::PETK::gui::tab3ActionsValidateButton] && [winfo exists $::PETK::gui::tab3ActionsValidateButton]} {
        $::PETK::gui::tab3ActionsValidateButton configure -state normal
    }

    if {[info procs ::PETK::gui::updateTab3ScrollRegion] ne ""} {
        after idle ::PETK::gui::updateTab3ScrollRegion
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
    
    set calc_mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set calc_mode $::PETK::gui::semCalculationMode
    }

    # Input Configuration
    $widget insert end "INPUT CONFIGURATION:\n" {header}
    # Prefer the SEM-ready (centered) path; fall back to the raw analytePDB
    # loaded in tab 2 if the centering workflow has not been run yet, and to
    # an explicit "(not set)" otherwise. This keeps the summary readable even
    # before validation has run.
    set sem_pdb_for_summary "(not set)"
    if {[info exists ::PETK::gui::semAnalytePDB] && $::PETK::gui::semAnalytePDB ne ""} {
        set sem_pdb_for_summary $::PETK::gui::semAnalytePDB
    } elseif {[info exists ::PETK::gui::analytePDB] && $::PETK::gui::analytePDB ne ""} {
        set sem_pdb_for_summary "$::PETK::gui::analytePDB  (uncentered — run Center Analyte in tab 2)"
    }
    $widget insert end "  Analyte PDB: $sem_pdb_for_summary\n"
    set conda_env_for_summary "(not set)"
    if {[info exists ::PETK::gui::condaEnvironment] && $::PETK::gui::condaEnvironment ne ""} {
        set conda_env_for_summary $::PETK::gui::condaEnvironment
    }
    $widget insert end "  Python Environment: $conda_env_for_summary\n"
    $widget insert end "\n"

    $widget insert end "CALCULATION MODE:\n" {header}
    set seed_text "42"
    if {[info exists ::PETK::gui::rotationSeed]} {
        set seed_text $::PETK::gui::rotationSeed
    }
    if {$calc_mode eq "rotation_scan"} {
        $widget insert end "  Mode: rotation_scan (rotation-based sampling)\n"
        $widget insert end "  Samples: $::PETK::gui::rotationSamples\n"
        $widget insert end "  Random seed: $seed_text\n"
    } elseif {$calc_mode eq "hybrid"} {
        $widget insert end "  Mode: hybrid (translocation + rotation)\n"
        $widget insert end "  Samples: $::PETK::gui::rotationSamples\n"
        $widget insert end "  Random seed: $seed_text\n"
        $widget insert end "  Z-sweep applied at every rotation.\n"
    } else {
        $widget insert end "  Mode: run (vertical passage)\n"
    }
    $widget insert end "\n"

    # Pore Geometry
    $widget insert end "PORE GEOMETRY:\n" {header}
    
    # Determine pore type correctly based on the two-level selection
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "biological"} {
        set pore_type "biological"
    } elseif {[info exists ::PETK::gui::membraneType]} {
        set pore_type $::PETK::gui::membraneType
    } elseif {[info exists ::PETK::gui::currentPoreType]} {
        set pore_type $::PETK::gui::currentPoreType
    } else {
        set pore_type "unknown"
    }
    
    $widget insert end "  Pore Type: $pore_type\n"
    
    if {$pore_type eq "cylindrical"} {
        $widget insert end "  Pore Diameter: $::PETK::gui::cylindricalDiameter Å\n"
        if {[info exists ::PETK::gui::cornerRadius]} {
            $widget insert end "  Corner Radius: $::PETK::gui::cornerRadius Å\n"
        }
        
    } elseif {$pore_type eq "double_cone" || $pore_type eq "doublecone"} {
        $widget insert end "  Inner Diameter: $::PETK::gui::innerDiameter Å\n"
        $widget insert end "  Outer Diameter: $::PETK::gui::outerDiameter Å\n"

    } elseif {$pore_type eq "conical"} {
        $widget insert end "  Top Diameter: $::PETK::gui::topDiameter Å\n"
        $widget insert end "  Bottom Diameter: $::PETK::gui::bottomDiameter Å\n"

    } elseif {$pore_type eq "biological"} {
        if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
            $widget insert end "  Biological Pore PDB: [file tail $::PETK::gui::selectedBioPore]\n"
            $widget insert end "  Full Path: $::PETK::gui::selectedBioPore\n"
        } elseif {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
            $widget insert end "  Biological Pore PDB: [file tail $::PETK::gui::customBioPorePath] (custom)\n"
            $widget insert end "  Full Path: $::PETK::gui::customBioPorePath\n"
        } else {
            $widget insert end "  Biological Pore PDB: Not selected\n" {warning}
        }
        
        if {[info exists ::PETK::gui::membraneZOffset]} {
            $widget insert end "  Membrane Z Offset: $::PETK::gui::membraneZOffset Å\n"
        } else {
            $widget insert end "  Membrane Z Offset: 0.0 Å (default)\n"
        }
    }
    
    $widget insert end "  Membrane Thickness: $::PETK::gui::nanoporeThickness Å\n"
    $widget insert end "\n"
    
    # Box Dimensions
    $widget insert end "BOX DIMENSIONS:\n" {header}
    $widget insert end "  Distance Padding: $::PETK::gui::sysPadding Å\n"

    set use_auto 0
    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        set use_auto 1
    }
    if {$calc_mode eq "rotation_scan"} {
        set use_auto 0
    }

    if {$use_auto} {
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
    
    # Add membrane conductivity for biological pores
    if {$pore_type eq "biological"} {
        if {[info exists ::PETK::gui::membraneConductivity]} {
            $widget insert end "  Membrane Conductivity: $::PETK::gui::membraneConductivity S/m\n"
        } else {
            $widget insert end "  Membrane Conductivity: 0.0000001 S/m (default)\n"
        }
    }
    
    $widget insert end "  Grid Resolution: $::PETK::gui::gridResolution Å\n"
    set prep_mode "pdb2pqr (CHARMM)"
    $widget insert end "  Structure Prep: $prep_mode\n"
    $widget insert end "  Default Radius: $::PETK::gui::semDefaultRadius Å\n"

    # Mesh mode (uniform / graded) and graded-mesh parameters.
    set mesh_mode_disp "uniform"
    if {[info exists ::PETK::gui::meshMode]} {
        set mesh_mode_disp $::PETK::gui::meshMode
    }
    $widget insert end "  Mesh Mode: $mesh_mode_disp\n"
    if {$mesh_mode_disp eq "multigrid"} {
        $widget insert end "    Fine cell size:   $::PETK::gui::gmshFineSize Å\n"
        $widget insert end "    Coarse cell size: $::PETK::gui::gmshCoarseSize Å\n"
        set fb_mode_disp "auto"
        if {[info exists ::PETK::gui::gmshFineBoxMode]} {
            set fb_mode_disp $::PETK::gui::gmshFineBoxMode
        }
        if {$fb_mode_disp eq "manual"} {
            $widget insert end \
                "    Fine box (Å):     $::PETK::gui::gmshFineBoxX × $::PETK::gui::gmshFineBoxY × $::PETK::gui::gmshFineBoxZ  (explicit)\n"
        } else {
            $widget insert end \
                "    Fine box auto-margin: $::PETK::gui::gmshAutoMargin Å (computed per pore)\n"
        }
    }
    $widget insert end "\n"
    
    if {$calc_mode eq "run" || $calc_mode eq "hybrid"} {
        # Movement Parameters (used in run and hybrid)
        $widget insert end "MOVEMENT PARAMETERS:\n" {header}
        $widget insert end "  Z Start: $::PETK::gui::zStartRange Å\n"
        $widget insert end "  Z End: $::PETK::gui::zEndRange Å\n"
        $widget insert end "  Z Step: $::PETK::gui::zStep Å\n"

        if {[string is double $::PETK::gui::zStartRange] && [string is double $::PETK::gui::zEndRange] && [string is double $::PETK::gui::zStep]} {
            set range [expr {abs($::PETK::gui::zEndRange - $::PETK::gui::zStartRange)}]
            set num_steps [expr {int($range / $::PETK::gui::zStep) + 1}]
            $widget insert end "  Number of Steps: $num_steps\n"
            $widget insert end "  Total Range: $range Å\n"
        }
        $widget insert end "\n"
    }
    if {$calc_mode eq "rotation_scan" || $calc_mode eq "hybrid"} {
        $widget insert end "ROTATION PARAMETERS:\n" {header}
        $widget insert end "  Samples: $::PETK::gui::rotationSamples\n"
        if {[info exists ::PETK::gui::rotationSeed]} {
            $widget insert end "  Random seed: $::PETK::gui::rotationSeed\n"
        }
        $widget insert end "\n"
    }
    
    # Output Settings
    $widget insert end "OUTPUT SETTINGS:\n" {header}
    $widget insert end "  Output Directory: $::PETK::gui::workdir\n"
    $widget insert end "  Output Prefix: $::PETK::gui::outputPrefix\n"
    $widget insert end "  Preview Frames: $::PETK::gui::semPreviewFrames\n"
    if {[info exists ::PETK::gui::semArbdEnabled] && $::PETK::gui::semArbdEnabled} {
        $widget insert end "  ARBD Export: ENABLED (ions: $::PETK::gui::semArbdIons, stride: $::PETK::gui::semArbdStride, wall: $::PETK::gui::semArbdWallHeight kcal/mol, T: $::PETK::gui::semArbdTemperature K)\n"
    } else {
        $widget insert end "  ARBD Export: disabled\n"
    }
    
    # Configure text tags for formatting
    $widget tag configure title -font {TkDefaultFont 12 bold}
    $widget tag configure header -font {TkDefaultFont 10 bold}
    $widget tag configure warning -foreground red
    
    $widget configure -state disabled
}

proc ::PETK::gui::syncMovementRangeToBox {} {
    if {[info exists ::PETK::gui::semCalculationMode] && $::PETK::gui::semCalculationMode ne "run"} {
        return
    }

    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        if {[info procs ::PETK::gui::calculateBoxDimensions] ne ""} {
            ::PETK::gui::calculateBoxDimensions
        }
    }

    set min_z ""
    set max_z ""

    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        if {[info exists ::PETK::gui::autoBoxZ] && $::PETK::gui::autoBoxZ ne ""} {
            if {[regexp {(-?\d+(?:\.\d+)?)\s+to\s+(-?\d+(?:\.\d+)?)\s+\(([-\d\.]+)\)} $::PETK::gui::autoBoxZ -> lo hi dim]} {
                set min_z $lo
                set max_z $hi
            }
        }
    }

    if {$min_z eq "" || $max_z eq ""} {
        if {[info exists ::PETK::gui::boxSizeZ] && [string is double -strict $::PETK::gui::boxSizeZ]} {
            set half [expr {double($::PETK::gui::boxSizeZ) / 2.0}]
            set min_z [expr {floor(-$half)}]
            set max_z [expr {ceil($half)}]
        } elseif {[info exists ::PETK::gui::calculatedBoxSizeZ] && [string is double -strict $::PETK::gui::calculatedBoxSizeZ]} {
            set half [expr {double($::PETK::gui::calculatedBoxSizeZ) / 2.0}]
            set min_z [expr {floor(-$half)}]
            set max_z [expr {ceil($half)}]
        }
    }

    if {$min_z eq "" || $max_z eq ""} {
        puts "Warning: Unable to determine Z bounds from current box dimensions."
        return
    }

    if {$min_z > $max_z} {
        set tmp $min_z
        set min_z $max_z
        set max_z $tmp
    }

    if {[info procs ::PETK::gui::updateMovementRangeFromBounds] ne ""} {
        ::PETK::gui::updateMovementRangeFromBounds $min_z $max_z
    } else {
        set ::PETK::gui::zStartRange [::PETK::gui::formatMovementValue $max_z]
        set ::PETK::gui::zEndRange [::PETK::gui::formatMovementValue $min_z]
    }

    ::PETK::gui::adjustMovementRangeForAnalyte $min_z $max_z

    puts [format "Analyte movement range updated to %s → %s Å" $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
}

proc ::PETK::gui::getAnalyteZSpan {} {
    if {![info exists ::PETK::gui::analyteMol] || $::PETK::gui::analyteMol eq ""} {
        return [list 0.0 0.0 0.0]
    }

    set molid $::PETK::gui::analyteMol
    if {[catch {molinfo $molid get numatoms} num_atoms] || $num_atoms <= 0} {
        return [list 0.0 0.0 0.0]
    }

    set selection "all"
    if {[info exists ::PETK::gui::analyteSelection] && $::PETK::gui::analyteSelection ne ""} {
        set selection $::PETK::gui::analyteSelection
    }

    if {[catch {atomselect $molid $selection} sel]} {
        return [list 0.0 0.0 0.0]
    }

    set span 0.0
    set z_min 0.0
    set z_max 0.0
    if {[catch {
        set minmax [measure minmax $sel]
        set min_coord [lindex $minmax 0]
        set max_coord [lindex $minmax 1]
        set z_min [lindex $min_coord 2]
        set z_max [lindex $max_coord 2]
        set span [expr {$z_max - $z_min}]
    } err]} {
        puts "Warning: Unable to determine analyte Z span: $err"
        set span 0.0
    }

    catch {$sel delete}

    if {$span < 0} {
        set span [expr {abs($span)}]
    }

    return [list $span $z_min $z_max]
}

proc ::PETK::gui::adjustMovementRangeForAnalyte {box_min_z box_max_z} {
    set span_info [::PETK::gui::getAnalyteZSpan]
    set analyte_span [lindex $span_info 0]

    if {$analyte_span <= 0} {
        set ::PETK::gui::zStartRange [::PETK::gui::formatMovementValue 0]
        set ::PETK::gui::zEndRange [::PETK::gui::formatMovementValue 0]
        puts "No analyte geometry available; movement range clamped to 0 Å."
        return
    }

    set half_span [expr {$analyte_span / 2.0}]
    set safe_start [expr {$box_max_z - $half_span}]
    set safe_end [expr {$box_min_z + $half_span}]

    if {$safe_start < $safe_end} {
        puts [format "Warning: Analyte height (%.2f Å) exceeds box capacity. Keeping original bounds." $analyte_span]
        set safe_start $box_max_z
        set safe_end $box_min_z
    }

    set ::PETK::gui::zStartRange [::PETK::gui::formatMovementValue $safe_start]
    set ::PETK::gui::zEndRange [::PETK::gui::formatMovementValue $safe_end]

    puts [format "Adjusted movement for analyte span %.2f Å: %s → %s Å" \
        $analyte_span $::PETK::gui::zEndRange $::PETK::gui::zStartRange]
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
    
    ::PETK::gui::updateStatusDisplay "Testing DOLFINx environment..."
    update
    
    # Test basic Python
    if {[catch {exec $::PETK::gui::pythonExecutable --version} python_version]} {
        ::PETK::gui::updateStatusDisplay "❌ Python test failed:\n$python_version"
        return
    }
    
    ::PETK::gui::appendStatusDisplay "✓ Python version: $python_version"
    
    # Test required packages based on setup.py
    set test_script "
import sys
import subprocess

# Core dependencies from setup.py install_requires
core_packages = \['numpy', 'scipy', 'MDAnalysis', 'matplotlib', 'pdb2pqr'\]

# DOLFINx ecosystem (should be installed via conda)
dolfinx_packages = \['dolfinx', 'mpi4py', 'petsc4py'\]

# SEM package itself
sem_packages = \['sem'\]

all_packages = core_packages + dolfinx_packages + sem_packages
missing = \[\]
available = \[\]
warnings = \[\]

for package in all_packages:
    try:
        if package == 'sem':
            # Test SEM package and its entry points
            import sem
            available.append('sem')
            print(f'{package} OK')
            
            # Test console scripts
            try:
                result = subprocess.run(\['sem', '--help'\], capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    print('sem CLI OK')
                else:
                    warnings.append('sem CLI may have issues')
            except Exception as e:
                warnings.append(f'sem CLI test failed: {e}')
                
        elif package == 'dolfinx':
            # Test DOLFINx specifically
            import dolfinx
            import dolfinx.fem
            import dolfinx.mesh
            available.append('dolfinx')
            print(f'{package} OK - version: {dolfinx.__version__}')
        elif package == 'pdb2pqr':
            # Test pdb2pqr module and CLI
            import pdb2pqr
            available.append('pdb2pqr')
            print(f'{package} OK')
            try:
                result = subprocess.run(\['pdb2pqr', '--help'\], capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    print('pdb2pqr CLI OK')
                else:
                    warnings.append('pdb2pqr CLI may have issues')
            except Exception as e:
                warnings.append(f'pdb2pqr CLI test failed: {e}')
            
        elif package == 'matplotlib':
            # Test matplotlib (might be matplotlib-base in conda)
            try:
                import matplotlib
                import matplotlib.pyplot
                available.append('matplotlib')
                print(f'{package} OK')
            except ImportError:
                try:
                    import matplotlib.pyplot
                    available.append('matplotlib')
                    print(f'{package} OK (base)')
                except ImportError:
                    missing.append(package)
                    print(f'{package} MISSING')
                    continue
                    
        else:
            __import__(package)
            available.append(package)
            print(f'{package} OK')
            
    except ImportError as e:
        missing.append(package)
        print(f'{package} MISSING: {e}')
    except Exception as e:
        missing.append(package)
        print(f'{package} FAILED: {e}')

print(f'AVAILABLE: {available}')
if missing:
    print(f'MISSING: {missing}')
if warnings:
    print(f'WARNINGS: {warnings}')
    
if not missing:
    print('ALL_OK')
else:
    print('SOME_MISSING')
"
    
    set exit_code [catch {exec conda run -n $::PETK::gui::condaEnvironment python -c $test_script 2>@1} result]
    
    # Parse the results regardless of warnings
    if {[string match "*ALL_OK*" $result]} {
        ::PETK::gui::appendStatusDisplay "✅ All required packages available!"
        ::PETK::gui::appendStatusDisplay "DOLFINx SEM environment ready for simulation."
        
        # Check for DOLFINx version info
        if {[regexp {dolfinx OK - version: ([^\r\n]+)} $result match dolfinx_version]} {
            ::PETK::gui::appendStatusDisplay "DOLFINx version: $dolfinx_version"
        }
        
    } elseif {[string match "*SOME_MISSING*" $result]} {
        if {[regexp {MISSING: \[(.*?)\]} $result match missing_list]} {
            # Clean up the missing list format
            set clean_missing [string map {"'" "" "," ", "} $missing_list]
            ::PETK::gui::appendStatusDisplay "⚠️ Missing packages: $clean_missing"
            
            # Provide installation hints for missing packages
            if {[string match "*dolfinx*" $clean_missing]} {
                ::PETK::gui::appendStatusDisplay "To install DOLFINx: conda install conda-forge::fenics-dolfinx"
            }
            if {[string match "*sem*" $clean_missing]} {
                ::PETK::gui::appendStatusDisplay "To install SEM package: pip install -e . (in SEM directory)"
            }
        }
    } else {
        ::PETK::gui::appendStatusDisplay "❌ Environment test failed:\n$result"
        return
    }
    
    # Show available packages
    if {[regexp {AVAILABLE: \[(.*?)\]} $result match available_list]} {
        set clean_available [string map {"'" "" "," ", "} $available_list]
        ::PETK::gui::appendStatusDisplay "Available packages: $clean_available"
    }
    
    # Show warnings if any
    if {[regexp {WARNINGS: \[(.*?)\]} $result match warnings_list]} {
        set clean_warnings [string map {"'" "" "," ", "} $warnings_list]
        ::PETK::gui::appendStatusDisplay "⚠️ Warnings: $clean_warnings"
    }
}

proc ::PETK::gui::createSEMEnvironment {} {
    set result [tk_messageBox -icon question -type yesno -title "Create DOLFINx Environment" \
        -message "This will create a new conda environment called 'sem-dolfinx' with DOLFINx.\n\nThis may take several minutes. Continue?"]
    
    if {$result eq "no"} {
        return
    }
    
    ::PETK::gui::updateStatusDisplay "Creating DOLFINx environment 'sem-dolfinx'..."
    ::PETK::gui::appendStatusDisplay "This may take several minutes, please wait..."
    update
    
    # Step 1: Create environment with Python 3.10
    set create_cmd [list conda create -n sem-dolfinx -y python=3.10]
    
    if {[catch {exec {*}$create_cmd} create_result]} {
        ::PETK::gui::updateStatusDisplay "❌ Environment creation failed"
        ::PETK::gui::appendStatusDisplay "Error details:\n$create_result"
        tk_messageBox -icon error -title "Environment Creation Error" \
            -message "Failed to create DOLFINx environment:\n$create_result"
        return
    }
    
    ::PETK::gui::appendStatusDisplay "✅ Base environment created successfully!"
    ::PETK::gui::appendStatusDisplay "Installing DOLFINx package..."
    update
    
    # Step 2: Install fenics-dolfinx from conda-forge
    set install_cmd [list conda install -n sem-dolfinx -y conda-forge::fenics-dolfinx]
    
    if {[catch {exec {*}$install_cmd} install_result]} {
        ::PETK::gui::updateStatusDisplay "❌ DOLFINx installation failed"
        ::PETK::gui::appendStatusDisplay "Error details:\n$install_result"
        tk_messageBox -icon error -title "DOLFINx Installation Error" \
            -message "Failed to install DOLFINx:\n$install_result"
        return
    }
    
    ::PETK::gui::appendStatusDisplay "✅ DOLFINx installed successfully!"
    ::PETK::gui::appendStatusDisplay "Installing pdb2pqr..."
    update

    # Step 3: Install pdb2pqr from conda-forge
    set pdb2pqr_cmd [list conda install -n sem-dolfinx -y -c conda-forge pdb2pqr]

    if {[catch {exec {*}$pdb2pqr_cmd} pdb2pqr_result]} {
        ::PETK::gui::updateStatusDisplay "❌ pdb2pqr installation failed"
        ::PETK::gui::appendStatusDisplay "Error details:\n$pdb2pqr_result"
        tk_messageBox -icon error -title "pdb2pqr Installation Error" \
            -message "Failed to install pdb2pqr:\n$pdb2pqr_result"
        return
    }

    ::PETK::gui::appendStatusDisplay "✅ pdb2pqr installed successfully!"
    ::PETK::gui::appendStatusDisplay "Installing local package in development mode..."
    update
    
    # Step 4: Install local package in development mode
    # Determine repository root (parent directory of the petk module dir)
    set module_dir [::PETK::gui::resourcePath]
    set project_root [file dirname $module_dir]
    set setup_path [file join $project_root setup.py]
    
    if {![file exists $setup_path]} {
        ::PETK::gui::updateStatusDisplay "⚠️ Warning: setup.py not found"
        ::PETK::gui::appendStatusDisplay "Expected setup.py at: $setup_path"
        ::PETK::gui::appendStatusDisplay "Please ensure PETK is installed in editable mode manually."
    } else {
        set pip_cmd [list conda run -n sem-dolfinx pip install -e $project_root]
        
        if {[catch {exec {*}$pip_cmd} pip_result]} {
            ::PETK::gui::updateStatusDisplay "⚠️ Warning: Local package installation failed"
            ::PETK::gui::appendStatusDisplay "Error details:\n$pip_result"
            ::PETK::gui::appendStatusDisplay "You may need to run the install manually:\nconda run -n sem-dolfinx pip install -e $project_root"
            # Don't return here - the environment is still usable
        } else {
            ::PETK::gui::appendStatusDisplay "✅ Local package installed successfully from $project_root!"
        }
    }
    
    # Final setup
    ::PETK::gui::updateStatusDisplay "✅ DOLFINx environment 'sem-dolfinx' created successfully!"
    ::PETK::gui::appendStatusDisplay "Environment is ready for use with DOLFINx."
    ::PETK::gui::appendStatusDisplay "To activate: conda activate sem-dolfinx"
    
    set ::PETK::gui::condaEnvironment "sem-dolfinx"
    ::PETK::gui::refreshCondaEnvironments
    ::PETK::gui::testPythonEnvironment
    
    tk_messageBox -icon info -title "Environment Created" \
        -message "DOLFINx environment 'sem-dolfinx' created successfully!\n\nTo use the environment:\nconda activate sem-dolfinx"
}

####################################################
# Action button function
####################################################
proc ::PETK::gui::validateSEMSetup {} {
    set ::PETK::gui::semCurrentStatus "Validating setup..."
    update
    
    set errors {}
    set warnings {}
    set info_messages {}
    set calc_mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set calc_mode $::PETK::gui::semCalculationMode
    }
    
    # ========== FILE MANAGEMENT SECTION ==========
    # Ensure workdir exists
    if {![file exists $::PETK::gui::workdir]} {
        if {[catch {file mkdir $::PETK::gui::workdir} mkdir_error]} {
            lappend errors "Cannot create output directory: $::PETK::gui::workdir ($mkdir_error)"
        } else {
            lappend warnings "Created output directory: $::PETK::gui::workdir"
        }
    } elseif {![file isdirectory $::PETK::gui::workdir]} {
        lappend errors "Output path exists but is not a directory: $::PETK::gui::workdir"
    } elseif {![file writable $::PETK::gui::workdir]} {
        lappend errors "Output directory is not writable: $::PETK::gui::workdir"
    }
    
    # Check and copy analyte PDB file to workdir
    set analyte_in_workdir ""
    set analyte_source ""
    
    if {$::PETK::gui::semAnalytePDB eq ""} {
        lappend errors "No analyte PDB file specified"
    } else {
        # Try to locate the analyte PDB file
        if {[file exists $::PETK::gui::semAnalytePDB]} {
            # Full path exists
            set analyte_source $::PETK::gui::semAnalytePDB
        } elseif {[file exists [file join [pwd] $::PETK::gui::semAnalytePDB]]} {
            # Found in current directory
            set analyte_source [file join [pwd] $::PETK::gui::semAnalytePDB]
        } elseif {[file exists [file join $::PETK::gui::workdir $::PETK::gui::semAnalytePDB]]} {
            # Found in workdir
            set analyte_source [file join $::PETK::gui::workdir $::PETK::gui::semAnalytePDB]
        }
        
        if {$analyte_source eq ""} {
            lappend errors "Analyte PDB file not found: $::PETK::gui::semAnalytePDB (searched: current dir, workdir)"
        }
    }
    
    if {$analyte_source ne ""} {
        # Copy analyte PDB to workdir if not already there
        set analyte_filename [file tail $analyte_source]
        set analyte_in_workdir [file join $::PETK::gui::workdir $analyte_filename]
        
        if {![file exists $analyte_in_workdir] || [file mtime $analyte_source] > [file mtime $analyte_in_workdir]} {
            if {[catch {file copy -force $analyte_source $analyte_in_workdir} copy_error]} {
                lappend errors "Failed to copy analyte PDB to workdir: $copy_error"
                set analyte_in_workdir ""
            } else {
                lappend info_messages "Copied analyte PDB to workdir: $analyte_filename"
            }
        }
    }
    
    # Check biological pore file if pore type is biological
    set bio_pore_in_workdir ""
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "biological"} {
        # First, try to find biological pore file in pwd/bio_pore
        set bio_pore_source ""
        set bio_pore_filename ""
        
        if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
            # Delegate to the canonical resolver from tab 1, which checks (in order):
            #   1. selectedBioPore as an absolute path
            #   2. ::PETK::gui::bioPorePathMap (built when tab 1 scans bio_pore dirs)
            #   3. ::PETK::gui::getBioPoreSearchDirs — both [pwd]/bio_pore AND the
            #      install-asset path returned by ::PETK::gui::resourcePath bio_pore
            # Fall back to the legacy two-path search if (somehow) the helper is
            # missing — e.g. tab 1 hasn't been sourced or initialised.
            set resolved ""
            if {[info procs ::PETK::gui::resolveBioPoreFile] ne ""} {
                set resolved [::PETK::gui::resolveBioPoreFile $::PETK::gui::selectedBioPore]
            }
            if {$resolved ne "" && [file exists $resolved]} {
                set bio_pore_source $resolved
                set bio_pore_filename [file tail $resolved]
            } elseif {[file exists $::PETK::gui::selectedBioPore]} {
                set bio_pore_source $::PETK::gui::selectedBioPore
                set bio_pore_filename [file tail $::PETK::gui::selectedBioPore]
            } else {
                set pwd_bio_pore [file join [pwd] "bio_pore" $::PETK::gui::selectedBioPore]
                if {[file exists $pwd_bio_pore]} {
                    set bio_pore_source $pwd_bio_pore
                    set bio_pore_filename $::PETK::gui::selectedBioPore
                }
            }
        } elseif {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
            # Custom biological pore selected from Tab 1
            if {[file exists $::PETK::gui::customBioPorePath]} {
                set bio_pore_source $::PETK::gui::customBioPorePath
                set bio_pore_filename [file tail $::PETK::gui::customBioPorePath]
            } else {
                set custom_candidates [list \
                    [file join [pwd] $::PETK::gui::customBioPorePath] \
                    [file join [pwd] "bio_pore" [file tail $::PETK::gui::customBioPorePath]] \
                    [file join $::PETK::gui::workdir [file tail $::PETK::gui::customBioPorePath]]]
                foreach candidate $custom_candidates {
                    if {[file exists $candidate]} {
                        set bio_pore_source $candidate
                        set bio_pore_filename [file tail $candidate]
                        break
                    }
                }
            }
        }
        
        if {$bio_pore_source eq ""} {
            lappend errors "No biological pore PDB file found. Selected: $::PETK::gui::selectedBioPore, Custom: $::PETK::gui::customBioPorePath"
        } else {
            # Copy biological pore to workdir
            set bio_pore_in_workdir [file join $::PETK::gui::workdir $bio_pore_filename]
            
            if {![file exists $bio_pore_in_workdir] || [file mtime $bio_pore_source] > [file mtime $bio_pore_in_workdir]} {
                if {[catch {file copy -force $bio_pore_source $bio_pore_in_workdir} copy_error]} {
                    lappend errors "Failed to copy biological pore PDB to workdir: $copy_error"
                    set bio_pore_in_workdir ""
                } else {
                    lappend info_messages "Copied biological pore PDB to workdir: $bio_pore_filename"
                    # Update the selectedBioPore to point to workdir version
                    if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
                        set ::PETK::gui::selectedBioPore $bio_pore_in_workdir
                    }
                    if {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
                        set ::PETK::gui::customBioPorePath $bio_pore_in_workdir
                    }
                }
            } else {
                # File already exists and is up to date
                if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
                    set ::PETK::gui::selectedBioPore $bio_pore_in_workdir
                }
                if {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
                    set ::PETK::gui::customBioPorePath $bio_pore_in_workdir
                }
            }
        }
        
        # Validate biological pore file in workdir
        if {$bio_pore_in_workdir ne "" && [file exists $bio_pore_in_workdir]} {
            # Basic PDB validation
            if {[file size $bio_pore_in_workdir] < 100} {
                lappend warnings "Biological pore PDB file is very small (< 100 bytes)"
            }
        } else {
            lappend errors "Biological pore PDB file not available in workdir"
        }
        
        # Validate membrane Z offset for biological pores
        if {[info exists ::PETK::gui::membraneZOffset]} {
            if {![string is double $::PETK::gui::membraneZOffset]} {
                lappend errors "Membrane Z offset is not a valid number: $::PETK::gui::membraneZOffset"
            }
        } else {
            lappend warnings "Membrane Z offset not specified, will use default"
        }
    }
    
    # Update semAnalytePDB to point to workdir version
    if {$analyte_in_workdir ne "" && [file exists $analyte_in_workdir]} {
        set ::PETK::gui::semAnalytePDB $analyte_in_workdir
    }
    
    # Verify both files are now in workdir (for biological pore setup)
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption eq "biological"} {
        set files_in_workdir 0
        if {$analyte_in_workdir ne "" && [file exists $analyte_in_workdir]} {
            incr files_in_workdir
        }
        if {$bio_pore_in_workdir ne "" && [file exists $bio_pore_in_workdir]} {
            incr files_in_workdir
        }
        
        if {$files_in_workdir == 2} {
            lappend info_messages "✓ Both analyte and biological pore files are now in workdir"
        } elseif {$files_in_workdir == 1} {
            lappend warnings "⚠ Only one file successfully copied to workdir"
        } else {
            lappend errors "✗ Neither file is available in workdir"
        }
    }
    
    # ========== PARAMETER VALIDATION SECTION ==========
    # Check numeric parameters that must be positive
    set positive_params {
        {nanoporeThickness "Membrane thickness"}
        {appliedVoltage "Applied voltage"}
        {bulkConductivity "Bulk conductivity"}
        {sysPadding "Distance putoff"}
        {gridResolution "Grid resolution"}
    }
    if {$calc_mode eq "run" || $calc_mode eq "hybrid"} {
        lappend positive_params {zStep "Z step size"}
    }
    # Multigrid-mesh parameters must be positive when multigrid mesh is selected.
    set mesh_mode_validate "uniform"
    if {[info exists ::PETK::gui::meshMode]} {
        set mesh_mode_validate $::PETK::gui::meshMode
    }
    if {$mesh_mode_validate eq "multigrid"} {
        lappend positive_params {gmshFineSize "Multigrid fine cell size"}
        lappend positive_params {gmshCoarseSize "Multigrid coarse cell size"}
        set fb_mode_validate "auto"
        if {[info exists ::PETK::gui::gmshFineBoxMode]} {
            set fb_mode_validate $::PETK::gui::gmshFineBoxMode
        }
        if {$fb_mode_validate eq "manual"} {
            lappend positive_params {gmshFineBoxX "Multigrid fine box X"}
            lappend positive_params {gmshFineBoxY "Multigrid fine box Y"}
            lappend positive_params {gmshFineBoxZ "Multigrid fine box Z"}
        } else {
            lappend positive_params {gmshAutoMargin "Multigrid auto-margin"}
        }
        # Sanity: fine should not exceed coarse.
        if {[info exists ::PETK::gui::gmshFineSize] && \
                [info exists ::PETK::gui::gmshCoarseSize] && \
                [string is double -strict $::PETK::gui::gmshFineSize] && \
                [string is double -strict $::PETK::gui::gmshCoarseSize]} {
            if {$::PETK::gui::gmshFineSize > $::PETK::gui::gmshCoarseSize} {
                lappend errors \
                    "Multigrid fine cell size ($::PETK::gui::gmshFineSize Å) is larger than coarse cell size ($::PETK::gui::gmshCoarseSize Å); the relationship should be reversed."
            }
        }
    }
    
    # Add pore-specific positive parameters based on pore type
    if {[info exists ::PETK::gui::currentPoreType]} {
        if {$::PETK::gui::currentPoreType eq "cylindrical" || $::PETK::gui::currentPoreType eq "Cylindrical"} {
            lappend positive_params {cylindricalDiameter "Pore diameter"}
        } elseif {$::PETK::gui::currentPoreType eq "double_cone" || $::PETK::gui::currentPoreType eq "Double Cone"} {
            lappend positive_params {innerDiameter "Inner diameter"}
            lappend positive_params {outerDiameter "Outer diameter"}
        } elseif {$::PETK::gui::currentPoreType eq "conical" || $::PETK::gui::currentPoreType eq "Conical"} {
            lappend positive_params {topDiameter "Top diameter"}
            lappend positive_params {bottomDiameter "Bottom diameter"}
        }
        # For biological pores, diameter parameters are not needed as they come from the PDB
    }
    
    foreach param_info $positive_params {
        set var_name [lindex $param_info 0]
        set display_name [lindex $param_info 1]
        
        if {[info exists ::PETK::gui::$var_name]} {
            set var_value [set ::PETK::gui::$var_name]
            
            if {![string is double $var_value]} {
                lappend errors "$display_name is not a valid number: $var_value"
            } elseif {$var_value < 0} {
                lappend errors "$display_name must be positive: $var_value"
            }
        } else {
            lappend warnings "$display_name variable not found"
        }
    }

    # Check Z coordinate parameters (can be negative, just need to be valid numbers).
    # Required in any mode that performs translocation: "run" and "hybrid".
    if {$calc_mode eq "run" || $calc_mode eq "hybrid"} {
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

        if {[string is double $::PETK::gui::zStartRange] && [string is double $::PETK::gui::zEndRange] && [string is double $::PETK::gui::zStep]} {
            set z_range [expr {abs($::PETK::gui::zStartRange - $::PETK::gui::zEndRange)}]
            if {$z_range <= $::PETK::gui::zStep} {
                lappend errors "Z range (|Zend - Zstart| = $z_range) must be greater than Z step size ($::PETK::gui::zStep)"
            }
        }
    }
    # Rotation samples (and seed) required in any mode that samples rotations:
    # "rotation_scan" and "hybrid".
    if {$calc_mode eq "rotation_scan" || $calc_mode eq "hybrid"} {
        if {![string is integer -strict $::PETK::gui::rotationSamples] || $::PETK::gui::rotationSamples <= 0} {
            lappend errors "Rotation samples must be a positive integer"
        }
        if {[info exists ::PETK::gui::rotationSeed] && \
                ![string is integer -strict $::PETK::gui::rotationSeed]} {
            lappend errors "Random seed must be an integer (got: $::PETK::gui::rotationSeed)"
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
    
    # ========== RESULTS DISPLAY ==========
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

    if {[info exists info_messages] && [llength $info_messages] > 0} {
        append message "INFO:\n"
        foreach info_msg $info_messages {
            append message "• $info_msg\n"
        }
        append message "\n"
    }
    
    if {[llength $errors] == 0} {
        append message "✓ Setup validation PASSED\n"
        append message "Ready to run SEM calculations\n"
        set ::PETK::gui::semValidationPassed 1
        set ::PETK::gui::semCurrentStatus "Validation passed - ready to run"
        set icon "info"
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
    set use_vdw_json "false"
    set use_pdb2pqr_json "true"
    set sem_force_field "CHARMM"
    
    # Get relative path for PDB file (relative to output directory)
    set output_dir [file dirname $output_file]
    set pdb_path ""
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
    set config_calc_mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set config_calc_mode $::PETK::gui::semCalculationMode
    }
    append json_content "  \"metadata\": {\n"
    append json_content "    \"generated_by\": \"PETK GUI\",\n"
    append json_content "    \"timestamp\": \"[clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]\",\n"
    append json_content "    \"calculation_mode\": \"$config_calc_mode\",\n"
    append json_content "    \"version\": \"1.0\"\n"
    append json_content "  },\n"
    
    # Input section
    append json_content "  \"input\": {\n"
    append json_content "    \"moving_pdb\": \"$pdb_path\"\n"
    append json_content "  },\n"
    
    # Pore geometry section
    append json_content "  \"pore_geometry\": {\n"
    
    # Determine pore type - check both variables for consistency
    set pore_type_raw "solid-state"
    if {[info exists ::PETK::gui::poreOption] && $::PETK::gui::poreOption ne ""} {
        set pore_type_raw $::PETK::gui::poreOption
    }
    if {$pore_type_raw eq "solid-state" && [info exists ::PETK::gui::membraneType] && $::PETK::gui::membraneType ne ""} {
        set pore_type_raw $::PETK::gui::membraneType
    }
    set pore_type [string tolower $pore_type_raw]
    switch -- $pore_type {
        "doublecone" { set pore_type "double_cone" }
        "double_cone" { set pore_type "double_cone" }
        "conical" { set pore_type "conical" }
        "biological" - "bin_file" - "cylindrical" {}
        default { set pore_type "cylindrical" }
    }

    append json_content "    \"pore_type\": \"$pore_type\",\n"
    
    # Add pore-specific parameters
    if {$pore_type eq "cylindrical"} {
        set cylinder_diameter $::PETK::gui::cylindricalDiameter
        if {![string is double -strict $cylinder_diameter]} {
            set cylinder_diameter 100.0
        }
        set pore_radius [expr {$cylinder_diameter / 2.0}]
        append json_content "    \"pore_radius\": $pore_radius,\n"
        if {[info exists ::PETK::gui::cornerRadius] && [string is double -strict $::PETK::gui::cornerRadius]} {
            append json_content "    \"corner_radius\": $::PETK::gui::cornerRadius,\n"
        }
        append json_content "    \"membrane_thickness\": $::PETK::gui::nanoporeThickness\n"
        
    } elseif {$pore_type eq "double_cone"} {
        set inner_diameter $::PETK::gui::innerDiameter
        set outer_diameter $::PETK::gui::outerDiameter
        if {![string is double -strict $inner_diameter]} {
            set inner_diameter 100.0
        }
        if {![string is double -strict $outer_diameter]} {
            set outer_diameter 300.0
        }
        set inner_radius [expr {$inner_diameter / 2.0}]
        set outer_radius [expr {$outer_diameter / 2.0}]
        append json_content "    \"pore_radius\": $inner_radius,\n"
        append json_content "    \"outer_radius\": $outer_radius,\n"
        append json_content "    \"membrane_thickness\": $::PETK::gui::nanoporeThickness\n"

    } elseif {$pore_type eq "conical"} {
        set top_diameter $::PETK::gui::topDiameter
        set bottom_diameter $::PETK::gui::bottomDiameter
        if {![string is double -strict $top_diameter]} {
            set top_diameter 240.0
        }
        if {![string is double -strict $bottom_diameter]} {
            set bottom_diameter 120.0
        }
        set top_radius [expr {$top_diameter / 2.0}]
        set bottom_radius [expr {$bottom_diameter / 2.0}]
        append json_content "    \"top_radius\": $top_radius,\n"
        append json_content "    \"bottom_radius\": $bottom_radius,\n"
        append json_content "    \"membrane_thickness\": $::PETK::gui::nanoporeThickness\n"

    } elseif {$pore_type eq "biological"} {
        # Handle biological pore PDB path with robust file location logic
        set bio_pore_path ""
        
        set selected_path ""
        if {[info exists ::PETK::gui::selectedBioPore] && $::PETK::gui::selectedBioPore ne ""} {
            set selected_path $::PETK::gui::selectedBioPore
        } elseif {[info exists ::PETK::gui::customBioPorePath] && $::PETK::gui::customBioPorePath ne ""} {
            set selected_path $::PETK::gui::customBioPorePath
        }

        if {$selected_path ne ""} {
            # Check if selected_path points to a file in workdir (after validation)
            if {[file exists $selected_path]} {
                set bio_pore_filename [file tail $selected_path]
                set bio_pore_workdir_path [file join $::PETK::gui::workdir $bio_pore_filename]
                
                # If the selected_path is already in workdir, use relative path
                if {[file normalize $selected_path] eq [file normalize $bio_pore_workdir_path]} {
                    set bio_pore_path "../$bio_pore_filename"
                } else {
                    # If it's elsewhere, still use the filename (assume validation copied it)
                    set bio_pore_path "../$bio_pore_filename"
                }
            } else {
                # selected_path doesn't exist, try to find it in expected locations
                set bio_pore_filename [file tail $selected_path]
                
                # Check if it exists in workdir
                set workdir_path [file join $::PETK::gui::workdir $bio_pore_filename]
                if {[file exists $workdir_path]} {
                    set bio_pore_path "../$bio_pore_filename"
                } else {
                    # Check if it exists in pwd/bio_pore
                    set pwd_bio_pore [file join [pwd] "bio_pore" $bio_pore_filename]
                    if {[file exists $pwd_bio_pore]} {
                        set bio_pore_path "../bio_pore/$bio_pore_filename"
                    } else {
                        # Use the selected name as-is (validation should have caught missing files)
                        set bio_pore_path "../$bio_pore_filename"
                    }
                }
            }
        }
        
        append json_content "    \"biological_pore_pdb\": \"$bio_pore_path\",\n"
        append json_content "    \"membrane_thickness\": $::PETK::gui::nanoporeThickness,\n"
        
        # Add membrane Z offset with validation
        if {[info exists ::PETK::gui::membraneZOffset] && $::PETK::gui::membraneZOffset ne ""} {
            # Validate that it's a number (similar to validateSEMSetup)
            if {[string is double $::PETK::gui::membraneZOffset]} {
                append json_content "    \"membrane_z_offset\": $::PETK::gui::membraneZOffset\n"
            } else {
                # Use default if invalid
                append json_content "    \"membrane_z_offset\": 0.0\n"
            }
        } else {
            append json_content "    \"membrane_z_offset\": 0.0\n"
        }
    }
    
    append json_content "  },\n"
    
    # Simulation section
    append json_content "  \"simulation\": {\n"
    append json_content "    \"voltage\": $::PETK::gui::appliedVoltage,\n"
    append json_content "    \"bulk_conductivity\": $::PETK::gui::bulkConductivity,\n"
    append json_content "    \"grid_resolution\": $::PETK::gui::gridResolution,\n"
    append json_content "    \"use_vdw_radii\": $use_vdw_json,\n"
    append json_content "    \"use_pdb2pqr\": $use_pdb2pqr_json,\n"
    append json_content "    \"force_field\": \"$sem_force_field\",\n"
    append json_content "    \"default_radius\": $::PETK::gui::semDefaultRadius,\n"

    # If graded mesh is selected, append mesh_engine and gmsh_* fields. The
    # fine box is either explicit X/Y/Z or auto-computed per pore widest +
    # 2 × margin in XY and membrane_thickness + 2 × margin in Z.
    set mesh_mode_local "uniform"
    if {[info exists ::PETK::gui::meshMode]} {
        set mesh_mode_local $::PETK::gui::meshMode
    }
    if {$mesh_mode_local eq "multigrid"} {
        append json_content "    \"membrane_conductivity\": $::PETK::gui::membraneConductivity,\n"
        append json_content "    \"mesh_engine\": \"gmsh\",\n"
        append json_content "    \"gmsh_fine_size\": $::PETK::gui::gmshFineSize,\n"
        append json_content "    \"gmsh_coarse_size\": $::PETK::gui::gmshCoarseSize,\n"

        set fb_mode "auto"
        if {[info exists ::PETK::gui::gmshFineBoxMode]} {
            set fb_mode $::PETK::gui::gmshFineBoxMode
        }
        if {$fb_mode eq "manual"} {
            set fb_x $::PETK::gui::gmshFineBoxX
            set fb_y $::PETK::gui::gmshFineBoxY
            set fb_z $::PETK::gui::gmshFineBoxZ
        } else {
            # Auto-compute from pore widest + 2*margin and membrane + 2*margin.
            set margin $::PETK::gui::gmshAutoMargin
            set widest_diam_A 0.0
            if {[info exists ::PETK::gui::membraneType]} {
                if {$::PETK::gui::membraneType eq "cylindrical" && \
                        [info exists ::PETK::gui::cylindricalDiameter] && \
                        [string is double -strict $::PETK::gui::cylindricalDiameter]} {
                    set widest_diam_A $::PETK::gui::cylindricalDiameter
                } elseif {$::PETK::gui::membraneType eq "doublecone" && \
                        [info exists ::PETK::gui::outerDiameter] && \
                        [string is double -strict $::PETK::gui::outerDiameter]} {
                    set widest_diam_A $::PETK::gui::outerDiameter
                } elseif {$::PETK::gui::membraneType eq "conical"} {
                    set top_d 0.0
                    set bot_d 0.0
                    if {[info exists ::PETK::gui::topDiameter] && \
                            [string is double -strict $::PETK::gui::topDiameter]} {
                        set top_d $::PETK::gui::topDiameter
                    }
                    if {[info exists ::PETK::gui::bottomDiameter] && \
                            [string is double -strict $::PETK::gui::bottomDiameter]} {
                        set bot_d $::PETK::gui::bottomDiameter
                    }
                    set widest_diam_A [expr {max($top_d, $bot_d)}]
                }
            }
            set fb_x [expr {$widest_diam_A + 2.0 * $margin}]
            set fb_y $fb_x
            set membrane_A 0.0
            if {[info exists ::PETK::gui::nanoporeThickness] && \
                    [string is double -strict $::PETK::gui::nanoporeThickness]} {
                set membrane_A $::PETK::gui::nanoporeThickness
            }
            set fb_z [expr {$membrane_A + 2.0 * $margin}]
        }
        append json_content "    \"gmsh_fine_box\": \[$fb_x, $fb_y, $fb_z\]\n"
    } else {
        append json_content "    \"membrane_conductivity\": $::PETK::gui::membraneConductivity\n"
    }

    append json_content "  },\n"
    
    # Movement section
    set movement_z_start $::PETK::gui::zStartRange
    set movement_z_end $::PETK::gui::zEndRange
    if {[info exists ::PETK::gui::semCalculationMode] && $::PETK::gui::semCalculationMode eq "rotation_scan"} {
        # Pure rotation scan evaluates at a single z (set both bounds to 0); hybrid
        # mode keeps the user-configured z-sweep so each rotation gets a full pass.
        set movement_z_start 0.0
        set movement_z_end 0.0
    }
    append json_content "  \"movement\": {\n"
    append json_content "    \"z_start\": $movement_z_start,\n"
    append json_content "    \"z_end\": $movement_z_end,\n"
    append json_content "    \"z_step\": $::PETK::gui::zStep\n"
    append json_content "  },\n"
    
    # Output section
    append json_content "  \"output\": {\n"
    append json_content "    \"output_prefix\": \"$::PETK::gui::outputPrefix\",\n"

    # ARBD export block (optional). Emitted only when semArbdEnabled is 1.
    set arbd_enabled 0
    if {[info exists ::PETK::gui::semArbdEnabled] && \
            [string is true -strict $::PETK::gui::semArbdEnabled]} {
        set arbd_enabled 1
    } elseif {[info exists ::PETK::gui::semArbdEnabled] && \
              [string is integer -strict $::PETK::gui::semArbdEnabled] && \
              $::PETK::gui::semArbdEnabled > 0} {
        set arbd_enabled 1
    }
    if {$arbd_enabled} {
        # Parse the ions string "POT:1, CLA:-1" into JSON [["POT",1],["CLA",-1]]
        set ions_json_parts [list]
        foreach pair [split $::PETK::gui::semArbdIons ","] {
            set pair [string trim $pair]
            if {$pair eq ""} continue
            set kv [split $pair ":"]
            if {[llength $kv] != 2} continue
            set name [string trim [lindex $kv 0]]
            set valence [string trim [lindex $kv 1]]
            if {$name eq "" || ![string is double -strict $valence]} continue
            lappend ions_json_parts "\[\"$name\", $valence\]"
        }
        set ions_json "\[[join $ions_json_parts ", "]\]"

        set arbd_stride $::PETK::gui::semArbdStride
        if {![string is integer -strict $arbd_stride]} { set arbd_stride 0 }
        set arbd_wall $::PETK::gui::semArbdWallHeight
        if {![string is double -strict $arbd_wall]} { set arbd_wall 100.0 }
        set arbd_temp $::PETK::gui::semArbdTemperature
        if {![string is double -strict $arbd_temp]} { set arbd_temp 295.0 }

        append json_content "    \"preview_frames\": $::PETK::gui::semPreviewFrames,\n"
        append json_content "    \"arbd_export\": {\n"
        append json_content "      \"ions\": $ions_json,\n"
        append json_content "      \"stride\": $arbd_stride,\n"
        append json_content "      \"wall_height\": $arbd_wall,\n"
        append json_content "      \"temperature_K\": $arbd_temp\n"
        append json_content "    }\n"
    } else {
        append json_content "    \"preview_frames\": $::PETK::gui::semPreviewFrames\n"
    }
    append json_content "  },\n"
    
    # Box dimensions section
    append json_content "  \"box_dimensions\": {\n"
    
    set use_auto 0
    if {[info exists ::PETK::gui::autoCalculateBoxDimensions] && $::PETK::gui::autoCalculateBoxDimensions} {
        set use_auto 1
    }
    
    if {$use_auto} {
        # Auto-calculated mode - parse range strings
        set x_range [parseBoxRange $::PETK::gui::autoBoxX]
        set y_range [parseBoxRange $::PETK::gui::autoBoxY] 
        set z_range [parseBoxRange $::PETK::gui::autoBoxZ]
        
        append json_content "    \"x\": $x_range,\n"
        append json_content "    \"y\": $y_range,\n"
        append json_content "    \"z\": $z_range\n"
    } else {
        # Manual mode - convert sizes to ranges (assume centered at origin)
        set box_x $::PETK::gui::boxSizeX
        set box_y $::PETK::gui::boxSizeY
        set box_z $::PETK::gui::boxSizeZ
        if {![string is double -strict $box_x]} {set box_x 150.0}
        if {![string is double -strict $box_y]} {set box_y 150.0}
        if {![string is double -strict $box_z]} {set box_z 150.0}
        set x_half [expr {$box_x / 2.0}]
        set y_half [expr {$box_y / 2.0}]
        set z_half [expr {$box_z / 2.0}]
        
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
    set mode "run"
    if {[info exists ::PETK::gui::semCalculationMode] && $::PETK::gui::semCalculationMode ne ""} {
        set mode $::PETK::gui::semCalculationMode
    }
    
    switch -- $mode {
        "run" {
            ::PETK::gui::runVerticalPreview
        }
        "rotation_scan" {
            ::PETK::gui::runRotationPreview
        }
        "hybrid" {
            # Hybrid combines a Z-sweep with rotational sampling. The most
            # informative single preview is the vertical trajectory; the
            # rotational sampling is exercised in the actual run.
            ::PETK::gui::runVerticalPreview
        }
        default {
            tk_messageBox -icon info -title "Preview Unavailable" \
                -message "Preview simulation is not available for calculation mode \"$mode\"."
        }
    }
}

proc ::PETK::gui::runVerticalPreview {} {
    # Run a preview simulation by generating a DCD trajectory showing
    # the analyte moving from zStartRange to zEndRange.
    # Uses corrected VMD animate syntax.
    
    puts "=== SEM PREVIEW SIMULATION (VERTICAL) ==="
    
    # Use the stored analyte molecule ID
    set analyte_molid [::PETK::gui::getAnalyteMoleculeID]
    if {$analyte_molid < 0} {
        tk_messageBox -icon error -title "Preview Error" \
            -message "Unable to locate the analyte molecule for preview."
        return
    }
    
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

proc ::PETK::gui::runRotationPreview {} {
    puts "=== SEM PREVIEW SIMULATION (ROTATION) ==="
    
    set analyte_molid [::PETK::gui::getAnalyteMoleculeID]
    if {$analyte_molid < 0} {
        tk_messageBox -icon error -title "Preview Error" \
            -message "Unable to locate the analyte molecule for preview."
        return
    }
    
    set total_samples 10
    if {[info exists ::PETK::gui::rotationSamples] && \
            [string is integer -strict $::PETK::gui::rotationSamples] && \
            $::PETK::gui::rotationSamples > 0} {
        set total_samples [expr {int($::PETK::gui::rotationSamples)}]
    }
    
    set preview_cap 8
    if {[info exists ::PETK::gui::semPreviewFrames] && \
            [string is integer -strict $::PETK::gui::semPreviewFrames] && \
            $::PETK::gui::semPreviewFrames > 0} {
        set preview_cap [expr {int($::PETK::gui::semPreviewFrames)}]
    }
    
    set num_frames [expr {min($total_samples, $preview_cap)}]
    if {$num_frames < 3} {
        set num_frames [expr {min($total_samples, 3)}]
    }
    if {$num_frames < 1} {
        tk_messageBox -icon error -title "Preview Error" \
            -message "Rotation preview requires at least one sample. Increase the rotation samples value."
        return
    }
    
    puts "Generating rotation preview using $num_frames orientation frames (of $total_samples requested samples)"
    
    if {[::PETK::gui::generateRotationTrajectory $analyte_molid $num_frames]} {
        ::PETK::gui::loadAndPlaySEMTrajectory $analyte_molid
        mol rename $analyte_molid "Analyte: [file tail $::PETK::gui::analytePDB] (Rotation Preview)"
        
        tk_messageBox -type ok -icon info -title "Rotation Preview Ready" \
            -message "Rotation preview trajectory generated successfully!\n\nFrames: $num_frames\nSamples requested: $total_samples\n\nPlayback now covers a smooth sweep of the sampled orientations."
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

proc ::PETK::gui::generateRotationTrajectory {molid num_frames} {
    # Generate a trajectory that sweeps the analyte through several orientations.
    #
    # Args:
    #   molid: VMD molecule ID
    #   num_frames: Number of preview frames/orientations to generate
    
    if {$num_frames < 1} {
        return 0
    }
    
    set success 0
    if {[catch {
        set output_dir [file join $::PETK::gui::workdir "sem_preview"]
        file mkdir $output_dir
        
        set dcd_file [file join $output_dir "rotation_preview.dcd"]
        set pdb_file [file join $output_dir "rotation_preview.pdb"]
        
        puts "Rotation preview files:"
        puts "  PDB: $pdb_file"
        puts "  DCD: $dcd_file"
        
        set all_atoms [atomselect $molid "all"]
        set original_coords [$all_atoms get {x y z}]
        set com [measure center $all_atoms]
        set comx [lindex $com 0]
        set comy [lindex $com 1]
        set comz [lindex $com 2]
        set num_atoms [$all_atoms num]
        
        puts "Analyte center of mass: ([format %.3f $comx], [format %.3f $comy], [format %.3f $comz])"
        puts "Atoms in selection: $num_atoms"
        
        set pi [expr {acos(-1)}]
        set denom [expr {$num_frames > 1 ? double($num_frames - 1) : 1.0}]
        
        set current_frames [molinfo $molid get numframes]
        if {$current_frames > 1} {
            animate delete beg 1 end -1 $molid
        }
        
        for {set frame 0} {$frame < $num_frames} {incr frame} {
            if {$frame > 0} {
                animate dup $molid
            }
            animate goto $frame
            
            set fraction [expr {$frame / $denom}]
            set yaw [expr {$fraction * 2.0 * $pi}]
            set pitch [expr {0.5 * sin($fraction * 2.0 * $pi)}]
            set roll [expr {0.7 * cos($fraction * 2.0 * $pi)}]
            
            set cos_yaw [expr {cos($yaw)}]
            set sin_yaw [expr {sin($yaw)}]
            set cos_pitch [expr {cos($pitch)}]
            set sin_pitch [expr {sin($pitch)}]
            set cos_roll [expr {cos($roll)}]
            set sin_roll [expr {sin($roll)}]
            
            set new_coords {}
            foreach coord $original_coords {
                set x [lindex $coord 0]
                set y [lindex $coord 1]
                set z [lindex $coord 2]
                
                # Translate to COM
                set dx [expr {$x - $comx}]
                set dy [expr {$y - $comy}]
                set dz [expr {$z - $comz}]
                
                # Apply yaw (Z rotation)
                set x1 [expr {$dx * $cos_yaw - $dy * $sin_yaw}]
                set y1 [expr {$dx * $sin_yaw + $dy * $cos_yaw}]
                set z1 $dz
                
                # Apply pitch (X rotation)
                set y2 [expr {$y1 * $cos_pitch - $z1 * $sin_pitch}]
                set z2 [expr {$y1 * $sin_pitch + $z1 * $cos_pitch}]
                set x2 $x1
                
                # Apply roll (Y rotation)
                set x3 [expr {$x2 * $cos_roll + $z2 * $sin_roll}]
                set z3 [expr {$z2 * $cos_roll - $x2 * $sin_roll}]
                set y3 $y2
                
                lappend new_coords [list [expr {$comx + $x3}] [expr {$comy + $y3}] [expr {$comz + $z3}]]
            }
            
            set frame_atoms [atomselect $molid "all" frame $frame]
            $frame_atoms set {x y z} $new_coords
            $frame_atoms delete
            
            set yaw_deg [format %.1f [expr {$yaw * 180.0 / $pi}]]
            set pitch_deg [format %.1f [expr {$pitch * 180.0 / $pi}]]
            set roll_deg [format %.1f [expr {$roll * 180.0 / $pi}]]
            puts "  Frame [expr {$frame + 1}]/$num_frames => yaw ${yaw_deg}°, pitch ${pitch_deg}°, roll ${roll_deg}°"
        }
        
        puts "Writing rotation preview trajectory..."
        animate goto 0
        $all_atoms writepdb $pdb_file
        animate write dcd $dcd_file beg 0 end -1 $molid
        $all_atoms delete
        
        set ::PETK::gui::currentTrajectoryDCD $dcd_file
        set ::PETK::gui::currentTrajectoryPDB $pdb_file
        set ::PETK::gui::currentTrajectoryFrames $num_frames
        
        set success 1
        
    } error]} {
        puts "ERROR generating rotation trajectory: $error"
        tk_messageBox -type ok -icon error -title "Rotation Preview Failed" \
            -message "Failed to generate rotation preview trajectory:\n\n$error"
        set success 0
    }
    
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
    
    set calc_mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set calc_mode $::PETK::gui::semCalculationMode
    }

    # Run calculations
    set config_file [file join $::PETK::gui::workdir "config.json"]
    ::PETK::gui::run_sem_calculation $config_file
}

# Main procedure for running SEM calculations
proc ::PETK::gui::run_sem_calculation {config_file} {
    # Ensure validation passed before running
    if {![info exists ::PETK::gui::semValidationPassed] || $::PETK::gui::semValidationPassed != 1} {
        set error_msg "Cannot run calculations: SEM setup validation has not passed.\nPlease run validation first and fix any errors."
        puts $error_msg
        tk_messageBox -icon error -title "Calculation Error" -message $error_msg
        return -code error $error_msg
    }

    # Prepare results directory
    set results_dir [file join $::PETK::gui::workdir "results"]
    if {![file exists $results_dir]} {
        file mkdir $results_dir
        puts "Created results directory: $results_dir"
    } else {
        puts "Results directory already exists: $results_dir"
    }

    # Copy config into results directory for the CLI run
    set config_basename [file tail $config_file]
    set local_config [file join $results_dir $config_basename]
    if {![file exists $local_config] || [file mtime $config_file] > [file mtime $local_config]} {
        file copy -force $config_file $local_config
        puts "Copied config file to results directory"
    }

    # Build command
    set calc_mode "run"
    if {[info exists ::PETK::gui::semCalculationMode]} {
        set calc_mode $::PETK::gui::semCalculationMode
    }
    # Prepare base sem command
    set env_name $::PETK::gui::condaEnvironment
    set base_cmd [list conda run -n $::PETK::gui::condaEnvironment python -m sem]
    if {$calc_mode eq "rotation_scan" || $calc_mode eq "hybrid"} {
        set scan_label [expr {$calc_mode eq "hybrid" ? "hybrid" : "rotation scan"}]
        if {![string is integer -strict $::PETK::gui::rotationSamples] || $::PETK::gui::rotationSamples <= 0} {
            set msg "Rotation samples must be a positive integer before running ${scan_label}."
            puts $msg
            tk_messageBox -icon error -title "Run Error" -message $msg
            return -code error $msg
        }
        set rotation_samples [expr {int($::PETK::gui::rotationSamples)}]
        set rotation_seed 42
        if {[info exists ::PETK::gui::rotationSeed] && \
                [string is integer -strict $::PETK::gui::rotationSeed]} {
            set rotation_seed [expr {int($::PETK::gui::rotationSeed)}]
        }
        set final_cmd [concat $base_cmd [list rotation_scan $config_basename \
            --samples $rotation_samples --seed $rotation_seed --mode run]]
        if {[info exists ::PETK::gui::reuseSemMesh] && $::PETK::gui::reuseSemMesh} {
            lappend final_cmd --reuse-mesh
        }
    } else {
        set final_cmd [concat $base_cmd [list run $config_basename]]
    }

    ::PETK::gui::updateStatusDisplay "Running SEM command... this may take a while. See sem_run.log for details."

    # Run from inside results directory and capture output
    set log_file [file join $results_dir "sem_run.log"]
    if {[file exists $log_file]} {
        file delete -force $log_file
    }

    set original_dir [pwd]
    cd $results_dir

    set log_chan ""
    if {[catch {
        set log_chan [open $log_file w]
        fconfigure $log_chan -encoding utf-8 -translation lf
        puts $log_chan "Command: [join $final_cmd { }]"
        puts $log_chan "Started: [clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]"
        flush $log_chan
        exec {*}$final_cmd >@ $log_chan 2>@ $log_chan
        puts $log_chan "Finished: [clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]"
        close $log_chan
    } run_error]} {
        catch {close $log_chan}
        cd $original_dir
        set message "Error running SEM calculation:\n$run_error\n\nCheck $log_file for details."
        puts $message
        tk_messageBox -icon error -title "Calculation Error" -message $message
        return -code error $run_error
    }

    cd $original_dir

    # Report generated files
    set pattern_types {"*.png" "*.dat" "*.txt" "*.csv" "*.out" "*.h5" "*.xdmf"}
    set result_files {}
    foreach pattern $pattern_types {
        foreach file [glob -nocomplain -directory $results_dir $pattern] {
            lappend result_files $file
        }
    }
    if {[llength $result_files] > 0} {
        puts "Generated files:"
        foreach file $result_files {
            puts "  - $file"
        }
    } else {
        puts "No result files found"
    }

    set message "SEM calculation completed successfully.\nResults saved in:\n$results_dir\n\nLog file: $log_file"
    tk_messageBox -icon info -title "Calculation Complete" -message $message
    return $results_dir
}

# ===========================
# Plotting enhancements (conductance difference & rotation histogram)
# ===========================

proc ::PETK::gui::applyPlotEnhancements {} {
    if {![namespace exists ::PETK::gui::plot]} {
        after 100 ::PETK::gui::applyPlotEnhancements
        return
    }

    proc ::PETK::gui::plot::readVoltageFromConfig {config_path} {
        if {$config_path eq "" || ![file exists $config_path]} {
            return ""
        }
        set fh [open $config_path r]
        set content [read $fh]
        close $fh
        if {[regexp {"voltage"\s*:\s*([-0-9.eE\+]+)} $content -> val]} {
            if {[catch {expr {double($val) / 1000.0}} conv]} {
                return ""
            }
            return $conv
        }
        return ""
    }

    proc ::PETK::gui::plot::parseResultsSummary {results_file} {
        set summary [dict create]
        if {$results_file eq "" || ![file exists $results_file]} {
            return $summary
        }

        set fh [open $results_file r]
        set currents {}
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
            set fields [split $line]
            if {[llength $fields] >= 2} {
                if {![catch {expr {double([lindex $fields 1])}} cur]} {
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
        return $summary
    }

    proc ::PETK::gui::plot::readOutputPrefixFromConfig {config_path} {
        if {$config_path eq "" || ![file exists $config_path]} {
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
        set metadata [dict create open_current "" voltage_v "" results_file "" config_path "" output_prefix $prefix]

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
            if {$prefix eq "" || $prefix eq "vertical_movement"} {
                set derived [::PETK::gui::plot::readOutputPrefixFromConfig $config_path]
                if {$derived ne ""} {
                    dict set metadata output_prefix $derived
                }
            }
        }

        return $metadata
    }

    proc ::PETK::gui::plot::parseFile {filename open_current voltage_v} {
        if {![file exists $filename]} {
            return {}
        }

        set fh [open $filename r]
        set metadata {}

        while {[gets $fh line] >= 0} {
            set line [string trim $line]
            if {[string match "#*" $line]} {
                if {[regexp {# Position (\d+)/(\d+)} $line -> pos total]} {
                    dict set metadata position $pos
                    dict set metadata total_positions $total
                } elseif {[regexp {# Z_position: ([\-0-9.eE\+]+)} $line -> zpos]} {
                    dict set metadata z_position $zpos
                } elseif {[regexp {# Current: ([\-0-9.eE\+]+)} $line -> current]} {
                    dict set metadata current $current
                } elseif {[regexp {# Blockage: ([\-0-9.eE\+]+)} $line -> blockage]} {
                    dict set metadata blockage $blockage
                }
            } elseif {[llength [split $line]] >= 2} {
                set values [split $line]
                dict set metadata z_position [lindex $values 0]
                dict set metadata current [lindex $values 1]
                if {[llength $values] >= 3} {
                    dict set metadata blockage [lindex $values 2]
                }
                break
            }
        }
        close $fh

        if {[dict exists $metadata current] && $open_current ne "" && $voltage_v ne ""} {
            if {![catch {expr {double($voltage_v)}}]} {
                set current_val [expr {double([dict get $metadata current])}]
                if {$voltage_v != 0.0} {
                    # Data files report currents in nA; convert to nS by dividing by V
                    set delta_g_ns [expr {(double($open_current) - $current_val) / double($voltage_v)}]
                    dict set metadata delta_g $delta_g_ns
                }
            }
        }

        return $metadata
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

    proc ::PETK::gui::plot::parseDirectory {directory prefix} {
        set alldata {}
        set pattern [file join $directory "${prefix}_*.dat"]
        set files [lsort -dictionary [glob -nocomplain $pattern]]

        set meta [::PETK::gui::plot::collectRunMetadata $directory $prefix]
        set open_current [dict get $meta open_current]
        set voltage_v [dict get $meta voltage_v]
        set results_file ""
        if {[dict exists $meta results_file]} {
            set results_file [dict get $meta results_file]
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

    proc ::PETK::gui::plot::extractPlotArrays {alldata} {
        set zpos_array {}
        set current_array {}
        set delta_g_array {}

        foreach entry $alldata {
            if {[dict exists $entry z_position]} {
                lappend zpos_array [dict get $entry z_position]
            }
            if {[dict exists $entry current]} {
                lappend current_array [dict get $entry current]
            }
            if {[dict exists $entry delta_g]} {
                lappend delta_g_array [dict get $entry delta_g]
            }
        }
        return [list $zpos_array $current_array $delta_g_array]
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

        set data_arrays [::PETK::gui::plot::extractPlotArrays $alldata]
        set zpos_array [lindex $data_arrays 0]
        set current_array [lindex $data_arrays 1]
        set delta_g_array [lindex $data_arrays 2]

        if {[llength $zpos_array] == 0} {
            tk_messageBox -icon info -title "Plot Data" -message "Parsed files did not contain Z-position data."
            return ""
        }

        set entry0 [lindex $alldata 0]
        set open_current_meta ""
        if {[dict exists $entry0 open_current]} {
            set open_current_meta [dict get $entry0 open_current]
        }

        if {$openPoreCurrentMode eq "auto" && $open_current_meta ne ""} {
            set openPoreCurrent $open_current_meta
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
            if {![dict exists $summary open_current] || ![dict exists $summary currents]} {
                continue
            }

            set open_current [dict get $summary open_current]
            set currents [dict get $summary currents]
            if {[llength $currents] == 0} {
                continue
            }
            set measured_current [lindex $currents end]
            set delta_ns [expr {(double($open_current) - double($measured_current)) / double($voltage_v)}]
            set entry [dict create label [file tail $dir] delta $delta_ns results_file $results_path voltage_v $voltage_v open_current $open_current config_path $config_path]
            lappend entries $entry
        }
        return $entries
    }

    # Render hybrid_currents.csv (full z × rotation trace) as a single
    # Normalized-Current-vs-z plot. At each z, the mean over rotations is the
    # marker; asymmetric error bars span [min, max] across rotations (i.e.
    # lower = mean - min, upper = max - mean). Rotations dropped by the overlap
    # check contribute NaNs and are excluded via np.nanmean / np.nanmin /
    # np.nanmax. Also writes hybrid_summary.csv with per-z statistics (mean,
    # std, min, max, n_rotations) on normalized current for downstream analysis.
    proc ::PETK::gui::plot::plotHybridCurrentsCsv {csv_path} {
        if {![file exists $csv_path]} {
            tk_messageBox -icon error -title "Plot Error" \
                -message "Hybrid CSV not found: $csv_path"
            return
        }
        set out_dir [file dirname $csv_path]
        set output_path [file join $out_dir "hybrid_normcurrent_minmax.png"]
        set summary_path [file join $out_dir "hybrid_summary.csv"]
        set output_path_py  [string map {"\\" "\\\\"} $output_path]
        set summary_path_py [string map {"\\" "\\\\"} $summary_path]
        set csv_path_py     [string map {"\\" "\\\\"} $csv_path]
        set python_script_path [file join $out_dir "plot_hybrid_currents.py"]

        set python_code [format {import csv
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

csv_path = r"%s"
output_path = r"%s"
summary_path = r"%s"

# ---- Load (idx, z) -> normalized_current. Group by rotation index. ----------
by_rot = defaultdict(dict)   # idx -> {z: normalized_current}
rot_meta = {}                # idx -> (rx, ry, rz)
all_z = set()
with open(csv_path, "r", newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        try:
            idx = int(row["index"])
            z = float(row["z_A"])
            norm = float(row["normalized_current"])
        except (KeyError, ValueError):
            continue
        by_rot[idx][z] = norm
        all_z.add(z)
        rot_meta[idx] = (row.get("rx", ""), row.get("ry", ""), row.get("rz", ""))

if not by_rot or not all_z:
    print("No usable rows in", csv_path)
else:
    sorted_idx = sorted(by_rot.keys())
    z_grid = np.array(sorted(all_z))
    n_rot = len(sorted_idx)
    n_z = len(z_grid)

    # ---- Build (n_rot, n_z) matrix; NaN where a (z, rot) cell is missing. ----
    mat = np.full((n_rot, n_z), np.nan)
    for i, idx in enumerate(sorted_idx):
        for j, z in enumerate(z_grid):
            if z in by_rot[idx]:
                mat[i, j] = by_rot[idx][z]

    # NaN-safe stats across rotations at each z (skipped rotations contribute
    # NaN and are excluded).
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_norm = np.nanmean(mat, axis=0)
        std_norm  = np.nanstd(mat,  axis=0)
        min_norm  = np.nanmin(mat, axis=0)
        max_norm  = np.nanmax(mat, axis=0)
        n_valid   = np.sum(~np.isnan(mat), axis=0)

    # ---- Errorbar plot: mean marker + asymmetric min/max error bars. --------
    # yerr lower = mean - min, upper = max - mean. Clamp tiny negatives that
    # can arise from float jitter when min == mean.
    lower = mean_norm - min_norm
    upper = max_norm - mean_norm
    lower = np.where(lower < 0, 0.0, lower)
    upper = np.where(upper < 0, 0.0, upper)
    yerr = np.vstack([lower, upper])
    n_used = int(np.nanmax(n_valid)) if np.any(n_valid > 0) else 0

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.errorbar(
        z_grid, mean_norm, yerr=yerr,
        fmt="o-",
        capsize=3,
        elinewidth=1.5,
        markersize=4,
        color="C0",
        label=f"mean (n={n_used} rotations)",
    )
    ax.set_xlabel("Z position (A)", fontsize=14, weight="bold")
    ax.set_ylabel("Normalized Current", fontsize=14, weight="bold")
    ax.set_title(
        f"Normalized Current with min-max error bars over {n_used} rotations",
        fontsize=14, weight="bold",
    )
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    print(f"Saved {output_path}")

    # ---- Per-z summary CSV for downstream analysis. -------------------------
    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z_A", "mean_norm_current", "std_norm_current",
                    "min_norm_current", "max_norm_current", "n_rotations"])
        for k, z in enumerate(z_grid):
            row = [
                f"{z:.4f}",
                f"{mean_norm[k]:.6e}" if np.isfinite(mean_norm[k]) else "",
                f"{std_norm[k]:.6e}"  if np.isfinite(std_norm[k])  else "",
                f"{min_norm[k]:.6e}"  if np.isfinite(min_norm[k])  else "",
                f"{max_norm[k]:.6e}"  if np.isfinite(max_norm[k])  else "",
                int(n_valid[k]),
            ]
            w.writerow(row)
    print(f"Saved {summary_path}")

    plt.show()
} $csv_path_py $output_path_py $summary_path_py]

        if {[catch {
            set fh [open $python_script_path w]
            puts $fh $python_code
            close $fh
        } write_error]} {
            tk_messageBox -icon error -title "Plot Error" \
                -message "Failed to write hybrid plot script:\n$write_error"
            return
        }

        set python_exec "python3"
        if {[info exists ::PETK::gui::pythonExecutable] && $::PETK::gui::pythonExecutable ne ""} {
            set python_exec $::PETK::gui::pythonExecutable
        }
        if {[catch {set plot_output [exec $python_exec $python_script_path]} err]} {
            tk_messageBox -icon error -title "Plot Error" \
                -message "Hybrid plot failed:\n$err"
        } else {
            puts $plot_output
        }
    }

    proc ::PETK::gui::plot::plotRotationHistogram {} {
    set candidate_dirs {}
    if {[info exists ::PETK::gui::workdir] && $::PETK::gui::workdir ne ""} {
        set workdir [file normalize $::PETK::gui::workdir]
    } else {
        set workdir [pwd]
    }

    set results_root [file join $workdir "results"]
    lappend candidate_dirs $results_root
    lappend candidate_dirs [file join $results_root "rotations"]
    lappend candidate_dirs [file join $workdir "rotation" "results"]
    lappend candidate_dirs [file join $workdir "rotation" "results" "rotations"]

    set install_root [::PETK::gui::resourcePath rotation results]
    lappend candidate_dirs $install_root
    lappend candidate_dirs [file join $install_root "rotations"]

    # Normalize and deduplicate candidates
    set norm_candidates {}
    foreach dir $candidate_dirs {
        if {$dir eq ""} {continue}
        set norm_dir [file normalize $dir]
        if {[lsearch -exact $norm_candidates $norm_dir] == -1} {
            lappend norm_candidates $norm_dir
        }
    }

    # First pass: detect hybrid output (full z × rotation trace). If found, render
    # a current-vs-z overlay (one curve per rotation) and return early — this
    # plot supersedes the histogram for hybrid runs.
    foreach dir $norm_candidates {
        if {![file isdirectory $dir]} { continue }
        set hybrid_csv [file join $dir "hybrid_currents.csv"]
        if {[file exists $hybrid_csv]} {
            ::PETK::gui::plot::plotHybridCurrentsCsv $hybrid_csv
            return
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
        set searched [join $norm_candidates "\n  • "]
        tk_messageBox -icon info -title "Plot Rotation Scan" \
            -message "No rotation scan data found in any of the expected locations:\n  • $searched\nRun a rotation scan or hybrid run first, or copy the results into one of these directories."
        return
    }

    set deltas {}
    set labels {}
    set ref_current ""
    set voltage_v ""
    foreach entry $entries {
        lappend deltas [dict get $entry delta]
        lappend labels [dict get $entry label]
        if {$ref_current eq "" && [dict exists $entry open_current]} {
            set ref_current [dict get $entry open_current]
        }
        if {$voltage_v eq "" && [dict exists $entry voltage_v]} {
            set voltage_v [dict get $entry voltage_v]
        }
    }

    if {$voltage_v eq "" || $voltage_v == 0.0} {
        tk_messageBox -icon error -title "Plot Error" -message "Could not determine voltage from rotation configurations."
        return
    }

    set count [llength $deltas]
    set sum 0.0
    foreach delta $deltas { set sum [expr {$sum + $delta}] }
    set mean [expr {$sum / $count}]
    set variance 0.0
    foreach delta $deltas { set variance [expr {$variance + pow(($delta - $mean), 2)}] }
    set stddev [expr {sqrt($variance / $count)}]
    set min_delta [::tcl::mathfunc::min {*}$deltas]
    set max_delta [::tcl::mathfunc::max {*}$deltas]

    set csv_path [file join $base_dir "rotation_delta_g.csv"]
    if {[catch {
        set csv_fh [open $csv_path w]
        puts $csv_fh "label,delta_G_nS"
        foreach entry $entries {
            puts $csv_fh "[dict get $entry label],[format %.6f [dict get $entry delta]]"
        }
        close $csv_fh
    } csv_error]} {
        puts "Warning: Could not write $csv_path: $csv_error"
    }

    set labels_formatted {}
    foreach label $labels {
        set escaped [string map {"'" "\\'"} $label]
        lappend labels_formatted "'$escaped'"
    }
    set labels_py [join $labels_formatted ", "]

    set deltas_formatted {}
    foreach delta $deltas {
        lappend deltas_formatted [format "%.6f" $delta]
    }
    set deltas_py [join $deltas_formatted ", "]

    set ref_current_val [format "%.6f" $ref_current]
    set voltage_val [format "%.6f" $voltage_v]
    set plot_title [format "Histogram of Conductance Differences (%s)" [file tail $::PETK::gui::workdir]]
    set plot_title_escaped [string map {"\"" "\\\""} $plot_title]
    set output_path [file join $base_dir "conductance_diff_histogram.png"]
    set output_path_py [string map {"\\" "\\\\"} $output_path]
    set python_script_path [file join $base_dir "plot_rotation_hist.py"]

    set python_code [format {import os
import matplotlib.pyplot as plt
import numpy as np

labels = [%s]
deltas = [%s]
ref_current = %s
voltage_v = %s
output_path = r"%s"
title = "%s"

if not deltas:
    print("No data found.")
else:
    print("Reference current (nA):", ref_current)
    print("Voltage (V):", voltage_v)
    for lbl, delta in zip(labels, deltas):
        print(f"{lbl}: ΔG = {delta:.6f} nS")
    min_delta = min(deltas)
    max_delta = max(deltas)
    spread = max_delta - min_delta
    arr = np.array(deltas)
    iqr = np.subtract(*np.percentile(arr, [75, 25]))
    if len(arr) > 1:
        bin_width = 2.0 * iqr / (len(arr) ** (1.0 / 3.0)) if iqr > 0 else 0
        if not np.isfinite(bin_width) or bin_width <= 0:
            std = arr.std(ddof=1)
            bin_width = 3.5 * std / (len(arr) ** (1.0 / 3.0)) if std > 0 else 0
    else:
        bin_width = 0
    if bin_width <= 0:
        bin_width = spread / max(1, int(np.sqrt(len(arr)))) if spread > 0 else 0.01
    num_bins = max(1, int(np.ceil(spread / bin_width)))
    plt.figure(figsize=(8, 6))
    plt.hist(deltas, bins=num_bins, edgecolor='black')
    plt.xlabel(r"$\Delta G$ (nS)", fontsize=16, weight='bold')
    plt.ylabel('Count', fontsize=16, weight='bold')
    plt.title(title, fontsize=16, weight='bold')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
    print(f"Histogram saved as {output_path}")
} $labels_py $deltas_py $ref_current_val $voltage_val $output_path_py $plot_title_escaped]

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
        set stats_msg [format "Rotation scans: %d\nMean ΔG: %.6f nS\nStd. dev: %.6f nS\nMin ΔG: %.6f nS\nMax ΔG: %.6f nS\nSummary saved to: %s\nPlot saved to: %s" \
            $count $mean $stddev $min_delta $max_delta $csv_path $output_path]
        tk_messageBox -icon info -title "Rotation ΔG Summary" -message $stats_msg
    }

    catch {file delete -force $python_script_path}
}

}

if {![info exists ::PETK::gui::plotEnhancementsScheduled]} {
    set ::PETK::gui::plotEnhancementsScheduled 1
    after 100 ::PETK::gui::applyPlotEnhancements
}
