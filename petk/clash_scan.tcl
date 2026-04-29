proc petk_scan_clashes {molid_ana molid_pore buffer {frame_all 1} {overlap_distance_threshold ""}} {
    global petk_clash_ana petk_clash_pore
    array unset petk_clash_ana
    array unset petk_clash_pore

    set sel_ana  [atomselect $molid_ana  "all"]
    set sel_pore [atomselect $molid_pore "all"]

    set ana_idx [$sel_ana get index]
    set pore_idx [$sel_pore get index]
    set ana_rad [$sel_ana get radius]
    set pore_rad [$sel_pore get radius]

    set max_a 0.0
    foreach r $ana_rad { if {$r > $max_a} { set max_a $r } }
    set max_p 0.0
    foreach r $pore_rad { if {$r > $max_p} { set max_p $r } }

    if {$overlap_distance_threshold ne ""} {
        set cutoff [expr {$overlap_distance_threshold + $buffer}]
    } else {
        set cutoff [expr {$max_a + $max_p + $buffer}]
    }

    set nframes [molinfo $molid_ana get numframes]
    set start 0
    set end [expr {$frame_all ? ($nframes - 1) : 0}]

    for {set frame $start} {$frame <= $end} {incr frame} {
        $sel_ana frame $frame
        $sel_ana update
        $sel_pore frame $frame
        $sel_pore update

        set ana_pos [$sel_ana get {x y z}]
        set pore_pos [$sel_pore get {x y z}]

        set res [measure contacts $cutoff $sel_ana $sel_pore]
        set ana_sel_idx  [lindex $res 0]
        set pore_sel_idx [lindex $res 1]

        set clash_ana {}
        set clash_pore {}

        set n [llength $ana_sel_idx]
        for {set k 0} {$k < $n} {incr k} {
            set ai [lindex $ana_sel_idx  $k]
            set pi [lindex $pore_sel_idx $k]

            set apos [lindex $ana_pos  $ai]
            set ppos [lindex $pore_pos $pi]

            set dx [expr {[lindex $apos 0] - [lindex $ppos 0]}]
            set dy [expr {[lindex $apos 1] - [lindex $ppos 1]}]
            set dz [expr {[lindex $apos 2] - [lindex $ppos 2]}]
            set dist [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]

            if {$overlap_distance_threshold ne ""} {
                if {$dist <= ($overlap_distance_threshold + $buffer)} {
                    lappend clash_ana  [lindex $ana_idx  $ai]
                    lappend clash_pore [lindex $pore_idx $pi]
                }
            } else {
                set ra [lindex $ana_rad  $ai]
                set rp [lindex $pore_rad $pi]
                if {$dist <= ($ra + $rp + $buffer)} {
                    lappend clash_ana  [lindex $ana_idx  $ai]
                    lappend clash_pore [lindex $pore_idx $pi]
                }
            }
        }

        set clash_ana  [lsort -unique $clash_ana]
        set clash_pore [lsort -unique $clash_pore]

        if {[llength $clash_ana] > 0} {
            set petk_clash_ana($frame)  $clash_ana
            set petk_clash_pore($frame) $clash_pore
            puts "frame $frame: clashes ana [llength $clash_ana], pore [llength $clash_pore]"
        }
    }

    $sel_ana delete
    $sel_pore delete
    puts "Scan done."
}

proc petk_show_clashes {molid_ana molid_pore frame} {
    global petk_clash_ana petk_clash_pore petk_rep_ana petk_rep_pore

    if {![info exists petk_clash_ana($frame)]} {
        puts "No clashes at frame $frame"
        return
    }

    if {![info exists petk_rep_ana]} { set petk_rep_ana -1 }
    if {![info exists petk_rep_pore]} { set petk_rep_pore -1 }

    if {$petk_rep_ana < 0} {
        set petk_rep_ana [molinfo $molid_ana get numreps]
        mol addrep $molid_ana
        mol modstyle $petk_rep_ana $molid_ana VDW
        mol modcolor $petk_rep_ana $molid_ana ColorID 1
    }
    if {$petk_rep_pore < 0} {
        set petk_rep_pore [molinfo $molid_pore get numreps]
        mol addrep $molid_pore
        mol modstyle $petk_rep_pore $molid_pore VDW
        mol modcolor $petk_rep_pore $molid_pore ColorID 3
    }

    mol modselect $petk_rep_ana  $molid_ana  "index $petk_clash_ana($frame)"
    mol modselect $petk_rep_pore $molid_pore "index $petk_clash_pore($frame)"
}

proc petk_play_clashes {molid_ana molid_pore} {
    global petk_clash_ana
    foreach frame [lsort -integer [array names petk_clash_ana]] {
        animate goto $frame
        petk_show_clashes $molid_ana $molid_pore $frame
        display update
        after 150
    }
}

proc petk_print_radii_by_element {molid label} {
    set sel [atomselect $molid "all"]
    set elems [$sel get element]
    set rads  [$sel get radius]

    array set elem_rads {}
    set n [llength $elems]
    for {set i 0} {$i < $n} {incr i} {
        set e [string toupper [string trim [lindex $elems $i]]]
        if {$e eq ""} { set e "UNKNOWN" }
        set r [format "%.3f" [lindex $rads $i]]
        if {![info exists elem_rads($e)]} { set elem_rads($e) {} }
        if {[lsearch -exact $elem_rads($e) $r] < 0} {
            lappend elem_rads($e) $r
        }
    }

    puts "=== Radii by element: $label ==="
    foreach e [lsort [array names elem_rads]] {
        puts [format "%-6s : %s" $e [lsort -unique $elem_rads($e)]]
    }
    $sel delete
}