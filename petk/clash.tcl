proc petk_radius_clash {molid_ana molid_pore buffer {frame 0} {cutoff ""}} {
    set sel_ana  [atomselect $molid_ana  "all" frame $frame]
    set sel_pore [atomselect $molid_pore "all" frame $frame]

    set ana_idx  [$sel_ana get index]
    set pore_idx [$sel_pore get index]
    set ana_pos  [$sel_ana get {x y z}]
    set pore_pos [$sel_pore get {x y z}]
    set ana_rad  [$sel_ana get radius]
    set pore_rad [$sel_pore get radius]

    if {$cutoff eq ""} {
        set max_a 0.0
        foreach r $ana_rad  { if {$r > $max_a} { set max_a $r } }
        set max_p 0.0
        foreach r $pore_rad { if {$r > $max_p} { set max_p $r } }
        set cutoff [expr {$max_a + $max_p + $buffer}]
    }

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

        set ra [lindex $ana_rad  $ai]
        set rp [lindex $pore_rad $pi]

        if {$dist <= ($ra + $rp + $buffer)} {
            lappend clash_ana  [lindex $ana_idx  $ai]
            lappend clash_pore [lindex $pore_idx $pi]
        }
    }

    set clash_ana  [lsort -unique $clash_ana]
    set clash_pore [lsort -unique $clash_pore]

    if {[llength $clash_ana] == 0} {
        puts "No radius-sum clashes found."
        return
    }

    # highlight analyte clash atoms
    set rep_ana [molinfo $molid_ana get numreps]
    mol addrep $molid_ana
    mol modselect $rep_ana $molid_ana "index $clash_ana"
    mol modstyle  $rep_ana $molid_ana VDW
    mol modcolor  $rep_ana $molid_ana ColorID 1

    # highlight pore clash atoms
    set rep_pore [molinfo $molid_pore get numreps]
    mol addrep $molid_pore
    mol modselect $rep_pore $molid_pore "index $clash_pore"
    mol modstyle  $rep_pore $molid_pore VDW
    mol modcolor  $rep_pore $molid_pore ColorID 3

    puts "Clashes: analyte [llength $clash_ana], pore [llength $clash_pore]"
}
