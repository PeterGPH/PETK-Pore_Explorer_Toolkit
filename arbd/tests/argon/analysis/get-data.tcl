
set ID [mol new ../BrownDyn.out.0.dcd first 50 last 100 waitfor all]
set sel [atomselect $ID all]
set out rho.dat

set rmax 10
set dr 0.1

## center on first atom
set size [expr 57*0.7017544]
set sel [atomselect $ID "index 1"]
set all [atomselect $ID "all"]

set dim [expr 1.5117018*57]
set dim "$dim $dim $dim"


set last [expr [molinfo $ID get numframes]-1]

for {set f 0} {$f <= $last} {incr f} {
    animate goto $f
    molinfo $ID set {a b c} $dim
    molinfo $ID set {alpha beta gamma} {90 90 90}
}

# set result [measure gofr $all $all delta 0.2 usepbc 1 first 0 last $last]
set result [measure gofr $all $all delta 0.1 usepbc 1 first 0 last $last step 5]
lassign $result rs gofrs int hist frames

set ch [open $out w]
foreach r $rs g $gofrs {
    puts $ch "$r $g"
}
