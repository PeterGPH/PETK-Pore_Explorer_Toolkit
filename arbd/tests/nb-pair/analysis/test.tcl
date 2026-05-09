
set ID [mol new ../BrownDyn.out.0.dcd waitfor all]
set sel [atomselect $ID all]
set out rho.dat

set rmax 10
set dr 0.1

array set rs ""
array set count ""


## center on first atom
set size [expr 57*0.7017544]
set sel [atomselect $ID "index 1"]
set all [atomselect $ID "all"]

for {set f 0} {$f < [molinfo $ID get numframes]} {incr f} {
    animate goto $f
    foreach dir {a b c} {
	molinfo $ID set $dir $size
    }
    foreach dir {alpha beta gamma} {
	molinfo $ID set $dir 90
    }
    $sel frame $f
    $all frame $f
    $all moveby [vecinvert [join [$sel get {x y z}]]]
}

package require pbctools
pbc wrap -sel all -center origin -all

for {set i 0} {$i*$dr < $rmax} {incr i} {
    set rs($i) 0
}
set tot 0

foreach rl [measure bond {0 1} frame all] {
    set i [expr {int(double($rl)/$dr)}]
    incr rs($i)
    incr tot
}

