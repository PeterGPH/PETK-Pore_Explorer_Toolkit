# Tcl package index for the PETK GUI. Lets VMD load it as a package:
#   lappend auto_path /path/to/PETK-Pore_Explorer_Toolkit/petk
#   package require petk_gui
#   ::PETK::gui::petk_gui
package ifneeded petk_gui 1.0 [list source [file join $dir petk_gui.tcl]]
