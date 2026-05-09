all:
	$(MAKE) $(patsubst %.plt.py,%.pdf, $(wildcard *.plt.py))

boltzmann-check.pdf: boltzmann-check.plt.py rho.dat plotTools.py
	./$< $@

rho.dat:  ../out.0.dcd get-data.tcl
	vmd -dispdev text -args $< $@ < get-data.tcl

.SECONDARY:
