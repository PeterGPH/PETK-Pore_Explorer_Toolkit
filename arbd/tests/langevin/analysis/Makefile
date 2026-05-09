all:
	$(MAKE) rho.dat
	$(MAKE) $(patsubst %.plt.py,%.pdf, $(wildcard *.plt.py))

plot.pdf: plot.plt.py rho.dat plotTools.py
	./$< $@

rho.dat: get-data.tcl ../BrownDyn.out.0.dcd
	vmd -dispdev text < $< 

.SECONDARY:
