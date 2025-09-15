#-*-makefile-*-
bld/plasimmod_tl.o: bld/resmod_tl.o
bld/oceanmod_tl.o: bld/resmod_tl.o
bld/miscmod_tl.o: bld/plasimmod_tl.o
bld/landmod_tl.o: bld/plasimmod_tl.o
bld/simba_tl.o: bld/landmod_tl.o
bld/surfmod_tl.o: bld/plasimmod_tl.o bld/simba_tl.o bld/landmod_tl.o bld/icemod_tl.o
bld/radmod_tl.o: bld/plasimmod_tl.o
bld/tafmod_tl.o: bld/plasimmod_tl.o bld/oceanmod_tl.o bld/radmod_tl.o bld/miscmod_tl.o bld/landmod_tl.o bld/surfmod_tl.o
bld/icemod_tl.o: bld/resmod_tl.o
bld/seamod_tl.o: bld/plasimmod_tl.o
bld/rainmod_tl.o: bld/plasimmod_tl.o
bld/fluxmod_tl.o: bld/plasimmod_tl.o
bld/wrapper_tl.o: bld/plasimmod_tl.o bld/tafmod_tl.o bld/icemod_tl.o bld/landmod_tl.o bld/oceanmod_tl.o bld/radmod_tl.o bld/rainmod_tl.o bld/seamod_tl.o bld/simba_tl.o
bld/mpimod_stub_tl.o: bld/plasimmod_tl.o
bld/fftmod_flow_directives_tl.o: bld/plasimmod_tl.o
bld/fftmod_tl.o: bld/plasimmod_tl.o
bld/cpl_stub_tl.o: bld/plasimmod_tl.o
bld/guimod_stub_tl.o: bld/plasimmod_tl.o
bld/p_earth_tl.o: bld/radmod_tl.o
bld/plasim_tl.o: bld/plasimmod_tl.o bld/seamod_tl.o bld/radmod_tl.o bld/oceanmod_tl.o bld/icemod_tl.o bld/landmod_tl.o bld/rainmod_tl.o bld/tafmod_tl.o
bld/legmod_tl.o: bld/plasimmod_tl.o
bld/outmod_tl.o: bld/plasimmod_tl.o
