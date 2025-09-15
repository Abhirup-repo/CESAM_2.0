#-*-makefile-*-
bld/plasimmod.o: bld/resmod.o
bld/oceanmod.o: bld/resmod.o
bld/miscmod.o: bld/plasimmod.o
bld/landmod.o: bld/plasimmod.o
bld/simba.o: bld/landmod.o
bld/surfmod.o: bld/plasimmod.o
bld/radmod.o: bld/plasimmod.o
bld/tafmod.o: bld/plasimmod.o bld/oceanmod.o bld/radmod.o bld/miscmod.o bld/landmod.o bld/surfmod.o
bld/icemod.o: bld/resmod.o
bld/seamod.o: bld/plasimmod.o
bld/rainmod.o: bld/plasimmod.o
bld/fluxmod.o: bld/plasimmod.o
bld/wrapper.o: bld/plasimmod.o bld/tafmod.o
bld/mpimod_stub.o: bld/plasimmod.o
bld/fftmod_flow_directives.o: bld/plasimmod.o
bld/fftmod.o: bld/plasimmod.o
bld/cpl_stub.o: bld/plasimmod.o
bld/guimod_stub.o: bld/plasimmod.o
bld/p_earth.o: bld/radmod.o
bld/plasim.o: bld/plasimmod.o bld/oceanmod.o bld/seamod.o bld/icemod.o bld/landmod.o bld/radmod.o bld/rainmod.o bld/tafmod.o
bld/legmod.o: bld/plasimmod.o
bld/outmod.o: bld/plasimmod.o
