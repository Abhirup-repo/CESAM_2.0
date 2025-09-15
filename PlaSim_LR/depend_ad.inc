#-*-makefile-*-
bld/plasimmod_ad.o: bld/resmod_ad.o
bld/oceanmod_ad.o: bld/resmod_ad.o
bld/miscmod_ad.o: bld/plasimmod_ad.o
bld/landmod_ad.o: bld/plasimmod_ad.o
bld/simba_ad.o: bld/landmod_ad.o bld/taf_ad.o
bld/surfmod_ad.o: bld/plasimmod_ad.o
bld/radmod_ad.o: bld/plasimmod_ad.o bld/taf_ad.o
bld/srvio_ad.o: bld/plasimmod_ad.o 
bld/observation_ad.o: bld/plasimmod_ad.o bld/srvio_ad.o
bld/tafmod_ad.o: bld/plasimmod_ad.o bld/oceanmod_ad.o bld/radmod_ad.o bld/miscmod_ad.o bld/landmod_ad.o bld/surfmod_ad.o
bld/icemod_ad.o: bld/resmod_ad.o
bld/seamod_ad.o: bld/plasimmod_ad.o bld/icemod_ad.o bld/oceanmod_ad.o
bld/rainmod_ad.o: bld/plasimmod_ad.o bld/taf_ad.o
bld/fluxmod_ad.o: bld/plasimmod_ad.o
bld/wrapper_ad.o: bld/plasimmod_ad.o bld/tafmod_ad.o bld/taf_ad.o bld/plasim_ad.o bld/simba_ad.o bld/radmod_ad.o bld/rainmod_ad.o bld/icemod_ad.o bld/seamod_ad.o bld/landmod_ad.o bld/oceanmod_ad.o bld/surfmod_ad.o bld/fluxmod_ad.o bld/miscmod_ad.o
bld/mpimod_stub_ad.o: bld/plasimmod_ad.o
bld/fftmod_flow_directives_ad.o: bld/plasimmod_ad.o
bld/fftmod_ad.o: bld/plasimmod_ad.o
bld/cpl_stub_ad.o: bld/plasimmod_ad.o
bld/guimod_stub_ad.o: bld/plasimmod_ad.o
bld/p_earth_ad.o: bld/radmod_ad.o
bld/plasim_ad.o: bld/plasimmod_ad.o bld/seamod_ad.o bld/radmod_ad.o bld/taf_ad.o bld/oceanmod_ad.o bld/icemod_ad.o bld/landmod_ad.o bld/rainmod_ad.o bld/tafmod_ad.o
bld/legmod_ad.o: bld/plasimmod_ad.o
bld/outmod_ad.o: bld/plasimmod_ad.o
