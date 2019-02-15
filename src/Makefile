CC = gcc
CFLAGS = -I../include/superdarn -I../include/base -I../include/general -I/usr/include -I/usr/include/superdarn -D_SIGMA_MAX
LDFLAGS = -lm -L../lib -lz -lfit -ldmap -laacgm -lradar -lrcnv -lcfit -lraw -lrscan -lrtime -lgfortran
FFLAGS = -freal-4-real-8
AR = ar
ARFLAGS = cr
#LDFLAGS = -lm -Llib -lz
all:	gs_lib.a
gs_lib.a: get_fof2_iri.o calc_Distance.o calc_sigma.o zenith_angle.o Solution.o genetic_search.o claster_analysis.o spherical.o FixParams.o
	$(AR) $(ARFLAGS) gs_lib.a get_fof2_iri.o calc_Distance.o calc_sigma.o zenith_angle.o Solution.o genetic_search.o claster_analysis.o spherical.o FixParams.o

get_fof2_iri.o: get_fof2_iri.c get_fof2_iri.h
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o get_fof2_iri.o get_fof2_iri.c 

genetic_search.o: genetic_search.c genetic_search.h
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o genetic_search.o genetic_search.c 
claster_analysis.o: claster_analysis.c claster_analysis.h
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o claster_analysis.o claster_analysis.c 
# iri-src/iri2007.a:
# 	cd ./iri-src && make
# 	rm -f d
# 	ln -s ./iri-src/d d

calc_sigma.o: 	calc_sigma.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o calc_sigma.o calc_sigma.c 

zenith_angle.o: zenith_angle.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o zenith_angle.o zenith_angle.c 

invmag.o: 	invmag.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o invmag.o invmag.c 
Solution.o: 	Solution.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o Solution.o Solution.c 
calc_Distance.o: 	calc_Distance.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o calc_Distance.o calc_Distance.c 
spherical.o: 	spherical.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o spherical.o spherical.c 
FixParams.o: 	FixParams.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c -o FixParams.o FixParams.c 

distclean:
	rm -f *.o

clean:
	rm -f process_noise_and_range
	rm -f *.o *.a
	rm -f d
#	cd ./iri-src && make clean
	
