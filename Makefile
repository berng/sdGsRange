CC = gcc
CFLAGS = -Iinclude/superdarn -Iinclude/base -Iinclude/general -I/usr/include -Iusr/include/superdarn -I./src -D_SIGMA_MAX 
# -DDEBUG
LDFLAGS = -lm -Llib -lz -lfit.1 -ldmap.1 -laacgm.1 -lradar.1 -lrcnv.1 -lcfit.1 -lraw.1 -lrscan.1 -lrtime.1 -lgfortran
FFLAGS = -freal-4-real-8
#LDFLAGS = -lm -Llib -lz
all:	process_noise_and_range
process_noise_and_range: process_noise_and_range.c src/gs_lib.a   iri-src/iri2007.a
# invmag.o calc_Distance.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o process_noise_and_range process_noise_and_range.c  src/gs_lib.a  	iri-src/iri2007.a 
	rm -f *.o

iri-src/iri2007.a:
	cd ./iri-src && make
	rm -f d
	ln -s ./iri-src/d d
src/gs_lib.a:
	cd ./src && make all distclean

clean:
	rm -f process_noise_and_range
	rm -f *.o
	rm -f d
	cd ./iri-src && make clean
	cd ./src && make clean
	
