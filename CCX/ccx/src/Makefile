#OPT = -O
# DO NOT SET AN OPT LEVEL ON CCX, ARPACK OR SPOOLES!!!
OPT = 
CFLAGS = -Wall $(OPT) -m64 -std=c99 -I ../../../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1
FFLAGS = -Wall $(OPT) -m64

CC=gcc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.8p2.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../../SPOOLES.2.2
LIBDIR=/mingw/x86_64-w64-mingw32/lib

LIBS = \
       $(DIR)/MT/src/spoolesMT.a \
       $(DIR)/spooles.a \
       ../../../ARPACK/libarpack_INTEL64.a \
       $(LIBDIR)/libpthread.a \
       $(LIBDIR)/libm.a

ccx_2.8p2_64MT.exe: $(OCCXMAIN) ccx_2.8p2_64MT.a  $(LIBS)
	./date.pl
	$(CC) $(CFLAGS) -c ccx_2.8p2.c
	$(FC) -static-libgcc -static-libgfortran -Wall $(OPT) -o $@ $(OCCXMAIN) ccx_2.8p2_64MT.a $(LIBS)

ccx_2.8p2_64MT.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
                                                                               
clean:
	rm *.o
