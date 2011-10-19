
########## Setup File Lists ##########

CFILES = Timer_test1.cc ../Timer.cc
 
OFILES = Timer_test1.o ../Timer.o 

.SUFFIXES:	.cc .o

# TIMER -- turns Timers on/off
# PAPI  -- turns PAPI counters on/off (won't work unless u have PAPI)
# COUNTER - specifies what counters to use with PAPI

C++	= g++  -O2 -DTIMER #-DPAPI -DTIMER_COUNTER=1#-DNDEBUG 
C++FLAGS = # -DUSE_MPI -I/usr/include

INCLUDES = -I/usr/local/papi/include -I/usr/local/include
#-I/usr/local/apps/mpi/include

LIBS = -lm -L.. -lUtilities #/usr/local/lib/libpapi.a #-lmpi -lm
#LIBS = -lm -L.. -lUtilities /usr/local/papi/lib/libpapi.a #-lmpi -lm
# LIBS =  -lm -L/usr/local/apps/mpi/lib -lmpich -L/usr/local/lib -lgs
# LIBS = -L/usr/lpp/ppe.poe/lib -lmpi -lm
# LIBS = -L/home/noel/src/mpich/build/LINUX/ch_shmem/lib -lmpich
#-lrpc -L/usr/local/lib -lgs

##########  creating executable  ##########

t: $(OFILES)
	$(C++) -o $@ $(OFILES) $(LIBS)

.cc.o: 
	${C++} ${C++FLAGS} ${INCLUDES} -o $*.o -c $*.cc

depend:
	makedepend ${CFILES}

clean:
	rm -f *~ core* *.o t

clobber:
	rm -f *~ core* *.o *.ii ii_files/*.ii t


# DO NOT DELETE THIS LINE -- make depend depends on it.

