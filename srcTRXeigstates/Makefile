CC=icpc
FC=ifort
LAPACK_PATH= /Users/jtcantin/addLibs/lapack-3.5.0/
LAPACKE_PATH= /Users/jtcantin/addLibs/lapack-3.5.0/lapacke/include/

SUP_PATH1= ./

INCLUDES= -I $(SUP_PATH1) -I${MKLROOT}/include -I/home/jtcantin/eigen 
#-I/home/tianruix/code/3dexc/eigen-eigen-6b38706d90a9
#-I/global/software/Eigen3
#
CFLAGS=-c -D LINUX_MKL $(INCLUDES)
FCFLAGS=-c -D LINUX_MKL $(INCLUDES) 

#LIB_SOURCES= $(SUP_PATH1)quantumPercolation.cpp $(SUP_PATH1)lapackCustInterface.cpp
LIB_SOURCES= $(SUP_PATH1)main.cpp $(SUP_PATH1)lapackCustInterface.cpp $(SUP_PATH1)matrix_set.cpp

LDFLAGS= -mkl


SOURCES= $(LIB_SOURCES)
SOURCES2 = $(SOURCES:.cpp=.o)
OBJECTS= $(SOURCES2:.f=.o)
EXECUTABLE=QmPercEigstates_LINUX_TRX

#print_vars:
#	echo $(SOURCES)
#	echo $(SOURCES2)
#	echo $(OBJECTS)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.f.o:
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
