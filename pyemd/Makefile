OCTAVEDIR = D:/Octave/Octave-8.2.0/mingw64
CXX = gcc
CFLAGS = -m64 -Wall -O3 -fPIC -I$(OCTAVEDIR)/include/octave-8.2.0/octave
COMMON = emd.c
HEADER = emd.h
OBJECT = emd.o
MAIN = c_emd.c

# MEX = $(OCTAVEDIR)/bin/
# MEX_OPTION = CFLAGS="$(CFLAGS)" CC=$(CXX)
# # comment the following for 32-bit systems
# MEX_OPTION += -largeArrayDims

all: octave

octave: c_emd.mex

c_emd.mex: libEMD.a $(MAIN)
	mkoctfile --mex libEMD.a $(MAIN) -o c_emd.mex 

libEMD.a: emd.o
	gcc-ar ruv libEMD.a emd.o

# c_emd.o: $(MAIN)
# 	$(CXX) $(CFLAGS) \
# 	-I"$(OCTAVEDIR)/include/octave-8.2.0/octave/" \
# 	-L"$(OCTAVEDIR)/lib/" \
# 	-c $(MAIN) -DMX_COMPAT_32 \
# 		-fno-common -fexceptions \
# 		-O2 -fwrapv -flto

emd.o: $(COMMON) $(HEADER)
	$(CXX) $(CFLAGS) -c emd.c -DMX_COMPAT_32 \
		-fno-common -fexceptions \
		-O2 -fwrapv -flto

clean:
	del -f *~ *.a *.o *.oct *.obj
