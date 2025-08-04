CXX = g++
PATHBASE = /opt/homebrew
CFLAGS = -I$(PATHBASE)/include
LINK = -v -L$(PATHBASE)/lib -L. -lgmpxx -lgmp -lm
ARCHIVE = /usr/bin/ar cr
OBJECTS = Posit32.o PositBase.o bfloat16.o gmp_helpers.o

# all: Test2 matrix

all: libPosit.a

libPosit.a: PositBase.o Posit32.o Quire.o bfloat16.o FFT.o gmp_helpers.o helpers.o
	$(ARCHIVE) libPosit.a PositBase.o Posit32.o helpers.o gmp_helpers.o FFT.o bfloat16.o Quire.o

Test.o: Test.cc PositBase.hh Posit32.hh gmp_helpers.hh
	$(CXX) $(CFLAGS) -c Test.cc

Test: Test.o libPosit.a
	$(CXX) $(CFLAGS) -o Test Test.o $(LINK) -lPosit

Test2: Test2.o libPosit.a
	$(CXX) $(CFLAGS) -o Test2 Test2.o -lPosit $(LINK)
	
Test2.o: Test2.cc
	$(CXX) $(CFLAGS) -c Test2.cc

PositBase.o: PositBase.cc PositBase.hh gmp_helpers.hh
	$(CXX) $(CFLAGS) -c PositBase.cc

# note that Posit32nc is still present, but not needed
Posit32.o: Posit32.cc Posit32.hh PositBase.hh
	$(CXX) $(CFLAGS) -c Posit32.cc
	
helpers.o: helpers.cc helpers.hh mmio.hh Posit32.hh half.hh bfloat16.hh
	$(CXX) $(CFLAGS) -c helpers.cc

HuffmanCoding: HuffmanCoding.cc
	$(CXX) $(CFLAGS) -o HuffmanCoding HuffmanCoding.cc $(LINK)

bfloat16.o: bfloat16.cc bfloat16.hh
	$(CXX) $(CFLAGS) -c bfloat16.cc

Matrix.o: Matrix.cc Matrix.hh mmio.hh helpers.hh
	$(CXX) $(CFLAGS) -c Matrix.cc

FFT.o: FFT.cc FFT.hh mmio.hh helpers.hh bfloat16.hh Posit32.hh Matrix.hh
	$(CXX) $(CFLAGS) -c FFT.cc

Benchmarks: Benchmarks.o libPosit.a
	$(CXX) $(CFLAGS) -c Benchmarks.cc -lPosit $(LINK)

Benchmarks.o: Benchmarks.cc Matrix.hh mmio.hh helpers.hh FFT.hh Quire.hh
	$(CXX) $(CFLAGS) -c Benchmarks.cc

gmp_helpers.o: gmp_helpers.cc gmp_helpers.hh
	$(CXX) $(CFLAGS) -c gmp_helpers.cc

Quire.o: Quire.cc Quire.hh bfloat16.hh helpers.hh Matrix.hh Posit32.hh
	$(CXX) $(CFLAGS) -c Quire.cc

clean:
	rm -f *.o

distclean: clean
	rm -f Test Test2 HuffmanCoding Benchmarks libPosit.a

