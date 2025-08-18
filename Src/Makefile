CXX = c++
PATHBASE = /opt/homebrew
CFLAGS = -std=c++17 -I$(PATHBASE)/include -I/usr/local/include -I/opt/X11/include
LINK = -v -L$(PATHBASE)/lib -L/usr/local -L. -L/opt/X11/lib -lX11 -lgmpxx -lgmp -lboost_filesystem -lm
ARCHIVE = /usr/bin/ar cr
OBJECTS = Posit32.o PositBase.o bfloat16.o gmp_helpers.o

# all: Test2 matrix

all: Benchmarks NewBenchmarks

libPosit.a: PositBase.o Posit32.o Quire.o bfloat16.o FFT.o gmp_helpers.o helpers.o mmio.o
	$(ARCHIVE) libPosit.a PositBase.o Posit32.o helpers.o gmp_helpers.o FFT.o bfloat16.o Quire.o mmio.o

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

Posit32.o: Posit32.cc Posit32.hh PositBase.hh
	$(CXX) $(CFLAGS) -c Posit32.cc

# note that Posit32nc is still present, but not needed
Posit32nc.o: Posit32nc.cc Posit32nc.hh PositBase.hh
	$(CXX) $(CFLAGS) -c Posit32nc.cc
	
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
	$(CXX) $(CFLAGS) -o Benchmarks Benchmarks.o -lPosit $(LINK)

Benchmarks.o: Benchmarks.cc Matrix.hh mmio.hh helpers.hh FFT.hh Quire.hh
	$(CXX) $(CFLAGS) -c Benchmarks.cc

NewBenchmarks: NewBenchmarks.o libPosit.a
	$(CXX) $(CFLAGS) -o NewBenchmarks NewBenchmarks.o -lPosit $(LINK)

NewBenchmarks.o: NewBenchmarks.cc Matrix.hh mmio.hh helpers.hh FFT.hh Quire.hh
	$(CXX) $(CFLAGS) -c NewBenchmarks.cc

gmp_helpers.o: gmp_helpers.cc gmp_helpers.hh
	$(CXX) $(CFLAGS) -c gmp_helpers.cc

Quire.o: Quire.cc Quire.hh bfloat16.hh helpers.hh Matrix.hh Posit32.hh
	$(CXX) $(CFLAGS) -c Quire.cc

clean:
	rm -f *.o

distclean: clean
	rm -f Test Test2 HuffmanCoding Benchmarks libPosit.a

