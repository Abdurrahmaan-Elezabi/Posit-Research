CXX = g++
PATHBASE = /opt/homebrew
CFLAGS = -I$(PATHBASE)/include
LINK = -v -L$(PATHBASE)/lib -L. -lgmpxx -lgmp -lm
ARCHIVE = /usr/bin/ar cr
OBJECTS = Posit32.o PositBase.o bfloat16.o gmp_helpers.o

# all: Test2 matrix

all: libPosit.a

libPosit.a: PositBase.o Posit32.o
	$(ARCHIVE) libPosit.a PositBase.o Posit32.o

Test.o: Test.cc PositBase.hh Posit32.hh gmp_helpers.hh
	$(CXX) $(CFLAGS) -c Test.cc

Test: Test.o libPosit.a
	$(CXX) $(CFLAGS) -o Test Test.o $(LINK) -lPosit

Test2: Test2.o SoftPosit $(OBJECTS)
	$(CXX) $(CFLAGS) -o Test2 Test2.o $(OBJECTS) $(LINK)
	
Test2.o: Test2.cc
	$(CXX) $(CFLAGS) -c Test2.cc

PositBase.o: PositBase.cc PositBase.hh
	$(CXX) $(CFLAGS) -c PositBase.cc

Posit32.o: Posit32.cc Posit32.hh
	$(CXX) $(CFLAGS) -c Posit32.cc

bfloat16.o: bfloat16.cc bfloat16.hh
	$(CXX) $(CFLAGS) -c bfloat16.cc

gmp_helpers.o: gmp_helpers.cc gmp_helpers.hh
	$(CXX) $(CFLAGS) -c gmp_helpers.cc
	
clean:
	rm -f *.o
	cd MatrixAlgebra && rm -f *.o

distclean: clean
	rm Test Test2 libPosit.a

