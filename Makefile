
DEBUG=0
PROFILE=0
OPTIMIZE=1

CXX=g++

ifeq ($(DEBUG),1)
DEBUGFLAGS=-g -fno-inline
OPTIMIZE=0
else
ifeq ($(DEBUG),2)
DEBUGFLAGS=-g -fno-inline -DDEBUG
OPTIMIZE=0
else
DEBUGFLAGS=
endif
endif

ifeq ($(PROFILE),1)
PROFILEFLAGS=-pg
OPTIMIZE=0
else
PROFILEFLAGS=
endif

ifeq ($(OPTIMIZE),1)
OPTIMIZEFLAGS=-O2
else
OPTIMIZEFLAGS=
endif

CXXFLAGS=-Wall -Wno-deprecated $(DEBUGFLAGS) $(PROFILEFLAGS) $(OPTIMIZEFLAGS)
CFLAGS=-Wall $(DEBUGFLAGS) $(PROFILEFLAGS) $(OPTIMIZEFLAGS)

ALL=tc tt_BBEd_floats

all: $(ALL)

tc: tc.o chronometer.o Ed.o Samples.o LocalDistances.o debug.o Useful.o Tt.o Aesa.o Elem.o k.o Heap.o Bubu.o KBestVector.o
	$(CXX) $(CXXFLAGS) -o tc $^

tc.o: tc.cc
	$(CXX) $(CXXFLAGS) -c $<

tt_BBEd_floats: tt_BBEd_floats.o Ed.o Samples.o LocalDistances.o debug.o Useful.o
	$(CXX) $(CXXFLAGS) -o tt_BBEd_floats $^

tt_BBEd_floats.o: tt_BBEd_floats.cc
	$(CXX) $(CXXFLAGS) -c $<

debug.o: debug.cc debug.hh
	$(CXX) $(CXXFLAGS) -c $<

Ed.o: Ed.cc Ed.hh
	$(CXX) $(CXXFLAGS) -c $< 

KBestVector.o: KBestVector.cc KBestVector.hh
	$(CXX) $(CXXFLAGS) -c $< 

Bubu.o: Bubu.cc Bubu.hh
	$(CXX) $(CXXFLAGS) -c $< 

Samples.o: Samples.cc Samples.hh
	$(CXX) $(CXXFLAGS) -c $<

LocalDistances.o: LocalDistances.cc LocalDistances.hh
	$(CXX) $(CXXFLAGS) -c $<

Useful.o: Useful.cc Useful.hh
	$(CXX) $(CXXFLAGS) -c $<

Aesa.o: Aesa.cc Aesa.hh
	$(CXX) $(CXXFLAGS) -c $<

chronometer.o: chronometer.c chronometer.h
	gcc $(CFLAGS) -c $<

Elem.o: Elem.cc Elem.hh
	$(CXX) $(CXXFLAGS) -c $<

k.o: k.cc k.hh
	$(CXX) $(CXXFLAGS) -c $<

Heap.o: Heap.cc
	$(CXX) $(CXXFLAGS) -c $<

clean: 
	rm -rf *.o $(ALL)
