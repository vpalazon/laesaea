/*===================================================================================

  File: Heap.hh

  ===================================================================================

    Copyright (C) 2016 Vicente Palazón-González

    This program is free software; you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this
    program (file LICENSE.txt); if not, write to the Free Software Foundation, Inc.,
    675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the author at:

    palazon@uji.es
    
  =================================================================================== */

#ifndef _HEAP_HH
#define _HEAP_HH

#include <vector>
#include <utility>
#include <map>

using namespace std;

class MinHeap {
public:
    vector<float> buffer;
    unsigned size;

    MinHeap(unsigned n, vector<float> &buffer);
    void heapify(unsigned i);
    float min();
    float extract_min();
    void insert(float item);
    void decrease(unsigned i, float value);
    unsigned is_empty();
    void printBuffer();
};

class IndexedMinHeap {
public:
    vector< pair<float, unsigned> > buffer;
    map<unsigned, unsigned> index;
    unsigned size;
    IndexedMinHeap() {};
    IndexedMinHeap(unsigned n, vector< pair<unsigned, float> > &buffer);
    void init(unsigned n);
    void heapify(unsigned i);
    void min(unsigned &key, float &score);
    void extract_min(unsigned &key, float &score);
    void insert(unsigned key, float score);
    void decrease(unsigned key, float score);
    unsigned contains(unsigned key);
    unsigned is_empty();
    float getitem(unsigned key);
    void printBuffer();
};


#endif
