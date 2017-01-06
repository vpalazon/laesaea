/*===================================================================================

  File: Heap.cc

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

#include <vector>
#include <iostream>
#include <algorithm>
#include "Heap.hh"

using namespace std;

MinHeap::MinHeap(unsigned n, vector<float> &buffer) {
    this->buffer.push_back(NULL);

    for (unsigned i=0; i<buffer.size(); i++)
        this->buffer.push_back(buffer[i]);

    for (unsigned i=buffer.size(); i<n+1; i++)
        this->buffer.push_back(NULL);

    this->size = buffer.size();

    for (unsigned i=this->size/2; i>0; i--)
        this->heapify(i);
}

void MinHeap::heapify(unsigned i) {
    unsigned l, r, smallest;
    float aux;
    
    while (1) {
        l = 2*i;
        r = 2*i+1;
        if (l <= this->size and this->buffer[l] < this->buffer[i])
            smallest = l;
        else
            smallest = i;
        if (r <= this->size and this->buffer[r] < this->buffer[smallest])
            smallest = r;
        if (smallest == i) 
            break;
        aux = this->buffer[i];
        this->buffer[i] = this->buffer[smallest];
        this->buffer[smallest] = aux;
        i = smallest;
    }
}

float MinHeap::min() {
    if (this->size == 0) {
        cerr << "Error: The heap is empty." << endl;
        exit(1);
    }
    return this->buffer[1];
}

float MinHeap::extract_min() {
    if (this->size == 0) {
        cerr << "Error: The heap is empty." << endl;
        exit(1);
    }
    float m = this->buffer[1];
    this->buffer[1] = this->buffer[this->size];
    this->size -= 1;
    this->heapify(1);
    return m;
}

void MinHeap::insert(float item) {
    float aux;
    this->size += 1;
    unsigned i = this->size;
    this->buffer[i] = item;
    while (i > 1 and this->buffer[i/2] > this->buffer[i]) {
        aux = this->buffer[i];
        this->buffer[i] = this->buffer[i/2];
        this->buffer[i/2] = aux;
        i /= 2;
    }
}

void MinHeap::decrease(unsigned i, float value) {
    float aux;
    if (value < this->buffer[i]) {
        this->buffer[i] = value;
        while (i > 1 and this->buffer[i/2] > this->buffer[i]) {
            aux = this->buffer[i];
            this->buffer[i] = this->buffer[i/2];
            this->buffer[i/2] = aux;
            i /= 2;
        }
    }
}
    
unsigned MinHeap::is_empty() {
    return this->size == 0;
}

// Debugging
void MinHeap::printBuffer() {
    for (unsigned i=0; i<this->buffer.size(); i++)
        cout << this->buffer[i] << " ";
    cout << endl;
}

///////////////////////////////////////////////////////////////////////////

IndexedMinHeap::IndexedMinHeap(unsigned n, vector< pair<unsigned, float> > &buffer) {
    
    for (unsigned i=0; i<buffer.size(); i++)
        this->index[buffer[i].first] = i+1;

    this->buffer.push_back(pair<float, unsigned>(NULL, NULL));
    
    for (unsigned i=0; i<buffer.size(); i++)
        this->buffer.push_back(pair<float, unsigned>(buffer[i].second, 
                                                     buffer[i].first));
    for (unsigned i=buffer.size(); i<n+1; i++)
        this->buffer.push_back(pair<float, unsigned>(NULL, NULL));
    
    this->size = buffer.size();

    for (unsigned i=this->size/2; i>0; i--)
        this->heapify(i);
}

void IndexedMinHeap::init(unsigned n) {
    for (unsigned i=0; i<n+1; i++)
        this->buffer.push_back(pair<float, unsigned>(NULL, NULL));
    this->size = 0;
}

void IndexedMinHeap::heapify(unsigned i) {
    unsigned l, r, smallest;
    
    while (1) {
        l = 2*i;
        r = 2*i+1;
        if (l <= this->size and this->buffer[l].first < this->buffer[i].first)
            smallest = l;
        else
            smallest = i;
        if (r <= this->size and this->buffer[r].first < this->buffer[smallest].first)
            smallest = r;
        if (smallest == i)
            break;
        this->index[this->buffer[i].second] = smallest;
        this->index[this->buffer[smallest].second] = i;
        swap(this->buffer[i].first, this->buffer[smallest].first);
        swap(this->buffer[i].second, this->buffer[smallest].second);

        i = smallest;
    }
}
    
void IndexedMinHeap::min(unsigned &key, float &score) {
    if (this->size == 0) {
        cerr << "Error: The heap is empty." << endl;
        exit(1);
    }
    key = this->buffer[1].second;
    score = this->buffer[1].first;
}
 
void IndexedMinHeap::extract_min(unsigned &key, float &score) {
    if (this->size == 0) {
        cerr << "Error: The heap is empty." << endl;
        exit(1);
    }
    float m_first = this->buffer[1].first;
    float m_second = this->buffer[1].second;
    
    this->index.erase(m_second);
    
    this->buffer[1] = this->buffer[this->size];
    this->index[this->buffer[this->size].second] = 1;
    this->size -= 1;
    this->heapify(1);
    key = m_second;
    score = m_first;
}

void IndexedMinHeap::insert(unsigned key, float score) {
    this->size += 1;
    unsigned i = this->size;

    this->buffer[i].first = score;
    this->buffer[i].second = key;
    this->index[key] = i;

    while (i>1 and this->buffer[i/2].first > this->buffer[i].first) {
        this->index[this->buffer[i].second] = i/2;
        this->index[this->buffer[i/2].second] = i;
        swap(this->buffer[i].first, this->buffer[i/2].first);
        swap(this->buffer[i].second, this->buffer[i/2].second);
        i /= 2;
    }
}

void IndexedMinHeap::decrease(unsigned key, float score) {
    unsigned i = this->index[key];
    if (score >= this->buffer[i].first)
        return;
    this->buffer[i].first = score;
    this->buffer[i].second = key;
    while (i>1 and this->buffer[i/2].first > this->buffer[i].first) {
        this->index[this->buffer[i].second] = i/2;
        this->index[this->buffer[i/2].second] = i;
        swap(this->buffer[i].first, this->buffer[i/2].first);
        swap(this->buffer[i].second, this->buffer[i/2].second);
        i /= 2;
    }
}
 
unsigned IndexedMinHeap::contains(unsigned key) {
    return this->index.find(key) != this->index.end();
}

unsigned IndexedMinHeap::is_empty() {
    return this->size == 0;
}

float IndexedMinHeap::getitem(unsigned key) {
    return this->buffer[this->index[key]].first;
}

// Debugging
void IndexedMinHeap::printBuffer() {
    for (unsigned i=0; i<this->size+1; i++)
        cout << "[" << this->buffer[i].first << ", " 
             << this->buffer[i].second << "] ";
    cout << endl;
}
