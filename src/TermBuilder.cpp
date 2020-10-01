//
// Created by caduser on 09.09.19.
//

#include "TermBuilder.h"
#include "InputData.h"
#include "model.h"

#include <iostream>
#include <iomanip>

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

TermBuilder::TermBuilder(std::string _name, int _edgesize, InputData *_inputData) {
    hyperedges.resize(_edgesize);
    chyperedges.resize(_edgesize);
    name = _name;
    inputData = _inputData;
}

void TermBuilder::checkpoint() {
    backup_hyperedges = hyperedges;
}

void TermBuilder::restore() {
    hyperedges = backup_hyperedges;
    copy();
}

void TermBuilder::dumpMetatdata(std::ostream &ostream) {
    for (int i = 0; i < labels.size(); i++) {
        ostream << labels[i].type << ",";
    }
    ostream << "\n";
}

std::string TermBuilder::getName() { return name; }

void TermBuilder::dumpNodes(std::ostream &ostream) {
    assert(hyperedges.size() != 0);
    for (int j = 0; j < hyperedges[0].size(); j++) {
        ostream << name << std::setfill('0') <<
                "_" << std::setw(8) << j << ";" <<
                name << std::setfill('0') <<
                "_" << std::setw(8) << j << "\n";
        for (int i = 0; i < hyperedges.size(); i++) {
            ostream << labels[i].type
                    << std::setfill('0') << "_" << std::setw(8) << hyperedges[i][j]
                    << ";" << labels[i].type
                    << std::setfill('0') << "_" << std::setw(8) << hyperedges[i][j] << "\n";

        }
    }
}

void TermBuilder::dumpEdges(std::ostream &ostream) {
    assert(hyperedges.size() != 0);
    for (int j = 0; j < hyperedges[0].size(); j++) {
        for (int i = 0; i < hyperedges.size(); i++) {
            ostream << labels[i].type
                    << std::setfill('0') << "_" << std::setw(8) << hyperedges[i][j]
                    << ";" << name << std::setfill('0') <<
                    "_" << std::setw(8) << j <<
                    ";\"" << name << "\"\n";

        }
    }
}


void TermBuilder::add(const std::vector <BarrierIndex> edge) {
    assert(hyperedges.size() == edge.size());
    for (int i = 0; i < hyperedges.size(); i++) {
        if (inputData != NULL) {
            bool verified = inputData->verifyIndex(edge[i]);
            assert(verified);
        }
        hyperedges[i].push_back(edge[i]());
    }
    if (labels.size() == 0) {
        for (int i = 0; i < hyperedges.size(); i++) {
            labels.push_back(edge[i]);
        }
    }
}

void TermBuilder::set(int edge, int offset, BarrierIndex value) {
    assert(value.sharesType(labels[edge]));
    hyperedges[edge][offset] = value();
}

int TermBuilder::get(int edge, int offset) {
    return hyperedges[edge][offset];
}

BarrierIndex TermBuilder::getTyped(int edge, int offset) {
    return BarrierIndex(labels[edge].getType(),hyperedges[edge][offset]);
}


void TermBuilder::allocate() {

    for (int i = 0; i < chyperedges.size(); i++) {
        gpuErrchk(cudaMalloc(&chyperedges[i], hyperedges[i].size() * sizeof(int)));
    }
}

void TermBuilder::release() {
    for (int i = 0; i < chyperedges.size(); i++) {
        cudaFree(chyperedges[i]);
    }
}

void TermBuilder::copy() {
    for (int i = 0; i < chyperedges.size(); i++) {
        cudaMemcpy(chyperedges[i], hyperedges[i].data(), hyperedges[i].size() * sizeof(int), cudaMemcpyHostToDevice);
    }
}

int TermBuilder::assign(void **result) {
    size = count();
    result[0] = &size;
    for (int i = 0; i < chyperedges.size(); i++) {
        result[i + 1] = chyperedges[i];
    }
    return chyperedges.size();
}

int TermBuilder::count() {
    return hyperedges[0].size();
}

int TermBuilder::countIdxArrays() {
    return hyperedges.size();
}

#undef gpuErrchk
