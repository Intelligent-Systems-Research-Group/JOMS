#pragma once

#include <vector>
#include "types.h"
#include "InputData.h"
#include "model.h"

class TermBuilder {
    std::vector<std::vector<int>> hyperedges;
    std::vector<BarrierIndex> labels;
    std::vector<int*> chyperedges;

    std::vector<std::vector<int>> backup_hyperedges;
    int size;
    std::string name;
    InputData* inputData;

public:
    TermBuilder(std::string _name, int _edgesize, InputData* _inputData = NULL);
    void checkpoint();

    void restore();

    void dumpMetatdata(std::ostream& ostream);

    std::string getName();

    void dumpNodes(std::ostream& ostream);

    void dumpEdges(std::ostream& ostream);


    void add(const std::vector<BarrierIndex> edge);

    void set(int edge, int offset, BarrierIndex value);
    int get(int edge, int offset);
	BarrierIndex getTyped(int edge, int offset);

    void allocate();
    void release();

    void copy();

    int assign(void** result);

    int count();

    int countIdxArrays();
};

