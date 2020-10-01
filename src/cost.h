#pragma once
#include "types.h"
#include "InputData.h"
#include "TermBuilder.h"
#include "model.h"

class Cost
{
private:
	InputData* inputData;
	std::vector<TermBuilder>* terms;
	Scalar w_fitSqrt;
	Scalar w_surface;

	Scalar catmulSimple(const std::vector<BarrierIndex>& args);
public:
	Cost(InputData* _inputData, std::vector<TermBuilder>* _terms);
	Scalar compute(Scalar _w_fitSqrt, Scalar _w_surface);
};
