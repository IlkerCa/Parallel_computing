#pragma once
#include "utilities/OBJLoader.hpp"
#include <vector>

union sse_float4 {
    float __attribute__ ((vector_size (16))) vector;
    float elements[4];
	float4 help;
};

void sse_test(Mesh &);