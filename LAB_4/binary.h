#pragma once
#include "functions.h"

#include <random>
#include <algorithm>

RunResult run_single_real(int bits_per_dim, std::mt19937 &rng, int func_id);
