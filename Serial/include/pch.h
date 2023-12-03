// Precompiled header

// std includes
#include <cmath>
#include <string>
#include <chrono>

#include <iostream>
#include <fstream>
#include <exception>

// Load CUDA library
#ifdef CUDA

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#endif
