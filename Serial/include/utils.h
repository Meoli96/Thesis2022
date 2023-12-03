#include "pch.h"

void readCSV(std::string, double *, size_t, uint8_t);

size_t ceil2pow(uint16_t);

size_t ij2idx(int i, int j);

void loadRadiusGMAT(double *&, std::string, size_t& );

void printMatrix(double *mat, uint16_t degree, std::string file);