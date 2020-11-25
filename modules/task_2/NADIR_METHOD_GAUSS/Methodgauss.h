// Copyright 2020 Nadir Mohammed
#ifndef MODULES_TASK_2_NADIR_METHOD_GAUSS_METHODGAUSS_H_
#define MODULES_TASK_2_NADIR_METHOD_GAUSS_METHODGAUSS_H_

std::vector<std::vector<double>> methodGauss(std::vector<std::vector<double>>array, int row);
double* convertVector(std::vector<std::vector<double>> vec, double* solution2);
double* methodGaussParallel(const double* array,  int row,int col);
 
#endif // MODULES_TASK_2_NADIR_METHOD_GAUSS_METHODGAUSS_H_



