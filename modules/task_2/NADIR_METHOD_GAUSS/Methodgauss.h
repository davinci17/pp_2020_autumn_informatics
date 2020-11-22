#ifndef MODULES_TASK_2_NADIR_M_HORIZONTAL_METHOD_GAUSS_HORIZONTAL_METHOD_H_
#define MODULES_TASK_2_NADIR_M_HORIZONTAL_METHOD_GAUSS_HORIZONTAL_METHOD_H_


#include<vector>


void displayVector(std::vector<std::vector<double>> vec);
std::vector<std::vector<double>> GenerateRandom();
std::vector<std::vector<double>> methodGauss(std::vector<std::vector<double>> vec);
std::vector<std::vector<double>> methodGaussParallel(std::vector<std::vector<double>> vec);

#endif // MODULES_TASK_2_NADIR_M_HORIZONTAL_METHOD_GAUSS_HORIZONTAL_METHOD_H_

