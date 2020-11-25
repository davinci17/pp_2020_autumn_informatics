// Copyright 2020 Nadir Mohammed
#include "../../../modules/task_2/NADIR_METHOD_GAUSS/Methodgauss.h"
#include<iostream>
#include<iomanip>
#include<mpi.h>
#include <stddef.h>
#include <math.h>
#include <ctime>





    
void methodGauss(const double* array, double* result, int size) {
    double* temp_array = new double[size * (size + 1)];
    for (int i = 0; i < size * (size + 1); i++) {
        temp_array[i] = array[i];
    }


    int p_Row;
    int* p_Pla = new int[size];
    int* p_Ind = new int[size];

    for (int i = 0; i < size; i++) {
        p_Ind[i] = -1;
        p_Pla[i] = 0;
    }

    int Ind;  

    for (Ind = 0; Ind < size; Ind++) {
       
        double n_max = 0;
        for (int i = 0; i < size; i++) {
            if ((p_Ind[i] == -1) && ((fabs(temp_array[i * (size + 1) + Ind]) > n_max))) {
                p_Row = i;
                n_max = temp_array[i * (size + 1) + Ind];
            }
        }

        p_Pla[Ind] = p_Row;
        p_Ind[p_Row] = Ind;

      
        double p_val, p_fra;
        p_val = temp_array[p_Row * (size + 1) + Ind];
        for (int i = 0; i < size; i++) {
            if (p_Ind[i] == -1) {
                p_fra = temp_array[i * (size + 1) + Ind] / p_val;
                for (int j = Ind; j < size; j++) {
                    temp_array[i * (size + 1) + j] -= p_fra * temp_array[p_Row * (size + 1) + j];
                }

                temp_array[i * (size + 1) + size] -= p_fra * temp_array[p_Row * (size + 1) + size];
            }
        }
    }

  
    int rowIdx, row;
    for (int i = (size - 1); i >= 0; i--) {
        rowIdx = p_Pla[i];
        result[i] = temp_array[rowIdx * (size + 1) + size] / temp_array[rowIdx * (size + 1) + i];
        temp_array[rowIdx * (size + 1) + i] = 1;

        for (int j = 0; j < i; j++) {
            row = p_Pla[j];
            temp_array[row * (size + 1) + size] -= temp_array[row * (size + 1) + i] * result[i];
            temp_array[row * (size + 1) + i] = 0;
        }
    }

    delete[] temp_array;
    delete[] p_Pla;
    delete[] p_Ind;
}

void methodGaussParallel(const double* array, double* solution, int col, int row) {
    double* temp_array = new double[row * col];
    for (int i = 0; i < row * col; i++) {
        temp_array[i] = array[i];
    }

    int p_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &p_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int numberRow = row / p_size;
  
    const int restRow = row % p_size;

    int stan_R;
    if (rank < restRow) {
        stan_R = numberRow + 1;
    }
    else {
        stan_R = numberRow;
    }

    double* sub_array = new double[col * stan_R];
    int* n_Element = new int[p_size];
    int* deplac = new int[p_size];

    deplac[0] = 0;
    for (int i = 0; i < p_size; i++) {
        if (i < restRow) {
            n_Element[i] = (numberRow + 1) * col;
        }
        else {
            n_Element[i] = numberRow * col;
        }
        if (i > 0) {
            deplac[i] = deplac[i - 1] + n_Element[i - 1];
        }
    }

    
    MPI_Scatterv(temp_array, n_Element, deplac, MPI_DOUBLE, sub_array, stan_R * col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* p_row = new double[col];

    for (int i = 0; i < deplac[rank] / col; i++) {
        int rac = 0;
        int sum = 0;
        for (int j = 0; j < p_size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) {
                rac = j;
                break;
            }
            rac++;
        }

        MPI_Bcast(p_row, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);

        for (int j = 0; j < n_Element[rank] / col; j++) {
            double s = p_row[i] / sub_array[j * col + i];
            for (int z = i; z < col; z++) {
                sub_array[j * col + z] = s * sub_array[j * col + z] - p_row[z];
            }
        }
    }

    for (int i = 0; i < n_Element[rank] / col; i++) {
        for (int j = 0; j < col; j++) {
            p_row[j] = sub_array[i * col + j];
        }
        MPI_Bcast(p_row, col, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        for (int j = i + 1; j < n_Element[rank] / col; j++) {
            double s = p_row[deplac[rank] / col + i] / sub_array[j * col + i + deplac[rank] / col];
            for (int z = i + deplac[rank] / col; z < col; z++) {
                sub_array[j * col + z] = s * sub_array[j * col + z] - p_row[z];
            }
        }
    }

    int rac = 0;
    for (int i = deplac[rank] / col + n_Element[rank] / col; i < row; i++) {
        int sum = 0;
        for (int j = 0; j < p_size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) {
                rac = j;
                break;
            }
            rac++;
        }
        MPI_Bcast(p_row, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);
    }

    MPI_Gatherv(sub_array, stan_R * col, MPI_DOUBLE, temp_array, n_Element, deplac, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = row - 1; i >= 0; i--) {
            double b = temp_array[i * col + col - 1];
            for (int j = row - 1; j >= i + 1; j--) {
                b -= temp_array[i * col + j] * solution[j];
            }
            solution[i] = b / temp_array[i * col + i];
        }
    }

    delete[] temp_array;
    delete[] sub_array;
    delete[] n_Element;
    delete[] deplac;
    delete[] p_row;
}
