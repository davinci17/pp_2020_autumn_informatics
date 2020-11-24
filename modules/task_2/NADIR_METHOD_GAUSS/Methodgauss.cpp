// Copyright 2020 Nadir Mohammed
#include<iostream>
#include<iomanip>
#include<mpi.h>
#include <stddef.h>
#include <math.h>
#include <ctime>
#include "../../../modules/task_2/nadir_method_gauss/MethodGauss.h"

void methodGauss(const double* array, double* solution, int row) {
    double* temp_array = new double[row *(row+1)];
    for (int i = 0; i < row * (row+1); i++) {
        temp_array[i] = array[i];
    }
    int p_Row, marge, r_ind, row;

    double p_Value, p_Factor;

    int* p_Pos = new int[row];

    int* p_Iter = new int[row];

    for (int i = 0; i < row; i++) {p_Iter[i] = -1; p_Pos[i] = 0;}

    for (marge = 0; marge < row; marge++) 

    {double maxValue = 0;
        for (int i = 0; i < row; i++) {
            if ((p_Iter[i] == -1) && ((fabs(temp_array[i * (row + 1) + marge]) > maxValue))) 
            {
        p_Row = i; maxValue = temp_array[i * (row + 1) + marge];
            }
        }
        p_Pos[marge] = p_Row;
        p_Iter[p_Row] = marge;

        p_Value = temp_array[p_Row * (row + 1) + marge];
        for (int i = 0; i < row; i++) {
            if (p_Iter[i] == -1) {
                p_Factor = temp_array[i * (row + 1) + marge] / p_Value;
                for (int j = marge; j < row; j++) {
                    temp_array[i * (row + 1) + j] -= p_Factor * temp_array[p_Row * (row + 1) + j];
                }
                temp_array[i * (row + 1) + row] -= p_Factor * temp_array[p_Row * (row + 1) + row];}
        }
    }

    for (int i = (row - 1); i >= 0; i--) {
        r_ind = p_Pos[i];
        solution[i] = temp_array[r_ind * (row + 1) + row] / temp_array[r_ind * (row + 1) + i];
        temp_array[r_ind * (row + 1) + i] = 1;

        for (int j = 0; j < i; j++) {row = p_Pos[j];
            temp_array[row * (row + 1) + row] -= temp_array[row * (row + 1) + i] * solution[i];
            temp_array[row * (row + 1) + i] = 0;}
    }
    delete[] temp_array;
    delete[] p_Pos;
    delete[] p_Iter;
}

void methodGaussParallel(const double* array, double* solution, int row,int col) {
    int size, rank;
    double* array_temp = new double[row * col];
    for (int i = 0; i < row * col; i++) { array_temp[i] = array[i]; }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n_Row = row / size;
    int rest_Row = row % size;

    int stan_Row;
    if (rank < rest_Row) { stan_Row = n_Row + 1; }
    else { stan_Row = n_Row; }

    double* sub_array = new double[col * stan_Row];
    int* n_Element = new int[size];
    int* n_displac = new int[size];

    n_displac[0] = 0;
    for (int i = 0; i < size; i++) {

        if (i < rest_Row) { n_Element[i] = (n_Row + 1) * col; }
        else { 
            n_Element[i] = n_Row * col; 
             }
        if (i > 0) { 
            n_displac[i] = n_displac[i - 1] + n_Element[i - 1]; 
        }
    }
    MPI_Scatterv(array_temp, n_Element, n_displac, MPI_DOUBLE, sub_array, stan_Row * col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* part_trans = new double[col];

    for (int i = 0; i < n_displac[rank] / col; i++)
    {
        int rac = 0, sum = 0;
        for (int j = 0; j < size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) { 
                rac = j; break;
            }  
            rac++;
        }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);

        for (int j = 0; j < n_Element[rank] / col; j++) {
            double s = part_trans[i] / sub_array[j * col + i];
            for (int z = i; z < col; z++) { 
                sub_array[j * col + z] = s * sub_array[j * col + z] - part_trans[z]; 
            }
        }
    }

    for (int i = 0; i < n_Element[rank] / col; i++) {
        for (int j = 0; j < col; j++) {
            part_trans[j] = sub_array[i * col + j];
        }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        for (int j = i + 1; j < n_Element[rank] / col; j++) {
            double s = part_trans[n_displac[rank] / col + i] / sub_array[j * col + i + n_displac[rank] / col];
            for (int z = i + n_displac[rank] / col; z < col; z++) 
            {
                sub_array[j * col + z] = s * sub_array[j * col + z] - part_trans[z];
            }
        }
    }

    int rac = 0;
    for (int i = n_displac[rank] / col + n_Element[rank] / col; i < row; i++) {
        int sum = 0;
        for (int j = 0; j < size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) { rac = j; break; 
            }
            rac++;
        }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);
    }

    MPI_Gatherv(sub_array, stan_Row * col, MPI_DOUBLE, array_temp, n_Element, n_displac, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = row - 1; i >= 0; i--) 
        {
            double b = array_temp[i * col + col - 1];
            for (int j = row - 1; j >= i + 1; j--)
            {
                b -= array_temp[i * col + j] * solution[j];
            }
            solution[i] = b / array_temp[i * col + i];
        }
    }
   
  
    delete[] array_temp;
    delete[] sub_array;
    delete[] n_Element;
    delete[] n_displac;
    delete[] part_trans;
   
}  
  
