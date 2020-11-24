#include<iostream>
#include<vector>
#include<random>
#include<iomanip>
#include<mpi.h>
#include "../../../modules/task_2/nadir_method_gauss/methodGauss.h"


int row, col;
std::vector<std::vector<double>> methodGauss(std::vector<std::vector<double>>array) {

    double* x = new double[row];
    double temp, s;
    for (int j = 0; j < row - 1; j++)
    {
        for (int i = j + 1; i < row; i++) {
            temp = array[i][j] / array[j][j];
            for (int k = 0; k < row + 1; k++) { array[i][k] -= array[j][k] * temp; }
        }
    }
    std::cout << "\n --------------------------------------------------------------------------\n";
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < row + 1; j++)
        {
            std::cout << std::setw(8) << std::setprecision(4) << std::fixed << array[i][j];
        } std::cout << std::endl;
    }
    std::cout << "\n --------------------------------------------------------------------------\n";
    for (int i = row - 1; i >= 0; i--)
    {
        s = 0;
        for (int j = i + 1; j < row; j++) {
            s += array[i][j] * x[j];
            x[i] = (array[i][row] - s) / array[i][i];
        }
    }
    return array;
}

double* methodGaussParallel(const double* array, double* solution) {
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
        else { n_Element[i] = n_Row * col; }
        if (i > 0) { n_displac[i] = n_displac[i - 1] + n_Element[i - 1]; }
    }
    MPI_Scatterv(array_temp, n_Element, n_displac, MPI_DOUBLE, sub_array, stan_Row * col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* part_trans = new double[col];

    for (int i = 0; i < n_displac[rank] / col; i++)
    {
        int rac = 0, sum = 0;
        for (int j = 0; j < size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) { rac = j; break; }  rac++;
        }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);

        for (int j = 0; j < n_Element[rank] / col; j++) {
            double s = part_trans[i] / sub_array[j * col + i];
            for (int z = i; z < col; z++) { sub_array[j * col + z] = s * sub_array[j * col + z] - part_trans[z]; }
        }
    }

    for (int i = 0; i < n_Element[rank] / col; i++) {
        for (int j = 0; j < col; j++) { part_trans[j] = sub_array[i * col + j]; }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        for (int j = i + 1; j < n_Element[rank] / col; j++) {
            double s = part_trans[n_displac[rank] / col + i] / sub_array[j * col + i + n_displac[rank] / col];
            for (int z = i + n_displac[rank] / col; z < col; z++) {
                sub_array[j * col + z] = s * sub_array[j * col + z] - part_trans[z];
            }
        }
    }

    int rac = 0;
    for (int i = n_displac[rank] / col + n_Element[rank] / col; i < row; i++) {
        int sum = 0;
        for (int j = 0; j < size; j++) {
            sum += n_Element[j] / col;
            if (i < sum) { rac = j; break; }
            rac++;
        }

        MPI_Bcast(part_trans, col, MPI_DOUBLE, rac, MPI_COMM_WORLD);
    }

    MPI_Gatherv(sub_array, stan_Row * col, MPI_DOUBLE, array_temp, n_Element, n_displac, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = row - 1; i >= 0; i--) {
            double b = array_temp[i * col + col - 1];
            for (int j = row - 1; j >= i + 1; j--) {
                b -= array_temp[i * col + j] * solution[j];
            }
            solution[i] = b / array_temp[i * col + i];
        }
    }
    std::cout << "\n --------------------------------------------------------------------------\n";
    for (int i = 0; i < row * col; i++)
    {
        std::cout << std::setprecision(4) << std::setw(8) << std::fixed << array_temp[i];
        if ((i + 1) % col == 0) { std::cout << std::endl; }
    }
    std::cout << "\n --------------------------------------------------------------------------\n";
  
    delete[] array_temp;
    delete[] sub_array;
    delete[] n_Element;
    delete[] n_displac;
    delete[] part_trans;
    return array_temp;
}
