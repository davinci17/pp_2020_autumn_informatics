#include "Methodgauss.h"
#include<iostream>
#include<vector>
#include<random>
#include<iomanip>
#include<mpi.h>


void displayVector(std::vector<std::vector<double>> vec) {
    std::cout << std::setprecision(2) << std::fixed;
    for (int i = 0; i < vec.size(); i++) {
        for (int j = 0; j < vec.size(); j++) {
            std::cout << vec[i][j] << "\t";

        }
        std::cout << std::endl;
    }
}
std::vector<std::vector<double>> GenerateRandom()
{
    std::cout << std::setprecision(2) << std::fixed << std::endl;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distr1(2, 5);
    int size = distr1(gen);
    std::vector<std::vector<double>> vec(size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::uniform_real_distribution<double> distr(0.0, 9.0);
            vec[i].push_back(distr(gen));
        }
        std::cout << std::endl;
    }
    // displayVector(vec);
    std::cout << "********************************" << std::endl;
    return vec;
}
std::vector<std::vector<double>> methodGauss(std::vector<std::vector<double>> vec)
{

    vec = GenerateRandom();
    // vector<std::vector<double>> vec{ {1,2} ,{3,4} };
    double result = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i][i] != 0.0) {
            for (int j = 0; j < vec[i].size(); j++)
            {
                if (i != j)
                {
                    result = vec[j][i] / vec[i][i];

                    for (int k = 0; k < vec[i].size(); k++)
                    {
                        vec[j][k] = vec[j][k] - result * vec[i][k];
                    }
                }
            }
        }
    }
    //displayVector(vec);
    return vec;
}

std::vector<std::vector<double>> methodGaussParallel(std::vector<std::vector<double>> vec)
{
    vec = GenerateRandom();

    //displayVector(vec);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int numberRow = vec.size() / size;
    int rest_nRow = vec.size() % size;
    int n_Element = numberRow * vec.size();
    std::vector< std::vector<double>> sub_vec(n_Element);
    for (int i = 0; i < n_Element; i++)
    {
        for (int j = 0; j < vec[i].size(); j++)
        {
            sub_vec[i][j] = vec[i][j];

        }
    }
    std::vector<std::vector<double>> rest_vec(rest_nRow);
    for (int i = numberRow; i < vec.size(); i++)
    {
        for (int j = 0; j < vec[i].size(); j++)
        {
            rest_vec[i][j] = vec[i][j];

        }
    }
    MPI_Scatter(&vec, numberRow * vec.size(), MPI_DOUBLE, &sub_vec, n_Element, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        if (rest_nRow != 0)
            MPI_Bcast(&rest_nRow, rest_nRow * vec.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    std::vector<std::vector<double>> new_vec = methodGauss(vec);

    for (int i = 0; i < n_Element; i++)
    {
        for (int j = 0; j < n_Element; j++)
        {
            sub_vec[i][j] = new_vec[i][j];
        }

    }
    std::vector<std::vector<double>> new_subVec(rest_nRow);
    for (int i = numberRow; i < vec.size(); i++)
    {
        for (int j = 0; j < vec[i].size(); j++)
        {
            rest_vec[i][j] = new_vec[i][j];

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(&sub_vec, vec.size() * numberRow, MPI_DOUBLE, &vec, vec.size() * numberRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        if (rest_nRow != 0)
            MPI_Gather(&rest_vec, vec.size() * rest_nRow, MPI_DOUBLE, &vec[numberRow], vec.size() * numberRow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // displayVector(vec);
    return vec;
}
