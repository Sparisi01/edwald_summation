#include <math.h>

#ifndef STATISTIC_H
#define STATISTIC_H

double mean(const double *array, size_t start_index, size_t end_index)
{
    double sum = 0;
    for (size_t i = start_index; i < end_index; i++)
    {
        sum += array[i];
    }

    return sum / (end_index - start_index);
}

double var(const double *array, size_t start_index, size_t end_index)
{
    double m = mean(array, start_index, end_index);
    double sum = 0;
    for (size_t i = start_index; i < end_index; i++)
    {
        sum += (array[i] - m) * (array[i] - m);
    }
    return sum / (end_index - start_index);
}

double stddev(const double *array, size_t start_index, size_t end_index)
{
    return sqrt(var(array, start_index, end_index));
}

double randGauss(double sigma)
{
    double x = rand() / (RAND_MAX + 1.);
    double y = rand() / (RAND_MAX + 1.);

    return sigma * sqrt(-2 * log(1 - x)) * cos(2 * PI * y);
}

#endif