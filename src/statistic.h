#include <math.h>

#ifndef STATISTIC_H
#define STATISTIC_H

double mean(double *array, size_t start_index, size_t end_index)
{
    double sum = 0;
    for (size_t i = start_index; i <= end_index; i++)
    {
        sum += array[i];
    }

    return sum / (end_index - start_index);
}

double var(double *array, size_t start_index, size_t end_index)
{
    double m = mean(array, start_index, end_index);
    double sum = 0;
    for (size_t i = start_index; i <= end_index; i++)
    {
        sum += (array[i] - m) * (array[i] - m);
    }
    return sum / (end_index - start_index);
}

double stddev(double *array, size_t start_index, size_t end_index)
{
    return sqrt(var(array, start_index, end_index));
}

// Return a random value gaussian distributed
double randGauss(double mean, double sigma)
{
    double x = rand() / (RAND_MAX + 1.);
    double y = rand() / (RAND_MAX + 1.);

    return sigma * sqrt(-2 * log(1 - x)) * cos(2 * PI * y) + mean;
}

#endif