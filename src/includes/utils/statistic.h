#include "../constants.h"
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

double randUnif(double min, double max)
{
    return min + rand() / (RAND_MAX + 1.) * (max - min);
}

/**
 * @brief Return a random double with a gaussian distribution
 *
 * @param mean
 * @param sigma
 * @return
 */
double randGauss(double mean, double sigma)
{
    double x = randUnif(0, 1);
    double y = randUnif(0, 1);

    return mean + sigma * sqrt(-2 * log(1 - x)) * cos(2 * PI * y);
}

#endif