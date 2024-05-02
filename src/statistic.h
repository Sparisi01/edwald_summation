#include <math.h>

double mean(const double* array, size_t start_index, size_t end_index) {
    double sum = 0;
    for (size_t i = start_index; i < end_index; i++) {
        sum += array[i];
    }

    return sum / (end_index - start_index);
}

double var(const double* array, size_t start_index, size_t end_index) {
    double m = mean(array, start_index, end_index);
    double sum = 0;
    for (size_t i = start_index; i < end_index; i++) {
        sum += (array[i] - m) * (array[i] - m);
    }
    return sum / (end_index - start_index);
}

double stddev(const double* array, size_t start_index, size_t end_index) {
    return sqrt(var(array, start_index, end_index));
}
