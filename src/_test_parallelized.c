#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constants.h"

typedef struct Partial_for_info
{
    int first_index;
    int last_index;

} Partial_for_info;

void *partial_for(void *args)
{
    Partial_for_info arg = *(Partial_for_info *)args;
    for (size_t i = arg.first_index; i < arg.last_index; i++)
    {
        for (size_t j = 0; j < N_PARTICLES; j++)
        {
            sin(j * i);
        }
    }

    return NULL;
}

int main(int argc, char const *argv[])
{

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    const int n_thread = 1;
    pthread_t *thread_ids = (pthread_t *)malloc(n_thread * sizeof(pthread_t));
    if (!thread_ids)
    {
        perror("Error allocating threads");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n_thread; i++)
    {
        Partial_for_info info_tmp = {.first_index = N_PARTICLES / n_thread * i, N_PARTICLES / n_thread * (i + 1)};
        pthread_create(&thread_ids[i], NULL, partial_for, &info_tmp);
    }

    for (int i = 0; i < n_thread; i++)
        pthread_join(thread_ids[i], NULL);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("%lf", cpu_time_used);

    return 0;
}
