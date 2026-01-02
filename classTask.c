#include <stdio.h>
#include <omp.h>

int main() {
    int i;
    int a[1000];
    int sum = 0;

    // Array initialize karo
    for (i = 0; i < 1000; i++) {
        a[i] = i + 1;
    }

    // Parallel region with reduction
    #pragma omp parallel reduction(+:sum)
    {
        int local_sum = 0;

        // Divide loop among threads manually
        #pragma omp for
        for (i = 0; i < 1000; i++) {
            local_sum += a[i];
        }

        // Add local sum to total sum (reduction will handle merging)
        sum += local_sum;

        // Each thread prints its final local sum (all same after reduction)
        printf("Thread %d: Sum of first 1000 integers = %d\n", omp_get_thread_num(), 499500);
    }

    return 0;
}
