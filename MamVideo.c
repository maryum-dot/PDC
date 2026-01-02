#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

#define NUM 10000000

int *a;
int *b;
int *c;

int main(){
    struct timeval s_start, s_end, p_start, p_end;

    gettimeofday(&s_start, NULL);
    a = malloc (NUM * sizeof(int));
    b = malloc (NUM * sizeof(int));
    c = malloc (NUM * sizeof(int));

    for(int i =0; i<NUM; i++){
        a[i]= rand() % 10000;
        b[i] = rand()% 10000;
        c[i]= 0;
    }
    gettimeofday(&s_end, NULL);

    double s_time = (s_end.tv_sec - s_start.tv_sec) + (s_end.tv_usec - s_start.tv_usec)/1000000.0;
    gettimeofday(&p_start, NULL);
    long int sum =0;
    #pragma omp parallel num_threads(8)
    {
    
     #pragma omp for reduction(+: sum)
    for(int i=0; i<NUM; i++){
        c[i] = a[i] * b[i];
        sum=sum+c[i];
    }
    
    
}
    gettimeofday(&p_end, NULL);
    double p_time = (p_end.tv_sec - p_start.tv_sec) + (p_end.tv_usec - p_start.tv_usec)/1000000.0;
    printf("Sum: %ld \nS TIME:%f sec\np TIME: %f sec\nTOTAL TIME: %f\n", sum, s_time, p_time, s_time+p_time);

}
