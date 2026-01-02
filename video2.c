#include<stdio.h>
#include<omp.h>

int main(){

    #pragma omp parallel num_threads(4)
    {
        #pragma omp for schedule (static, 3)
        for(int i=0; i<20; i++){

                int tid = omp_get_thread_num();
                // int tid = 1;
                printf("Thread %d is working on iteration %d\n", tid, i);
            
        }
    }
    return 0;
}