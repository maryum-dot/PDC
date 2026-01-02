#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

#define ARRAY_SIZE 10000

int arr[ARRAY_SIZE];
long long partial_sums[20]; // To store thread results (max 20 threads)
int numThreads;

// Function to generate random numbers in array
void generateArray() {
    for (int i = 0; i < ARRAY_SIZE; i++) {
        arr[i] = rand() % 100 + 1; // Numbers between 1 and 100
    }
}

// Serial sum
long long serialSum() {
    long long sum = 0;
    for (int i = 0; i < ARRAY_SIZE; i++) {
        sum += arr[i];
    }
    return sum;
}

// Structure to pass data to threads
typedef struct {
    int start;
    int end;
    int index; // To store result in partial_sums[index]
} ThreadData;

void* partialSum(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    long long sum = 0;
    for (int i = data->start; i < data->end; i++) {
        sum += arr[i];
    }
    partial_sums[data->index] = sum;
    pthread_exit(NULL);
}

long long parallelSum(int threads) {
    pthread_t tid[threads];
    ThreadData tdata[threads];
    int blockSize = ARRAY_SIZE / threads;

    for (int i = 0; i < threads; i++) {
        tdata[i].start = i * blockSize;
        tdata[i].end = (i == threads - 1) ? ARRAY_SIZE : tdata[i].start + blockSize;
        tdata[i].index = i;
        pthread_create(&tid[i], NULL, partialSum, &tdata[i]);
    }

    for (int i = 0; i < threads; i++) {
        pthread_join(tid[i], NULL);
    }

    long long total = 0;
    for (int i = 0; i < threads; i++) {
        total += partial_sums[i];
    }
    return total;
}

double getTimeInMs() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000.0) + (tv.tv_usec / 1000.0);
}

int main() {
    srand(time(NULL));
    generateArray();

    double start, end;

    // Serial
    start = getTimeInMs();
    long long serialResult = serialSum();
    end = getTimeInMs();
    double serialTime = end - start;

    printf("\nSerial Sum: %lld\n", serialResult);
    printf("Serial Time: %.4f ms\n\n", serialTime);

    int threadCounts[] = {4, 5, 10};

    for (int i = 0; i < 3; i++) {
        int t = threadCounts[i];
        start = getTimeInMs();
        long long parallelResult = parallelSum(t);
        end = getTimeInMs();
        double parallelTime = end - start;

        printf("Threads: %d\n", t);
        printf("Parallel Sum: %lld\n", parallelResult);
        printf("Parallel Time: %.4f ms\n", parallelTime);
        printf("Speedup: %.4f\n\n", serialTime / parallelTime);
    }

    return 0;
}
