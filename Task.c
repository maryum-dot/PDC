#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

#define NUM_SPOTS 5     // Total parking spots
#define NUM_CARS 10     // Total cars trying to park

sem_t parkingLot;       // Counting semaphore

void* car(void* arg) {
    int id = *(int*)arg;

    printf("Car %d is trying to enter the parking lot...\n", id);

    // Car waits if no spot is available
    sem_wait(&parkingLot);
    printf("Car %d has PARKED.\n", id);

    // Stay parked for a while
    sleep(2);

    // Leave and free the spot
    printf("Car %d is LEAVING the parking lot.\n", id);
    sem_post(&parkingLot);

    return NULL;
}

int main() {
    pthread_t cars[NUM_CARS];
    int car_ids[NUM_CARS];

    // Initialize semaphore with number of spots
    sem_init(&parkingLot, 0, NUM_SPOTS);

    // Create car threads
    for (int i = 0; i < NUM_CARS; i++) {
        car_ids[i] = i + 1;
        pthread_create(&cars[i], NULL, car, &car_ids[i]);
        sleep(1); // stagger arrivals for clarity
    }

    // Join car threads
    for (int i = 0; i < NUM_CARS; i++) {
        pthread_join(cars[i], NULL);
    }

    sem_destroy(&parkingLot);
    return 0;
}
