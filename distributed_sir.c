/* distributed_sir.c
   Distributed SIR epidemic simulation across N cities using MPI.
   Compile: mpicc -O2 -o distributed_sir distributed_sir.c
   Run example: mpirun -np 4 ./distributed_sir 100 0.3 0.1 0.01
   args: <timesteps> <beta> <gamma> <travel_fraction>
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Default parameters (can be overridden by command-line)
    int timesteps = 200;
    double beta = 0.3;   // infection rate
    double gamma = 0.1;  // recovery rate
    double travel_frac = 0.01; // fraction of city population that may travel from i to each neighbor per timestep

    if (argc >= 2) timesteps = atoi(argv[1]);
    if (argc >= 3) beta = atof(argv[2]);
    if (argc >= 4) gamma = atof(argv[3]);
    if (argc >= 5) travel_frac = atof(argv[4]);

    // Each city initial population (for demo: all equal, but easy to change)
    int local_pop = 10000 + (rank * 100); // slight variation by rank
    double S = (double)(local_pop - 1); // susceptible
    double I = 1.0 * (rank == 0 ? 10 : 0); // start with 10 infected in city 0 (rank 0)
    double R = 0.0;

    // Normalize if rank 0 had extra I
    if (rank == 0) S = local_pop - I;

    // Build travel matrix (broadcast from root)
    // We'll use a simple ring + small fraction to all others:
    // travel_rate[i][j] fraction of population of i that moves to j each timestep.
    double *travel_flat = NULL; // flattened matrix of size nprocs*nprocs

    if (rank == 0) {
        travel_flat = (double*) malloc(sizeof(double) * nprocs * nprocs);
        // Initialize
        for (int i = 0; i < nprocs; ++i) {
            double base = 0.0;
            for (int j = 0; j < nprocs; ++j) {
                travel_flat[i * nprocs + j] = 0.0;
            }
            // simple pattern:
            // send travel_frac to the next and previous city (ring)
            int next = (i + 1) % nprocs;
            int prev = (i - 1 + nprocs) % nprocs;
            travel_flat[i * nprocs + next] = travel_frac;
            travel_flat[i * nprocs + prev] = travel_frac;
            // small uniform trickle to all other cities (excluding self) but keep totals reasonable
            double trickle = travel_frac / 10.0;
            for (int j = 0; j < nprocs; ++j) {
                if (j != i && j != next && j != prev) {
                    travel_flat[i * nprocs + j] = trickle;
                }
            }
            // ensure self stays positive fraction (people who stay)
            double out_sum = 0.0;
            for (int j = 0; j < nprocs; ++j) out_sum += travel_flat[i * nprocs + j];
            if (out_sum >= 1.0) {
                // scale down all outgoing fractions proportionally
                double scale = 0.5 / out_sum;
                for (int j = 0; j < nprocs; ++j) travel_flat[i * nprocs + j] *= scale;
                out_sum = 0.0;
                for (int j = 0; j < nprocs; ++j) out_sum += travel_flat[i * nprocs + j];
            }
            travel_flat[i * nprocs + i] = 1.0 - out_sum;
        }
    } else {
        // allocate local buffer to receive broadcast
        travel_flat = (double*) malloc(sizeof(double) * nprocs * nprocs);
    }

    // Broadcast the flattened travel matrix to all processes
    MPI_Bcast(travel_flat, nprocs * nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Buffers for gathering S,I,R from all processes
    double *S_all = (double*) malloc(sizeof(double) * nprocs);
    double *I_all = (double*) malloc(sizeof(double) * nprocs);
    double *R_all = (double*) malloc(sizeof(double) * nprocs);

    // Print header from root
    if (rank == 0) {
        printf("Distributed SIR simulation: processes=%d timesteps=%d beta=%.3f gamma=%.3f travel_frac=%.4f\n",
               nprocs, timesteps, beta, gamma, travel_frac);
        printf("Columns: timestep total_S total_I total_R\n");
    }

    // Main time-stepping loop
    for (int t = 0; t <= timesteps; ++t) {
        // Gather S,I,R from all processes
        MPI_Allgather(&S, 1, MPI_DOUBLE, S_all, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(&I, 1, MPI_DOUBLE, I_all, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(&R, 1, MPI_DOUBLE, R_all, 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // Compute global totals (for reporting)
        double local_total_S = S;
        double local_total_I = I;
        double local_total_R = R;
        double total_S = 0.0, total_I = 0.0, total_R = 0.0;
        MPI_Reduce(&local_total_S, &total_S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_total_I, &total_I, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_total_R, &total_R, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("%4d %12.1f %12.1f %12.1f\n", t, total_S, total_I, total_R);
            fflush(stdout);
        }

        if (t == timesteps) break;

        // --- Movement step ---
        // After gathering S_all/I_all/R_all, compute incoming populations for this city (rank)
        double S_in = 0.0, I_in = 0.0, R_in = 0.0;
        for (int src = 0; src < nprocs; ++src) {
            double frac = travel_flat[src * nprocs + rank];
            S_in += S_all[src] * frac;
            I_in += I_all[src] * frac;
            R_in += R_all[src] * frac;
        }

        // Update local S,I,R to incoming values (movement done)
        S = S_in;
        I = I_in;
        R = R_in;

        // --- SIR local transitions (after movement) ---
        double local_pop_double = S + I + R;
        if (local_pop_double <= 0.0) continue;

        double new_infections = beta * S * I / local_pop_double;
        double new_recoveries = gamma * I;

        // Clamp to available counts
        if (new_infections > S) new_infections = S;
        if (new_recoveries > I + new_infections) new_recoveries = I + new_infections; // avoid negative

        S -= new_infections;
        I += new_infections - new_recoveries;
        R += new_recoveries;

        // Numerical stability: avoid tiny negatives due to fp errors
        if (S < 1e-9) S = 0.0;
        if (I < 1e-9) I = 0.0;
        if (R < 1e-9) R = 0.0;
    }

    free(travel_flat);
    free(S_all);
    free(I_all);
    free(R_all);

    MPI_Finalize();
    return 0;
}
