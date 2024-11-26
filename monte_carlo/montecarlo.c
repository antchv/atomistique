#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
const double boltzmann = 1.38064852e-23;
const double L = 5.0;//size of the box

// Function declarations
void Initialization(int num_particles, double** POSITION_X, double** POSITION_Y);
double Potential(double r);
double Energy(int num_particles, double* POSITION_X, double* POSITION_Y);
void PBC(int num_particles, double* POSITION_X, double* POSITION_Y);

int main() {
    // Define parameters
    const int NUM_PART = 20;
    double T = 200.0;
    int sim_step = 10;
    int cycle = 10;
    FILE* file; // File pointer
    double current_energy;
    double energy_change;
    int particle;

    // Initialize the position of particles and use PBC to ensure particles are in the box
    double* POSITION_X;
    double* POSITION_Y;
    Initialization(NUM_PART, &POSITION_X, &POSITION_Y);
    PBC(NUM_PART, POSITION_X, POSITION_Y);

    // Calculate the initial energy of the system
    double initial_energy = Energy(NUM_PART, POSITION_X, POSITION_Y);

    // Monte Carlo Loop
    for (int i = 0; i < sim_step; i++) {
        for (int j = 0; j < cycle ; j++) {
            // Select a particle
            particle = rand() % NUM_PART;

            // Propose a move
            double dx = ((double)rand() / RAND_MAX) * 0.1; // Adjust the step size as needed
            double dy = ((double)rand() / RAND_MAX) * 0.1;

            // Update positions
            POSITION_X[particle] += dx;
            POSITION_Y[particle] += dy;

            // Apply PBC
            PBC(NUM_PART, POSITION_X, POSITION_Y);

            // Calculate the energy change
            current_energy = Energy(NUM_PART, POSITION_X, POSITION_Y);
            energy_change = current_energy - initial_energy;

            // Accept or reject the move
            if (energy_change < 0) { // When energy change is negative, the system will always accept the move
                initial_energy += energy_change;
            } 
            else { // When energy change is positive, the system will accept the move with a probability of exp(-energy_change/(boltzmann * T))
                double p = ((double)rand() / RAND_MAX);
                if (p < exp(-energy_change / (boltzmann * T))) {
                    initial_energy += energy_change;
                } 
                else {
                    // Reject the move and revert the position change
                    POSITION_X[particle] -= dx;
                    POSITION_Y[particle] -= dy;
                }
            }
        }
            // Open the file for appending data
            file = fopen("output.txt", "a");

            // Check if the file was opened successfully
            if (file == NULL) {
                perror("Unable to open the file");
                return 1;
            }

            // Write data to the file
            fprintf(file, "%d, %lf, %lf, %lf\n", i, POSITION_X[particle], POSITION_Y[particle], current_energy);

            // Close the file
            fclose(file);
    }

    // Free dynamically allocated memory
    free(POSITION_X);
    free(POSITION_Y);

    return 0;
}

void Initialization(int num_particles, double** POSITION_X, double** POSITION_Y) {
    *POSITION_X = (double*)malloc(num_particles * sizeof(double));
    *POSITION_Y = (double*)malloc(num_particles * sizeof(double));

    int num_rows = (int)sqrt(num_particles); // Number of rows in the square grid
    double spacing = L / num_rows; // Spacing between particles

    int particle_index = 0;
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_rows; j++) {
            (*POSITION_X)[particle_index] = (i + 0.5) * spacing; // Add 0.5 for centering
            (*POSITION_Y)[particle_index] = (j + 0.5) * spacing;
            particle_index++;
        }
    }
}


double Potential(double r) {
    // Lennard-Jones potential
    double epsilon = 1.0; // Adjust epsilon as needed
    double potential = 4 * epsilon * (pow(1 / r, 12) - pow(1 / r, 6));
    return potential;
}

double Energy(int num_particles, double* POSITION_X, double* POSITION_Y) {
    double energy = 0;
    double r = 0;
    double Rc = 2.5;
    for (int i = 0; i < num_particles; i++) {
        for (int j = 0; j < num_particles; j++) {
            if (i != j) {
                r = sqrt(pow(POSITION_X[i] - POSITION_X[j], 2) + pow(POSITION_Y[i] - POSITION_Y[j], 2));
                if (r < Rc) {
                    energy += Potential(r);
                }
            }
        }
    }
    return energy;
}

void PBC(int num_particles, double* POSITION_X, double* POSITION_Y) {

    for (int i = 0; i < num_particles; i++) {
        if (POSITION_X[i] > L) {
            POSITION_X[i] -= L;
        }
        if (POSITION_X[i] < 0) {
            POSITION_X[i] += L;
        }
        if (POSITION_Y[i] > L) {
            POSITION_Y[i] -= L;
        }
        if (POSITION_Y[i] < 0) {
            POSITION_Y[i] += L;
        }
    }
}
