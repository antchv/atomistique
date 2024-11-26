#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
const double boltzmann = 1.38064852e-23;
double L = 5.0; // size of the box

// Function declarations
void Initialization(int num_particles, double **POSITION_X, double **POSITION_Y);
double Potential(double r);
double Energy(int num_particles, double *POSITION_X, double *POSITION_Y);
void Forces(int num_particles, double *POSITION_X, double *POSITION_Y, double *FORCE_X, double *FORCE_Y);
double Virial(int num_particles, double *POSITION_X, double *POSITION_Y, double *FORCE_X, double *FORCE_Y);
double Pressure(int num_particles, double virial);
void PBC(int num_particles, double *POSITION_X, double *POSITION_Y);

int main()
{
    // Define parameters
    const int NUM_PART = 128;
    double T = 200.0;
    int sim_steps = 10;
    int cycles = 10;
    FILE *file; // File pointer

    // Variable to store the pressure
    double pressure;

    // Monte Carlo Loop
    for (double current_L = L; current_L <= 8.0; current_L += 0.1)
    {
        // Set the current size of the box
        L = current_L;

        // Initialize the position of particles and use PBC to ensure particles are in the box
        double *POSITION_X;
        double *POSITION_Y;
        Initialization(NUM_PART, &POSITION_X, &POSITION_Y);
        PBC(NUM_PART, POSITION_X, POSITION_Y);

        // Allocate memory for forces
        double *FORCE_X;
        double *FORCE_Y;
        FORCE_X = (double *)malloc(NUM_PART * sizeof(double));
        FORCE_Y = (double *)malloc(NUM_PART * sizeof(double));

        // Calculate the initial energy of the system
        double initial_energy = Energy(NUM_PART, POSITION_X, POSITION_Y);

        // Perform all simulation steps for the current volume
        for (int i = 0; i < sim_steps; i++)
        {
            for (int j = 0; j < cycles; j++)
            {
                // Select a particle
                int particle = rand() % NUM_PART;

                // Propose a move
                double dx = ((double)rand() / RAND_MAX) * 0.1; // Adjust the step size as needed
                double dy = ((double)rand() / RAND_MAX) * 0.1;

                // Update positions
                POSITION_X[particle] += dx;
                POSITION_Y[particle] += dy;

                // Apply PBC
                PBC(NUM_PART, POSITION_X, POSITION_Y);

                // Calculate the energy change
                double current_energy = Energy(NUM_PART, POSITION_X, POSITION_Y);
                double energy_change = current_energy - initial_energy;

                // Accept or reject the move
                if (energy_change < 0)
                { // When energy change is negative, the system will always accept the move
                    initial_energy += energy_change;
                }
                else
                { // When energy change is positive, the system will accept the move with a probability of exp(-energy_change/(boltzmann * T))
                    double p = ((double)rand() / RAND_MAX);
                    if (p < exp(-energy_change / (boltzmann * T)))
                    {
                        initial_energy += energy_change;
                    }
                    else
                    {
                        // Reject the move and revert the position change
                        POSITION_X[particle] -= dx;
                        POSITION_Y[particle] -= dy;
                    }
                }
            }

            // Calculate the forces
            Forces(NUM_PART, POSITION_X, POSITION_Y, FORCE_X, FORCE_Y);

            // Calculate the virial
            double virial = Virial(NUM_PART, POSITION_X, POSITION_Y, FORCE_X, FORCE_Y);

            // Calculate the pressure
            pressure = Pressure(NUM_PART, virial);

            // Free dynamically allocated memory for forces after each cycle
            free(FORCE_X);
            free(FORCE_Y);

            // Reallocate memory for forces for the next cycle
            FORCE_X = (double *)malloc(NUM_PART * sizeof(double));
            FORCE_Y = (double *)malloc(NUM_PART * sizeof(double));
        }

        // Open the file for appending data
        file = fopen("output.txt", "a");

        // Check if the file was opened successfully
        if (file == NULL)
        {
            perror("Unable to open the file");
            return 1;
        }

        // Write data to the file (only pressure and volume)
        fprintf(file, "%lf, %lf\n", current_L, pressure);

        // Close the file
        fclose(file);

        // Free dynamically allocated memory for each volume
        free(POSITION_X);
        free(POSITION_Y);
    }

    return 0;
}

void Initialization(int num_particles, double **POSITION_X, double **POSITION_Y)
{
    *POSITION_X = (double *)malloc(num_particles * sizeof(double));
    *POSITION_Y = (double *)malloc(num_particles * sizeof(double));

    int num_rows = (int)sqrt(num_particles); // Number of rows in the square grid
    double spacing = L / num_rows;            // Spacing between particles

    int particle_index = 0;
    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_rows; j++)
        {
            (*POSITION_X)[particle_index] = (i + 0.5) * spacing; // Add 0.5 for centering
            (*POSITION_Y)[particle_index] = (j + 0.5) * spacing;
            particle_index++;
        }
    }
}

double Potential(double r)
{
    // Lennard-Jones potential
    double epsilon = 1.0; // Adjust epsilon as needed
    double potential = 4 * epsilon * (pow(1 / r, 12) - pow(1 / r, 6));
    return potential;
}

double Energy(int num_particles, double *POSITION_X, double *POSITION_Y)
{
    double energy = 0;
    double r = 0;
    double Rc = 2.5;
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < num_particles; j++)
        {
            if (i != j)
            {
                r = sqrt(pow(POSITION_X[i] - POSITION_X[j], 2) + pow(POSITION_Y[i] - POSITION_Y[j], 2));
                if (r < Rc)
                {
                    energy += Potential(r);
                }
            }
        }
    }
    return energy;
}

void PBC(int num_particles, double *POSITION_X, double *POSITION_Y)
{
    for (int i = 0; i < num_particles; i++)
    {
        if (POSITION_X[i] > L)
        {
            POSITION_X[i] -= L;
        }
        if (POSITION_X[i] < 0)
        {
            POSITION_X[i] += L;
        }
        if (POSITION_Y[i] > L)
        {
            POSITION_Y[i] -= L;
        }
        if (POSITION_Y[i] < 0)
        {
            POSITION_Y[i] += L;
        }
    }
}

void Forces(int num_particles, double *POSITION_X, double *POSITION_Y, double *FORCE_X, double *FORCE_Y)
{
    // Initialize forces to zero
    for (int i = 0; i < num_particles; i++)
    {
        FORCE_X[i] = 0.0;
        FORCE_Y[i] = 0.0;
    }

    // Calculate forces
    double epsilon = 1.0; // Adjust epsilon as needed
    double sigma = 1.0;   // Adjust sigma as needed
    double Rc = 2.5;
    double r, r_inv, r_inv2, r_inv6, r_inv12;
    double dx, dy;

    for (int i = 0; i < num_particles; i++)
    {
        for (int j = i + 1; j < num_particles; j++)
        {
            dx = POSITION_X[j] - POSITION_X[i];
            dy = POSITION_Y[j] - POSITION_Y[i];

            // Apply PBC
            if (dx > 0.5 * L)
                dx -= L;
            if (dx < -0.5 * L)
                dx += L;
            if (dy > 0.5 * L)
                dy -= L;
            if (dy < -0.5 * L)
                dy += L;

            r = sqrt(dx * dx + dy * dy);

            if (r < Rc)
            {
                const double epsilon_numerical = 1e-12;
                r_inv = 1.0 / fmax(r, epsilon_numerical);
                //r_inv = 1.0 / r;
                r_inv2 = r_inv * r_inv;
                r_inv6 = r_inv2 * r_inv2 * r_inv2;
                r_inv12 = r_inv6 * r_inv6;

                double force_scalar = 24.0 * epsilon * (2.0 * sigma * sigma * r_inv12 - sigma * sigma * r_inv6);

                FORCE_X[i] += force_scalar * dx;
                FORCE_Y[i] += force_scalar * dy;

                FORCE_X[j] -= force_scalar * dx;
                FORCE_Y[j] -= force_scalar * dy;
            }
        }
    }
}

double Virial(int num_particles, double *POSITION_X, double *POSITION_Y, double *FORCE_X, double *FORCE_Y)
{
    double virial = 0.0;

    for (int i = 0; i < num_particles; i++)
    {
        virial += FORCE_X[i] * POSITION_X[i] + FORCE_Y[i] * POSITION_Y[i];
    }

    return virial;
}

double Pressure(int num_particles, double virial)
{
    return virial / (2.0 * L * L);
}
