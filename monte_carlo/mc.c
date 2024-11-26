#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define N 1280
#define SIZE 10.0
#define STEPS 1000
#define BETA 0.01

typedef struct {
    double x;
    double y;
} Position;

Position positions[N];
double energie_actuelle;
double trajectoires_x[STEPS][N];
double trajectoires_y[STEPS][N];

// Initialisation des positions avec conditions périodiques
void initialiser_positions() {
    for (int i = 0; i < N; i++) {
        positions[i].x = rand() / (double)RAND_MAX * SIZE;
        positions[i].y = rand() / (double)RAND_MAX * SIZE;
    }
}

// Déplacement avec conditions périodiques
void deplacer_particule(Position *position, double delta_x, double delta_y) {
    position->x += delta_x;
    position->y += delta_y;

    // Appliquer les conditions aux limites périodiques
    position->x = fmod(position->x + SIZE, SIZE);
    position->y = fmod(position->y + SIZE, SIZE);
}

double calculer_energie_potentielle() {
    double energie = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = positions[i].x - positions[j].x;
            double dy = positions[i].y - positions[j].y;
            double distance = sqrt(dx * dx + dy * dy);

            // Potentiel de Lennard-Jones
            double r_over_sigma = 1 / distance;
            double r_over_sigma_sqr = r_over_sigma * r_over_sigma;
            double r_over_sigma_sixth = r_over_sigma_sqr * r_over_sigma_sqr * r_over_sigma_sqr;
            double lj_potential = 4.0 * 1 * (r_over_sigma_sixth * r_over_sigma_sixth - r_over_sigma_sqr);

            energie += lj_potential;
        }
    }

    return energie;
}

// Sauvegarde des positions dans des fichiers CSV distincts pour x et y
void sauvegarder_trajectoires_csv() {
    FILE *fichier_x = fopen("trajectoires_x.csv", "w");
    FILE *fichier_y = fopen("trajectoires_y.csv", "w");

    if (fichier_x == NULL || fichier_y == NULL) {
        perror("Erreur lors de l'ouverture des fichiers");
        exit(EXIT_FAILURE);
    }

    for (int step = 0; step < STEPS; step++) {
        for (int i = 0; i < N; i++) {
            fprintf(fichier_x, "%lf,", trajectoires_x[step][i]);
            fprintf(fichier_y, "%lf,", trajectoires_y[step][i]);
        }
        fprintf(fichier_x, "\n");
        fprintf(fichier_y, "\n");
    }

    fclose(fichier_x);
    fclose(fichier_y);
}

int main() {
    srand(42); // Initialisation du générateur de nombres aléatoires

    initialiser_positions();
    energie_actuelle = calculer_energie_potentielle();
    clock_t start_time = clock();  // Start measuring time
    // Boucle Monte Carlo avec conditions périodiques
    for (int step = 0; step < STEPS; step++) {
        int particule = rand() % N;

        double delta_x = (rand() % 3) - 1;
        double delta_y = (rand() % 3) - 1;

        // Sauvegarde des positions pour les trajectoires
        for (int i = 0; i < N; i++) {
            trajectoires_x[step][i] = positions[i].x;
            trajectoires_y[step][i] = positions[i].y;
        }

        deplacer_particule(&positions[particule], delta_x, delta_y);

        double nouvelle_energie = calculer_energie_potentielle();

        double delta_energie = nouvelle_energie - energie_actuelle;

        if (delta_energie < 0 || rand() / (double)RAND_MAX < exp(-delta_energie * BETA)) {
            energie_actuelle = nouvelle_energie;
        } else {
            // Restauration des positions précédentes si le déplacement est rejeté
            positions[particule].x -= delta_x;
            positions[particule].y -= delta_y;
        }
    }
    clock_t end_time = clock();  // Stop measuring time
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", cpu_time_used);
    // Sauvegarde des trajectoires dans des fichiers CSV distincts
    sauvegarder_trajectoires_csv();

    // Le reste du code pour l'analyse des résultats...

    return 0;
}
