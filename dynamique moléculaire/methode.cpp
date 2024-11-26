#include <iostream>
#include <vector>
#include <cmath>
#include "methode.h"
#include "fcts.h"


// Méthode d'Euler explicite (2D)
void euler(fonctionInteraction f, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Vx, std::vector<double>& Vy, std::vector<double>& FXant, std::vector<double>& FYant) {
    int N = X.size();

    // Actualisation des positions
    for(int i = 0; i < N; i++) {
        X[i] += Vx[i] * dt;
        Y[i] += Vy[i] * dt;
    }

    // Actualisation des vitesses
    for(int i = 0; i < N; i++) {
        Vx[i] += FXant[i] / m * dt;
        Vy[i] += FYant[i] / m * dt;
    }

    // Actualisation de la force
    calculForces(f,X,Y,FXant,FYant);
}

// Méthode de "Velocity-Verlet" (2D)
void velocityVerlet(fonctionInteraction f, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Vx, std::vector<double>& Vy, std::vector<double>& FXant, std::vector<double>& FYant) {
    int N = X.size();
    double tm = dt / 2 / m;
    
    // Actualisation des positions
    for(int i = 0; i < N; i++) {
        X[i] += Vx[i] * dt + FXant[i] * tm * dt;
        Y[i] += Vy[i] * dt + FYant[i] * tm * dt;
    }

    // Création de la force postérieure
    std::vector<double> FXpost(N,0), FYpost(N,0);
    calculForces(f,X,Y,FXpost,FYpost);

    // Actualisation des vitesses
    for(int i = 0; i < N; i++) {
        Vx[i] += (FXpost[i] + FXpost[i]) * tm;
        Vy[i] += (FYpost[i] + FYpost[i]) * tm;
    }

    // Actualisation des forces
    FXant = FXpost;
    FYant = FYpost;
}