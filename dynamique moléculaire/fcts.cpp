#include <iostream>
#include <vector>
#include <cmath>
#include "fcts.h"


// Interaction de Lennard-Jones (2D)
void lennardJones(const double& x1, const double& y1, const double& x2, const double& y2, double& f12x, double& f12y, double& u12) {
    double r2 = 1./((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2));
    double r6 = r2*r2*r2;
    double r12 = r6*r6;
    double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    double sigma12 = sigma6 * sigma6;

    f12x = 24 * epsilon * (2 * sigma12 * r12 - sigma6 * r6) * (x2 - x1) * r2;
    f12y = 24 * epsilon * (2 * sigma12 * r12 - sigma6 * r6) * (y2 - y1) * r2;
    u12 = 4 * epsilon * (sigma12 * r12 - sigma6 * r6);
}

// Actualisation des Forces (2D)
void calculForces(fonctionInteraction f, const std::vector<double>& X, const std::vector<double>& Y, std::vector<double>& Fx, std::vector<double>& Fy) {
    int N = X.size();
    double f12x, f12y, u12;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(j != i) {
                f(X[i],Y[i],X[j],Y[j],f12x,f12y,u12);
                Fx[i] += f12x;
                Fy[i] += f12y;
            }
        }
    }
}

// Energie potentielle (2D)
double Ep(fonctionInteraction f, const std::vector<double>& X, const std::vector<double>& Y) {
    double f12x, f12y, u12, energie = 0.;
    int N = X.size();

    for(int i = 0; i < N; i++) {
        for(int j = i+1; j < N; j++) {
            f(X[i],Y[i],X[j],Y[j],f12x,f12y,u12);
            energie += u12;
        }
    }
    return energie;
}

// Energie cinétique (2D)
double Ec(const std::vector<double>& Vx, const std::vector<double>& Vy) {
    double energie = 0.;
    int N = Vx.size();

    for(int i = 0; i < N; i++) {
        energie += 0.5 * m * (Vx[i] * Vx[i] + Vy[i] * Vy[i]);
    }
    return energie;
}