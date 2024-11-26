#include <iostream>
#include <cmath>
#include "methode.h"


// Méthode d'Euler pour déterminer la position à l'instant d'après
double eulerQ(const double q, const double qp, const double dt) {
    return q + qp * dt;
}

// Méthode d'Euler pour déterminer la vitesse à l'instant d'après
double eulerQp(const double qp, const double qpp, const double dt) {
    return qp + qpp * dt;
}

// Méthode de velocity-Verlet pour déterminer la position à l'instant d'après
double velocityVerletQ(const double q, const double qp, const double qpp, const double dt) {
    return q + qp * dt + 0.5 * qpp * dt * dt;
}

// Méthode de velocity-Verlet pour déterminer la vitesse à l'instant d'après
double velocityVerletQp(const double qp, const double qppAv, const double qppAp, const double dt) {
    return qp + 0.5 * (qppAv + qppAp) * dt;
}

// Méthode pour remettre la particule sortie de la boite dans la boite
void boite(double& q, double L) {
    q = fmod(q, L);
    if(q < -L/2) {q += L;}
    if(q > L/2) {q -= L;}
}