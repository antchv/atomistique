#include <iostream>
#include <vector>
#include "fcts.h"


// Energie cinétique (2D)
double Ec(const std::vector<double>& Xp, const std::vector<double>& Yp) {
    double energie = 0.;
    int N = Xp.size();

    for(int i = 0; i < N; i++) {
        energie += 0.5 * (Xp[i] * Xp[i] + Yp[i] * Yp[i]);
    }
    return energie;
}

// Energie potentielle de Lennard-Jones entre 2 particules 1 et 2 (2D)
double u12(const double x1, const double y1, const double x2, const double y2) {
    double r2 = 1. / ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    if(r2 < rc2) {
        return 0.;
    }
    else {
        double r6 = r2 * r2 * r2;
        double rc6 = rc2 * rc2 * rc2;
        return 4 * (r6 * r6 - r6) - 4 * (rc6 * rc6 - rc6);
    }
}

// Energie potentielle de Lennard-Jones (2D)
double EpLJ(const std::vector<double>& X, const std::vector<double>& Y) {
    double energie = 0.;
    int N = X.size();

    for(int i = 0; i < N; i++) {
        for(int j = i+1; j < N; j++) {
            energie += u12(X[j]-L, Y[j]-L, X[i], Y[i]);
            energie += u12(X[j]-L, Y[j], X[i], Y[i]);
            energie += u12(X[j]-L, Y[j]+L, X[i], Y[i]);
            energie += u12(X[j], Y[j]-L, X[i], Y[i]);
            energie += u12(X[j], Y[j], X[i], Y[i]);
            energie += u12(X[j], Y[j]+L, X[i], Y[i]);
            energie += u12(X[j]+L, Y[j]-L, X[i], Y[i]);
            energie += u12(X[j]+L, Y[j], X[i], Y[i]);
            energie += u12(X[j]+L, Y[j]+L, X[i], Y[i]);
        }
    }
    return energie;
}

// Energie mécanique de Lennard-Jones (2D)
double EmLJ(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Xp, const std::vector<double>& Yp) {
    return Ec(Xp,Yp) + EpLJ(X,Y);
}

// Force, suivant x, de Lennard-Jones entre 2 particules 1 et 2 (2D)
double f12x(const double x1, const double y1, const double x2, const double y2) {
    double r2 = 1. / ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    if(r2 < rc2) {
        return 0.;
    }
    else {
        double r6 = r2 * r2 * r2;
        return 24 * (2 * r6 * r6 - r6) * (x2 - x1) * r2;
    }
}

// Force, suivant y, de Lennard-Jones entre 2 particules 1 et 2 (2D)
double f12y(const double x1, const double y1, const double x2, const double y2) {
    double r2 = 1. / ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    if(r2 < rc2) {
        return 0.;
    }
    else {
        double r6 = r2 * r2 * r2;
        return 24 * (2 * r6 * r6 - r6) * (y2 - y1) * r2;
    }
}

// Force de lennard-Jones totale, suivant x, qu'exerce toute les particules sur la particule i
double FLJX(int i, const std::vector<double>& X, const std::vector<double>& Y) {
    double Fx = 0.;
    int N = X.size();
    for(int j = 0; j < N; j++) {
        if(i!=j) {
            Fx += f12x(X[j]-L, Y[j]-L, X[i], Y[i]);
            Fx += f12x(X[j]-L, Y[j], X[i], Y[i]);
            Fx += f12x(X[j]-L, Y[j]+L, X[i], Y[i]);
            Fx += f12x(X[j], Y[j]-L, X[i], Y[i]);
            Fx += f12x(X[j], Y[j], X[i], Y[i]);
            Fx += f12x(X[j], Y[j]+L, X[i], Y[i]);
            Fx += f12x(X[j]+L, Y[j]-L, X[i], Y[i]);
            Fx += f12x(X[j]+L, Y[j], X[i], Y[i]);
            Fx += f12x(X[j]+L, Y[j]+L, X[i], Y[i]);
        }
    }
    return Fx;
}

// Force de lennard-Jones totale, suivant y, qu'exerce toute les particules sur la particule i
double FLJY(int i, const std::vector<double>& X, const std::vector<double>& Y) {
    double Fy = 0.;
    int N = X.size();
    for(int j = 0; j < N; j++) {
        if(i!=j) {
            Fy += f12y(X[j]-L, Y[j]-L, X[i], Y[i]);
            Fy += f12y(X[j]-L, Y[j], X[i], Y[i]);
            Fy += f12y(X[j]-L, Y[j]+L, X[i], Y[i]);
            Fy += f12y(X[j], Y[j]-L, X[i], Y[i]);
            Fy += f12y(X[j], Y[j], X[i], Y[i]);
            Fy += f12y(X[j], Y[j]+L, X[i], Y[i]);
            Fy += f12y(X[j]+L, Y[j]-L, X[i], Y[i]);
            Fy += f12y(X[j]+L, Y[j], X[i], Y[i]);
            Fy += f12y(X[j]+L, Y[j]+L, X[i], Y[i]);
        }
    }
    return Fy;
}