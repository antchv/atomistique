#ifndef FCTS_H
#define FCTS_H

#include <vector>

extern const double epsilon;
extern const double sigma;
extern const double m;

// Création du type de donnée correspondant à la fonction qui calcule la force
typedef void (*fonctionInteraction)(const double& x1, const double& y1, const double& x2, const double& y2, double& f12x, double& f12y, double& u12);

void lennardJones(const double& x1, const double& y1, const double& x2, const double& y2, double& f12x, double& f12y, double& u12);
void calculForces(fonctionInteraction f, const std::vector<double>& X, const std::vector<double>& Y, std::vector<double>& Fx, std::vector<double>& Fy);
double Ep(fonctionInteraction f, const std::vector<double>& X, const std::vector<double>& Y);
double Ec(const std::vector<double>& Vx, const std::vector<double>& Vy);

#endif