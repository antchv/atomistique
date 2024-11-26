#ifndef METHODE_H
#define METHODE_H

#include <vector>

extern const double dt;
extern const double m;

// Création du type de donnée correspondant à la fonction qui calcule la force
typedef void (*fonctionInteraction)(const double& x1, const double& y1, const double& x2, const double& y2, double& f12x, double& f12y, double& u12);

void euler(fonctionInteraction f, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Vx, std::vector<double>& Vy, std::vector<double>& FXant, std::vector<double>& FYant);
void velocityVerlet(fonctionInteraction f, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Vx, std::vector<double>& Vy, std::vector<double>& FXant, std::vector<double>& FYant);

#endif