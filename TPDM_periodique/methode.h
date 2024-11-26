#ifndef METHODE_H
#define METHODE_H

double eulerQ(const double q, const double qp, const double dt);
double eulerQp(const double qp, const double qpp, const double dt);
double velocityVerletQ(const double q, const double qp, const double qpp, const double dt);
double velocityVerletQp(const double qp, const double qppAv, const double qppAp, const double dt);
void boite(double& q, double L);

#endif