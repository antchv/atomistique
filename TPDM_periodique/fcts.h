#ifndef FCTS_H
#define FCTS_H

extern const double rc2;
extern const double L;

double Ec(const std::vector<double>& Xp, const std::vector<double>& Yp);
double u12(const double x1, const double y1, const double x2, const double y2);
double EpLJ(const std::vector<double>& X, const std::vector<double>& Y);
double EmLJ(const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Xp, const std::vector<double>& Yp);
double f12x(const double x1, const double y1, const double x2, const double y2);
double f12y(const double x1, const double y1, const double x2, const double y2);
double FLJX(int i, const std::vector<double>& X, const std::vector<double>& Y);
double FLJY(int i, const std::vector<double>& X, const std::vector<double>& Y);

#endif