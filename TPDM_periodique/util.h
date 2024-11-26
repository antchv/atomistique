#ifndef UTIL_H
#define UTIL_H

extern const double tf;

int compteLigneFichier(const std::string& nomFichier);
std::vector<double> extraireColonne(const std::string& nomFichier, int colonne);
void ecrireEntete(std::ofstream& fichier, double dt, int N, std::string methode);
void ecrirePos(const std::string& nomFichier, const std::vector<std::vector<double>>& Q, int nbPoint, double dt, int N, std::string methode);
void ecrireEnNu(const std::string& nomFichier, const std::vector<double>& x, const std::vector<double>& y, double dt, int N, std::string methode);

#endif