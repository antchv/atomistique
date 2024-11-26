#ifndef UTIL_H
#define UTIL_H

std::vector<double> extraireColonne(const std::string& nomFichier, int colonne);
int nombreDeLignesDansCSV(const std::string& nomFichier);
void remplirCC(std::string& nom);
std::vector<std::string> avoirFichiersDossier(const std::string& cheminDossier);
void ecriturePos(const std::string& nomFichier, const std::vector<double>& X, const std::vector<double>& Y);
void ecritureEn(const std::string& nomFichier, const double& temps, const double& energie);
bool suppFichier(const std::string& nomFichier);

#endif