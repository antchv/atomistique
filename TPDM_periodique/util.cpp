#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include "util.h"


// Retourne le nombre de lignes d'un fichier
int compteLigneFichier(const std::string& nomFichier) {
    std::ifstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << " ! Impossible d'ouvrir le fichier " << nomFichier << std::endl;
        return -1; // Une valeur négative indique une erreur
    }

    int nbLigne = 0;
    std::string line;
    while (std::getline(fichier, line)) {
        nbLigne++;
    }
    fichier.close();

    return nbLigne;
}

// Retourne la ieme colonne du fichier
std::vector<double> extraireColonne(const std::string& nomFichier, int colonne) {
    std::vector<double> colonne_ieme;

    // Ouvrir le fichier en mode lecture
    std::ifstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << " ! Impossible d'ouvrir le fichier " << nomFichier << std::endl;
        return colonne_ieme; // Retourne une colonne vide en cas d'erreur
    }

    std::string ligne;
    while (std::getline(fichier, ligne)) {
        std::istringstream stream_ligne(ligne);
        double valeur;
        int colonne_actuelle = -1;

        while (stream_ligne >> valeur) {
            colonne_actuelle++;

            // Si la colonne actuelle correspond à la colonne souhaitée, ajouter la valeur à la colonne_ieme
            if (colonne_actuelle == colonne) {
                colonne_ieme.push_back(valeur);
                break; // Passer à la ligne suivante
            }
        }
    }
    fichier.close();

    return colonne_ieme;
}

// écriture de l'entête
void ecrireEntete(std::ofstream& fichier, double dt, int N, std::string methode) {
    fichier << "# précision  : " << dt << " s\n";
    fichier << "# durée      : " << tf << " s\n";
    fichier << "# particules : " << N << std::endl;
    fichier << "# méthode    : " << methode << std::endl;
}

// Pour écrire les positions
void ecrirePos(const std::string& nomFichier, const std::vector<std::vector<double>>& Q, int nbPoint, double dt, int N, std::string methode) {
    int nbLigne = Q.size();
    
    // Ouvrir le fichier en écriture
    std::ofstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << " ! Impossible d'ouvrir le fichier " << nomFichier << std::endl;
        return;
    }

    // Ecriture de l'entete
    ecrireEntete(fichier, dt, N, methode);

    // Écrire la ligne de format
    fichier << std::setiosflags(std::ios::scientific) << std::setprecision(6);

    // Écrire les lignes de la matrice Q
    int sautLigne = nbLigne / nbPoint;
    for (int i = 0; i < nbLigne; i += sautLigne) {
        for (size_t j = 0; j < Q[i].size(); j++) {
            fichier << Q[i][j];
            if (j < Q[i].size() - 1) {
                fichier << ",";
            }
        }
        fichier << "\n"; // Retour à la ligne
    }

    std::cout << " > Le fichier " << nomFichier << " a bien été créé" << std::endl;

    // Fermer le fichier
    fichier.close();
}

// Pour écrire l'énergie ou nu
void ecrireEnNu(const std::string& nomFichier, const std::vector<double>& x, const std::vector<double>& y, double dt, int N, std::string methode) {
    // Vérifier que les vecteurs x et y ont la même taille
    if (x.size() != y.size()) {
        std::cerr << " ! Les vecteur x et y n'ont pas la même taille" << std::endl;
        return;
    }

    // Ouvrir le fichier en mode écriture
    std::ofstream fichier(nomFichier);

    // Vérifier si le fichier est ouvert avec succès
    if (!fichier.is_open()) {
        std::cerr << " ! Impossible d'ouvrir le fichier " << nomFichier << std::endl;
        return;
    }

    // Écrire les données
    ecrireEntete(fichier,dt,N,methode);
    for (size_t i = 0; i < x.size(); ++i) {
        fichier << x[i] << "," << y[i] << std::endl;
    }

    std::cout << " > Le fichier " << nomFichier << " a bien été créé" << std::endl;

    // Fermer le fichier
    fichier.close();
}