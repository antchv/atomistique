#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNIN
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <experimental/filesystem> // Header file for pre-standard implementation
using namespace std::experimental::filesystem::v1;
#include "util.h"


bool suppFichier(const std::string& nomFichier) {
    std::ifstream fichier(nomFichier);
    if (fichier.good()) {
        fichier.close();
        if (std::remove(nomFichier.c_str()) == 0) {
            std::cout << "Le fichier '" << nomFichier << "' a été détruit avec succès." << std::endl;
            return true;
        } else {
            std::cerr << "Erreur lors de la suppression du fichier '" << nomFichier << "'." << std::endl;
            return false;
        }
    } else {
        std::cout << "Le fichier '" << nomFichier << "' n'existe pas." << std::endl;
        return false;
    }
}

std::vector<double> extraireColonne(const std::string& nomFichier, int colonne) {
    std::vector<double> colonne_ieme;

    // Ouvrir le fichier en mode lecture
    std::ifstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier : " << nomFichier << std::endl;
        return colonne_ieme; // Retourne une colonne vide en cas d'erreur
    }

    std::string ligne;
    while (std::getline(fichier, ligne)) {
        std::istringstream stream_ligne(ligne);
        double valeur;
        int colonne_actuelle = 0;

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

int nombreDeLignesDansCSV(const std::string& nomFichier) {
    int nombreLignes = 0;

    // Ouvrir le fichier en mode lecture
    std::ifstream fichier(nomFichier);
    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier : " << nomFichier << std::endl;
        return nombreLignes; // Retourne 0 en cas d'erreur
    }

    std::string ligne;
    while (std::getline(fichier, ligne)) {
        nombreLignes++;
    }

    fichier.close();
    return nombreLignes;
}

std::vector<std::string> avoirFichiersDossier(const std::string& cheminDossier) {
    std::vector<std::string> nomFichiers;

    try {
        for (const auto& entree : std::experimental::filesystem::directory_iterator(cheminDossier)) {
            if (std::experimental::filesystem::is_regular_file(entree)) {
                nomFichiers.push_back(entree.path().filename().string());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }

    return nomFichiers;
}

void remplirCC(std::string& nom) {
    std::vector<std::string> fichiers;
    fichiers = avoirFichiersDossier("./condition_initiale");
    int taille = fichiers.size();
    int index;

    std::cout << "************************" << std::endl;

    for(int i = 0; i < taille; i++) {
        std::cout << i << " -> " << fichiers[i] << std::endl;
    }

    std::cout << "************************\nchoix : ";

    std::cin >> index;
    if(index >=0 && index < taille) {
        nom = fichiers[index];
    }
    else {
        std::cout << "Erreur dans le choix du fichier" << std::endl;
    }
}

void ecriturePos(const std::string& nomFichier,
                 const std::vector<double>& X,
                 const std::vector<double>& Y) {
    // Ouvrir un fichier en mode append
    std::ofstream fichier(nomFichier, std::ios::app);

    // Vérifier si l'ouverture a réussi
    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << nomFichier << std::endl;
    }
    else {
        // Écrire dans le fichier
        fichier << std::setiosflags(std::ios::scientific) << std::setprecision(6);
        
        for (size_t i = 0; i < X.size(); i++) {
            fichier << X[i];
            if (i < X.size() - 1) {
                fichier << ",";
            }
        }

        fichier << std::endl;

        for (size_t i = 0; i < Y.size(); i++) {
            fichier << Y[i];
            if (i < Y.size() - 1) {
                fichier << ",";
            }
        }
        
        fichier << std::endl;

        // Fermer le fichier
        fichier.close();
    }
}

void ecritureEn(const std::string& nomFichier, const double& temps, const double& energie) {
    // Ouvrir un fichier en mode append
    std::ofstream fichier(nomFichier, std::ios::app);

    // Vérifier si l'ouverture a réussi
    if (!fichier.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << nomFichier << std::endl;
    }
    else {
        // Écrire dans le fichier
        fichier << std::setiosflags(std::ios::scientific) << std::setprecision(6);
        fichier << temps << "," << energie << std::endl;

        // Fermer le fichier
        fichier.close();
    }
}