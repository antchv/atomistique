#include <iostream>
#include <vector>
#include <cmath>
#include "fcts.h"
#include "methode.h"
#include "util.h"


// Nom des fichiers
const std::string nomSim = choixSimu("./condition_initiale");        // Nom de la simulation
const std::string nomCI = "./condition_initiale/" + nomSim + ".csv";     // Nom du fichier à étudier
const std::string nomPosX = "./resultats/" + nomSim + "_position_x.csv"; // Nom du fichier des positions suivant X
const std::string nomPosY = "./resultats/" + nomSim + "_position_y.csv"; // Nom du fichier des positions suivant Y
const std::string nomEn = "./resultats/" + nomSim + "_energie.csv";      // Nom du fichier des énergies
const std::string nomNu = "./resultats/" + nomSim + "_nu.csv";           // Nom du fichier pour calculer nu

// Paramètres de la simulation
const double tf = 30.;                    // Temps final
const int nbFic = 1e3;                    // # de pts des fichiers (Pos,En), doir être supérieur à nbSim
const int nbNu = 5e0;                     // # de pts du fichier (Nu)
const int N = compteLigneFichier(nomCI);  // # de particule(s)
const int nbSimMin = 1e3;
const int nbSimMax = 1e4;

int main() {
	// Allocation mémoire des matrices position, vitesse et force (2D)
	// colonne = particule, ligne = instant donné
	std::vector<std::vector<double>> X(nbSimMax, std::vector<double>(N, 0));
	std::vector<std::vector<double>> Xp(nbSimMax, std::vector<double>(N, 0));
	std::vector<std::vector<double>> Xpp(nbSimMax, std::vector<double>(N, 0));
	std::vector<std::vector<double>> Y(nbSimMax, std::vector<double>(N, 0));
	std::vector<std::vector<double>> Yp(nbSimMax, std::vector<double>(N, 0));
	std::vector<std::vector<double>> Ypp(nbSimMax, std::vector<double>(N, 0));

	// Implémentation des Conditions Initiales
	X[0] = extraireColonne(nomCI, 0);
	Xp[0] = extraireColonne(nomCI, 2);
	Y[0] = extraireColonne(nomCI, 1);
	Yp[0] = extraireColonne(nomCI, 3);

	// Allocation mémoire des quantités qui permettent de calculer nu et l'énergie
	std::vector<double> deltaE(nbNu, -1. * EmLJ(X[0], Y[0], Xp[0], Yp[0])), dt(nbNu), nbSim(nbNu);
	std::vector<double> Em(nbFic), t(nbFic);

	// Boucle pour calculer nbSim
	for (int i = 0; i < nbNu; i++) {
		nbSim[i] = nbSimMin + i * (nbSimMax - nbSimMin) / (nbNu - 1);
	}

	// Boucle principale -> calcul des positions, vitesses et forces
	// for(int k = 0; k < nbNu; k++) {
	// 	dt[k] = tf / nbSim[k];

	// 	for(int i = 0; i < nbSim[k]-1; i++) {
	// 		for(int j = 0; j < N; j++) {
	// 			// Actualisation des forces
	// 			Xpp[i][j] = FLJX(j,X[i],Y[i]);
	// 			Ypp[i][j] = FLJY(j,X[i],Y[i]);

	// 			// Actualisation des positions
	// 			X[i+1][j] = eulerQ(X[i][j],Xp[i][j],dt[k]);
	// 			Y[i+1][j] = eulerQ(Y[i][j],Yp[i][j],dt[k]);

	// 			// Actualisation des vitesses
	// 			Xp[i+1][j] = eulerQp(Xp[i][j],Xpp[i][j],dt[k]);
	// 			Yp[i+1][j] = eulerQp(Yp[i][j],Ypp[i][j],dt[k]);
	// 		}
	// 	}
	// 	deltaE[k] += EmLJ(X[nbSim[k]-1],Y[nbSim[k]-1],Xp[nbSim[k]-1],Yp[nbSim[k]-1]);
	// }

	// Boucle principale -> calcul des positions, vitesses et forces
	for (int k = 0; k < nbNu; k++) {
		dt[k] = tf / nbSim[k];

		for (int i = 0; i < nbSim[k] - 1; i++) {
			for (int j = 0; j < N; j++) {
				// Actualisation des forces
				Xpp[i][j] = FLJX(j, X[i], Y[i]);
				Ypp[i][j] = FLJY(j, X[i], Y[i]);

				// Actualisation des positions
				X[i + 1][j] = velocityVerletQ(X[i][j], Xp[i][j], Xpp[i][j], dt[k]);
				Y[i + 1][j] = velocityVerletQ(Y[i][j], Yp[i][j], Ypp[i][j], dt[k]);
			}
			for (int j = 0; j < N; j++) {
				// Actualisation des forces
				Xpp[i + 1][j] = FLJX(j, X[i + 1], Y[i + 1]);
				Ypp[i + 1][j] = FLJY(j, X[i + 1], Y[i + 1]);

				// Actualisation des vitesses
				Xp[i + 1][j] = velocityVerletQp(Xp[i][j], Xpp[i][j], Xpp[i + 1][j], dt[k]);
				Yp[i + 1][j] = velocityVerletQp(Yp[i][j], Ypp[i][j], Ypp[i + 1][j], dt[k]);
			}
		}
		deltaE[k] += EmLJ(X[nbSim[k] - 1], Y[nbSim[k] - 1], Xp[nbSim[k] - 1], Yp[nbSim[k] - 1]);
	}

	// Boucle pour créer les données sur l'énergie
	int valeurPasEnergie = nbSim[nbNu - 1] / nbFic;
	for (int i = 0; i < nbFic; i++) {
		t[i] = i * valeurPasEnergie * dt[nbNu - 1];
		Em[i] = EmLJ(X[i * valeurPasEnergie], Y[i * valeurPasEnergie],
			Xp[i * valeurPasEnergie], Yp[i * valeurPasEnergie]);
	}

	// Sauvegarder les resultats
	ecrirePos(nomPosX, X, nbFic);
	ecrirePos(nomPosY, Y, nbFic);
	ecrireEnNu(nomEn, t, Em);
	ecrireEnNu(nomNu, dt, deltaE);

	return 0;
}