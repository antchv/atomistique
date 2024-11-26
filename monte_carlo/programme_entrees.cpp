#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

void creation_fichier(int N, const vector<float>& X, const vector<float>& Y, const vector<float>& Z)
{
    ofstream fichier("parametres.txt");

    fichier << std::setiosflags(std::ios::scientific) << std::setprecision(2);

    for(int i = 0; i < N; i++)
    {
        fichier << X[i] <<" ";
        fichier << Y[i] <<" ";
        fichier << Z[i] <<" \n";
    }
}

int main()
{
    int N;
    int type_distrib;

    cout << "Veuillez entrer le nombre de particules : ";
    cin >> N;

    vector<float> X(N, 0);
    vector<float> Y(N, 0);
    vector<float> Z(N, 0);

    int rows, cols, depth;

    cout << "Entrez le nombre de lignes du quadrillage : ";
    cin >> rows;

    cout << "Entrez le nombre de colonnes du quadrillage : ";
    cin >> cols;

    cout << "Entrez le nombre de profondeur du quadrillage : ";
    cin >> depth;

    float spacing_x = pow(2., 1. / 6.);
    float spacing_y = pow(2., 1. / 6.);
    float spacing_z = pow(2., 1. / 6.);

    float grid_width = (cols - 1) * spacing_x;
    float grid_height = (rows - 1) * spacing_y;
    float grid_depth = (depth - 1) * spacing_z;

    int particle_count = 0;

    for (int k = 0; k < depth; k++) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                X[particle_count] = (j * spacing_x) - (grid_width / 2.0);
                Y[particle_count] = (i * spacing_y) - (grid_height / 2.0);
                Z[particle_count] = (k * spacing_z) - (grid_depth / 2.0);
                particle_count++;
            }
        }
    }

    creation_fichier(N, X, Y, Z);

    return 0;
}
