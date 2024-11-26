import matplotlib.pyplot as plt
import numpy as np

# Charger les données depuis le fichier
data = np.genfromtxt('pressure_output.txt', delimiter=',', skip_header=1)

# Extraire les colonnes du fichier de données

volume = data[:, 0]
pressure = data[:, 1]

# Tracer le graphique P(V)
plt.plot(volume, pressure, marker='o', linestyle='-', color='b')
plt.title('Pression en fonction du Volume (P(V))')
plt.xlabel('Volume')
plt.ylabel('Pression')
plt.grid(True)
plt.show()
