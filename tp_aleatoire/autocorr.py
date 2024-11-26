import matplotlib.pyplot as plt

# Lire les données depuis le premier fichier
with open('autocorrelation.dat', 'r') as file:
    lines = file.readlines()

# Extraire les décalages et les valeurs d'autocorrélation du premier fichier
lag_values = []
autocorrelation_values = []

for line in lines:
    lag, autocorrelation = map(float, line.split())
    lag_values.append(lag)
    autocorrelation_values.append(autocorrelation)

# Lire les données depuis le deuxième fichier
with open('autocorrelation2.dat', 'r') as file:
    lines2 = file.readlines()

# Extraire les décalages et les valeurs d'autocorrélation du deuxième fichier
autocorrelation_values2 = []

for line in lines2:
    lag, autocorrelation = map(float, line.split())
    autocorrelation_values2.append(autocorrelation)

# Créer un subplot avec deux graphiques côte à côte
fig, axs = plt.subplots(1, 2, figsize=(12, 4))

# Premier graphique
axs[0].scatter(lag_values, autocorrelation_values, s=1)
axs[0].set_title('rng classique')
axs[0].set_xlabel('Décalage')
axs[0].set_ylabel('Autocorrélation')
axs[0].set_ylim([-0.02, 0.02])
axs[0].grid(True)

# Deuxième graphique
axs[1].scatter(lag_values, autocorrelation_values2, s=1)
axs[1].set_title('mersene twister')
axs[1].set_xlabel('Décalage')
axs[1].set_ylabel('Autocorrélation')
axs[1].set_ylim([-0.02, 0.02])
axs[1].grid(True)

# Ajuster la disposition pour éviter les chevauchements
plt.tight_layout()

# Afficher le subplot
plt.show()
