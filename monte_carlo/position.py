import pygame
import pandas as pd
import sys


# Chargement des données de position X et Y depuis les fichiers CSV
nom_fichier_x = "./resultats/cristal100_position_x.csv"
nom_fichier_y = "./resultats/cristal100_position_y.csv"

# On charge les valeurs avec pandas
dfx = pd.read_csv(nom_fichier_x, header=None, comment='#')
dfy = pd.read_csv(nom_fichier_y, header=None, comment='#')

# on stocke les données dans des variables
x_data = dfx.values.tolist()
y_data = dfy.values.tolist()
x_min = min(x_data[0])
x_max = max(x_data[0])
y_min = min(y_data[0])
y_max = max(y_data[0])



frame_x = -10.
frame_y = -10.
frame_width = 30
frame_height = 30

# Initialisation de pygame
pygame.init()

# Paramètres de la fenêtre
width, height = 1000, 600
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Animation de particules")

# Couleurs
white = (255, 255, 255)
black = (0, 0, 0)
taillePart = 10

# Durée de la vidéo en secondes
desired_duration = 30.

# Calcul du nombre total de frames et de la durée par frame
num_frames = min(len(x_data), len(y_data))
frame_duration = desired_duration / num_frames

# Initialisation de pygame
pygame.init()

# Fonction pour convertir les coordonnées du monde en coordonnées d'écran
def world_to_screen(x, y):
    screen_x = float((x - frame_x) / frame_width * width)
    screen_y = float(height - (y - frame_y) / frame_height * height)
    return screen_x, screen_y

# Fonction pour afficher les particules à une certaine frame
def draw_particles(frame_num):
    screen.fill(white)
    for i in range(len(x_data[frame_num])):
        x = float(x_data[frame_num][i])
        y = float(y_data[frame_num][i])
        screen_x, screen_y = world_to_screen(x, y)
        pygame.draw.circle(screen, black, (screen_x, screen_y), taillePart)



# Boucle principale
running = True
frame_num = 0
clock = pygame.time.Clock()

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    draw_particles(frame_num)

    pygame.display.flip()

    frame_num += 1
    if frame_num >= num_frames:
        running = False

    clock.tick(1 / frame_duration)

# Fermeture de pygame
pygame.quit()