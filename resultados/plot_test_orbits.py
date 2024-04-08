#RUN THIS IN SAME FOLDER AS 'pos_global.npy' which holds a 3-dim array (iterations, bodies, axes)

import matplotlib.pyplot as plt
# import matplotlib.animation as anim
import numpy as np

def external_plot(center = 0):
    pos = np.load('resultados/pos_global.npy')
    orbitplot(center, pos)

def orbitplot(center : int, pos):
    if center > 8 or center <= -1:
        print('error: el argumento debe estar entre 0 y 8')
        return
    
    if center != 0:
        pos -= pos[:, center, np.newaxis]
    
    fig, axs = plt.subplots(ncols=2, figsize=(20, 20))

    colors = ['gold', 'lightsteelblue', 'orange', 'olivedrab', 'indianred', 'peru', 'moccasin', 'lightcyan', 'royalblue']
    planet_names = ['Sol', 'Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter','Saturno', 'Urano', 'Neptuno']
    for i in range(8,-1,-1):
        if i == center:
            continue
        if i <= 4:
            axs[0].plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names[i]}', color=colors[i])
        else:
            axs[1].plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names[i]}', color=colors[i])

    axs[0].set_title(f'Órbita de los planetas interiores respecto a {planet_names[center]}')
    axs[1].set_title(f'Órbita de los planetas exteriores respecto a {planet_names[center]}')

    for ax in axs:
        ax.scatter(0, 0, color=colors[center], label=planet_names[center])
        ax.set_aspect('equal')
        ax.set_xlabel('Posición X (UA)')
        ax.set_ylabel('Posición Y (UA)')
        ax.legend(loc = 'upper right')

    plt.show()
    
if __name__ == '__main__':
    pos = np.load('pos_global.npy')
    orbitplot(3, pos)
