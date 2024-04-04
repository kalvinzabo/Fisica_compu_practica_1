#RUN THIS IN SAME FOLDER AS 'pos_global.npy' which holds a 3-dim array (iterations, bodies, axes)

import matplotlib.pyplot as plt
# import matplotlib.animation as anim
import numpy as np

def myplot():
    pos = np.load('resultados/pos_global.npy')

    # vel = np.load('vel_global')
    # a = np.load('ac_global')

    fig, axs = plt.subplots(ncols=2, figsize=(10, 10))

    planet_names = ['Sol', 'Mercurio', 'Venus', 'Tierra', 'Marte']
    for i in range(0, 5):
        axs[0].plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names[i]}')
        
    planet_names_ext = ['Júpiter','Saturno', 'Urano', 'Neptuno']
    for i in range(5, 9):
        axs[1].plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names_ext[i-5]}')

    axs[0].set_title('Órbita de los planetas interiores respecto al Sol')
    axs[1].set_title('Órbita de los planetas exteriores respecto al Sol')

    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xlabel('Posición X (UA)')
        ax.set_ylabel('Posición Y (UA)')
        ax.legend()

    plt.show()
