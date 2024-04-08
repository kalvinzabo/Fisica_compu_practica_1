#RUN THIS IN SAME FOLDER AS 'pos_global.npy' which holds a 3-dim array (iterations, bodies, axes)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridsp
import numpy as np

plt.rcParams['text.usetex'] = True

def external_plot(center = 0):
    pos = np.load('resultados/pos_global.npy')
    energias_individuales = np.load('resultados/energias_cuerpos.npy')
    energia_total = np.load('resultados/E_global.npy')
    L_individuales = np.load('resultados/L_cuerpos.npy')
    L_total = np.load('resultados/L_global.npy')
    t = np.load('resultados/t_global.npy')
    orbitplot(center, pos, energia_total, energias_individuales, L_total, L_individuales, t)

def orbitplot(center : int, pos, e_total, e_ind, l_total, l_ind, t):
    if center > 8 or center <= -1:
        print('error: el argumento debe estar entre 0 y 8')
        return
    
    if center != 0:
        pos -= pos[:, center, np.newaxis]
    
    fig = plt.figure(figsize=(20, 20))

    gs = gridsp.GridSpec(3, 2, height_ratios=[3, 1, 1], width_ratios=[1, 1])

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[1, :])
    ax4 = plt.subplot(gs[2, :])
    axs = [ax1,ax2,ax3,ax4]

    plt.subplots_adjust(hspace = 0.35)

    # Now ax1 and ax2 are square subplots in the first row,
    # and ax3 and ax4 are rectangular subplots in the second and third rows.
    # fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15, 15))

    colors = ['gold', 'lightsteelblue', 'orange', 'olivedrab', 'indianred', 'peru', 'moccasin', 'lightcyan', 'royalblue']
    planet_names = ['Sol', 'Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter','Saturno', 'Urano', 'Neptuno']
    l_modulos = np.linalg.norm(l_ind, axis=2)
    l_modulos_total = np.linalg.norm(l_total, axis=1)

    for i in range(8,-1,-1):
        
        ax3.plot(t, e_ind[:, i], label=f'{planet_names[i]}', color=colors[i])
        ax4.plot(t, l_modulos[:, i], label=f'{planet_names[i]}', color=colors[i])
        
        if i == center:
            continue
        if i <= 4:
            ax1.plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names[i]}', color=colors[i])
        else:
            ax2.plot(pos[:, i, 0], pos[:, i, 1], label=f'{planet_names[i]}', color=colors[i])
    
    ax3.plot(t, e_total, label='Total del sistema', color='purple')
    ax4.plot(t, l_modulos_total, label='Total del sistema', color='purple')
        

    ax1.set_title(f'Órbita de los planetas interiores respecto a {planet_names[center]}')
    ax2.set_title(f'Órbita de los planetas exteriores respecto a {planet_names[center]}')
    ax3.set_title(f'Conservación de la energía total')
    ax4.set_title(f'Conservación del momento angular total')

    for i in range(2):
        axs[i].scatter(0, 0, color=colors[center], label=planet_names[center])
        axs[i].set_xlabel('Posición X (UA)')
        axs[i].set_ylabel('Posición Y (UA)')
        axs[i].set_aspect('equal')
        axs[i].legend(loc = 'upper left', bbox_to_anchor=(1,1))

    # axs[2].set_xlabel(r'Tiempo (1 $\approx$ 58.1 dias)')
    axs[2].set_ylabel(r'Energía \newline $(M_{Sol}\cdot UA^2/(58,1 \: dias)^2)$')
    axs[2].legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    
    axs[2].set_xlim(min(t), max(t))
    axs[3].set_xlim(min(t), max(t))

    axs[3].set_xlabel(r'Tiempo (1 $\approx$ 58.1 dias)')
    axs[3].set_ylabel(r'Momento angular \newline $(M_{Sol}\cdot UA^2/(58,1 dias))$')
    

    plt.show()
    
if __name__ == '__main__':
    pos = np.load('pos_global.npy')
    energias_individuales = np.load('energias_cuerpos.npy')
    energia_total = np.load('E_global.npy')
    L_individuales = np.load('L_cuerpos.npy')
    L_total = np.load('L_global.npy')
    t = np.load('t_global.npy')
    orbitplot(0, pos, energia_total, energias_individuales, L_total, L_individuales, t)
