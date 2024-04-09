#RUN THIS IN SAME FOLDER AS 'pos_global.npy' which holds a 3-dim array (iterations, bodies, axes)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridsp
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch

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

    # Create an inset axes for the energy and the angular momentum plot
    ax3_inset = inset_axes(axs[2], width=1, height=0.75, loc='upper left', bbox_to_anchor=(1.01, 2), bbox_transform=axs[2].transAxes, borderpad=0)
    ax4_inset = inset_axes(axs[3], width=1, height=0.75, loc='upper left', bbox_to_anchor=(1.01, -0.08), bbox_transform=axs[3].transAxes, borderpad=0)

    for i in range(8,-1,-1):
        
        ax3.plot(t, e_ind[:, i], label=f'{planet_names[i]}', color=colors[i])
        ax4.plot(t, l_modulos[:, i], label=f'{planet_names[i]}', color=colors[i])
        ax3_inset.plot(t, e_ind[:, i], color=colors[i])
        ax4_inset.plot(t, l_modulos[:, i], label=f'{planet_names[i]}', color=colors[i])

        
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

    ax3_inset.set_xlim(1000, 1100)  # specify the limits for x-axis
    ax3_inset.set_ylim(-0.0000025, 0.0000005)  # specify the limits for y-axis

    ax4_inset.set_xlim(1000, 1100)  # specify the limits for x-axis
    ax4_inset.set_ylim(-0.0000015, 0.000005)  # specify the limits for y-axis
    
    # Draw zoom lines for the energy plot
    xyA = (1000, -0.0000025)  # coordinates of the lower left corner of the zoomed area in the main plot
    xyB = (0, 1)  # coordinates of the lower left corner of the inset plot
    con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA="data", coordsB="axes fraction", axesA=axs[2], axesB=ax3_inset, color="gainsboro", linewidth =2)
    axs[2].add_artist(con)

    xyA = (1100, 0.0000005)  # coordinates of the upper right corner of the zoomed area in the main plot
    xyB = (1, 0)  # coordinates of the upper right corner of the inset plot
    con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA="data", coordsB="axes fraction", axesA=axs[2], axesB=ax3_inset, color="gainsboro", linewidth =2)
    axs[2].add_artist(con)

    # Draw zoom lines for the angular momentum plot
    xyA = (1000, -0.0000015)  # coordinates of the lower left corner of the zoomed area in the main plot
    xyB = (0, 0)  # coordinates of the lower left corner of the inset plot
    con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA="data", coordsB="axes fraction", axesA=axs[3], axesB=ax4_inset, color="gainsboro", linewidth =2)
    axs[3].add_artist(con)

    xyA = (1100, 0.000005)  # coordinates of the upper right corner of the zoomed area in the main plot
    xyB = (1, 1)  # coordinates of the upper right corner of the inset plot
    con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA="data", coordsB="axes fraction", axesA=axs[3], axesB=ax4_inset, color="gainsboro", linewidth =2)
    axs[3].add_artist(con)

    plt.show()
    
if __name__ == '__main__':
    pos = np.load('pos_global.npy')
    energias_individuales = np.load('energias_cuerpos.npy')
    energia_total = np.load('E_global.npy')
    L_individuales = np.load('L_cuerpos.npy')
    L_total = np.load('L_global.npy')
    t = np.load('t_global.npy')
    orbitplot(0, pos, energia_total, energias_individuales, L_total, L_individuales, t)
