#Diseñar y escribir el programa que simula el comportamiento del Sistema Solar. Usando datos reales para la
#excentricidad de las órbitas, las distancias y las masas de los diferentes planetas, podemos dar una condición
#inicial realista para el Sistema Solar y así estudiar su dinámica y estabilidad frente a perturbaciones. Incluir
#todas las interacciones gravitatorias entre los planetas.
#Comprobar la validez de la simulación midiendo los períodos de rotación de los diferentes planetas y
#comparándolos con los reales. Mostrar la órbita relativa respecto a la Tierra del resto de planetas, en función del
#tiempo.

#Librerias
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from resultados import plot_test_orbits
import timeit

# plt.rcParams['text.usetex'] = True

#Condiciones iniciales
t0 = 0

    #Constantes físicas
G = 6.67e-11
c = 149.598e9

    #Masas (kg)
M_sol = 1.9885e30
m_mercurio = 3.301e23
m_venus = 4.8673e24
m_tierra = 5.9722e24
m_marte = 6.4169e23
m_jupiter = 1.89813e27
m_saturno = 5.6832e26
m_urano = 8.6811e25
m_neptuno = 1.02409e26
masas = np.array([M_sol, m_mercurio, m_venus, m_tierra, m_marte, m_jupiter, m_saturno, m_urano, m_neptuno])
#Reescalo la masa con reespecto a la masa del sol siendo las nuevas masas m'=m/M_sol
masas_r = masas/M_sol

    #Posiciones iniciales (m)
perihelio_sol = np.array([0.,0.,0.])
perihelio_mercurio = np.array([4.6e10, 0., 0.])
perihelio_venus = np.array([1.0748e11, 0., 0.])
perihelio_tierra = np.array([1.47095e11, 0., 0.])
perihelio_marte = np.array([2.0665e11, 0., 0.])
perihelio_jupiter = np.array([7.40595e11, 0., 0.])
perihelio_saturno = np.array([1.357554e12, 0., 0.])
perihelio_urano = np.array([2.732696e12, 0., 0.])
perihelio_neptuno = np.array([4.47105e12, 0., 0.])
perihelios = np.array([perihelio_sol, perihelio_mercurio, perihelio_venus, perihelio_tierra, perihelio_marte, perihelio_jupiter, perihelio_saturno, perihelio_urano, perihelio_neptuno])
#Reescalo las distancias con respecto la distancia sol-tierra siendo las nuevas distancias r'=r/r_sol_tierra
perihelios_r = perihelios/c

    #Velocidades (m/s)
v_max_sol = np.array([0.,0.,0.])
v_max_mercurio = np.array([0.,58.97e3, 0.])
v_max_venus = np.array([0.,35.26e3, 0.])
v_max_tierra = np.array([0.,30.29e3, 0.])
v_max_marte = np.array([0.,26.5e3, 0.])
v_max_jupiter = np.array([0.,13.72e3, 0.])
v_max_saturno = np.array([0.,10.14e3, 0.])
v_max_urano = np.array([0.,7.13e3, 0.])
v_max_neptuno = np.array([0.,5.47e3, 0.])
velocidades_max = np.array([v_max_sol, v_max_mercurio, v_max_venus, v_max_tierra, v_max_marte, v_max_jupiter, v_max_saturno, v_max_urano, v_max_neptuno])
#Reescalo las velocidades teniendo en cuenta que las distancias las reescalamos como r'=r/r_sol_tierra y los tiempos t'=t*(G*M_sol/d_sol_tierra**3)**0.5
velocidades_max_r = (velocidades_max*(c**(1.5)))/(c*np.sqrt(G*M_sol))

    #Excentricidades
e_sol = 0.
e_mercurio = 0.206
e_venus = 0.007
e_tierra = 0.017
e_marte = 0.09341233
e_jupiter = 0.04839266
e_saturno = 0.05415060
e_urano = 0.047
e_neptuno = 0.00858587
epsilons_teo = np.array([e_sol, e_mercurio, e_venus, e_tierra, e_marte, e_jupiter, e_saturno, e_urano, e_neptuno])

    #Periodos (días)
T_sol = 0.
T_mercurio = 88.0
T_venus = 224.7
T_tierra = 365.2
T_marte = 687.0
T_jupiter = 4331
T_saturno = 10747
T_urano = 30589
T_neptuno = 59800
T_teo = np.array([T_sol, T_mercurio, T_venus, T_tierra, T_marte, T_jupiter, T_saturno, T_urano, T_neptuno])

#Variables necesarias para el algoritmo y el ploteo
h = 0.01
t_limit = 1100
iterations=int(t_limit//h)

cuerpos = len(masas)
coordenadas = len(perihelios[0])

time = np.linspace(t0,t_limit,iterations)
pos = np.zeros((iterations, cuerpos, coordenadas))
vel = np.zeros((iterations, cuerpos, coordenadas))
a = np.zeros((iterations, cuerpos, coordenadas))
energia_total = np.zeros((iterations))
energias_cuerpos = np.zeros((iterations,cuerpos))
momento_angular_total = np.zeros((iterations, coordenadas))
momento_angular_cuerpos = np.zeros((iterations,cuerpos,coordenadas))
periodos = np.zeros(cuerpos)
is_r_disminuyendo = [False]*9    #list of flags to check if each planet has started approaching its starting point (used to calculate the orbit period)
r_anteriores = [0.]*9     #this one holds the last distance so we can check if it is increasing or decreasing
maximas_distancias_r = np.zeros(cuerpos)
epsilons = np.zeros(cuerpos)
epsilons_energia = np.zeros(cuerpos)
epsilons_energia_0 = np.zeros(cuerpos)
epsilons_fit = np.zeros(cuerpos)
var_epsilons_fit = np.zeros(cuerpos)

#Bloque de funciones.
#Función de cálculo de aceleraciones. Las calculamos de manera más sencilla debido al reescalamiento que hemos hecho a las masas, distancias y tiempos
def calculo_aceleracion (mass=masas_r, dis=perihelios_r):
    aceleraciones = np.zeros((9,3))
    len_mass = len(mass)
    for i in range(1, len_mass):
        for j in range(len_mass):
            if j == i:
                continue
            aceleraciones[i] += -1*(mass[j]*(dis[i]-dis[j]))/np.linalg.norm(dis[i]-dis[j])**3
    return aceleraciones

#Función de cálculo de energía. Calculamos la E cinética y la E potencial y las sumamos
def calculo_energia_total(iteracion, pos, vel, masas=masas_r, energias_cuerpos=energias_cuerpos):
    energia_total = 0
    for i in range(len(masas)):
        energia_cinetica = 0.5 * masas[i] * np.linalg.norm(vel[i])**2
        energia_potencial = 0
        for j in range(len(masas)):
            if i == j:
                continue
            r = np.linalg.norm(pos[i] - pos[j])
            energia_potencial -= masas[i] * masas[j] / r
        energia_cuerpo = energia_cinetica + energia_potencial
        energia_total += energia_cuerpo
        energias_cuerpos[iteracion, i] = energia_cuerpo
    return energia_total

#Función de cálculo del momento angular total. Calculamos L = m (r x v)
def calculo_momento_angular_total(iteracion, pos, vel, masas=masas_r, mom_ang_cuerpos=momento_angular_cuerpos):
    momento_angular_total = np.zeros(3)
    for i in range(len(masas)):
        mom_ang_cuerpo = masas[i] * np.cross(pos[i], vel[i])
        momento_angular_total += mom_ang_cuerpo
        mom_ang_cuerpos[iteracion,i] = mom_ang_cuerpo
    return momento_angular_total

#Función para calcular el afelio
def actualizar_maxima_distancia(num_iteracion : int, planeta:int):
    global maximas_distancias_r
    r = np.linalg.norm(pos[num_iteracion, planeta])
    if r > maximas_distancias_r[planeta]:
        maximas_distancias_r[planeta] = r

#Función de cálculo del periodo. Asigna un periodo a un cuerpo si su distancia con su posición inicial estaba disminuyendo y empieza de nuevo a aumentar.
def calculo_periodo_orbital(num_iteracion, periodos=periodos, pos_global=pos, step=h, posiciones_iniciales=perihelios_r):
    for j in range(1,len(posiciones_iniciales)):
        if periodos[j] != 0:
            continue

        actualizar_maxima_distancia(num_iteracion, j)

        r = np.linalg.norm(pos_global[num_iteracion][j] - posiciones_iniciales[j])
        
        if r_anteriores[j] > r:     #distance is decreasing
            is_r_disminuyendo[j] = True
        
        elif r_anteriores[j] == r:  #distance did not change (error)
            raise Exception(f"The distance to the starting point of planet {j} is not changing. The planet is not moving or the distance is being updated incorrectly.")
        
        elif is_r_disminuyendo[j]:   #distance is increasing (is not decreasing nor equal) and was decreasing until now.
            print(f'El cuerpo {j} ha completado una órbita en {num_iteracion} iteraciones (aprox. {num_iteracion*step*58.1} días)')
            periodos[j] = num_iteracion * step
            continue
        
        r_anteriores[j] = r     #after entering the first case (or none in the first half of the orbit), update last distance

#Función de la elipse para hacer curve_fit
def ellipse(x,a,b):
    return b * np.sqrt(1-((x)/a)**2)

#Introduzco las posiciones, velocidades, aceleraciones, la energía y el momento angular iniciales en los arrays que dependen del tiempo
pos[0] = perihelios_r
vel[0] = velocidades_max_r
a[0] = calculo_aceleracion()
energia_total[0] = calculo_energia_total(0, pos[0], vel[0])
momento_angular_total[0] = calculo_momento_angular_total(0, pos[0], vel[0])

start_time = timeit.default_timer()

#Algoritmo de Verlet
for t in range (iterations-1):
    w = vel[t] + h/2*a[t]
    pos[t+1] = pos[t] + h*w
    a[t+1] = calculo_aceleracion(dis=pos[t+1])
    vel[t+1] = w + h/2*a[t+1]

    #Calculo la energía y el momento angular del sistema en la nueva posición.
    energia_total[t+1] = calculo_energia_total(t+1, pos[t+1], vel[t+1])
    momento_angular_total[t+1] = calculo_momento_angular_total(t+1, pos[t+1], vel[t+1])

    #Compruebo si algún planeta ha completado su órbita.
    calculo_periodo_orbital(t+1)


    if t % 10000 == 0:
        print(f'Iteración número {t+1} de {iterations-1}')

end_time = timeit.default_timer()
execution_time = end_time - start_time

print(f'Finished simulating in {execution_time} time')

#Reescalamos los periodos, maximas distancias, energias y momentos angulares calculados a los reales
periodos_r = periodos*(c**3/(G*M_sol))**0.5*1.15741e-5
maximas_distancias = c*maximas_distancias_r

energias_cuerpos_medias = np.mean(energias_cuerpos, axis=0)
# energias_cuerpos_stderr = np.std(energias_cuerpos, axis=0)

momentos_angulares_medios = np.mean(momento_angular_cuerpos, axis=0)    #el truco es pensar que el eje que pasas (aqui iteraciones) es el que eliminas. ademas 0 es el eje de las filas (cada fila son 9 cuerpos) y te quedas con una unica fila (de 9 cuerpos)
modulos_momentos_angulares = np.linalg.norm(momentos_angulares_medios, axis=1)
modulos_momentos_angulares_0 = np.linalg.norm(momento_angular_cuerpos[0], axis=1)
# momentos_angulares_stderr = np.std(momento_angular_cuerpos, axis=0)

for i in range(1,cuerpos):

#Calculamos las excentricidades predichas según la forma de nuestra órbita.
    epsilons[i] = (maximas_distancias_r[i] - perihelios_r[i][0])/(maximas_distancias_r[i] + perihelios_r[i][0])

#Calculamos las excentricidades en relación con la energía
    energia_total_cuerpo = energias_cuerpos_medias[i] * M_sol * ((c**(1.5))/(c*np.sqrt(G*M_sol)))**(-2)     #esto está desnormalizando
    energia_0_cuerpo = energias_cuerpos[0,i] * M_sol * ((c**(1.5))/(c*np.sqrt(G*M_sol)))**(-2)     #esto está desnormalizando
#    print(f'energia sin normalizar era: {energias_cuerpos_medias[i]}, desnormalizada es: {energia_total_cuerpo}')

    momento_angular_total_cuerpo = modulos_momentos_angulares[i] * M_sol * ((c**(1.5))/(c*np.sqrt(G*M_sol)))**(-1) * c
    momento_angular_0_cuerpo = modulos_momentos_angulares_0[i] * M_sol * ((c**(1.5))/(c*np.sqrt(G*M_sol)))**(-1) * c
#    print(f'momento normalizado era: {momentos_angulares_medios[i]}, modulo {modulos_momentos_angulares[i]}, desnormalizado: {momento_angular_total_cuerpo}')
    
    factor = (2*energia_total_cuerpo*momento_angular_total_cuerpo**2)/(G**2*(M_sol)**2*(masas[i])**3)
    factor_0 = (2*energia_0_cuerpo*momento_angular_0_cuerpo**2)/(G**2*(M_sol)**2*(masas[i])**3)
    epsilons_energia[i] = np.sqrt(1+factor)
    epsilons_energia_0[i] = np.sqrt(1+factor_0)

#Calculamos las excentricidades haciendo curve_fit
    x = pos[:,i,0]
    y = pos[:,i,1]

    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)

    a_init = np.max(np.sqrt(x_centered**2 + y_centered**2))
    b_init = np.max(np.abs(y_centered))

    params, params_covariance = curve_fit(ellipse, x_centered, y_centered, p0=[a_init,b_init])
    a_opt, b_opt = params
    print(a_opt,b_opt)
    std_devs = np.sqrt(np.diag(params_covariance))
    std_dev_a, std_dev_b = std_devs

    epsilons_fit[i] = np.sqrt(1 - (b_opt/a_opt)**2)
    var_epsilons_fit[i] = np.sqrt((b_opt**2 / (a_opt**3 * epsilons_fit[i]) * std_dev_a)**2 + (-b_opt / (a_opt * epsilons_fit[i]) * std_dev_b)**2)

#Comparamos los periodos y las excentricidades calculadas
    print('\033[1m' + f'Para el cuerpo {i}:' + '\033[0m')
    print(f'El periodo medido es {periodos_r[i]} días. Podemos compararlo con su periodo real que es {T_teo[i]} días.')
    print(f'Su error relativo es {abs(periodos_r[i]-T_teo[i])/T_teo[i]*100}%.')
    print(f'La excentricidad real es {epsilons_teo[i]}.')
    print(f'La excentricidad medida según la órbita predicha es {epsilons[i]}. Su error relativo respecto a la teórica es {abs(epsilons[i]-epsilons_teo[i])/epsilons_teo[i]*100}%.')
    print(f'La excentricidad obtenida por el ajuste de mínimos cuadrados es {epsilons_fit[i]} \pm {var_epsilons_fit[i]}. Su error relativo respecto a la teórica es {abs(epsilons_fit[i]-epsilons_teo[i])/epsilons_teo[i]*100}%.')
    print(f'La excentricidad calculada en función de la E_0 total y el L_0 total es {epsilons_energia_0[i]}. Su error relativo respecto a la teórica es {abs(epsilons_energia_0[i]-epsilons_teo[i])/epsilons_teo[i]*100}%.')
    print(f'Podemos comparar las excentricidades medidas e inicialmente calculadas en la simulación, atendiendo al error relativo que es {abs(epsilons[i]-epsilons_energia_0[i])/epsilons_energia_0[i]*100}%.')
    print(f'El factor es: {factor}')
    if factor < -1:
        print(f'No se puede calcular para este cuerpo la excentricidad en función de la E total y el L total obtenidos en esta simulación.')
        print()
        continue 
    print(f'La excentricidad calculada en función de la E total y el L total predichos es {epsilons_energia[i]}. Su error relativo respecto a la teórica es {abs(epsilons_energia[i]-epsilons_teo[i])/epsilons_teo[i]*100}%.')
    print(f'Podemos comparar las excentricidades medidas y calculadas en la simulación, atendiendo al error relativo que es {abs(epsilons[i]-epsilons_energia[i])/epsilons[i]*100}%.')
    print()

np.save('resultados/pos_global.npy', pos)
np.save('resultados/vel_global.npy', vel)
np.save('resultados/ac_global.npy', a)
np.save('resultados/E_global.npy', energia_total)
np.save('resultados/energias_cuerpos.npy', energias_cuerpos)
np.save('resultados/L_global.npy', momento_angular_total)
np.save('resultados/L_cuerpos.npy', momento_angular_cuerpos)
np.save('resultados/t_global.npy', time)
np.savetxt('resultados/maximas_distancias.txt', maximas_distancias)
#Save maximum distances of each planet

plot_test_orbits.external_plot()

# for i in range(len(masas)):

#     np.save(f'resultados/planeta_{i}_pos.npy', pos[:,i,:])
#     np.save(f'resultados/planeta_{i}_vel.npy', vel[:,i,:])
#     np.save(f'resultados/planeta_{i}_ac.npy', a[:,i,:])
    
#     np.savetxt(f'resultados/planeta_{i}_pos.txt', pos[:,i,:], fmt='%.4f')
#     np.savetxt(f'resultados/planeta_{i}_vel.txt', vel[:,i,:], fmt='%.4f')
#     np.savetxt(f'resultados/planeta_{i}_ac.txt', a[:,i,:], fmt='%.4f')
    
    # np.savetxt(f'resultados/planeta_{i}_pos.txt', pos[:,i,:])
    # np.savetxt(f'resultados/planeta_{i}_vel.txt', vel[:,i,:])
    # np.savetxt(f'resultados/planeta_{i}_ac.txt', a[:,i,:])

# Plot de posiciones, velocidades y aceleraciones (comentado porque está en un archivo aparte que trabaja con los resultados)
# fig, ax = plt.subplots(figsize=(10, 10))

# for i in range(1, 9):
#     ax.plot(pos[:, i, 0], pos[:, i, 1], label=f'Planeta {i}')

# ax.set_aspect('equal')
# ax.set_xlabel('Posición X (UA)')
# ax.set_ylabel('Posición Y (UA)')
# ax.set_title('Órbita de los cuerpos respecto al Sol')
# ax.legend()

# plt.show()
