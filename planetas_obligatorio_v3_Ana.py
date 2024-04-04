#Diseñar y escribir el programa que simula el comportamiento del Sistema Solar. Usando datos reales para la
#excentricidad de las órbitas, las distancias y las masas de los diferentes planetas, podemos dar una condición
#inicial realista para el Sistema Solar y así estudiar su dinámica y estabilidad frente a perturbaciones. Incluir
#todas las interacciones gravitatorias entre los planetas.
#Comprobar la validez de la simulación midiendo los períodos de rotación de los diferentes planetas y
#comparándolos con los reales. Mostrar la órbita relativa respecto a la Tierra del resto de planetas, en función del
#tiempo.

#Librerias
import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams['text.usetex'] = True

#Condiciones iniciales
t0 = 0

    #Constantes físicas
G = 6.67e-11

    #Masas (kg)
M_sol = 1.99e30
m_mercurio = 3.30e23
m_venus = 4.87e24
m_tierra = 5.97e24
m_marte = 6.42e23
m_jupiter = 1.898e27
m_saturno = 5.68e26
m_urano = 8.68e25
m_neptuno = 1.02e26
masas = np.array([M_sol, m_mercurio, m_venus, m_tierra, m_marte, m_jupiter, m_saturno, m_urano, m_neptuno])
#Reescalo la masa con reespecto a la masa del sol siendo las nuevas masas m'=m/M_sol
masas_r = masas/M_sol

    #Posiciones iniciales (m)
d_sol = np.array([0.,0.,0.])
d_mercurio = np.array([5.79e10, 0., 0.])
d_venus = np.array([1.082e11, 0., 0.])
d_tierra = np.array([1.496e11, 0., 0.])
d_marte = np.array([2.280e11, 0., 0.])
d_jupiter = np.array([7.785e11, 0., 0.])
d_saturno = np.array([1.432e12, 0., 0.])
d_urano = np.array([2.867e12, 0., 0.])
d_neptuno = np.array([4.515e12, 0., 0.])
distancias = np.array([d_sol, d_mercurio, d_venus, d_tierra, d_marte, d_jupiter, d_saturno, d_urano, d_neptuno])
#Reescalo las distancias con respecto la distancia sol-tierra siendo las nuevas distancias r'=r/r_sol_tierra
distancias_r = distancias/d_tierra[0]

    #Velocidades (m/s)
v_sol = np.array([0.,0.,0.])
v_mercurio = np.array([0.,47.4e3, 0.])
v_venus = np.array([0.,35e3, 0.])
v_tierra = np.array([0.,29.8e3, 0.])
v_marte = np.array([0.,24.1e3, 0.])
v_jupiter = np.array([0.,13.1e3, 0.])
v_saturno = np.array([0.,9.7e3, 0.])
v_urano = np.array([0.,6.8e3, 0.])
v_neptuno = np.array([0.,5.4e3, 0.])
velocidades = np.array([v_sol, v_mercurio, v_venus, v_tierra, v_marte, v_jupiter, v_saturno, v_urano, v_neptuno])
#Reescalo las velocidades teniendo en cuenta que las distancias las reescalamos como r'=r/r_sol_tierra y los tiempos t'=t*(G*M_sol/d_sol_tierra**3)**0.5
velocidades_r = (velocidades*(d_tierra[0]**(1.5)))/(d_tierra[0]*np.sqrt(G*M_sol))

    #Excentricidades
e_sol = 0.
e_mercurio = 0.206
e_venus = 0.007
e_tierra = 0.017
e_marte = 0.094
e_jupiter = 0.049
e_saturno = 0.052
e_urano = 0.047
e_neptuno = 0.010
epsilons = np.array([e_sol, e_mercurio, e_venus, e_tierra, e_marte, e_jupiter, e_saturno, e_urano, e_neptuno])

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
coordenadas = len(distancias[0])

t = np.linspace(t0,t_limit,iterations)
pos = np.zeros((iterations, cuerpos, coordenadas))
vel = np.zeros((iterations, cuerpos, coordenadas))
a = np.zeros((iterations, cuerpos, coordenadas))
energia_total = np.zeros((iterations))
momento_angular_total = np.zeros((iterations, coordenadas))
periodos = np.zeros(cuerpos)
is_r_disminuyendo = [False]*9    #list of flags to check if each planet has started approaching its starting point (used to calculate the orbit period)
r_anteriores = [0.]*9     #this one holds the last distance so we can check if it is increasing or decreasing

#Bloque de funciones.
#Función de cálculo de aceleraciones. Las calculamos de manera más sencilla debido al reescalamiento que hemos hecho a las masas, distancias y tiempos
def calculo_aceleracion (mass=masas_r, dis=distancias_r):
    aceleraciones = np.zeros((9,3))
    len_mass = len(mass)
    for i in range(1, len_mass):
        for j in range(len_mass):
            if j == i:
                continue
            aceleraciones[i] += -1*(mass[j]*(dis[i]-dis[j]))/np.linalg.norm(dis[i]-dis[j])**3
    return aceleraciones

#Función de cálculo de energía. Calculamos la E cinética y la E potencial y las sumamos
def calculo_energia_total(pos, vel, masas=masas_r):
    energia_total = 0
    for i in range(len(masas)):
        energia_cinetica = 0.5 * masas[i] * np.linalg.norm(vel[i])**2
        energia_potencial = 0
        for j in range(len(masas)):
            if i == j:
                continue
            r = np.linalg.norm(pos[i] - pos[j])
            energia_potencial -= masas[i] * masas[j] / r
        energia_total += energia_cinetica + energia_potencial
    return energia_total

#Función de cálculo del momento angular total. Calculamos L = m (r x v)
def calculo_momento_angular_total(pos, vel, masas=masas_r):
    momento_angular_total = np.zeros(3)
    for i in range(len(masas)):
        momento_angular_total += masas[i] * np.cross(pos[i], vel[i])
    return momento_angular_total

#Función de cálculo del periodo. Asigna un periodo a un cuerpo si su distancia con su posición inicial estaba disminuyendo y empieza de nuevo a aumentar.
def calculo_periodo_orbital(num_iteracion, periodos=periodos, pos_global=pos, step=h, posiciones_iniciales=distancias_r):
    for j in range(1,len(posiciones_iniciales)):
        if periodos[j] != 0:
            continue

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

#Introduzco las posiciones, velocidades, aceleraciones, la energía y el momento angular iniciales en los arrays que dependen del tiempo
pos[0] = distancias_r
vel[0] = velocidades_r
a[0] = calculo_aceleracion()
energia_total[0] = calculo_energia_total(pos[0], vel[0])
momento_angular_total[0] = calculo_momento_angular_total(pos[0], vel[0])

#Algoritmo de Verlet
for t in range (iterations-1):
    w = vel[t] + h/2*a[t]
    pos[t+1] = pos[t] + h*w
    a[t+1] = calculo_aceleracion(dis=pos[t+1])
    vel[t+1] = w + h/2*a[t+1]

    #Calculo la energía y el momento angular del sistema en la nueva posición.
    energia_total[t+1] = calculo_energia_total(pos[t+1], vel[t+1])
    momento_angular_total[t+1] = calculo_momento_angular_total(pos[t+1], vel[t+1])

    #Compruebo si algún planeta ha completado su órbita.
    calculo_periodo_orbital(t+1)


    if t % 10000 == 0:
        print(f'Iteración número {t+1} de {iterations-1}')

#Reescalamos los periodos calculados a los 
periodos_r = periodos*(d_tierra[0]**3/(G*M_sol))**0.5*1.15741e-5

for i in range(cuerpos):
    print('Para el cuerpo', i, ', el periodo medido es', periodos_r[i],'días. Podemos compararlo con su periodo real que es', T_teo[i] ,'días. Su error relativo es', abs(periodos_r[i]-T_teo[i])/T_teo[i]*100, '%.')

np.save('resultados/pos_global.npy', pos)
np.save('resultados/vel_global.npy', vel)
np.save('resultados/ac_global.npy', a)
np.save('resultados/E_global.npy', energia_total)
np.save('resultados/L_global.npy', momento_angular_total)
np.save('resultados/t_global.npy', t)

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
