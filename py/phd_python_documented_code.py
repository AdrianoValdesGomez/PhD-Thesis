import numpy as np
import math
from matplotlib import pyplot as plt
def orto(x):
    """Given a vector x, this function generates another one that is orthogonal with respect
    to the natural scalar product on 3 dimensional euclidean space E3"""
    if np.dot(x,x) == 0:
        return 'No se puede: ese es el vector cero!'
    else:
        if 0 not in x:
            v1 = 1
            v2 = -(x[0]/x[1])
            v3 = 0
            #return np.array([v1,v2,v3])
        else:
            if x[0] == 0:
                if x[1] == 0:
                    v1 = 1
                    v2 = 0
                    v3 = 0
                else:
                    v1 = 0
                    v2 = 0
                    v3 = 1
            elif x[1] == 0:
                v1 = 0
                v2 = 1
                v3 = 0
            else:
                v1 = 0
                v2 = 0
                v3 = 1
        return np.array([v1,v2,v3])
    

def base_ort_nor(x):
    """This function generates a basis of vetors orthogonal, with respect to the natural scalar
    product on 3 dimensional euclidean space, to the vector x"""
    y = orto(x)
    v1 = y/np.linalg.norm(y)
    z = np.cross(x,v1)
    v2 = z/np.linalg.norm(z)
    return v1, v2


def vector_des(v1,v2):
    """This function generates a unit random vetor, with unifom distribution among all possible directions,
    in the plane spanned by the linearly independent vectors v1 and v2"""
    na = 2 * np.pi*np.random.rand()
    vn = v1*np.cos(na) + v2*np.sin(na)
    return vn/np.linalg.norm(vn)




def vector_q(x,s):
    """This function returns a vector in the direction of x of length tan(s)"""
    q = np.tan(s)
    return q*x


def nuevo_r(r, vector_q):
    """This function returns a point on the unit sphere, originaly in r, a geodesic distant s, in the direction
    of vector_q"""
    y = r + vector_q
    y = y/np.linalg.norm(y)
    return y


def actualiza(r,s):
    """This function updates the position of a particleasembles one a function that generates a 
    random vector in the plane orthogonal to r and of displacement magintude s"""
    v1, v2 = base_ort_nor(r)
    pre_q = vector_des(v1,v2)
    q = vector_q(pre_q, s)
    return nuevo_r(r, q)


def act_n(lista, D, delta_t):
    """This function updates all the elements in a list using the actualiza function on each element"""
    l = []
    for v in lista:
        s = ese(D,delta_t)
        l.append(actualiza(v,s))
    return l


def b_steps_(ri,rf, n):
    """ This function returns n points between ri and rf along the geodesic whichs pass through them"""
    l = [ri]
    r0 = ri
    lamb = (np.dot(ri,rf))/((np.linalg.norm(ri))*(np.linalg.norm(rf)))
    
    if abs(lamb) > 1:
        #print('Is fucked up: there was a rounding ')
        if lamb < 0:
            lamb = -1
        else:
            lamb = 1
    
    
    
    theta = np.arccos(lamb)
    #if theta < 1e17:
        #return l
    if theta == 0:
        return [ri,rf]
    
    else:

        normal = np.cross(ri, rf)/ np.linalg.norm(np.cross(ri,rf))
        for i in range(1,n + 1):
            #vi = rot_theta(r0, theta/n, normal)
            vi = rot_finita(r0, -normal, theta/n)
            l.append(vi)
            r0 = vi
        return l




def trans_s_c(r,theta, phi):
    """This function transform from spherical coordinates (r,theta,phi) into cartesian coordinates (x,y,z)"""
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)* np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z

# Version mas reciente (26 de marzo del 2018)

def trans_c_s(x,y,z):
    """This function transform from cartesian coordinates (x,y,z) into spherical coordinates (r,theta,phi)"""
    r = np.sqrt(x**2 + y**2 + z**2)
    #print r
    cociente = z/r
    if abs(cociente) > 1:
        #print "cociente mayor a 1"
        if cociente < 0:
            theta = np.arccos(-1)
            #print "cociente negativo", theta
        else:
            theta = np.arccos(1)
            #print theta
    else:
        
        theta = np.arccos(z/r)
        #print theta
    #Aqui hay un problema relevante: cada vez que y o x sean nulos, habra un problema
    #de indefinicion de operacion
    if x == 0.:
        #print "x=0"
        if y == 0.:
            #print "y=0"
            phi = 2*np.pi*np.random.rand()
        else:
            if y > 0:
                phi = np.pi/2
            else:
                phi = -np.pi/2
    else:
        
        if x < 0:
            #print "x<0"
            phi = np.arctan(y/x) + np.pi
        else:
            #print "x debería ser positivo"
            if y < 0:
                #print "y<0 tambien!"
                phi = np.arctan(y/x) + 2*np.pi
            else:
                #print "x<0 pero y>0"
                phi = np.arctan(y/x)
            
    return r, theta, phi


def r_uni(theta, phi):
    """This function returns the unit vector in the direction in which r increases at the point (theta,phi)
    of the unit sphere"""
    x = np.sin(theta)*np.cos(phi)
    y = np.cos(theta)*np.cos(phi)
    z = np.cos(theta)
    return np.array([x,y,z])


def theta_uni(theta, phi):
    """This function returns the unit vector in the direction in which theta increases at the point (theta,phi)
    of the unit sphere"""
    x = np.cos(theta)*np.cos(phi)
    y = np.cos(theta)*np.sin(phi)
    z = -np.sin(theta)
    return np.array([x,y,z])


def phi_uni(theta, phi):
    """This function returns the unit vector in the direction in which phi increases at the point (theta,phi)
    of the unit sphere"""
    x = -np.sin(phi)
    y = np.cos(phi)
    z = 0
    return np.array([x,y,z])


def nombre(s):
    """This funtion returns a string of length 4 in which s (the imput) is a number and it has the necessary
    zeros to the left in order of being of 4 digit length"""
    diferencia = 5 - len(str(s))
    ceros = '' 
    for i in range(diferencia):
        ceros = ceros + '0'
    variable = ceros + str(s)
    return variable


def var(D, delta_t):
    """This function returns the varianza of a brownian particle in the 2 dimensional infinite plane, of self
    difussion coeffient D, and in an interval of time of length delta_t"""
    return 4 * D * delta_t



def ese(D, delta_t):
    """This function returns the magnitude of the random displacement displacement of a brownian particle
    that has a self-diffusion coeffient D, in a time interval of length delta_t"""
    return abs(np.random.normal(loc = 0., scale = np.sqrt(var(D,delta_t)),size = None))


def rot_finita(r_ini, N, Phi):
    """This function rotates a vector r_ini around the direction given by the vector N, an angle Phi. It uses 
    algebraic operations only and three trigonometric evaluations."""
    n = N/np.linalg.norm(N)
    r_fin = np.cos(Phi)*r_ini + (np.dot(n,r_ini))*(1 - np.cos(Phi))*n + (np.sin(Phi))*(np.cross(r_ini,n))
    return r_fin


#Funcion que regresa una lista de n numpy arrays que son l
def Trayectoria(ri,rf,n):
    """This function returns a list of a sequence of n points equally angular spaced between r1 and rf"""
    l = [ri]
    r0 = ri
    theta = np.arccos((np.dot(ri,rf))/((np.linalg.norm(ri))*(np.linalg.norm(rf))))
    N = np.cross(ri, rf)
    
    for i in range(1,n + 1):
        vi = rot_finita(r0, N, theta/n)
        l.append(vi)
        r0 = vi
    return l




def plot_particles(lista, vpolar, vazim, numero, titulo):
    """This function plots a list of vectors given in lista on the surface of a unit sphere with a view 
    an angle vpolar of inclination with respect to the plane x-y and rotated an angle of vazim angular units
    in the direction around the z axis"""
    from mpl_toolkits.mplot3d import axes3d
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    from itertools import product, combinations
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection='3d')
    #ax.set_aspect("equal")
    ax._axis3don = False


    #draw sphere
    R = 1
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x=R*np.cos(u)*np.sin(v)
    y=R*np.sin(u)*np.sin(v)
    z=R*np.cos(v)
    #ax.plot_surface(x, y, z, color="r", alpha = 0.15)

    ax.plot_surface(x, y, z, cmap=cm.YlGnBu_r,rstride=1, cstride=1, alpha = 0.10, linewidth = 0.10)
    #ax.view_init(vpolar, vazim)
    
    
    #draw an arrow or a set of arrow
    #ax.quiver(0,0,1.5,0,0,1, length=0.5, arrow_length_ratio = .5, color = "b")
    #ax.quiver(1.5,0,0,1,0,0, length=0.5, arrow_length_ratio = .5, color ="g")
    #ax.quiver(0,1.5,0,0,1,0, length=0.5, arrow_length_ratio = .5, color ="r")
    #draw patch
    #u, v = np.mgrid[0:2*np.pi:50j, 0:(np.pi/7):50j]
    #x=R*np.cos(u)*np.sin(v)
    #y=R*np.sin(u)*np.sin(v)
    #z=R*np.cos(v)
    #ax.plot_surface(x, y, z, color="r", alpha = 0.25)    
    
    #draw points
    for p in lista:
        ax.scatter([p[0]],[p[1]],[p[2]],color="b",s=15, alpha = 0.25)
    
    ax.view_init(vpolar, vazim)
    fig.savefig('{}_Img_{}.png'.format(titulo,nombre(numero)))
    
    #plt.show()
    plt.close()

    

def polo_n(n, R):
    """This function generates a list of n vectors in the position [0,0,1]; the north pole in the unit sphere"""
    l = []
    for i in range(n):
        l.append(np.array([0,0,R]))
    return l



def obs_uniforme(N, R, size):
    
    list_obs = []
    omega = np.cos(size)
    while len(list_obs) < N:
        x, y, z = np.random.uniform(-1,1), np.random.uniform(-1,1), np.random.uniform(-1,1)
        v = np.array([x, y, z])
        norma = np.linalg.norm(v)
        if norma <= R:
            n = v/norma
            if not np.dot(n, np.array([0.,0.,1.]))/R > omega:
                list_obs.append(R*n)
    
    return list_obs    

#Possible Initial Distributions

#Distribución en el polo sur

def polo_s(n):
    """This function returns a list of n vectors all in the position [0,0,-1]"""
    l = []
    for i in range(n):
        l.append(np.array([0,0,-1]))
    return l

#Dsitribución en el ecuador

def dist_ecuador(n):
    """This function returns a list of n vectors all uniformly distributed on the equator of a unit sphere"""
    l = []
    for i in range(n):
        x, y, z = trans_s_c(1.,np.pi/2, np.random.uniform(0,2*np.pi))
        l.append(np.array([x,y,z]))
    return l




########################################################################################################################
# FIELDS



coff10 = np.sqrt((3./(4*np.pi)))

def field_y10(ri, v0):
    """This function generates a vector in the direction of a field proportional to the gradient of the 
    spherical harmonic Y10"""
    r, theta, phi = trans_c_s(ri[0],ri[1],ri[2])
    
    theta_field = coff10 * (np.sin(theta))
    
    phi_field = 0.
    
    field = v0*(theta_field)*theta_uni(theta, phi)
    
    return field


# FIELD partial_theta Y20

coff20 = np.sqrt((5./(16*np.pi)))

def field_y20(ri, v0):
    """This function generates a vector in the direction of a field proportional to the gradient of the 
    spherical harmonic Y20"""
    r, theta, phi = trans_c_s(ri[0], ri[1], ri[2])
    
    theta_field = coff20 * (-6 * np.cos(theta) * np.sin(theta))
    
    phi_field = 0.
    
    field = v0*(theta_field)*theta_uni(theta, phi)
    
    return field



coff30 = np.sqrt((7./(16*np.pi)))

def field_y30(ri, v0):
    """This function generates a vector in the direction of a field proportional to the gradient of the 
    spherical harmonic Y30"""
    r, theta, phi = trans_c_s(ri[0], ri[1], ri[2])
    
    theta_field = coff30 * (-15 * np.cos(theta)**2 * np.sin(theta) + 3*np.sin(theta))
    
    phi_field = 0.
    
    field = v0 * (theta_field) * theta_uni(theta, phi)
    
    return field


coff43 =  - 105 * np.sqrt((5./(math.factorial(7) * np.pi)))

def field_y43(ri, v0):
    """This function generates a vector in the direction of a field proportional to the gradient of the 
    spherical harmonic Y43"""
    
    r, theta, phi = trans_c_s(ri[0], ri[1], ri[2])
    
    theta_field = coff43 * (np.cos(3 * phi) * np.sin(theta)**2 * (- np.sin(theta)**2 + 3 * np.cos(theta)**2 ))
    
    phi_field = coff43 * (- 3 * np.cos(theta) * np.sin(theta)**3 * np.sin(3 * phi))
    
    
    field = v0 * ( theta_field * theta_uni(theta, phi) + phi_field * phi_uni(theta, phi) )
    
    return field



###########################################################################################################################
###########################################################################################################################

# TIME DEPENDENT FIELDS


def field_harmonic(ri,v0,t, omega):
    """This field is harmonic in time and has azimuthal symmetry"""
    
    r, theta, phi = trans_c_s(ri[0], ri[1], ri[2])
    
    theta_field = np.sin(8 * theta) * np.cos(omega * t)
    
    phi_field = 0.
    
    
    field = v0 * ( theta_field * theta_uni(theta, phi) + phi_field * phi_uni(theta, phi) )
    
    return field





def nuevo_r_field(r, vector_q, field):
    """This function updates the position r of a brownian particle urged by the white noise (vector_r)
    and an external field acting on the tangent plane of a unit sphere at r"""
    v_sum = field + vector_q
    xf = np.linalg.norm(v_sum)
    if xf == 0:
        nfield = 0
    else:
        fvecuni = (v_sum)/xf
        nfield = np.tan(xf)*fvecuni
    
    y = r + nfield
    y = y/np.linalg.norm(y)
    return y



def actualiza_field(r, s, v0, field):
    """This function generates the random white noise tangent to the unit sphere at r and applies
    the updating recipe to consider an external field also"""
    v1, v2 = base_ort_nor(r)
    pre_q = vector_des(v1,v2)
    q = vector_q(pre_q, s)
    return nuevo_r_field(r, q, field(r,v0))



def act_n_field(lista, v0, field, D, dt):
    """This funtion updates a list of positions of a set of brownian particles in which there is also
    an external field acting on them"""
    l = []
    for v in lista:
        s = ese(D,dt)
        l.append(actualiza_field(v, s, v0, field))
    return l



#def act_n2_field(arreglo, D, delta_t, v0, field):
#    
#    #s = ese(D,delta_t)

#    arreglo = np.apply_along_axis(actualiza_field, 1, arreglo, s, v0, field)
#    return arreglo




##############################################################################################################
##############################################################################################################

# Some Analysis functions



def mean_var_hist_theta(lista):
    """This function returns from a list of vectors a list of the theta angle of each vector
     along with the mean and the variance of this set"""
    thetas = []
    for r in lista:
        cociente = r[2]
        if abs(cociente) > 1.:
            if cociente < 0:
                theta = np.arccos(-1)
            else:
                theta = np.arccos(1)
        else:

            theta = np.arccos(r[2])
            
        thetas.append(theta)
        
    return thetas, np.mean(thetas), np.var(thetas)



def mean_var_hist_phis(lista):
    """This function returns from a list of vectors a list of the phis angles of each vector
     along with the mean and the variance of this set"""
    phis = []
    for r in lista:
        
        if r[0] == 0:
            if r[1] == 0:
                phi = 2*np.pi*np.random.rand()
        else:
            if r[0] < 0:
                
                phi = np.arctan(r[1]/r[0]) + np.pi
                    
            else:
                if r[1] < 0:
                    phi = np.arctan(r[1]/r[0]) + 2*np.pi
                else:
                
                    phi = np.arctan(r[1]/r[0])
                
        phis.append(phi)
    return phis, np.mean(phis), np.var(phis)





##############################################################################################################
##############################################################################################################
# Generalized Langevin Equation

def act_gLan(pos_ini, vel_ini, Ext_F, ud_rand_F, m, gamma, delta_t):
    """This function updates de position and velocity of a diffusion process restricted to move
    on the surface of the unit sphere"""
    
    vel_fin = vel_ini*(1 - (gamma/m) * delta_t) + ud_rand_F + Ext_F * delta_t
    
    if np.linalg.norm(vel_ini) > 0.:
        
        vel_fin = rot_finita(vel_fin, np.cross(vel_ini, pos_ini), np.linalg.norm(vel_ini * delta_t))

    s = np.linalg.norm(vel_ini * delta_t)
    
    #if s > 1e-8:
    if s > 0 and s < np.pi/2:
        
        uni_des = vel_ini * delta_t / s
    
        q = np.tan(s)
        pos_fin = pos_ini + q*uni_des
        pos_fin = pos_fin / np.linalg.norm(pos_fin)
    else:
        pos_fin = pos_ini
    #print(np.linalg.norm(pos_fin))
    return pos_fin, vel_fin


Kb = 13.8064852


def ese_V(T, m, delta_t, gamma):
    """This function returns the variance of the velocity of a brownian particle whose velocity is normally 
    distributed and is at absolute temperature T, has a mass m, the friction coefficient of the medium is gamma
    and the fluctuatuons is considered in a window of time delta_t"""
    sigma2 = 4 * Kb * T * gamma * delta_t / m**2
    
    return abs(np.random.normal(loc=0., scale= np.sqrt(sigma2), size=None))


def tangent_white_noise(x, delta_t, T, m, gamma):
    """This function returns a vector with a 2 dimensional gaussian distribution of variance 4 K T gamma/m**2
    on the ttangent plane to the unit sphere at x"""
    v1, v2 = base_ort_nor(x)
    s = ese_V(T, m, delta_t ,gamma )
    return s * vector_des(v1,v2)

def act_ensamble(l_pos, l_vel, Ext_F, U0, m, T, gamma, delta_t):
    """This function applies the updating recipe act_gLan to a list of initial positions and velocities
    l_pos, and l_vel and returns the updates of these lists"""
    l_npos = []
    l_nvel = []
    for i in range(len(l_pos)):
        pos, vel = l_pos[i], l_vel[i]
        random_wn = tangent_white_noise(pos, delta_t, T, m, gamma)
        #pos_fin, vel_fin = act_gLan(pos, vel, Ext_F(*args), random_wn, m, gamma, delta_t)
        pos_fin, vel_fin = act_gLan(pos, vel, Ext_F(pos, U0), random_wn, m, gamma, delta_t)
        #print(Ext_F(pos, U0))
        l_npos.append(pos_fin)
        l_nvel.append(vel_fin)
    return l_npos, l_nvel

def act_ensamble_td_field(l_pos, l_vel, Ext_F, U0, m, T, gamma, delta_t, t, omega):
    """This function applies the updating recipe act_gLan to a list of initial positions and velocities
    l_pos, and l_vel and returns the updates of these lists"""
    l_npos = []
    l_nvel = []
    for i in range(len(l_pos)):
        pos, vel = l_pos[i], l_vel[i]
        random_wn = tangent_white_noise(pos, delta_t, T, m, gamma)
        #pos_fin, vel_fin = act_gLan(pos, vel, Ext_F(*args), random_wn, m, gamma, delta_t)
        pos_fin, vel_fin = act_gLan(pos, vel, Ext_F(pos, U0, t, omega), random_wn, m, gamma, delta_t)
        #print(Ext_F(pos, U0))
        l_npos.append(pos_fin)
        l_nvel.append(vel_fin)
    return l_npos, l_nvel



def act_ensamble_arr(arr_pos, arr_vel, Ext_F, U0, m, T, gamma, delta_t):
    """This function updates the position and velocity of an array of initial positions and velocities
    It uses np.apply_along_axis instead of a list"""

    random_wn = tangent_white_noise(arr_pos, delta_t, T, m, gamma)
    Ext = Ext_pos(pos,U0)
    #s = ese(D,delta_t)
    #arreglo = np.apply_along_axis(actualiza_field, 1, arreglo, s, v0, field)
    #(pos, vel, Ext_F(pos, U0), random_wn, m,gamma,delta_t)
    pos_fin, vel_fin = np.apply_along_axis(act_gLan,1,arr_pos,arr_vel, Ext, random_wn, m, gamma, delta_t)

    return pos_fin, vel_fin


def n_vel_ini(x, n, delta_t, T, m, inten, gamma):
    """This functions returns a list of n vectors in the tangent plane to the unit sphere at x,
    normally distributed and multipplied by the factor inten"""
    l_vel = []
    for i in range(n):
        l_vel.append(inten * tangent_white_noise(x,delta_t,T,m,gamma))
    return l_vel





# Functions for the statistics

def theta_phi(l_pos):
    thetas = []
    phis = []
    for pos in l_pos:
        r,theta, phi = trans_c_s(pos[0],pos[1],pos[2])
        thetas.append(theta)
        phis.append(phi)
    return thetas, phis

def momentos_theta(thetas):
    promedio = np.mean(thetas)
    varianza = np.var(thetas)
    return promedio, varianza


##############################################################################################################
##############################################################################################################

# SIMULATIONS WITH OBSTACLES

def particion_esfera(ccero, Nphi):
    Ntheta = int(4*np.pi/(ccero*Nphi))
    print('Ntheta', Ntheta, 'Nphi', Nphi, 'Ntheta*Nphi', Ntheta*Nphi)
    sigmaPhi = 2*np.pi/Nphi
    deltaphi = 2*np.pi/Nphi
    thetas = []
    phis = [0]
    cociente = ccero/sigmaPhi
    for i in range(Ntheta + 1):
        theta = np.arccos(1 - (i)*cociente)
        thetas.append(theta)
    for j in range(Nphi):
        phis.append(phis[j] + deltaphi)
    return thetas, phis


def secuencia_part(tamini, Nfi, numero):
    l1, l2 = particion_esfera(4*np.pi/tamini, Nfi)
    particion = []
    for i in range(len(l2)):
        for j in range(len(l1)):
            x, y, z = trans_s_c(1, l1[j], l2[i])
            particion.append(np.array([x, y, z]))
            
    return plot_particles(particion, 45, 45, numero)



def coordenadas_centro(l1,l2):
    """Funcion que regresa las coordenadas del centro de dos arreglos para 
    las coordenadas theta y phi"""
    thetas_centro = []
    phis_centro = []
    for i in range(len(l1) - 1):
        theta_media = l1[i] + (l1[i + 1] - l1[i])/2.
        thetas_centro.append(theta_media)
    for j in range(len(l2) - 1):
        phi_media = l2[j] + (l2[j + 1] - l2[j])/2.
        phis_centro.append(phi_media)
    return thetas_centro, phis_centro


def secuencia_obs(N, Nfi, numero):
    l1_prima, l2_prima = particion_esfera(4*np.pi/N, Nfi)
    l1, l2 = coordenadas_centro(l1_prima, l2_prima)
    particion = []
    for i in range(len(l2)):
        for j in range(len(l1)):
            x, y, z = trans_s_c(1, l1[j], l2[i])
            particion.append(np.array([x, y, z]))
            
    print(len(particion))
    
    #return plot_particles(particion, 0, 0, numero)
    return particion


def plot_particle_traj_obs(lista_obstaculos, trayectoria,  vpolar, vazim, numero, title):
    from mpl_toolkits.mplot3d import axes3d
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    #import matplotlib.pyplot as plt
    #import numpy as np
    from itertools import product, combinations
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection='3d')
    #ax.set_aspect("equal")
    ax._axis3don = False
    ax.margins(.2,.2,.2)
    

    



    #draw sphere
    R = 1
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x=R*np.cos(u)*np.sin(v)
    y=R*np.sin(u)*np.sin(v)
    z=R*np.cos(v)
    #ax.plot_surface(x, y, z, color="r", alpha = 0.15)

    ax.plot_surface(x, y, z, cmap=cm.YlGnBu_r,rstride=1, cstride=1, alpha = 0.10, linewidth = 0.15)
    ax.view_init(vpolar, vazim)
    #draw patch
    #u, v = np.mgrid[0:2*np.pi:50j, 0:(np.pi/7):50j]
    #x=R*np.cos(u)*np.sin(v)
    #y=R*np.sin(u)*np.sin(v)
    #z=R*np.cos(v)
    #ax.plot_surface(x, y, z, color="r", alpha = 0.25)
    
    
    
    #draw obstacles
    
    for p in lista_obstaculos:
        ax.scatter([p[0]],[p[1]],[p[2]], color="b", s=10, alpha = 0.2)
    
    #draw trajectory
    for p in trayectoria:
        ax.scatter([p[0]],[p[1]],[p[2]], color="k",s=20, alpha = 0.7)
    
    #Plot the x positive direction
    
    #ax.quiver(1.5,0,0,1,0,0, length=0.5, arrow_length_ratio = .5)
    #ax.quiver(0,1.5,0,0,1,0, length=0.5, arrow_length_ratio = .5)
    #ax.quiver(0,0,1.5,0,0,1, length=0.5, arrow_length_ratio = .5)
    
    fig.savefig('{}-{}.png'.format(title,nombre(numero)))
    #ax.view_init(80, 30)
    plt.close()
    #plt.show()

def obs_uniforme(N, R, size):
    """This function generates a sequence of N obstacles uniformly distributed on the surface of a unit sphere
    outside a neiborhood of radious R of given angular aperture size"""
    list_obs = []
    omega = np.cos(size)
    while len(list_obs) < N:
        x, y, z = np.random.uniform(-1,1), np.random.uniform(-1,1), np.random.uniform(-1,1)
        v = np.array([x, y, z])
        norma = np.linalg.norm(v)
        if norma <= R:
            n = v/norma
            if not np.dot(n, np.array([0.,0.,1.]))/R > omega:
                list_obs.append(R*n)
    
    return list_obs

def obs_uniforme_S2(N, size):
    """This function generates a sequence of N obstacles uniformly distributed on the surface of a unit sphere
    outside a neiborhood of radious R of given angular aperture size"""
    list_obs = []
    omega = np.cos(size)
    while len(list_obs) < N:
        x, y, z = np.random.uniform(-1,1), np.random.uniform(-1,1), np.random.uniform(-1,1)
        v = np.array([x, y, z])
        norma = np.linalg.norm(v)
        if norma <= 1:
            n = v/norma
            if not np.dot(n, np.array([0.,0.,1.])) > omega:
                list_obs.append(n)
    
    return list_obs    


def puntos_obs_j(r_omega, theta_omega, n):
    r , theta, phi = trans_c_s(r_omega[0],r_omega[1],r_omega[2])
    rp = rot_finita(r_omega, phi_uni(theta, phi), theta_omega)
    puntos_obs_j = [rp]
    for i in range(1,n):
        x = rot_finita(rp, r_omega, 2*np.pi/n)
        puntos_obs_j.append(x)
        rp = x
    return puntos_obs_j



def puntos_obs(lista_obstaculos, size, n_countor):
    mis_obs = []
    for i in range(len(lista_obstaculos)):
        a = lista_obstaculos[i]
        b = size
        mis_obs = mis_obs + puntos_obs_j(a, b, n_countor)
    return mis_obs
#Collision_check es una función que, dada una trayectoria: una lista de vectores que
#pasan por puntos sucesivos de la trayectoria, verifica si alguna de estas posiciones
#interesecto a alguno de los obstáculos. En caso de que así sea, actualiza conforme una
#colision elastica. En caso de no intersectar a ningun obstaculo regresa una lista
#con dos vectores: posicion inicial y posicion final en ese orden.

def penetrate_obs(lista_vect, lista_obs, size):
    metiches = []
    for obs in lista_obs:
        r_omega, theta_omega = obs, size
        #frontera = .2
        #metiches = []
        for v in lista_vect:
            tamanho = np.cos(theta_omega)
            if np.dot(v,r_omega) >= tamanho:
                #print('Penetro el mother fucker obstacle')
                metiches.append(v)
                
            else:
                continue
    #print('no choco el mother fucker')
    #valor = False
    return metiches
    

def check_collision(lista_vect, lista_obs, size):
    for obs in lista_obs:
        r_omega, theta_omega = obs, size 
        for v in lista_vect:
            tamanho = np.cos(theta_omega)
            if np.dot(v,r_omega) > tamanho:
                return  True
            else:
                continue
    return False
    

def tangent_space(x1,x2,xo):
    np_prima = np.cross(np.cross(xo, x1), x1)
    nor_p = np_prima/np.linalg.norm(np_prima)
    up_prima = np.cross(np.cross(x1, x2), x1)
    up = up_prima/np.linalg.norm(up_prima)
    tp_prima = np.cross(x1, nor_p)
    tp = tp_prima/np.linalg.norm(tp_prima)
    y = (np.dot(up,tp))*tp - (np.dot(up, nor_p))*nor_p
    v_rot_prima = np.cross(x1, y)
    v_rot = v_rot_prima/np.linalg.norm(v_rot_prima)
    return v_rot

def compute_all(lista_vec):
    thetas = []
    phis = []
    thetaphis = []
    cos_thetas = [] 
    
    for v in lista_vec:
        r, theta, phi = trans_c_s(v[0],v[1],v[2])
        thetas.append(theta)
        phis.append(phi)
        thetaphis.append(theta*phi)
        cos_thetas.append(np.cos(theta))
    
    mean_theta, mean_phi = np.mean(thetas), np.mean(phis)
    
    var_theta, var_phi = np.var(thetas), np.var(phis)
    
    mean_theta_phi = np.mean(thetaphis)
    
    mean_cos_theta = np.mean(cos_thetas)
    
    return thetas, phis, mean_theta, mean_phi, var_theta, var_phi, mean_theta_phi, mean_cos_theta




def simulacion_obs(lista_obs, estructura, size_obs, initial_cond, sim_size, sensitiveness, D, delta_t):
    """This function simulates the diffusion of brownian particles on the unit two dimensional sphere in which
    there are certain prohibited regions, defined by its centers lista_obs and angular apertures size_obs"""
    to_update = initial_cond
    thetas_t = []
    phis_t = []
    m_theta_t = []
    m_phi_t = []
    v_theta_t = []
    v_phi_t = []
    cov_theta_phi_t = []
    m_cos_theta_t = []
    #plot_particle_traj_obs(estructura, to_update, 90,0,0)

 

    for i in range(sim_size):
        tentative_paths = []
        updated_pos_at_t = []
        tentative_pos = act_n(to_update, ese(D, delta_t))
        for j in range(len(tentative_pos)):
            tentative_paths.append(b_steps_(to_update[j], tentative_pos[j], sensitiveness))
        
        for path in tentative_paths:
            if check_collision(path, lista_obs, size_obs):
                for k in range(1,len(path)):
                    if check_collision([path[k]], lista_obs, size_obs):
                        updated_pos_at_t.append(path[k-1])
                        break
            else:
                updated_pos_at_t.append(path[-1])
        
        #plot_particle_traj_obs(estructura, updated_pos_at_t, 90,0,i + 1)
        thetas, phis, mean_theta, mean_phi, var_theta, var_phi, mean_theta_phi, mean_cos_theta = compute_all(updated_pos_at_t)
        thetas_t.append(thetas)
        phis_t.append(phis)
        m_theta_t.append(mean_theta)
        m_phi_t.append(mean_phi)
        v_theta_t.append(var_theta)
        v_phi_t.append(var_phi)
        cov_theta_phi_t.append(mean_theta_phi - mean_theta*mean_phi)
        m_cos_theta_t.append(mean_cos_theta)
        
        #print penetrate_obs(updated_pos_at_t, lista_obs , size_obs)
        to_update = updated_pos_at_t
    return thetas_t, phis_t, m_theta_t, m_phi_t, v_theta_t, v_phi_t, cov_theta_phi_t, m_cos_theta_t



def adapted_path(pos_ini, vel_ini, field, U0, m, T, gamma, delta_t, l_obs, size_obs, sensitiveness):
    to_update = pos_ini
    tentative_paths = []
    updated_pos_at_t = []
    # Here it would be desireble; more economic computational speaking, to update just the position
    # of the particles using the velocities, but leave the velocities untouch.
    tentative_pos, tentaive_vel = act_ensamble(pos_ini, vel_ini, field, U0, m, T, gamma, delta_t)

    # Generates the points between initial and final points for each particle called the path

    for j in range(len(tentative_pos)):
        tentative_paths.append(b_steps_(to_update[j], tentative_pos[j], sensitiveness))

    # In this part we check which paths intersect an obstacle

    for path in tentative_paths:
        if check_collision(path, l_obs, size_obs):
            for k in range(1,len(path)):
                if check_collision([path[k]], l_obs, size_obs):
                    updated_pos_at_t.append(path[k - 1])
                    break
        else:
            updated_pos_at_t.append(path[-1])
            
    return updated_pos_at_t


def adapted_path_tfield(pos_ini, vel_ini, field, U0, m, T, gamma, delta_t, l_obs, size_obs, sensitiveness, t, omega):
    to_update = pos_ini
    tentative_paths = []
    updated_pos_at_t = []
    # Here it would be desireble; more economic computational speaking, to update just the position
    # of the particles using the velocities, but leave the velocities untouch.
    tentative_pos, tentaive_vel = act_ensamble_td_field(pos_ini, vel_ini, field, U0, m, T, gamma, delta_t, t, omega)

    # Generates the points between initial and final points for each particle called the path

    for j in range(len(tentative_pos)):
        tentative_paths.append(b_steps_(to_update[j], tentative_pos[j], sensitiveness))

    # In this part we check which paths intersect an obstacle

    for path in tentative_paths:
        if check_collision(path, l_obs, size_obs):
            for k in range(1,len(path)):
                if check_collision([path[k]], l_obs, size_obs):
                    updated_pos_at_t.append(path[k - 1])
                    break
        else:
            updated_pos_at_t.append(path[-1])
            
    return updated_pos_at_t


def init_uni_dist_out_obs(n_part, l_obs, obs_size):
    """This function generates a uniform initial distribution of particles in which all of 
    them are outside of all obstacles in the list l_obs"""
    l_part = []
    omega = np.cos(obs_size)
    while len(l_part) < n_part:
        # We generate a possible location for particle n
        x, y, z = np.random.uniform(-1,1), np.random.uniform(-1,1), np.random.uniform(-1,1)
        v = np.array([x, y, z])
        norma = np.linalg.norm(v)
        if norma <= 1:
            n = v/norma
        # Then we check whether this particle lies inside one of the obstacles
            l_trues = []
            for obs in l_obs:
                if not np.dot(n, obs) > omega:
                    l_trues.append(True)
            if sum(l_trues) == len(l_obs):
                l_part.append(n)
                
    return l_part

def coeficiente(n,N,D):
    #N es el numero de pasos, por lo que el tiempo total transcurrido hasta ese momento es
    t = dt * N
    return ((2*n + 1.)/(4*np.pi)) * np.exp(-n*(n + 1)*D*t)


def distribucion(theta, N, orden,D):
    #lista = []
    c = np.zeros(orden)
    for n in range(orden):
        c[n] = coeficiente(n,N,D)
        #lista.append(coeficiente(n,N)*np.polynomial.legendre.legval(np.cos(theta), c))
    lista = np.polynomial.legendre.legval(np.cos(theta), c)
    #print lista
    #return sum(lista)
    return lista

def distribucion_weighted(theta, N, orden,D):
    #lista = []
    c = np.zeros(orden)
    for n in range(orden):
        c[n] = coeficiente(n,N,D)
        #lista.append(coeficiente(n,N)*np.polynomial.legendre.legval(np.cos(theta), c))
    lista = np.polynomial.legendre.legval(np.cos(theta), c)*np.sin(theta)*2*np.pi
    #print lista
    #return sum(lista)
    return lista


    def plot_free_diff(thetas, l_steps, title, dpi, save=False, show=True):
    """This function plots """
    theta = np.linspace(0.0,np.pi,1000)
    tiempos = []
    suma = 0.
    for i in range(len(thetas)):
        suma =+ i*dt
        tiempos.append(suma)
        
    mean_cos_thetas_t =[]
    for ensamble in thetas:
        ensamble = np.array(ensamble)
        cos_ens = np.cos(ensamble)
        mean_cos_thetas_t.append(np.mean(cos_ens))

    fig, ax = plt.subplots(figsize=(7,6))

    ax.spines['bottom'].set_position('zero')
    for i in l_steps:

        #plt.hist(thetas[i], bins=int(( np.array(thetas[i][:]).max() - np.array(thetas[i][:]).min() )*40),
        #         density=True,color='k', alpha=0.14, label='Numeric Solution', edgecolor='black', linewidth=.55);
        ax.hist(thetas[i], bins=int(( np.array(thetas[i][:]).max() - np.array(thetas[i][:]).min() )*40),
                 density=True,color='cadetblue', alpha=0.25, label='Numeric Solution', edgecolor='black', linewidth=.55);
        ax.set_ylim(-.5,6.2)

        if ((distribucion_weighted(theta,i,1000,D,dt)[np.argmax(distribucion_weighted(theta,i,1000,D,dt))]+ 0.35) < 6.5):
            ax.text(theta[np.argmax(distribucion_weighted(theta,i,1000,D,dt))], 
                     distribucion_weighted(theta,i,1000,D,dt)[np.argmax(distribucion_weighted(theta,i,1000,D,dt))]+ 0.2,
                     't={0:.2f}'.format(i*dt), fontsize=12)
        ax.plot(theta, distribucion_weighted(theta,i,1000,D,dt), color = 'k', alpha = .96,#
                         label=r"$P(\theta,t|0,0) = \sum P_{l}(\cos{\theta})\, \exp{[-l(l+1)D\,t]}$",linewidth=.3);
    ax.set_xlabel(r'$\theta$', size=19)
    ax.set_ylabel(r'$P(\theta,t|0,0)$', size=19)

    ax.set_xticks([i*.5 for i in range(7) ])
    ax.set_yticks([i for i in range(7)])

    ax.set_xticklabels([str(i*.5) for i in range(7) ], fontsize=13)
    ax.set_yticklabels([str(i) for i in range(7)] , fontsize=13)

    ins = inset_axes(ax, 
                        width="50%", # width = 30% of parent_bbox
                        height=2.0, # height : 1 inch
                        loc=1)
    #print(len(tiempos))
    ins.plot(tiempos, mean_cos_thetas_t, label= 'Numeric', alpha=.75, linewidth=2.5,
             color = 'cadetblue', linestyle=':')
    #plt.scatter(tiempos, mean_cos_thetas_t, label= 'Numeric', marker='o', s=5, alpha=.05)
    ins.plot(tiempos, np.exp(-(2*D*np.array(tiempos))), label='Analytic', color='k', linewidth=0.85, linestyle='-')
    ins.legend(fontsize=12)
    ins.set_xlabel(r"$t$", fontsize=15)
    ins.set_ylabel(r"$\langle \mathbf{n}(t)\cdot \mathbf{n}(0) \rangle$", fontsize=14);
    #ax.set_tight_layout();
    if save:
        fig.savefig(title, dpi=dpi)
    if show:
        plt.show()





def plot_particles_Ylm(lista, vpolar, vazim, numero):
    """This function plots an ensamble of particles on the surface of a unit sphere and a heatmap of the potential
    interaction given by the first Spherical Harmonic with azimuthal symmetry, a Legendre Polynomial. vpolar and vazim
    the angles for the perspective and numero is the number of a sequence to be appended to the saved exported file
    with the purpose of making animations"""

    #import matplotlib.pyplot as plt
    #import numpy as np
    from itertools import product, combinations
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection='3d')
    #ax.set_aspect("equal")
    ax._axis3don = False

    



    #draw sphere
    
    u, v = np.mgrid[0:2*np.pi:60j, 0:np.pi:60j]
    
    radio = np.sqrt(3./(4*np.pi))*np.cos(v)
    
    x=1.05*np.cos(u)*np.sin(v)
    y=1.05*np.sin(u)*np.sin(v)
    z=1.05*np.cos(v)
    #ax.plot_surface(x, y, z, color="r", alpha = 0.15)
    
    
    fcolors = np.sqrt(3./(4*np.pi))*np.cos(v)
    fmax, fmin = fcolors.max(), fcolors.min()
    fcolors = (fcolors - fmin)/(fmax - fmin)

    ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.viridis(fcolors), alpha = 0.38)

    #ax.plot_surface(x, y, z, cmap=cm.YlGnBu_r, rstride=1, cstride=1, alpha = 0.10, linewidth = 0.10)
    ax.view_init(vpolar, vazim)
    
    
      
    
    #draw particles
    for p in lista:
        ax.scatter([p[0]],[p[1]],[p[2]],color="cadetblue", s=20, alpha = 0.85,
                   linewidth=1.0, edgecolor='k')
    
    fig.savefig('SHY10_Field_X01_Img{}.png'.format(nombre(numero)))
    #ax.view_init(80, 30)
    #plt.show()
    plt.close()

