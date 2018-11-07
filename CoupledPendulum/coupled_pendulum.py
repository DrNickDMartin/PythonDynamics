from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def run():
    theta1_0 = np.deg2rad(40)
    theta2_0 = np.deg2rad(0)

    sol = solve_ivp(cp_pen, [0, 20], [0, theta1_0, 0, theta2_0], rtol=1e-6)

    plt.plot(sol.t, np.rad2deg(sol.y[1,:]),".-")
    plt.plot(sol.t, np.rad2deg(sol.y[3,:]),".-")
    plt.show()

def cp_pen(t, z):

    zdot = [0,0,0,0]

    thetadot1 = z[0]
    theta1 = z[1]
    thetadot2 = z[2]
    theta2 = z[3]

    x0 = 100e-3 # Neutral length of elastic band mm
    b = 550e-3 # Seperation distance mm
    k = 50 # Elastic band stiffness N/m
    a1 = 400e-3 #mm
    a2 = 450e-3 #mm
    l1 = 600e-3 #mm
    l2 = 600e-3 #mm
    m1 = 220e-3 #kg
    m2 = 220e-3 #kg
    g = 9.81 #m/s**2
    H1 = 130e-3 #mm
    W1 = 60e-3 #mm
    H2 = 130e-3 #mm
    W2 = 60e-3 #mm

    # Calculate the geometric variables
    d = a2*np.cos(theta2) - a1*np.cos(theta1)
    c = b + a1*np.sin(theta1) - a2*np.sin(theta2)
    beta = np.arctan(d/c)
    x = d/np.sin(beta)
    alpha1 = beta - theta1
    alpha2 = beta - theta2

    # Calculate the spring force
    Fk = k*(x0-x)
    Fg1 = m1*g
    Fg2 = m2*g

    # Calculate the mass moment of intertia
    Izz1 = (1/12)*m1*(H1**2 + W1**2) + m1*l1**2
    Izz2 = (1/12)*m2*(H2**2 + W2**2) + m2*l2**2

    zdot[0] = (-a1*Fk*np.cos(alpha1) - l1*Fg1*np.sin(theta1))/Izz1
    zdot[1] = thetadot1
    zdot[2] = (a2*Fk*np.cos(alpha2) - l2*Fg2*np.sin(theta2))/Izz2
    zdot[3] = thetadot2

    return zdot

if __name__=="__main__":
    run()