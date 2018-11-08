from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def run():
    theta0 = np.deg2rad(20)

    sol = solve_ivp(pend, [0, 20], [0, theta0], rtol=1e-6)

    plt.plot(sol.t, np.rad2deg(sol.y[1,:]),".-")
    plt.show()

def pend(t, z):

    zdot = [0,0]

    thetadot = z[0]
    theta = z[1]

    l = 600e-3 #mm
    m = 220e-3 #kg
    g = 9.81 #m/s**2
    H = 130e-3 #mm
    W = 60e-3 #mm
    Fg = m*g

    # Calculate the mass moment of intertia
    Izz = (1/12)*m*(H**2 + W**2) + m*l**2

    zdot[0] = (-l*Fg*np.sin(theta))/Izz
    zdot[1] = thetadot

    return zdot

if __name__=="__main__":
    run()