from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def run():
    theta0 = np.deg2rad(48)

    sol = solve_ivp(pend, [0, 20], [theta0, 0], rtol=1e-6, dense_output=True)

    theta = sol.y[0,:]
    thetadot = sol.y[1,:]
    tehtadotdot = (pend(sol.t, sol.y))[1]*0.5

    lin_exp = np.loadtxt('Pendulums/linearpen3.csv', delimiter=',', skiprows=1)
    ang_exp = np.loadtxt('Pendulums/angularpen1.csv', delimiter=',', skiprows=1)

    trim = np.logical_and(lin_exp[:,0] > 0.73,lin_exp[:,0] < 20.73)
    acc_x = lin_exp[trim,1]
    acc_y = lin_exp[trim,2]
    acc_z = lin_exp[trim,2]
    n = acc_x.shape[0]
    texp = np.linspace(0,20,n)

    trim = np.logical_and(ang_exp[:,0] > 6.6,ang_exp[:,0] < 26.6)
    aac_x = ang_exp[trim,1]
    aac_y = ang_exp[trim,2]
    aac_z = ang_exp[trim,3]
    texp2 = np.linspace(0,20,aac_x.shape[0])

    plt.scatter(texp[range(0,n,5)], acc_y[range(0,n,5)], facecolors='none', edgecolors='b')
    plt.plot(sol.t, tehtadotdot, 'r', linewidth=2)
    #plt.plot(texp, acc_y)
    #plt.plot(texp2, aac_z**2)
    plt.legend(['model','experiemnt'])
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (rad/s)')
    plt.show()

def pend(t, z):

    zdot = [0,0]

    thetadot = z[1]
    theta = z[0]

    l = 500e-3 #mm
    m = 199e-3 #kg
    g = 9.81 #m/s**2
    H = 150e-3 #mm
    W = 80e-3 #mm
    Fg = m*g
    c = 0.0014 #damping factor

    # Calculate the mass moment of intertia
    Izz = (1/12)*m*(H**2 + W**2) + m*l**2

    zdot[1] = (-l*Fg*np.sin(theta))/Izz - (c*thetadot)/Izz
    zdot[0] = thetadot

    return zdot

if __name__=="__main__":
    run()