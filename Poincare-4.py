import numpy as np
import matplotlib.pyplot as plt
import timeit
from scipy.integrate import odeint
from scipy.optimize import fsolve

start=timeit.default_timer()
lam = 1
E=1/12

def V(x,y):
    v = 1/2*(x**2+y**2)+lam*(x**2*y-y**3/3)
    return v

def solve_x0(x0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E
    return equation

def solve_xdot0(xdot0):
    equation = V(x0,y0)+1/2*(xdot0**2+ydot0**2)-E
    return equation


def f(r, t):
    x, xdot, y, ydot = r
    fx = xdot
    fy = ydot
    fxdot = -1 * (x + 2 * lam * x * y)
    fydot = -1 * (y + lam * (x ** 2 - y ** 2))
    return [fx, fxdot, fy, fydot]

def integration(t1,t2,r0):
    t = np.linspace(t1, t2, 3)
    r = odeint(f, r0, t)
    return [r[0,0], r[-1,0], r[0,1], r[-1,1], r[0,2], r[-1,2], r[0,3], r[-1,3], t[0], t[-1]]
#Συνολική Ενέργεια του συστήματος για την οποία θα υπολογίσουμε την τελευταία αρχικη συνθηκη




fig, ax = plt.subplots(1, figsize=(6, 4))
x0 = 0
y0 = 0
vy0 = np.arange(0, 0.45, 0.001)
for ydot0 in vy0:
    xdot0 = fsolve(solve_xdot0, [0.7, -0.7])
    yP = []  # hold data for plotting
    ydotP = []
    for vx0 in xdot0:

        r0 = [x0, vx0, y0, ydot0]
        j = 0
        max_sections = 100
        # Αρχική τιμή χρόνου
        t1 = 0.0

        # Χρονικό διάστημα
        dt = 0.01
        t2 = t1 + dt

        # Αρχικοποίηση Μεταβλητών
        ta = 0
        tb = 0
        x1 = 0
        x2 = 0
        y1 = 0
        y2 = 0
        ydot1 = 0
        ydot2 = 0
        xdot1 = 0
        xdot2 = 0

        while j < max_sections:
            [x1, x2, xdot1, xdot2, y1, y2, ydot1, ydot2, ta, tb] = integration(t1, t2, r0)
            if (x1 < 0 and x2 > 0):
                t_interpl = np.interp(0, [x1, x2], [ta, tb])
                [x1, x2, xdot1, xdot2, y1, y2, ydot1, ydot2, ta, tb] = integration(ta, t_interpl, r0)
                yP.append(y2)
                ydotP.append(ydot2)
                r0 = [x2, xdot2, y2, ydot2]

                t1 = tb
                t2 = tb + dt
                j += 1
                print(tb)
            else:
                r0 = [x2, xdot2, y2, ydot2]
                t1 = tb
                t2 = tb + dt


    ax.scatter(yP, ydotP,s=1.2)


stop=timeit.default_timer()
print(stop-start)

plt.xlabel('y')
plt.ylabel('ydot')
plt.title('Phase Space Plot of Henon-Heiles System')
# Αποθήκευση του διαγράμματος ως εικόνα
plt.savefig('Henon-Heiles Poincare Section_2.png')
plt.savefig("Henon-Heiles Poincare Section_2.pdf")
plt.show()

