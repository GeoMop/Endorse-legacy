import math
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from decimal import *


def take_closest(myList, myNumber):
    """
    Always returns last known value
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return (pos - 1)
    else:
        return (pos - 1)


# input section

v = 10  # actual velocity
d = 1e-1  # Diffusion coeeficient

rho = 2700  # rock density
kd = 0  # distribution coefficient of linear sorption
poro = 1  # porosity
direct_retardation = True  # retardation may be either computed or input directly
r = 1.0

lam = 0  # decay rate

l = 1000  # transport path length
t = 250  # simulation period
dt_nr = 2000  # number of time step (higher is better but slower)

file_in = 'conc.csv'  # file with transport BC
file_out = 'result.csv'  # results file

# input section end

# compute constants
if not direct_retardation:
    r = (1+kd*rho*(1-poro)/poro)
a = r/(2*math.sqrt(d*r))
b = v/(2*math.sqrt(d*r))
c = v/d
print(a, b, c)

if t < a*l/b:
    print('Simulation period is probably too short.')
    print('Consider prolonging it to value > ' + str(a*l/b))

# load input file
times = []
c_BCs = []
with open(file_in, 'r') as f_in:
    for line in f_in:
        line = line.rstrip().split(';')
        times.append(float(line[0]))
        c_BCs.append(float(line[1]))

# computation
c_BCs_discrete = []
c_BCs_discrete_decay = []
g_discrete = []
t_discrete = []
dt = t/dt_nr
for i in range(1, dt_nr+1):
    t_discrete.append(i*dt)
    c_BCs_discrete.append(c_BCs[take_closest(times, i*dt)])
    c_BCs_discrete_decay.append(c_BCs[take_closest(times, i*dt)]*math.exp(lam*i*dt))
    e_1 = -math.pow((a*l-b*(i*dt))/math.sqrt(i*dt), 2)
    e_2 = c*l-math.pow((a*l+b*(i*dt))/math.sqrt(i*dt), 2)
    exp_1 = Decimal(math.exp(1))**Decimal(e_1)
    exp_2 = Decimal(math.exp(1))**Decimal(e_2)
    mult_1 = a*l*math.pow(i*dt, -1.5) + b*math.pow(i*dt, -0.5)
    mult_2 = a*l*math.pow(i*dt, -1.5) - b*math.pow(i*dt, -0.5)
    g_discrete.append((1/(2*math.sqrt(math.pi))) * (float(exp_1) * mult_1 + float(exp_2) * mult_2))
    g_x_dt = [dt*j for j in g_discrete]
    print(i, exp_1, exp_2)

c_out = np.convolve(c_BCs_discrete, g_x_dt, mode='full')
c_out_decay = np.convolve(c_BCs_discrete_decay, g_x_dt, mode='full')

t_2 = []
c_out_RHS = []

with open(file_out, 'w') as f_out:
    f_out.write('Time; Concentration\n')
    for i in range(1, len(list(c_out))+1):
        t_2.append(i*dt)
        c_out_RHS.append(list(c_out_decay)[i-1]*math.exp(-i*dt*lam))
        f_out.write(str(i*dt) + ';' + str(list(c_out_decay)[i-1]*math.exp(-i*dt*lam)) + '\n')

#plt.plot(t_discrete, c_BCs_discrete, 'r--', t_2, list(c_out), 'g--', t_discrete, c_BCs_discrete_decay, 'b--')
#plt.plot(t_discrete, c_BCs_discrete, 'r--', t_2, list(c_out), 'g--', t_2, list(c_out_decay), 'b--')
plt.plot(t_discrete, c_BCs_discrete, 'r--', t_2, list(c_out), 'g--', t_2, c_out_RHS, 'b--')
plt.show()

# print(list(c_out))

print(sum(g_x_dt))
