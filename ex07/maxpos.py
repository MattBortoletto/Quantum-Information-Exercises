import numpy as np
import matplotlib.pyplot as plt

t_len = 1000
t = np.arange(0, t_len, step=1)
prob = np.loadtxt('prob_time_evol_0.10000.txt')

dx = 0.02
mean = np.zeros(t_len)
for i in range(t_len):
    mean[i] = sum(prob[:, i+1]*prob[:,0]*dx)


plt.plot(t, mean) 
plt.xlabel('t')
plt.grid()
plt.ylabel('Mean position')
plt.savefig('mean_evol.png')
plt.show()