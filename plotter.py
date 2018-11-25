import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import math 

sns.set_style('darkgrid')

probability = np.loadtxt('results/probability.txt')

x = np.linspace(0,1*math.pi,100)

plt.plot(x,np.sin(x),label='Integrating function')
plt.plot(x,probability,'o',c='orange',label='p(x)')
plt.xlabel('x axis',fontsize=40)
plt.ylabel('y axis',fontsize=40)
plt.legend(loc='best',fontsize=30)
plt.tick_params(axis='both', labelsize=20)
plt.show()
