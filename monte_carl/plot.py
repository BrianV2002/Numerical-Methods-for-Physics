import numpy as np
import matplotlib.pyplot as plt



initial_state = np.loadtxt('S_inicial.dat')
final_state = np.loadtxt('S_final.dat')


fig, ax = plt.subplots(1, 2, figsize=(10, 5))


ax[0].imshow(initial_state, cmap='Greys', interpolation='nearest')
ax[0].set_title("Estado Inicial")


ax[1].imshow(final_state, cmap='Greys', interpolation='nearest')
ax[1].set_title("Estado final")

plt.show()
