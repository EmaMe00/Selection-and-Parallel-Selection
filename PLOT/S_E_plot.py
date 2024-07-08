import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def compute_speedup(t1, tp):
    return [t1/tp[j] for j in range(len(tp))]

def compute_efficiency(t1, tp, p):
    return [t1/(tp[j]*p[j]) for j in range(len(tp))]


nomeFile = 'plotSORT.svg'
p = [1, 2, 3, 4, 5, 6, 7, 8]
t1 = [0.000161, 0.001486, 0.009013, 0.138889, 1.288136,17.746692]
tp = [
    [0.000161,0.000126,0.00011,0.000103,0.000126,0.000121,0.000525,0.001856],
    [0.001486,0.001299,0.001076,0.000857,0.000489,0.000475,0.000817,0.003928],
    [0.009013,0.006663,0.006258,0.005914,0.004684,0.00441,0.008495,0.013135],
    [0.138889,0.073733,0.059555,0.051007,0.043509,0.041047,0.058564,0.070927],
    [1.288136,0.85433,0.599473,0.546694,0.464436,0.43694,0.482981,0.519402],
    [17.746692,8.734381,6.181075,5.899794,4.880016,4.621592,4.430671,4.28057]
]

"""
nomeFile = 'plotNOSORT.svg'
t1 = [0.000365, 0.003564, 0.035986, 0.358161, 3.497114, 33.281954]
p = [1, 2, 3, 4, 5, 6, 7, 8]
tp = [
    [0.000365,0.000298,0.000236,0.000204,0.000194,0.000186,0.000335,0.000898],
    [0.003564,0.002823,0.00209,0.001595,0.001277,0.001168,0.001595,0.002406],
    [0.035986,0.026112,0.020344,0.016481,0.013136,0.011972,0.015744,0.018214],
    [0.358161,0.254193,0.190988,0.151802,0.121768,0.108266,0.114953,0.112982],
    [3.497114,2.513229,1.843699,1.478477,1.193963,1.043614,1.034724,1.036359],
    [33.281954,24.164159,18.67837,15.12063,12.942816,11.499078,11.121658,10.945536]
]
"""


# Calcolo di speedup ed efficienza
speedup = [compute_speedup(t1[i], tp[i]) for i in range(len(t1))]
efficiency = [compute_efficiency(t1[i], tp[i], p) for i in range(len(t1))]
execution_time = tp  # Definiamo il tempo di esecuzione come tp, visto che non è stato definito

# Funzione per interpolare i dati
def smooth_curve(x, y):
    interp_func = interp1d(x, y, kind='cubic')
    x_smooth = np.linspace(min(x), max(x), 300)
    y_smooth = interp_func(x_smooth)
    return x_smooth, y_smooth

# Plot di speedup ed efficienza
plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
# Plot delle curve lisce per speedup
for i in range(len(t1)):
    x_smooth, y_smooth = smooth_curve(p, speedup[i])
    plt.plot(x_smooth, y_smooth, label=f'(Problem size: $10^{{{4+i}}}$)')
# Aggiungi la bisettrice ideale
plt.plot(p, p, linestyle='--', color='k', label='Speedup ideale')
plt.title('Speedup', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Speedup', fontsize='small')
plt.legend(fontsize='small')  # Legenda con font più piccolo
plt.grid(True)

plt.subplot(1, 3, 2)
# Plot delle curve lisce per efficienza
for i in range(len(t1)):
    x_smooth, y_smooth = smooth_curve(p, efficiency[i])
    plt.plot(x_smooth, y_smooth, label=f'(Problem size: $10^{{{4+i}}}$)')
# Aggiungi la retta ideale y=1
plt.axhline(y=1, color='k', linestyle='--', label='Efficienza ideale')
plt.title('Efficienza', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Efficienza', fontsize='small')
plt.legend(fontsize='small')  # Legenda con font più piccolo
plt.grid(True)

# Plot del tempo di esecuzione con curve lisce
plt.subplot(1, 3, 3)
for i in range(len(t1)):
    x_smooth, y_smooth = smooth_curve(p, execution_time[i])
    plt.plot(x_smooth, y_smooth, label=f'(Problem size: $10^{{{4+i}}}$)')
plt.yscale('log')  # Scala logaritmica sull'asse y per il tempo di esecuzione
plt.title('Tempo di esecuzione - Numero di processori', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Tempo di esecuzione', fontsize='small')
plt.legend(fontsize='small')  # Legenda con font più piccolo
plt.grid(True)

# Salva il plot come SVG (immagine vettoriale)
plt.tight_layout()
plt.savefig(nomeFile, format='svg')

# Mostra i plot (opzionale)
plt.show()