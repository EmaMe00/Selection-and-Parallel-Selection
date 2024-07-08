import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def compute_speedup(t1, tp):
    return [t1/tp[j] for j in range(len(tp))]

def compute_efficiency(t1, tp, p):
    return [t1/(tp[j]*p[j]) for j in range(len(tp))]

nomeFile = 'plotMIX.svg'
# Definizione dei dati
p = [1, 2, 3, 4, 5, 6, 7, 8]

t1_Sort = [0.000161, 0.001486, 0.009013, 0.138889, 1.288136,17.746692]
tp_Sort = [
    [0.000161,0.000126,0.00011,0.000103,0.000126,0.000121,0.000525,0.001856],
    [0.001486,0.001299,0.001076,0.000857,0.000489,0.000475,0.000817,0.003928],
    [0.009013,0.006663,0.006258,0.005914,0.004684,0.00441,0.008495,0.013135],
    [0.138889,0.073733,0.059555,0.051007,0.043509,0.041047,0.058564,0.070927],
    [1.288136,0.85433,0.599473,0.546694,0.464436,0.43694,0.482981,0.519402],
    [17.746692,8.734381,6.181075,5.899794,4.880016,4.621592,4.430671,4.28057]
]

t1_noSort = [0.000365, 0.003564, 0.035986, 0.358161, 3.497114, 33.281954]
tp_noSort = [
    [0.000365,0.000298,0.000236,0.000204,0.000194,0.000186,0.000335,0.000898],
    [0.003564,0.002823,0.00209,0.001595,0.001277,0.001168,0.001595,0.002406],
    [0.035986,0.026112,0.020344,0.016481,0.013136,0.011972,0.015744,0.018214],
    [0.358161,0.254193,0.190988,0.151802,0.121768,0.108266,0.114953,0.112982],
    [3.497114,2.513229,1.843699,1.478477,1.193963,1.043614,1.034724,1.036359],
    [33.281954,24.164159,18.67837,15.12063,12.942816,11.499078,11.121658,10.945536]
]

# Calcolo di speedup ed efficienza
speedup_noSort = [compute_speedup(t1_noSort[i], tp_noSort[i]) for i in range(len(t1_noSort))]
efficiency_noSort = [compute_efficiency(t1_noSort[i], tp_noSort[i], p) for i in range(len(t1_noSort))]

speedup_Sort = [compute_speedup(t1_Sort[i], tp_Sort[i]) for i in range(len(t1_Sort))]
efficiency_Sort = [compute_efficiency(t1_Sort[i], tp_Sort[i], p) for i in range(len(t1_Sort))]

# Funzione per interpolare i dati
def smooth_curve(x, y):
    interp_func = interp1d(x, y, kind='cubic')
    x_smooth = np.linspace(min(x), max(x), 300)
    y_smooth = interp_func(x_smooth)
    return x_smooth, y_smooth

# Plot di speedup ed efficienza
plt.figure(figsize=(18, 6))

plt.subplot(1, 3, 1)
colors = ['b', 'g', 'r', 'c', 'm', 'y']
for i in range(len(t1_noSort)):
    x_smooth, y_smooth = smooth_curve(p, speedup_noSort[i])
    plt.plot(x_smooth, y_smooth, linestyle='-', color=colors[i])
    x_smooth, y_smooth = smooth_curve(p, speedup_Sort[i])
    plt.plot(x_smooth, y_smooth, linestyle='--', color=colors[i])
plt.plot(p, p, linestyle='--', color='k', label='Speedup ideale')
plt.title('Speedup', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Speedup', fontsize='small')
plt.legend([plt.Line2D([0], [0], color='k', linestyle='-'),
            plt.Line2D([0], [0], color='k', linestyle='--')],
           ['NoSort', 'Sort'], fontsize='small')
plt.grid(True)

plt.subplot(1, 3, 2)
for i in range(len(t1_noSort)):
    x_smooth, y_smooth = smooth_curve(p, efficiency_noSort[i])
    plt.plot(x_smooth, y_smooth, linestyle='-', color=colors[i])
    x_smooth, y_smooth = smooth_curve(p, efficiency_Sort[i])
    plt.plot(x_smooth, y_smooth, linestyle='--', color=colors[i])
plt.axhline(y=1, color='k', linestyle='--', label='Efficienza ideale')
plt.title('Efficienza', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Efficienza', fontsize='small')
plt.legend([plt.Line2D([0], [0], color='k', linestyle='-'),
            plt.Line2D([0], [0], color='k', linestyle='--')],
           ['NoSort', 'Sort'], fontsize='small')
plt.grid(True)

plt.subplot(1, 3, 3)
for i in range(len(t1_noSort)):
    x_smooth, y_smooth = smooth_curve(p, tp_noSort[i])
    plt.plot(x_smooth, y_smooth, linestyle='-', color=colors[i])
    x_smooth, y_smooth = smooth_curve(p, tp_Sort[i])
    plt.plot(x_smooth, y_smooth, linestyle='--', color=colors[i])
plt.yscale('log')
plt.title('Tempo di esecuzione - Numero di processori', fontsize='large')
plt.xlabel('Numero di processori', fontsize='small')
plt.ylabel('Tempo di esecuzione', fontsize='small')
plt.legend([plt.Line2D([0], [0], color='k', linestyle='-'),
            plt.Line2D([0], [0], color='k', linestyle='--')],
           ['NoSort', 'Sort'], fontsize='small')
plt.grid(True)

# Salva il plot come SVG (immagine vettoriale)
plt.tight_layout()
plt.savefig(nomeFile, format='svg')

# Mostra i plot (opzionale)
plt.show()
