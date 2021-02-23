import numpy as np
import matplotlib.pyplot as plt
import isr
import py_plot as pylt

def prob3():
    params = {
            'N_f': int(1e3), 'f_min': -2e6, 'f_max': 2e6, 'f': 430e6, 'n_e': 2e10, 
            'B': 3.5e-5, 'M': 16, 'T_e': 200, 'T_i': 200, 'aspect': 135
            }
    x, y = isr.isr(params)
    y = 10 * np.log10(y)

    plt.figure(figsize=(7,4))
    plt.plot(x, y)
    plt.xlabel('Frequency $f$')
    plt.ylabel('Echo power [dB]')
    plt.grid(alpha=0.4)
    # plt.savefig('gyrolines.pdf')
    plt.show()

def prob5():
    params = {
            'N_f': int(1e3), 'f_min': 3.5e6, 'f_max': 7.5e6, 'f': 933e6, 'n_e': 2e11,
            'B': 5e-5, 'M': 16, 'T_e': 2000, 'T_i': 2000, 'aspect': 180
            }
    
    T_E = [8000, 6000, 4000, 2000]
    data = []
    for T_e in T_E:
        params['T_e'] = T_e
        x, y = isr.isr(params)
        data.append((x, y))
        
    lab = [r'$T_\mathrm{e}=\,$' + f'{t_e}' + r'$\,\mathrm{K}$' for t_e in T_E]
    plt.figure('temps', figsize=(7,4))
    pylt.ridge_plot(data, 'squeeze', 'grid', xlabel='Frequency $f$', \
            ylabel='Echo power [dB]', labels=lab, figname='temps')
    # for d, l in zip(data, lab):
    #     plt.plot(d[0], d[1], label=l)
    # plt.legend()
    plt.savefig('temps.pdf')
    plt.show()

if __name__ == '__main__':
    prob5()
