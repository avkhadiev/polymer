import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import pi
import matplotlib.ticker as mticker

class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)
def colname(name, k):
    return "%s (k=%d)" % (name, k)
def mathenv(tex):
    return "$" +tex + "$"
def search(key, value):
    return key + "==" + "'" + value + "'"
def make_label(name, units):
    return mathenv(name) + ", " + mathenv("\\left[" + units + "\\right]")
def make_title(sim, yl, ell):
    return (sim + ": "
            + mathenv(yl + "\\left(" + "\\tau" + "\\right)") + ", "
            + mathenv("\\ell=" + "%5.3f" % ell))

sim_name = 'test'
sims = ['slerp','shortstep','md']
lbls = ['Slerp','Enhanced Slerp','MD']
fdir = "test/geodesic/"
outdir = "/Users/Arthur/stratt/lab_notebook/"
# xlabel, xunits
xl = "\\tau = \\ell^{(m)}/\\,\\ell^{\\mathrm{f}}"
xu = "1"
xr = [-0.01, 1.01]
fs = 14  # fontsize
observables = ['omega_proj', 'psi', 'delta_theta']
ylabels = ['\\hat{\\Omega} \\cdot \\hat{n}', '\\psi', '\\delta\\theta']
yunits = ['1', '\\mathrm{rad}', '\\mathrm{rad}']
yranges = [[-1.01, 1.01], [-0.2, pi + 0.2], [-pi/130., +pi/100.]]
# indices = [1, 2]                              # only plot first link
colors = ['b', 'r', 'k']
styles = ['--', '--', '-']
for (obs, yl, yu, yr) in zip(observables, ylabels, yunits, yranges):
    fig, ax = plt.subplots()
    # plt.title(make_title(sim, yl, ell), fontsize = fs)
    plt.xlabel(mathenv(xl), fontsize = fs + 2.5)
    plt.ylabel(make_label(yl, yu), fontsize = fs + 2.5)
    plt.xticks(fontsize=fs + 1.5)
    plt.yticks(fontsize=fs + 1.5)
    # plt.grid(True)
    for (sim, c, l, lbl) in zip(sims, colors, styles, lbls):
        data = pd.read_csv(fdir + sim_name + sim + "_inst_data.csv")
        meta = pd.read_csv(fdir + sim_name + sim + "_meta_data.csv")
        lengths = data['ell'].values
        ell = lengths[-1]
        tau = lengths / ell                     # affine parameter
        # select element in array closest to 0.1
        idx = (np.abs(tau-1.0)).argmin()
        plt.ylim(yr)
        if ((sim == 'slerp') and (obs == 'delta_theta')):
            plt.plot(tau[0:idx], np.zeros((tau[0:idx]).shape), color = c, linestyle = l, label = (lbl + " link"))
        else:
            col = colname(obs, 1)
            tex = meta.query(search("short_name", col))['tex_name'].item()  # wut label to use?
            var = data[col].values
            plt.plot(tau[0:idx], var[0:idx], color = c, linestyle = l, label = (lbl + " link"))
    plt.legend(prop={'size': fs})
    plt.tight_layout()
    plt.savefig(outdir + obs + ".png", dpi = 300)
    plt.close(fig)














