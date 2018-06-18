import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import pi
import matplotlib.ticker as mticker

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

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

def mathenv(tex):
    return "$" +tex + "$"
def search(key, value):
    return key + "==" + "'" + value + "'"
def make_label(label, coordinate, units):
    return (mathenv(label + "_{" + coordinate + "}") + ", "
        + mathenv("\\left[" + units + "\\right]"))
def make_name(observable, atom, coordinate):
    return observable + "_" + str(atom) + "_" + str(coordinate);

sim_names = ['test', 'test', 'test']
sims = ['md', 'shove', 'slerp']#, 'md']
sim_lbls = ['MD', 'Plerp', 'Slerp']
styles = ['-', '--', ':']                       # correspond to algorithm
fdir = "test/geodesic/"
outdir = "/Users/Arthur/stratt/lab_notebook/"
# observable to plot, which atoms, and how to plot them
observable = 'rpol_atom_1'
atoms = ['1', '7', '11']
colors = ['b', 'r', 'k']                        # correspond to atom
# unit, label, ranges, and other settings
unit = '\\mathrm{\\sigma}'
label = '\\mathrm{r}'
yr = (-1.05, 1.05)
xr = (-1.05, 1.05)
xcoor = 'y'
ycoor = 'z'
ms = 'o'                        # markerstyle for endpoints
fs = 14                         # fontsize
# -----------------------------------------------------------------------------
#                                      PLOT
# -----------------------------------------------------------------------------
# setup canvas
fig, ax = plt.subplots()
plt.xlabel(make_label(label, xcoor, unit), fontsize = fs + 2.5)
plt.ylabel(make_label(label, ycoor, unit), fontsize = fs + 2.5)
plt.xticks(fontsize = fs + 1.5)
plt.yticks(fontsize = fs + 1.5)
for sim_name, sim, sim_lbl, a_style in zip(sim_names, sims, sim_lbls, styles):
        # read in the data
        data = pd.read_csv(fdir + sim_name + sim + "_inst_data.csv")
        # determine index to select the right portion of the trajectories
        lengths = data['ell'].values
        ell = lengths[-1]
        tau = lengths / ell
        for atom, a_color in zip(atoms, colors):
            # get coordinates to plot
            x_val = data[make_name(observable, atom, xcoor)].values
            y_val = data[make_name(observable, atom, ycoor)].values
            # make a legend label
            a_label = (sim_lbl + " "
                + mathenv(label + "^{" + atom + "}_{" + xcoor + ycoor + "}"))
            # make the plot, mark endpoints; get the handle to mark the arrow
            line = plt.plot(x_val, y_val,
                ms, color = a_color, linestyle = a_style, label = a_label,
                markevery=[0,-1])[0]
            add_arrow(line)
# make final brushes like placing the legend, formatting the layout, etc.
# plt.xlim(xr)
# plt.ylim(yr)
plt.tight_layout()
# Shrink current axis
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
# Put a legend to the right of the current axis
plt.legend(prop={'size': fs - 2.0}, frameon = False,
    loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(outdir + observable + "_" + xcoor + ycoor + ".png", dpi = 300)
plt.close(fig)
