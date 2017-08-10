#pb.figure( self.figures )
#pb.plot(self.sampleTArray / self.T , self.tempArray / self.tempMax,
#        linestyle = '-', marker = 'o', color = 'r')
#pb.axhline(y = self.meanTemp() / self.tempMax,
#        linestyle = '--',  color='b')
#pb.axhline(y = self.iniTemp,
#        linestyle = '--',  color='k')
#pb.title("Instanteneous temperature (red) and Mean Temperature (blue)")
#pb.xlabel("t/T")
#pb.ylabel("Inst Temp / Max Temp")
#pb.grid(True)
import math
import pandas as pd
import matplotlib.pyplot as plt

class Plotter:
    nfigures = 0
    fontsize = 14
    markersize = 2
    def plot(self, dataframe, x_name, y_name):
        self.nfigures += 1
        plt.figure(self.nfigures)
        plt.plot(dataframe[x_name], dataframe[y_name], 'ro', markersize = self.markersize)
        #plt.axvline(x = pow(2.0, 1.0/6.0), linestyle = '--', color='b')
        plt.xlim((0.5, 2.0))
        plt.ylim((-4.0, 4.0))
        x_axisname = "$\mathbf{" + dataframe["axis_name_" + x_name].iloc[0] + "}$"
        y_axisname = "$\mathbf{" + dataframe["axis_name_" + y_name].iloc[0] + "}$"
        x_units = "$(" + dataframe["units_" + x_name].iloc[0] + ")$"
        y_units = "$(" + dataframe["units_" + y_name].iloc[0] + ")$"
        x_label = x_axisname + " " + x_units
        y_label = y_axisname + " " + y_units
        plt.xlabel(x_label, fontsize = self.fontsize)
        plt.ylabel(y_label, fontsize = self.fontsize)
        plt.grid(True)

if __name__ == '__main__':
    x_name = "distance"
    y_names = ["potential", "force", "virial"]
    plotter = Plotter()
    for y_name in y_names:
        x = pd.read_csv("test/ljpot_" + x_name + ".csv")
        y = pd.read_csv("test/ljpot_" + y_name + ".csv")
        xsuffix = "_" + x_name
        ysuffix = "_" + y_name
        xy = x.iloc[:,1:].join(y.iloc[:,1:], lsuffix=xsuffix, rsuffix=ysuffix)
        plt.rc('text', usetex=True)
        #plt.rc('font', family='serif')
        plotter.plot(xy, x_name, y_name)
    plt.show()
