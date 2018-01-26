import pandas as pd
import matplotlib.pyplot as plt

class Plotter:
    """Plotter Class"""
    def __init__(self, plot_directory, fontsize = 14):
        self.plot_directory = plot_directory
        self.fontsize = fontsize
        self.nfigures = 0
    def bind_to_data(self, sim_name, data_dir):
        self.sim_name = sim_name
        self.data_dir = data_dir
        self.meta_data = pd.read_csv(self.file_dir('meta_data'))
        self.inst_data = pd.read_csv(self.file_dir('inst_data'))
        self.run_summary = pd.read_csv(self.file_dir('run_summary'))
    def file_dir(self, name):
        return self.data_dir + self.sim_name + '_' + name + '.csv'
    def tex_name(self, short_name):
        return self.meta_data.query('short_name==' + '\"' + short_name + '\"')['tex_name'].item()
    def units(self, short_name):
        return self.meta_data.query('short_name==' + '\"' + short_name + '\"')['units'].item()
    def math_mode(self, string):
        return '$' + string + '$'
    def average(self, short_name):
        return "\\langle" + self.tex_name(short_name) + "\\rangle"
    def label(self, x_name, y_name):
        return self.math_mode(self.tex_name(y_name) + '(' + self.tex_name(x_name) + ')')
    def axis_label(self, axis_symbol, axis_unit):
        return self.math_mode('\\frac{' + axis_symbol + '}{' + axis_unit + '}')

class RelaxationTimePlotter(Plotter):
    def __init__(self, plot_directory, fontsize = 14):
        Plotter.__init__(self, plot_directory, fontsize)
        self.plot_title = 'Polymer Relaxation'
        self.xlim = [0.0, 100.0]
    def plot_inst_data(self):
        # plot kinetic energy
        plt.semilogx(self.inst_data['t'], self.inst_data['kep'], label = self.label('t', 'kep'), color = 'red', linestyle = '-')
        # plot potential energy
        plt.semilogx(self.inst_data['t'], self.inst_data['vljp'], label = self.label('t', 'vljp'), color = 'blue', linestyle = '-')
        # plot total energy
        tot_energy_label =  self.label('t', 'vljp') + self.math_mode('+') + self.label('t', 'kep')
        plt.semilogx(self.inst_data['t'], self.inst_data['vljp'] + self.inst_data['kep'], label = tot_energy_label, color = 'k', linestyle = '-')
    def plot_avg_data(self):
        # plot energy averages
        ke_avg_label = self.math_mode(self.average('kep'))
        pe_avg_label = self.math_mode(self.average('vljp'))
        plt.axhline(y = self.run_summary['kep_mean'].iloc[0], label = ke_avg_label, color = 'k', linestyle = '--')
        plt.axhline(y = self.run_summary['vljp_mean'].iloc[0], label = pe_avg_label, color = 'k', linestyle = '-.')
    def make_axes_and_title(self):
        plt.xlabel(self.axis_label(self.tex_name('t'), self.units('t')), fontsize = self.fontsize)
        plt.ylabel(self.axis_label('E', self.units('kep')), fontsize = self.fontsize)
        plt.title(self.plot_title)
    def make_plots(self, sim_name, data_dir):
        self.bind_to_data(sim_name, data_dir)
        self.make_energy_plot()
        self.make_temperature_plot()
    def make_energy_plot(self):
        self.nfigures += 1
        plt.figure(self.nfigures)
        self.plot_inst_data()
        self.plot_avg_data()
        self.make_axes_and_title()
        # configure x-axis range
        plt.xlim(self.xlim)
        plt.legend()
        plt.grid()
        directory = self.plot_directory + sim_name + "_relax_time_energy.png"
        plt.savefig(directory, dpi = 300)
    def make_temperature_plot(self):
        self.nfigures += 1
        plt.figure(self.nfigures)
        # instanteneous temperature
        plt.semilogx(self.inst_data['t'], self.inst_data['temp_kinp'], label = self.label('t', 'temp_kinp'), color = 'red', linestyle = '-')
        # average temperature
        temp_avg_label = self.math_mode(self.average('temp_kinp'))
        plt.axhline(y = self.run_summary['temp_kinp_mean'].iloc[0], label = temp_avg_label, color = 'k', linestyle = '--')
        # axes and title
        plt.xlabel(self.axis_label(self.tex_name('t'), self.units('t')), fontsize = self.fontsize)
        plt.ylabel(self.axis_label('k_{B}T', self.units('temp_kinp')), fontsize = self.fontsize)
        plt.title(self.plot_title)
        # configure x-axis range
        plt.xlim(self.xlim)
        plt.legend()
        plt.grid()
        directory = self.plot_directory + sim_name + "_relax_time_temperature.png"
        plt.savefig(directory, dpi = 300)

if __name__ == '__main__':
    import sys
    if (len(sys.argv) <= 2):
        print "Usage: plot_relaxation_time <sim_name> <data_dir> <image_dir>"
    else:
        sim_name = sys.argv[1]
        data_dir = sys.argv[2]
        image_dir = sys.argv[3]
        plotter = RelaxationTimePlotter(image_dir)
        plotter.make_plots(sim_name, data_dir)
