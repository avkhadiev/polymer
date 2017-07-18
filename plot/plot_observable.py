# -*- coding: utf-8 -*-
# option parser specifies the observable to plot and customizes the plot
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c", "config",
        action="store", type="string", dest="config_string",
        help="specify configuration via ':'-separated string")
(options, args) = parser.parse_args()
# TODO removes capitalization, replaces ' ' with '_' and '.' with '_'
def parse_name(string):
    return string
def parse_config(config_string)
    config = dict(item.split("=") for item in config_string.split(":"))
    return config
def make_file_name(simulation, observable):
    simulation = parse_name(simulation)
    observable = parse_name(observable)
    csv_file = "test/" + simulation + "_" + observable + ".csv"
    return csv_file
# if launched directly
if __name__ == '__main__':
    import sys
    import pandas as pd
    from plotter import Plotter
    if (len(sys.argv) <= 2):
        print "Usage: plot_observable <sim_name> <observable_name> <options>"
    else:
        simulation = sys.argv[1]
        observable = sys.argv[2]
        csv_file = make_file_name(simulation, observable)
        dataframe = pd.read_csv(csv_file)
        config = simulate.parse(options.config)
        plotter = Plotter(dataframe, config)
        plotter.make_canvas("time", observable)
        plotter.plot("time", observable)
        plotter.show_plot()
