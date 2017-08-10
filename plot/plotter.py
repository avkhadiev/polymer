# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from default_settings import default, parameteres
class Plotter:
    def __init__(self, dataframe, config):
        self.dataframe = dataframe
        self.parameters = {}
        self.__configure_parameters()
    def __configure_parameters(self):
        for parameter in parameteres:
            if parameter not in config.keys():
                self.parameters[parameter] = default[parameter]
    def make_axis_name(self):
        axis_name = "$\mathbf{" + self.dataframe["axis_name"].iloc[0] + "}$"
        return axis_name
    # make axis units given the name of the observable
    def make_axis_units(self):
        axis_units = "$(" + self.dataframe["units"].iloc[0] + ")$"
        return axis_units
    # make changes to the canvas according to saved settings
    def __prepare_canvas(self, x_name, y_name):
        # set y- axis label
        y_axis_name = self.make_axis_name()
        y_units = self.make_axis_units()
        x_label = "$\mathbf{t}$"
        y_label = y_axis_name + " " + y_units
        plt.xlabel(x_label, fontsize = self.parameters['fontsize'])
        plt.ylabel(y_label, fontsize = self.parameters['fontsize'])
    def make_canvas(self, x_name, y_name):
        plt.figure(1)
        self.__prepare_canvas(x_name, y_name)
    # looks up xmin and xmax;
    # finds the biggest and the smallest values of y between xmin and xmax
    # calculates the range based on the span of y
    def adjust_axes_range(x, y):
        
    def plot(self, x, y):
        plt.plot(self.dataframe[x], self.dataframe[y], 'ro', markersize= self.parameters['markersize'])
    def show_plot(self):
        plt.show()
