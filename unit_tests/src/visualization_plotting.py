#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class VisualizationPlotting:
    def __init__(self, yaml_file: str, results_dict: dict, visualization_df: pd.DataFrame, observable_df: pd.DataFrame, measurement_df: pd.DataFrame):
        """This class is designed to plot and save the visualization data.
        input:
            yaml_file: str - path to the YAML file
            result_dict: dict - dictionary containing the simulation results
            visualization_df: str - path to the visualization dataframe
            observable_df: str - path to the observable dataframe
            model: str - path to the SBML model
            measurement_df: str - path to the measurement dataframe
        """
        self.yaml_file = yaml_file
        self.results_dict = results_dict
        self.visualization_df = visualization_df
        self.observable_df = observable_df
        self.measurement_df = measurement_df

    def dynamic_plot(self):
        """Parses the information within the visualization dataframe and 
        plots simulation data accordingly. Information within the 
        visualization dataframe must conform to PEtab standards.

        output: plots the simulation data according to the visualization dataframe
        """

        results_dict = self.results_dict
        visualization_df = self.visualization_df

        unique_plots = visualization_df['plotId'].unique()

        num_rows, num_cols, plot_position_matrix = generate_position_matrix(len(unique_plots))

        unique_plt_positions = {}

        for i, plotId in enumerate(unique_plots):
            # Find the position from the matrix
            row, col = np.where(plot_position_matrix == i)
            row, col = row[0], col[0]
            unique_plt_positions[plotId] = (row, col)

        # # Create a subplot figure with the determined grid size
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 2 * num_rows))

        # Ensure axes is always a 2D array
        if num_rows == 1 and num_cols == 1:
            axes = np.array([[axes]])
        elif num_rows == 1 or num_cols == 1:
            axes = axes.reshape((max(num_rows, num_cols),1))

        for i, plotId in enumerate(visualization_df['plotId']):
            plot_type = visualization_df['plotTypeSimulation'][i]
            yValue = visualization_df['yValues'][i]
            xValue = visualization_df['xValues'][i]

            # Find the position from the matrix
            row, col = unique_plt_positions[plotId]

            if plot_type == 'ScatterPlot':
                condition = visualization_df['datasetId'][i]
                for cell in results_dict[condition]:
                    yvals = results_dict[condition][cell][yValue]['xoutS']
                    if xValue == 'time':
                        xvals = results_dict[condition][cell][yValue]['toutS']
                    else:
                        xvals = results_dict[condition][cell][xValue]['xoutS']
                    axes[row, col].scatter(yvals, xvals)

            elif plot_type == 'LinePlot':
                condition = visualization_df['datasetId'][i]
                for cell in results_dict[condition]:
                    yvals = results_dict[condition][cell][yValue]['xoutS']
                    if xValue == 'time':
                        xvals = results_dict[condition][cell][yValue]['toutS']
                    else:
                        xvals = results_dict[condition][cell][xValue]['xoutS']
                    axes[row, col].plot(yvals, xvals)

            elif plot_type == 'BarPlot':
                condition = visualization_df['datasetId'][i]
                for cell in results_dict[condition]:
                    yvals = results_dict[condition][cell][yValue]['xoutS']
                    if xValue == 'time':
                        xvals = results_dict[condition][cell][yValue]['toutS']
                    else:
                        xvals = results_dict[condition][cell][xValue]['xoutS']
                    axes[row, col].bar(yvals, xvals)


        # Remove unused subplots
        for i in range(len(unique_plots), num_rows * num_cols):
            row, col = np.where(plot_position_matrix == i)
            row, col = row[0], col[0]
                
            fig.delaxes(axes[row, col])
        # # Adjust layout
        plt.tight_layout()

        return fig

@staticmethod
def generate_position_matrix(num_plots):
    """Generates a position matrix based on the number of plots to be generated
    input:
        num_plots: int - number of plots to be generated

    output:
        num_rows: int - number of rows in the plot grid
        num_cols: int - number of columns in the plot grid
        position_matrix: np.array - matrix containing the position of each plot in the grid"""

    side_length = int(np.ceil(np.sqrt(num_plots)))
    position_matrix = np.arange(side_length ** 2).reshape((side_length, side_length)).T
    num_rows, num_cols = position_matrix.shape
    return num_rows, num_cols, position_matrix
    

