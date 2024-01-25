from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.path import Path

import matplotlib.patches as patches
import numpy as np
import pylab

from math import cos, sin, atan, pi
from utils.utils_08.ann_benchmarks import MAX_WEIGHT

class ANNPlotter(object) :
    """
    Plot a neural network, positions determination based on code from:
    https://stackoverflow.com/questions/29888233/how-to-visualize-a-neural-network
    """
    
    def __init__(self, net,  axes=None, 
                 vertical_distance_between_layers=10,
                 horizontal_distance_between_neurons=10, neuron_radius=0.9):
        self.net = net
        self.axes = axes
        if self.axes is None :
            self.axes = plt.figure('ANN Plotter').gca()
            
        self.vertical_distance_between_layers = vertical_distance_between_layers
        self.horizontal_distance_between_neurons =\
                horizontal_distance_between_neurons
        self.neuron_radius = neuron_radius
        
        self.number_of_neurons_in_widest_layer = np.max(
                [net.num_inputs,net.num_hidden,net.num_outputs])
        
        # now calculate positions of all neurons
        self.neuron_positions = [[]]
        x = self.__calculate_left_margin_so_layer_is_centered(
                                self.net.num_inputs)
        y = 0.
        for _ in range(self.net.num_inputs) :
            self.neuron_positions[-1].append((x,y))
            x += self.horizontal_distance_between_neurons
        
        if self.net.num_hidden > 0 :
            self.neuron_positions.append([])
            y += self.vertical_distance_between_layers
            x = self.__calculate_left_margin_so_layer_is_centered(
                                self.net.num_hidden)
            for _ in range(self.net.num_hidden) :
                self.neuron_positions[-1].append((x,y))
                x += horizontal_distance_between_neurons
        
        self.neuron_positions.append([])
        y += self.vertical_distance_between_layers
        x = self.__calculate_left_margin_so_layer_is_centered(
                            self.net.num_outputs)
        for _ in range(self.net.num_outputs) :
            self.neuron_positions[-1].append((x,y))
            x += horizontal_distance_between_neurons
        
        norm = colors.Normalize(vmin=-MAX_WEIGHT, vmax=MAX_WEIGHT)
        self.scalar_map = pylab.cm.ScalarMappable(norm=norm, cmap="jet")
        self.weights = []
    
    def __calculate_left_margin_so_layer_is_centered(self, number_of_neurons):
        return (self.horizontal_distance_between_neurons * 
                (self.number_of_neurons_in_widest_layer - number_of_neurons)
                ) / 2.0
    
    def __draw_synapse(self, neuron1, neuron2, weight, recurrent=False):
        color=self.scalar_map.to_rgba(weight)
        self.weights.append(weight)
        
        if not recurrent :
        
            angle = atan((neuron2[0] - neuron1[0]) / 
                         float(neuron2[1] - neuron1[1]))
            
            x_adjustment = self.neuron_radius * sin(angle)
            y_adjustment = self.neuron_radius * cos(angle)
            
            self.axes.arrow(neuron1[0], neuron1[1], 
                            neuron2[0] - x_adjustment - neuron1[0], 
                            neuron2[1] - y_adjustment - neuron1[1], 
                            head_width=0.5, head_length=0.5, fc=color, 
                            ec=color,length_includes_head=True,zorder=1)
            
            # show weight 1/3 of way up connection line
            
            label_x = neuron1[0] + (neuron2[0] - neuron1[0]) / 3.
            label_y = neuron1[1] + (neuron2[1] - neuron1[1]) / 3.
            
        else :
            
            if neuron2[1] == neuron1[1] : # same layer
                height_factor = 1.5
                width_factor = 2.
                if neuron1[0] > neuron2[0] :
                    height_factor = 4.
                    width_factor = 4.
                verts = [
                    neuron1,  
                    (neuron1[0] - self.horizontal_distance_between_neurons/width_factor, 
                     neuron1[1] + self.vertical_distance_between_layers/height_factor),
                    (neuron2[0] + self.horizontal_distance_between_neurons/width_factor, 
                     neuron2[1] + self.vertical_distance_between_layers/height_factor),
                    neuron2]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            else : # going to prior layer (actually doesn't exist right now)
                factor = -2.0
                if neuron2[0] > neuron1[0] : # left to right
                    factor = 2.0
                verts = [
                        neuron1,  
                        (neuron1[0] + self.horizontal_distance_between_neurons * factor, 
                         (neuron1[1] + neuron2[1])/2.0),
                        neuron2]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4]
            
            path = Path(verts, codes)
            
            patch = patches.PathPatch(path, ec=color, fc='none', zorder=1, aa=True)
            #patch = patches.FancyArrowPatch(path=path, ec=color, fc=color, zorder=1, aa=True, arrowstyle="-|>", mutation_scale=150**.5)
            patch_verts = patch.get_verts()
            self.axes.add_patch(patch)
            
            if neuron1 == neuron2 :
                vert = patch_verts[int(len(patch_verts)/2)]
            else :
                vert = patch_verts[int(len(patch_verts)/3)]
            
            label_x = vert[0]
            label_y = vert[1]  
            
            if patch_verts[-2][1] != neuron2[1] :
                angle = atan((patch_verts[-2][0] - neuron2[0]) /
                     float(patch_verts[-2][1] - neuron2[1]))
                dx = 0.5 * sin(angle)
                dy = 0.5 * cos(angle)
                x_adjustment = self.neuron_radius * sin(angle)
                y_adjustment = self.neuron_radius * cos(angle)
                
                self.axes.arrow(neuron2[0] + x_adjustment + dx,
                            neuron2[1] + y_adjustment + dy, -dx, -dy,
                            head_width=0.5, head_length=0.5, fc=color,
                            ec=color,length_includes_head=True,zorder=1,
                            aa=True)
                
            else:
                if neuron1[0]>neuron2[0]:
                    dx = 0.5
                else:
                    dx = -0.5
                
                self.axes.arrow((neuron1[0]+neuron2[0])/2-0.1,
                            patch_verts[int(len(patch_verts)/2-1)][1]+0.05, dx, 0,
                            head_width=0.5, head_length=0.5, fc=color,
                            ec=color,length_includes_head=True,zorder=1,
                            aa=True)
                
        # use default box of white with white border
        bbox = dict(boxstyle='round',
                    ec=(1.0, 1.0, 1.0),
                    fc=(1.0, 1.0, 1.0),
                    )
        self.axes.text(label_x, label_y,
                       "{:.2f}".format(weight),
                       color='black', size='x-small', 
                       verticalalignment='center',
                       horizontalalignment='center',bbox=bbox,zorder=2)
    
    def __draw_bias(self, layer):
        position = (self.number_of_neurons_in_widest_layer * 
                    self.horizontal_distance_between_neurons,
                    self.neuron_positions[layer][0][1] + 
                            self.vertical_distance_between_layers/2.)
        text = self.axes.text(position[0], position[1],
                        "Bias", color='black',horizontalalignment='center',
                        size='large')
        text.set_bbox(dict(edgecolor='black',facecolor='lightgray'))
        text.set_zorder(2)
        return position       
        
    def draw(self):
        for layer in self.neuron_positions :
            for neuron_position in layer :
                circle = plt.Circle(neuron_position, radius=self.neuron_radius,
                                    ec='black', fc='lightgray',aa=True)
                circle.set_zorder(2)
                self.axes.add_patch(circle)
        
        bias_positions = []
        bias_positions.append(self.__draw_bias(0))
        if (self.net.num_hidden > 0) :
            bias_positions.append(self.__draw_bias(1))        
        
        if self.net.num_hidden == 0 :
            for j in range(self.net.num_outputs) :
                for i in range(self.net.num_inputs) :
                    self.__draw_synapse(self.neuron_positions[0][i],
                                        self.neuron_positions[1][j],
                                        self.net.weights[0][i][j])
                self.__draw_synapse(bias_positions[0], 
                                    self.neuron_positions[1][j], 
                                    self.net.weights[0][-1][j])
            if self.net.recurrent :
                # for no hidden units, but recurrent, we feedback outputs to self
                for j in range(self.net.num_outputs) :
                    for i in range(self.net.num_outputs) :
                        self.__draw_synapse(self.neuron_positions[1][i],
                                            self.neuron_positions[1][j],
                                            self.net.weights[0][self.net.num_inputs + i][j], True)
        else :
            for j in range(self.net.num_hidden) :
                for i in range(self.net.num_inputs) :
                    self.__draw_synapse(self.neuron_positions[0][i],
                                        self.neuron_positions[1][j],
                                        self.net.weights[0][i][j]) 
                self.__draw_synapse(bias_positions[0], 
                                    self.neuron_positions[1][j], 
                                    self.net.weights[0][-1][j])
            for j in range(self.net.num_outputs) :  
                for i in range(self.net.num_hidden) :
                    self.__draw_synapse(self.neuron_positions[1][i],
                                        self.neuron_positions[2][j],
                                        self.net.weights[1][i][j])
                self.__draw_synapse(bias_positions[1], 
                                self.neuron_positions[2][j], 
                                self.net.weights[1][-1][j])
                
            if self.net.recurrent :
                # for hidden + recurrent, we feedback hidden to hidden
                for j in range(self.net.num_hidden) :
                    for i in range(self.net.num_hidden) :
                        self.__draw_synapse(self.neuron_positions[1][i],
                                            self.neuron_positions[1][j],
                                            self.net.weights[0][
                                                        self.net.num_inputs +
                                                        i][j], True)
        
        self.axes.relim()
        self.axes.autoscale(tight=False)
        xlim = self.axes.get_xlim()
        self.axes.set_xlim(xlim[0] - self.neuron_radius,
                           xlim[1] + self.horizontal_distance_between_neurons + self.neuron_radius)
        self.axes.set_aspect('equal')
        # turn off all x axis ticks        
        self.axes.get_xaxis().set_ticks([])
        self.axes.get_yaxis().set_ticks([layer[0][1] 
                                         for layer in self.neuron_positions])
        if self.net.num_hidden > 0 :
            self.axes.get_yaxis().set_ticklabels(["Input", "Hidden", "Output"])
        else :
            self.axes.get_yaxis().set_ticklabels(["Input", "Output"])
        
        self.scalar_map.set_array(np.asarray(self.weights))
        plt.colorbar(self.scalar_map, ax=self.axes)
