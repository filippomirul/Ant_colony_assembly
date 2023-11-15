# -*- coding: utf-8 -*-

from random import Random
from time import time
import math
import inspyred

import utils.utils_07.plot_utils

from matplotlib import *
from pylab import *

from matplotlib import collections  as mc

def readFileAsMatrix(file):
    with open(file) as f:
        lines = f.read().splitlines()
        matrix = []
        for line in lines:
            row = []
            for value in line.split():
                row.append(float(value))
            matrix.append(row)
        return matrix

def plotSolution(points, distances, solution, title):
    fig = figure(title)
    ax = fig.add_subplot(111)
    ax.scatter(*zip(*points))
    
    for i,p in enumerate(points):
        ax.annotate(str(i), p)
    
    # draw all possible path segments
    lines = []
    for i, p in enumerate(points):
        for j, q in enumerate(points):
            if distances[i][j] > 0 and i > j:
                lines.append((points[i], points[j]))
    lc = mc.LineCollection(lines, linewidths=.1)
    ax.add_collection(lc)

    # draw the solution
    lines = []
    for i in arange(len(solution)-1):
        lines.append((points[solution[i]], points[solution[i+1]]))
    lines.append((points[solution[0]], points[solution[-1]]))
    lc = mc.LineCollection(lines, linewidths=1, color="r")
    ax.add_collection(lc)

    #ax.set_title(title)
    ax.autoscale()
    ax.margins(0.1)
