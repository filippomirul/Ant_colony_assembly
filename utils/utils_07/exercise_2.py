# -*- coding: utf-8 -*-

from random import Random
from time import time
import math
import inspyred

import utils.utils_07.plot_utils

from matplotlib import *
from pylab import *

def readFileAsList(file):
    with open(file) as f:
        lines = f.read().splitlines()
        n = len(lines)
        array = np.empty(n,dtype=np.uint32)
        for i in range(n):
            array[i] = int(lines[i])
        return array
