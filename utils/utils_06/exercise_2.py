# -*- coding: utf-8 -*-

from pylab import *
import sys
from inspyred import ec
import plot_utils

from inspyred import benchmarks

import ga, es, pso
from inspyred_utils import NumpyRandomWrapper

"""
-------------------------------------------------------------------------
Edit this part to do the exercises
"""

num_vars = 100 # Number of dimensions of the search space

args = {}

# the problem class
args["problem_class"] = benchmarks.Sphere

# other problems to try,
# see  https://pythonhosted.org/inspyred/reference.html#module-inspyred.benchmarks

# unimodal
#benchmarks.Sphere
#benchmarks.Rosenbrock

# multimodal
#benchmarks.Griewank
#benchmarks.Ackley
#benchmarks.Rastrigin
#benchmarks.Schwefel

# common parameters
args["max_generations"] = 100 # Number of generations
args["pop_size"] = 50 # population size

# parameters for the GA
args["gaussian_stdev"] = 0.5 # Standard deviation of the Gaussian mutations
args["mutation_rate"] = 0.5 # fraction of loci to perform mutation on
args["tournament_size"] = 2
args["num_elites"] = 1 # number of elite individuals to maintain in each gen

# parameters for the ES
args["num_offspring"] = 100 #lambda
args["sigma"] = 1.0 # default standard deviation
args["strategy_mode"] = es.INDIVIDUAL # es.GLOBAL, es.INDIVIDUAL, None
args["mixing_number"] = 1 #rho

# parameters for the PSO
args["topology"] = pso.RING #pso.RING, pso.STAR
args["neighborhood_size"] = 5   #used only for the ring topology
args["inertia"] = 0.5
args["cognitive_rate"] = 2.1
args["social_rate"] = 2.1

"""
-------------------------------------------------------------------------
"""

display = True # Plot initial and final populations

if __name__ == "__main__":
    
    if len(sys.argv) > 1 :
        rng = NumpyRandomWrapper(int(sys.argv[1]))
    else :
        rng = NumpyRandomWrapper()

    # Run GA
    args["fig_title"] = "GA"
    best_individual, best_fitness, final_pop = ga.run_ga(rng,num_vars=num_vars,
                                              display=display,use_log_scale=True,
                                              **args)
    print("Best GA fitness:", best_fitness)

    # Run ES
    args["fig_title"] = "ES"
    best_individual, best_fitness, final_pop = es.run_es(rng,num_vars=num_vars,
                                                display=display,use_log_scale=True,
                                                **args)
    print("Best ES fitness:", best_fitness)

    # Run PSO
    args["fig_title"] = "PSO"
    best_individual, best_fitness, final_pop = pso.run_pso(rng,num_vars=num_vars,
                                                display=display,use_log_scale=True,
                                                **args)
    print("Best PSO fitness:", best_fitness)

    ioff()
    show()
