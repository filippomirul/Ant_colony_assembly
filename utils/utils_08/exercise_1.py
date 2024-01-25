from pylab import *
from random import Random
from ga import run_ga
import sys

from ann_benchmarks import Or, And, Xor
from ann_plotter import ANNPlotter

args = {}

"""   
-------------------------------------------------------------------------
Edit this part to do the exercises
"""

# problem
problem_class = Or

# parameters for the GA
args["num_hidden_units"] = 0 # Number of hidden units of the neural network
args["gaussian_stdev"] = 1.0 # Standard deviation of the Gaussian mutations
args["crossover_rate"]  = 0.8 # Crossover fraction
args["tournament_size"] = 2 # Tournament size
args["pop_size"] = 10 # Population size

args["num_elites"] = 1 # number of elite individuals to maintain in each gen
args["mutation_rate"] = 0.5 # fraction of loci to perform mutation on

# by default will use the problem's defined init_range
# uncomment the following line to use a specific range instead
#args["pop_init_range"] = [-500, 500] # Range for the initial population

args["use_bounder"] = True # use the problem's bounder to restrict values
# comment out the previously line to run unbounded

args["max_generations"] = 100 # Number of generations of the GA
display = True # Plot initial and final populations

"""
-------------------------------------------------------------------------
"""

args["fig_title"] = 'GA'

if __name__ == "__main__":
    if len(sys.argv) > 1 :
        rng = Random(int(sys.argv[1]))
    else :
        rng = Random()
        
    best_individual, best_fitness, final_pop = run_ga(rng, display=display,
                                                      problem_class=problem_class,**args)
    print("Best Individual", best_individual)
    print("Best Fitness", best_fitness)
    
    if display :
        net = problem_class(args["num_hidden_units"]).net
        net.set_params(best_individual)
        
        ann_plotter = ANNPlotter(net)
        ann_plotter.draw()
        
        ioff()
        show()
