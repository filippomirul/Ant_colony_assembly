from pylab import *
import utils.utils_06.plot_utils as plot_utils
from inspyred import ec, benchmarks

from inspyred.ec.observers import *

import matplotlib

STAR = 'star'
RING = 'ring'

def generator(random, args):
    return asarray([random.uniform(args["pop_init_range"][0],
                                   args["pop_init_range"][1]) 
                    for _ in range(args["num_vars"])])

def initial_pop_observer(population, num_generations, num_evaluations, 
                         args):
    if num_generations == 0 :
        args["initial_pop_storage"]["individuals"] = asarray([guy.candidate 
                                                 for guy in population]) 
        args["initial_pop_storage"]["fitnesses"] = asarray([guy.fitness 
                                          for guy in population]) 

def run_pso(random, display=False, num_vars=0, problem_class=benchmarks.Sphere, 
           maximize=False, use_bounder=True, **kwargs) :
    
    #create dictionaries to store data about initial population, and lines
    initial_pop_storage = {}
    
    algorithm = inspyred.swarm.PSO(random)
    algorithm.topology = inspyred.swarm.topologies.star_topology
    algorithm.terminator = ec.terminators.generation_termination
    
    if display :
        algorithm.observer = [plot_utils.plot_observer,initial_pop_observer]
    else :
        algorithm.observer = initial_pop_observer
    
    if "topology" in kwargs :
        if kwargs["topology"] is STAR:
            algorithm.topology = inspyred.swarm.topologies.star_topology
        elif kwargs["topology"] is RING:
            algorithm.topology = inspyred.swarm.topologies.ring_topology

    kwargs["num_selected"]=kwargs["pop_size"]
    problem = problem_class(num_vars)
    if use_bounder :
        kwargs["bounder"]=problem.bounder
    if "pop_init_range" in kwargs :
        kwargs["generator"]=generator
    else :
        kwargs["generator"]=problem.generator
    
    final_pop = algorithm.evolve(evaluator=problem.evaluator,  
                          maximize=problem.maximize,
                          initial_pop_storage=initial_pop_storage,
                          num_vars=num_vars,
                          **kwargs)
    
    final_pop_fitnesses = asarray([guy.fitness for guy in final_pop])
    final_pop_candidates = asarray([guy.candidate for guy in final_pop])

    sort_indexes = sorted(range(len(final_pop_fitnesses)), key=final_pop_fitnesses.__getitem__)
    final_pop_fitnesses = final_pop_fitnesses[sort_indexes]
    final_pop_candidates = final_pop_candidates[sort_indexes]
    
    best_guy = final_pop_candidates[0]
    best_fitness = final_pop_fitnesses[0]

    if display :
        # Plot the parent and the offspring on the fitness landscape 
        # (only for 1D or 2D functions)
        if num_vars == 1 :
            plot_utils.plot_results_1D(problem, initial_pop_storage["individuals"], 
                                  initial_pop_storage["fitnesses"], 
                                  final_pop_candidates, final_pop_fitnesses,
                                  'Initial Population', 'Final Population', kwargs)
    
        elif num_vars == 2 :
            plot_utils.plot_results_2D(problem, initial_pop_storage["individuals"], 
                                  final_pop_candidates, 'Initial Population',
                                  'Final Population', kwargs)

    return best_guy, best_fitness, final_pop
