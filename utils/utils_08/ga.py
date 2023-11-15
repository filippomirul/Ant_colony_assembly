from pylab import *
import utils.utils_08.plot_utils as plot_utils
from inspyred import ec, benchmarks
#from ann_benchmarks import NeuralNetworkBenchmark
from utils.utils_08.ann_benchmarks import NeuralNetworkBenchmark


def generate_offspring(random, x0, std_dev, num_offspring, display, kwargs) :
    
    x0 = asarray(x0, dtype=float64)
    
    problem = benchmarks.Sphere(len(x0))
    
    parent_fitness = problem.evaluator([x0], None)[0]

    algorithm = ec.EvolutionaryComputation(random)
    algorithm.terminator = ec.terminators.generation_termination
    algorithm.replacer = ec.replacers.generational_replacement    
    algorithm.variator = ec.variators.gaussian_mutation
    
    final_pop = algorithm.evolve(generator=(lambda random, args: x0.copy()),
                          evaluator=problem.evaluator,
                          pop_size=num_offspring,
                          maximize=False,
                          max_generations=1,
                          mutation_rate=1.0,
                          gaussian_stdev=std_dev)
    
    offspring_fitnesses = asarray([guy.fitness for guy in final_pop])
    offspring = asarray([guy.candidate for guy in final_pop])
    
    if display :        
        # Plot the parent and the offspring on the fitness landscape 
        # (only for 1D or 2D functions)
        if len(x0) == 1 :
            plot_utils.plot_results_1D(problem, x0, parent_fitness, offspring, 
                                  offspring_fitnesses, 'Parent', 'Offspring', kwargs)
    
        elif len(x0) == 2 :
            plot_utils.plot_results_2D(problem, asarray([x0]), offspring,
                                  'Parent', 'Offspring', kwargs)

    return (parent_fitness, offspring_fitnesses)

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

def run_ga(random,display=False, num_vars=0, problem_class=benchmarks.Sphere, 
           maximize=False, use_bounder=True, **kwargs) :
    
    #create dictionaries to store data about initial population, and lines
    initial_pop_storage = {}
    
    algorithm = ec.EvolutionaryComputation(random)
    algorithm.terminator = ec.terminators.generation_termination
    algorithm.replacer = ec.replacers.generational_replacement    
    algorithm.variator = [ec.variators.uniform_crossover,ec.variators.gaussian_mutation]
    algorithm.selector = ec.selectors.tournament_selection
    
    if display :
        algorithm.observer = [plot_utils.plot_observer,initial_pop_observer]
    else :
        algorithm.observer = initial_pop_observer 
    
    kwargs["num_selected"]=kwargs["pop_size"]  
    
    # special stuff for ANN Benchmarks    
    if issubclass(problem_class, NeuralNetworkBenchmark) :
        if "num_hidden_units" not in kwargs :
            kwargs["num_hidden_units"] = 0
        if "recurrent" not in kwargs :
            kwargs["recurrent"] = False
        problem = problem_class(kwargs["num_hidden_units"], kwargs["recurrent"])
        num_vars = problem.dimensions
    else :
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

    #best_guy = final_pop[0].candidate
    #best_fitness = final_pop[0].fitness
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
