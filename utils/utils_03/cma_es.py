from utils.utils_03.cma import *
from inspyred import benchmarks
from utils.utils_03.inspyred_utils import *
from utils.utils_03.plot_utils import *

from pylab import *

import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)

def run(random,display=False, num_vars=0, problem_class=benchmarks.Sphere, 
           **kwargs) :
    
    problem = problem_class(num_vars)
    if "pop_init_range" in kwargs :
        generator=generator
    else :
        generator=generator_wrapper(problem.generator)
    
    es = CMAEvolutionStrategy(generator(random, kwargs),
                                   kwargs["sigma"],
                                   {'popsize': kwargs["num_offspring"],
                                    'seed' : random.rand() * 100000,
                                    'CMA_mu' : kwargs["pop_size"]})
                                    #'CMA_elitist' : True})
    gen = 0
    while gen <= kwargs["max_generations"] :
        candidates = es.ask()    # get list of new solutions
        fitnesses = problem.evaluator(candidates, kwargs)
        
        if display :
            fitnesses_tmp = np.sort(fitnesses)
            average_fitness = mean(fitnesses_tmp)
            median_fitness = fitnesses[int(len(fitnesses_tmp)/2)]
            best_fitness = fitnesses_tmp[0]
            worst_fitness = fitnesses_tmp[-1]
            
            num_generations = gen
            num_evaluations = es.countevals
            
            if num_generations == 0 :
                initial_pop = asarray(candidates).copy()
                initial_fitnesses = asarray(fitnesses).copy()

                import matplotlib.pyplot as plt
                import numpy
                colors = ['black', 'blue', 'green', 'red']
                labels = ['average', 'median', 'best', 'worst']
                args = {}
                
                figure(kwargs["fig_title"] + ' (fitness trend)')
                plt.ion()
                data = [[num_evaluations], [average_fitness], [median_fitness], [best_fitness], [worst_fitness]]
                lines = []
                for i in range(4):
                    line, = plt.plot(data[0], data[i+1], color=colors[i], label=labels[i])
                    lines.append(line)
                args['plot_data'] = data
                args['plot_lines'] = lines
                plt.xlabel('Evaluations')
                plt.ylabel('Fitness')
            else :
                data = args['plot_data']
                data[0].append(num_evaluations)
                data[1].append(average_fitness)
                data[2].append(median_fitness)
                data[3].append(best_fitness)
                data[4].append(worst_fitness)
                lines = args['plot_lines']
                for i, line in enumerate(lines):
                    line.set_xdata(numpy.array(data[0]))
                    line.set_ydata(numpy.array(data[i+1]))
                args['plot_data'] = data
                args['plot_lines'] = lines
        
        es.tell(candidates, fitnesses)
        gen += 1
    
    final_pop = asarray(es.ask())
    final_pop_fitnesses = asarray(problem.evaluator(final_pop, kwargs))
    
    best_guy = es.best.x
    best_fitness = es.best.f
    
    if display :
        
        ymin = min([min(d) for d in data[1:]])
        ymax = max([max(d) for d in data[1:]])
        yrange = ymax - ymin
        plt.xlim((0, num_evaluations))
        plt.ylim((ymin - 0.1*yrange, ymax + 0.1*yrange))
        plt.draw()
        plt.legend()
        
        # Plot the parent and the offspring on the fitness landscape 
        # (only for 1D or 2D functions)
        if num_vars == 1 :
            plot_results_1D(problem, initial_pop, 
                                  initial_fitnesses, 
                                  final_pop, final_pop_fitnesses,
                                  'Initial Population', 'Final Population')
    
        elif num_vars == 2 :
            plot_results_2D(problem, initial_pop, 
                                  final_pop, 'Initial Population', 
                                  'Final Population')

    return best_guy, best_fitness
