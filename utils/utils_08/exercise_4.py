"""
3-input Boolean function
"""

from __future__ import print_function
import os
import neat
import visualize

from pylab import *

# 3-input Boolean function inputs and expected outputs.
inputs = [(0.0, 0.0, 0.0),
          (0.0, 0.0, 1.0),
          (0.0, 1.0, 0.0),
          (0.0, 1.0, 1.0),
          (1.0, 0.0, 0.0),
          (1.0, 0.0, 1.0),
          (1.0, 1.0, 0.0),
          (1.0, 1.0, 1.0)]
outputs = [(0.0,),(1.0,),(1.0,),(0.0,),(1.0,),(0.0,),(0.0,),(0.0,)]

num_generations = 100
num_runs = 1

config_files = ['config-feedforward-3input-function-nohidden.txt',
                'config-feedforward-3input-function-hidden.txt']

def eval_genomes(genomes, config):
    for genome_id, genome in genomes:
        genome.fitness = len(inputs)
        net = neat.nn.FeedForwardNetwork.create(genome, config)
        for xi, xo in zip(inputs, outputs):
            output = net.activate(xi)
            genome.fitness -= (output[0] - xo[0]) ** 2


if __name__ == '__main__':
    
    local_dir = os.path.dirname(__file__)
    
    if num_runs == 1:
    
        # Load configuration.
        config_file = os.path.join(local_dir, config_files[0])
        config = neat.Config(neat.DefaultGenome, neat.DefaultReproduction,
                             neat.DefaultSpeciesSet, neat.DefaultStagnation,
                             config_file)
                             
        # Create the population, which is the top-level object for a NEAT run.
        p = neat.Population(config)

        # Add a stdout reporter to show progress in the terminal.
        stats = neat.StatisticsReporter()
        p.add_reporter(neat.StdOutReporter(True))
        p.add_reporter(stats)

        # run NEAT for num_generations
        winner = p.run(eval_genomes, num_generations)

        # Display the winning genome.
        print('\nBest genome:\n{!s}'.format(winner))

        # Show output of the most fit genome against training data.
        print('\nOutput:')
        winner_net = neat.nn.FeedForwardNetwork.create(winner, config)
        for xi, xo in zip(inputs, outputs):
            output = winner_net.activate(xi)
            print("input {!r}, expected output {!r}, got {!r}".format(xi, xo, output))

        node_names = {-1:'A', -2: 'B', -3: 'C', 0:'f(A,B,C)'}
        visualize.draw_net(config, winner, filename='3-input Bool function', view=True, node_names=node_names)
        visualize.plot_stats(stats, ylog=False, view=True)
        visualize.plot_species(stats, view=True)

    else:
        
        results = []
        for file in config_files:
        
            # Load configuration.
            config_file = os.path.join(local_dir, file)
            config = neat.Config(neat.DefaultGenome, neat.DefaultReproduction,
                                 neat.DefaultSpeciesSet, neat.DefaultStagnation,
                                 config_file)

            best_fitnesses = []
            for i in range(num_runs):
                print('{0}/{1}'.format(i+1,num_runs))
                p = neat.Population(config)
                winner = p.run(eval_genomes, num_generations)
                best_fitnesses.append(winner.fitness)
            results.append(best_fitnesses)

        fig = figure('NEAT')
        ax = fig.gca()
        ax.boxplot(results)
        ax.set_xticklabels(['Without hidden nodes', 'With hidden nodes'])
        #ax.set_yscale('log')
        ax.set_xlabel('Condition')
        ax.set_ylabel('Best fitness')
        show()
