
from Bio import pairwise2
from Bio import SeqIO
import yaml
from math import modf
from Bio.Seq import Seq
from inspyred import swarm
from inspyred import benchmarks
from inspyred import ec
from inspyred.ec import selectors
from collections import deque
from collections import Counter
import numpy as np
import itertools
from random import Random
import random
import seaborn as sns
import os
import sys
import inspyred
import collections
collections.Iterable = collections.abc.Iterable
collections.Sequence = collections.abc.Sequence


def comstum_reads(seq: str, length_reads = 10, coverage = 5, verbose = False) -> list:
    
    """The function split the sequence in input in reads.
    The splitting is done using random numbers, the amount of reads is given by: (len(seq)/length_read)*coverage.
    """

    number_of_reads = int(len(seq)/length_reads) * coverage
    starting_pos = random.sample(range(0, len(seq)-length_reads+1), number_of_reads)
    reads = []

    for num in starting_pos:
        reads.append(seq[num:num+length_reads])

    if verbose == True:
        y = [0 for i in range(0,len(seq)+1)]
        for i in starting_pos:
            for j in range(i, i+length_reads+1):
                y[j] += 1 
        sns.set_theme(style="darkgrid")
        sns.lineplot(y)
        print(f"There are {y.count(0)} bases that have not been transcriped.")

    return reads


def eval_allign(reads:list, par = [2, -5, -22, -5]):
    
    """Function that evaulates the alignment

    reads: list of DNA sequences, each read is a string

    par: list of parameters to perform the alignment
    es (the examples represent the default parameters):
    match_score = 2,
    mismatch_penalty = -5,
    gap_open_penalty = -22,
    gap_extension_penalty = -5

    output:
    Matrix with the weigts (distances) between the reads (nodes)
    This matrix can host two kind of entires: alignment scores (integers) and lists. Each list reports: the index of the 
    beginning of the alignment with respect to the first read of the couple (read i), the index of the end of the alignment
    with respect to read i, the index of the beginning of the alignment with respect to the second read of the couple (j)
    and the index of the end of the alignment with respect to read j. All indexes are stored as included. 
    
    An entry ij is filled with the alignment score OR with the alignment list, with the other entry being saved in position
    ji. The decision is made on the basis of which of the two reads starts first in the pairwise alignment. 
    
    """
    
    length = len(reads)
    # initialization of the matrix
    weigth_matrix = [[0] * len(reads) for x in range(len(reads))]

    # The score of the allingment of read[1] to read[2] is the same of the opposite (read[2] to read[1])
    # So when the function found the diretionality of the allignment put the score in rigth spot and a 0 in the wrong one.
    visited = deque([j for j in range(length)])
    
    for i in range(length):

        for j in visited:

            if i == j:
                # the diagonal of the matrix has 0 score because we are allinging the same sequence
                continue
            else:
                # pairwise must return a positive score, if there is no one it return None
                allignment = pairwise2.align.localms(Seq(reads[i]), Seq(reads[j]), par[0], par[1], par[2], par[3])
                if len(allignment) == 0:
                    #   allignment = 1    to decide if put 1 or 0
                    continue
                else:
                    allignment = allignment[0]
                      
                start_r1 = allignment[3]                                              # Start of the alignment, with respect to the first read (i)
                end_r1 = allignment[4]                                                # End of the alignment,  with respect to the first read (i)
                start_r2 = allignment[3] - allignment[1].count('-')                   # Start of the alignment, with respect to the second read (j)
                end_r2 = start_r2 + (end_r1 - start_r1)                               # End of the alignment, with repsect to the second read (j)
                # return object [seqA, seqB, score, start(inc), end(ex)]

                if allignment[0][0] == "-":
                    # This means that the first reads in input has a gap at the beginning of the allignment.
                    # Therefore the first reads in input (i) is downstream,
                    # so I add the score in the matrix but in the position (j,i) instead of (i,j), where the alignment list is stored
                    weigth_matrix[j][i] = float(allignment[2])
                    weigth_matrix[i][j] = [start_r1, end_r1-1, start_r2, end_r2-1]    # Start and end included

                else:
                    # In the opposite case, where the i read is upstream, (i,j) has the score, while (j,i) hosts the alignment list                    
                    weigth_matrix[i][j] = float(allignment[2])
                    weigth_matrix[j][i] = [start_r1, end_r1-1, start_r2, end_r2-1]    # Start and end included

                    
        visited.popleft()

    return weigth_matrix 


def consensus_reconstructor(path:list, reads:list, positions:list) -> str:
    """Function to reconstruct the sequence
    path: is the candidate path of the graph, a list with the indexes of the reads
    reads: list of all the reads
    positions: is the matrix with the weights and the list with the positions 
    """
    
    # CONSENSUS MATRIX BUILDER: builds a matrix containing all the pairwise local alignments between the reads in the path; all the empty positions are filled with '-'.
    
    consensus_matrix = []
    consensus_matrix.append([i for i in reads[path[0].element[0]]])
    shift_required = 0                                                          # Initial shift is 0 (no shift, the first read starts where the consensus matrix starts)
    
    new_path = []
    for i in path:
        new_path.append(i.element)

    for i,j in new_path:
        consensus_matrix.append([i for i in reads[j]])
        
    iteration_number = -1
    # print(consensus_matrix)
    for i,j in new_path:
        iteration_number = iteration_number + 1
        nums = positions[j][i]
        start_r1 = nums[0]                                                      # Index of the begininng of the alignment, with respect to the read i, included
        # end_r1 = nums[1]                                                      # Index of the end of the alignment, with respect to the read i, included
        start_r2 = nums[2]                                                      # Index of the begininng of the alignment, with respect to the read j, included
        # end_r2 = nums[3]                                                      # Index of the end of the alignment, with respect to the read j, included
        shift_required = shift_required + (start_r1 - start_r2)                 # Amount of shifting required to translate the read j up to the correct position
        for i in range(shift_required):
            consensus_matrix[iteration_number+1].insert(0,"-")
        # consensus_matrix[iteration_number+1].insert(0, '-'*shift_required)    # Performing the shift
    # print(consensus_matrix)
    max_length = len(consensus_matrix[iteration_number+1])                      # Filling the empty positions with '-'
    for row in consensus_matrix:
        row.extend('-'*(max_length - len(row)))
    
    # CONSENSUS SEQUENCE BUILDER: goes through all columns in the consensus matrix and selects, for each column, the most frequent base. Then it elongates a string with the 
    # selected base. If two or more bases have the same frequency, the first one met is chosen.
    
    consensus_seq = ''
    for k in range(len(consensus_matrix[0])):                                        # All rows in the consensus matrix have the same length, so I choose the first one
        bases_in_current_column = []
        for row in consensus_matrix:                                                 # Appending all elements in the current consensus matrix column in the list 'bases_in_current_column' 
            if row[k] != '-':                                                        # If the element is a spacer, we do not add it
                bases_in_current_column.append(row[k])
        most_present_base = Counter(bases_in_current_column).most_common(1)[0][0]    # Most common base in 'bases_in_current_column' is found with the imported Counter class
        consensus_seq = consensus_seq + most_present_base                            # If two or more bases have the same frequency, Counter chooses by default the first base met.  
    
    return consensus_matrix, consensus_seq


class Assembly_problem():
    """Defines the de novo genome assembly problem.
    
    This class based on the Traveling Salesman problem defines the problem
    of assembling a new genome for which no reference is available (de novo assembly):
    given a set of genomic reads and their pairwise overlap score, find the
    path generating the longest consensus sequence. This problem assumes that 
    the ``weights`` parameter is an *n*-by-*n* matrix of pairwise 
    overlap among *n* reads. This problem is treated as a 
    maximization problem, so fitness values are determined to be the 
    proportional to the sum of the overlaps between each couple of reads
    (the weight of the edge) and the length of the final assembled sequence.
    
    Public Attributes:
    
    - *weights* -- the two-dimensional list of pairwise overlap 
    - *components* -- the set of ``TrailComponent`` objects constructed
        from the ``weights`` attribute, where the element is the ((source,
        destination), weight)
    - *bias* -- the bias in selecting the component of maximum desirability
        when constructing a candidate solution for ant colony optimization 
        (default 0.5)
    - *sigma* -- severity of the penality applied when the length of the 
        reconstructed sequence falls outside of the tolerated ± 4% symmetric 
        interval centered around the experimental lenght
    - *experimental_lenght* -- The algorithm needs an approximate length of the sequence to score
        properly the path. In fact the path return a candidate allignment, which is
        compared to a value experimentally derived.
    
    """
    def __init__(self, reads, experimental_length):
        self.sigma = 0.1
        self.experimental_length = experimental_length
        self.reads = reads
        self.weights = eval_allign(self.reads)
        self.components = [swarm.TrailComponent((i, j), value=(self.weights[i][j])) for i, j in itertools.permutations(range(len(self.weights)), 2) if (type(self.weights[i][j]) == float) and (self.weights[i][j] != 0)]
        self.bias = 0.65
        self.bounder = ec.DiscreteBounder([i for i in range(len(self.weights))])
        self.maximize = True
    
    def constructor(self, random, args):
        """Return a candidate solution for an ant colony optimization."""

        candidate = []
        feasible_components = [1]   #Fake initialization to allow while loop to start
        
        # We need to visit all the nodes that CAN be visited, the graph is directed and not complete, meaning we can have no more nodes to visit without visiting all the
        # nodes in the graph, thus, our termination condition is not visitin all the nodes but not having anymore feasible components
        while len(feasible_components) > 0:
            # At the start of the visit, all the components are feasible
            if len(candidate) == 0:
                feasible_components = self.components
            elif len(candidate) == len(self.weights) - 1: # All the nodes have been visited
                return candidate
            else:
                # Update feasible components and set of already visited nodes considering the node visited in the last iteration
                last = candidate[-1]
                already_visited = [c.element[0] for c in candidate]
                already_visited.extend([c.element[1] for c in candidate])
                already_visited = set(already_visited)
                feasible_components = [c for c in self.components if c.element[0] == last.element[1] and c.element[1] not in already_visited]
            # Choose a feasible component
            if len(feasible_components) == 0:
                return candidate
            if random.random() <= self.bias:
                next_component = max(feasible_components)
            else:
                next_component = selectors.fitness_proportionate_selection(random, feasible_components, {'num_selected': 1})[0]
            candidate.append(next_component)
        return candidate
    
    def evaluator(self, candidates, args):
            """Return the fitness values for the given candidates."""
            #P rappresenta la lunghezza stimata del nostro dna
            #Pm rappresenta P-0.04P
            #PM rappresenta P+0.04P
            #sigma = peso della penalità, che aumenta all'aumentare della distanza dal valore di lunghezza dal limite superiore o inferiore.
            fitness = []
            Pm = self.experimental_length - 0.04*self.experimental_length
            PM = self.experimental_length + 0.04*self.experimental_length
            for candidate in candidates:
                    total = 0
                    lencandidate = len(consensus_reconstructor(path = candidate, reads=self.reads, positions=self.weights)[1])
                    for c in candidate:
                        total += self.weights[c.element[0]][c.element[1]]
                    if lencandidate >= Pm:
                        if lencandidate <= PM :
                                fitness.append(total) 
                        else:
                                total=total-self.sigma*(lencandidate-PM)
                                fitness.append(total)
                    else:
                            total=total-self.sigma*(Pm-lencandidate)
                            fitness.append(total)
            return fitness


def run_simulation(sequence: str, pop:int, num_max_generations: int,
                    seed_r:int, evapor_rate:int, rate_of_learning: int, read_length:int, desired_coverage:int, display=True):
    """This is the function which runs the entire algorithm.
    sequence: complete sequence of the genome
    coverage: coverage of the allignment
    """
    # common parameters
    pop_size = pop
    max_generations = num_max_generations
    seed = seed_r
    prng = Random(seed)
    # ACS specific parameters
    evaporation_rate = evapor_rate
    learning_rate = rate_of_learning

    # generate the reads from the input sequence
    args = {}
    args["fig_title"] = "ACS"
    reads = comstum_reads(sequence, length_reads = read_length, verbose = True, coverage = desired_coverage)

    # run ACS
    expected_length = len(sequence)
    problem = Assembly_problem(reads = reads, experimental_length = expected_length)
    ac = inspyred.swarm.ACS(prng, problem.components)
    ac.observer = [plot_observer]
    ac.terminator = inspyred.ec.terminators.generation_termination
    final_pop = ac.evolve(generator=problem.constructor, 
                        evaluator=problem.evaluator, 
                        bounder=problem.bounder,
                        maximize=problem.maximize, 
                        pop_size=pop_size,
                        max_generations=max_generations,
                        evaporation_rate=evaporation_rate,
                        learning_rate=learning_rate,**args)
    best_ACS = max(ac.archive)

    # reconstruct sequence of best candidate
    reconstructed_seq = consensus_reconstructor(path = best_ACS.candidate, reads=reads, positions=eval_allign(reads))[1]

    # perform global alignment between reconstructed best candidate and reference sequence
    final_alignment = pairwise2.align.globalms(reconstructed_seq, ref, 3,-1,-3,-2)[0]

    # write results to "assembly_results" file
    with open("assembly_results", "w") as results_file:
        results_to_write = ["Reference sequence\n", ref, "\n", "\n","Sequence reconstructed by ACO assembler:\n", reconstructed_seq, "\n", "\n",
                        "Length of the reconstructed sequence:\n", str(len(reconstructed_seq)), "\n", "\n", "Score of the alignment:\n", str(final_alignment[2]), "\n", "\n", 
                        "Global alignment:\n", "\n"]
        results_file.writelines(results_to_write)
        global_ref = final_alignment[1]                                                                                            # the ref sequence, as in the global alignment
        global_reconstructed = final_alignment[0]                                                                                  # the reconstructed sequence, as in the global alignment
        splitted_global_ref = [global_ref[i:i + 100] for i in range(0, len(global_ref), 100)]                                      # the first, split every 100pb
        splitted_global_reconstructed = [global_reconstructed[i:i + 100] for i in range(0, len(global_reconstructed), 100)]        # the second, split every 100bp
        for i in range(len(splitted_global_ref)):
            results_file.write("reference")
            results_file.write("\n")
            results_file.write(splitted_global_ref[i])
            results_file.write("\n")
            results_file.write(splitted_global_reconstructed[i])
            results_file.write("\n")
            results_file.write("reconstructed")
            results_file.write("\n")
            results_file.write("\n")

    return


# input sequence management: by default it is not a genome, but a short 280 bp sequence of TP53, hosted in the Data folder together with longer genomic sequences
os.chdir(sys.path[0])                                               # changes working directory to directory of this python script 
support_list = []                                                   # holds the content of the file, line - wise
with open("./Data/TP53.txt", "r") as ref_handler:         # opens default file from the Data folder; the latter contains all test genomes + the short 280 bp TP53 sequence
    for line in ref_handler:
        support_list.append(line)
ref = "".join(support_list[1:]).replace("\n", "")                   # joins the content of the support list into a single string, skipping its first element (sequence header)
ref = ref[:500]                                                     # only a small nunber of bases considered for TP53 (in the report a run with 900bp was performed)

## Sequential call of the whole algorithm ##

# common parameters initialization
pop_size_default = 100
max_generations_default = 50
seed_default = 10
prng_default = Random(seed_default)
display_default = True
length_read_test = 160
coverage_default = 6

## ACS specific parameters initialization
evaporation_rate_default = 0.4
learning_rate_default = 0.1

# generate the reads from the input sequence
args = {}
args["fig_title"] = "ACS"
reads = comstum_reads(ref, length_reads = length_read_test, verbose = True, coverage = coverage_default)

from utils.utils_07.exercise_1 import *
from utils.utils_07.plot_utils import *

# run ACS
expected_length = len(ref)
problem = Assembly_problem(reads = reads, experimental_length = expected_length)
ac = inspyred.swarm.ACS(prng_default, problem.components)
ac.observer = [plot_observer]
ac.terminator = inspyred.ec.terminators.generation_termination
final_pop = ac.evolve(generator=problem.constructor, 
                    evaluator=problem.evaluator, 
                    bounder=problem.bounder,
                    maximize=problem.maximize, 
                    pop_size=pop_size_default,
                    max_generations=max_generations_default,
                    evaporation_rate=evaporation_rate_default,
                    learning_rate=learning_rate_default,**args)
best_ACS = max(ac.archive)

# reconstruct sequence of best candidate
reconstructed_seq = consensus_reconstructor(path = best_ACS.candidate, reads=reads, positions=eval_allign(reads))[1]

# perform local alignment between reconstructed best candidate and reference sequence
final_alignment = pairwise2.align.localms(reconstructed_seq, ref, 3,-1,-3,-2)[0]

# write results to "assembly_results" file
with open("assembly_results", "w") as results_file:
    results_to_write = ["Reference sequence\n", ref, "\n", "\n","Sequence reconstructed by ACO assembler:\n", reconstructed_seq, "\n", "\n",
                    "Length of the reconstructed sequence:\n", str(len(reconstructed_seq)), "\n", "\n", "Score of the alignment:\n", str(final_alignment[2]), "\n", "\n"]
    results_file.writelines(results_to_write)
    global_ref = final_alignment[1]                                                                                            # the ref sequence, as in the global alignment
    global_reconstructed = final_alignment[0]                                                                                  # the reconstructed sequence, as in the global alignment
    splitted_global_ref = [global_ref[i:i + 100] for i in range(0, len(global_ref), 100)]                                      # the first, split every 100pb
    splitted_global_reconstructed = [global_reconstructed[i:i + 100] for i in range(0, len(global_reconstructed), 100)]        # the second, split every 100bp

    mism_counter = 0
    for i in range(len(global_reconstructed)):                                                                                 # small script to count number of matches
        if global_reconstructed[i] == global_ref[i]:
            mism_counter = mism_counter + 1
    perc_of_matches = (mism_counter / len(global_ref))*100
    results_file.writelines(["Percentage of matches:\n", str(perc_of_matches), "\n", "\n"])

    results_file.writelines (["Local alignment:\n", "\n"])                                                                     # small script to print decently the local alignment
    for i in range(len(splitted_global_ref)):
        results_file.write("reference")
        results_file.write("\n")
        results_file.write(splitted_global_ref[i])
        results_file.write("\n")
        results_file.write(splitted_global_reconstructed[i])
        results_file.write("\n")
        results_file.write("reconstructed")
        results_file.write("\n")
        results_file.write("\n")

# results are saved in the assembly_results file