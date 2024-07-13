from functools import partial  
import random
import subprocess
import multiprocessing 
import time
import numpy as np
import math
import copy


class Solution:
    def __init__(sol,self):
        sol.index = [random.randint(1, max_index) for max_index in self.GeneMaxIndices]
    
    def combineSeq(sol,self):
        sequence_parts = [self.gene20ntInfo[gene][index - 1] for gene, index in zip(list(self.gene20ntInfo.keys()), sol.index)]
        sol.seq = self.RIBOJ + self.REPEAT + self.REPEAT.join(sequence_parts) + self.REPEAT


def run_RNAfold(sequence):
    cmd = ["RNAfold", "--noPS"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def evaluate_solution(sol):
    output = run_RNAfold(sol.seq)
    if output:
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
    else:
        energy = 0
    return energy

# Add the generated neighborhood solution
def generate_neighbor(self):
    newsolutionSet = []
    for sol in self.solutionSet:
        neighbor = copy.deepcopy(sol)
        index_to_modify = random.randint(0, len(neighbor.index) - 1)
        neighbor.index[index_to_modify] = random.randint(1,self.GeneMaxIndices[index_to_modify])
        newsolutionSet.append(neighbor)  
    return newsolutionSet


def metropolis_criterion(current_energy, new_energy, current_temperature):
    delta_energy = new_energy - current_energy
    if delta_energy < 0:
        return True
    else:
        if random.random() < math.exp(-delta_energy / current_temperature):
            return True
        else:
            return False

def get_uniqueMFE(population):
    unique_MFEs = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_MFEs.append(individual.energy)
    return unique_MFEs

def remove_duplicate_solution(solutionSet):
    unique_solutions = []
    seen_indices = set()
    for solution in solutionSet:
        index_tuple = tuple(solution.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_solutions.append(solution)
    return unique_solutions



class SA:
    def __init__(self,gene20ntInfo,insulatorSeq,repeatSeq,numSolutions,maxIterations,cooling_rate):
        # SA parameters
        self.numSolutions = numSolutions
        self.cooling_rate = cooling_rate
        
        self.generation = 1
        self.maxIterations = maxIterations
        
        # Gene sequence parameters
        self.gene20ntInfo = gene20ntInfo
        self.GeneMaxIndices = [len(spacers) for spacers in gene20ntInfo.values()]
        self.RIBOJ = insulatorSeq
        self.REPEAT = repeatSeq

        random.seed(42)

        # Collect the fitness value of each generation, 
        # as well as the average and minimum fitness value of each generation group, 
        # and the unique individuals of the last generation
        self.uniqueMFEs = []
        self.average_MFEs = []
        self.minimumMFEs = []
        self.final_solutions = []
        self.all_solutions  =[]
        
        self.converge_time = 0
        self.start_time = time.time()
        self.end_time = self.start_time
        
        # Create a pool of worker processes
        self.pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        self.solutionSet = [Solution(self,) for _ in range(self.numSolutions)]
        for sol in self.solutionSet:
            sol.combineSeq(self)

        # Parallelize the calculation of MFE for each individual in the population
        self.energies = self.pool.map(partial(evaluate_solution),self.solutionSet)
    
        for sol, energy in zip(self.solutionSet, self.energies):
            sol.energy = energy  
            
        #The initial temperature is related to the maximum target value difference
        self.initial_temperature = max(self.energies) - min(self.energies)
        self.temperature = self.initial_temperature

        self.uniqueMFEs.append(get_uniqueMFE(self.solutionSet))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.bestsol = min(self.solutionSet, key=lambda sol: sol.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.solutionSet)},vest solution index:{self.bestsol.index},best MFE:{self.bestsol.energy},average MFE:{self.average_MFEs[-1]}")
        

    def get_first_minMFE(self):
        self.bestsol = min(self.solutionSet, key=lambda sol: sol.energy)
        return  self.bestsol.energy
        
        
    def iter_once(self):
        # At the current temperature, iteratively searches the neighborhood of the solution
        self.newsolutionSet = generate_neighbor(self)
        for sol in self.newsolutionSet:
            sol.combineSeq(self)
            
        self.new_energies = self.pool.map(partial(evaluate_solution), self.newsolutionSet)
        for sol, new_energy in zip(self.newsolutionSet, self.new_energies):
            sol.energy = new_energy  
        
        # Determine whether to accept a new solution based on the Metropolis criterion
        for i in range(len(self.solutionSet)):
            if metropolis_criterion(self.energies[i], self.new_energies[i],  self.temperature):
                # Accept new solutions
                continue
            else:
                # Do not accept new solutions
                self.newsolutionSet[i] = self.solutionSet[i]
                self.new_energies[i] = self.energies[i]
        
        # Reduce the temperature
        self.temperature = self.temperature * self.cooling_rate
       
        self.best_newsol = min(self.newsolutionSet, key=lambda sol: sol.energy)
        self.best_newsol_index = self.newsolutionSet.index(self.best_newsol)

        # If the lowest value of the original solution set is less than the lowest value of the new solution set, 
        # replace the lowest individual in the new solution set with the lowest individual in the original solution set.
        if self.bestsol.energy < self.best_newsol.energy:
            self.newsolutionSet[self.best_newsol_index] = self.bestsol
            self.new_energies[self.best_newsol_index] = self.bestsol.energy

        # Replace the old generation with the offspring 
        self.solutionSet[:] = self.newsolutionSet
        self.energies[:] = self.new_energies

        self.uniqueMFEs.append(get_uniqueMFE(self.solutionSet))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        if self.minimumMFEs[-1] != min(self.uniqueMFEs[-1]):
            self.end_time = time.time() # Update convergence end time
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        if self.generation >= self.maxIterations-10:
            self.all_solutions.extend(sorted(self.solutionSet, key=lambda ind: ind.energy))
            self.all_solutions = remove_duplicate_solution(self.all_solutions)
        
        self.generation = self.generation + 1
        self.bestsol = min(self.solutionSet, key=lambda sol: sol.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.solutionSet)},best solution index:{self.bestsol.index},best MFE:{self.bestsol.energy},average MFE:{self.average_MFEs[-1]}")
        
        if self.generation == self.maxIterations:
            self.final_solutions = remove_duplicate_solution(self.solutionSet)
            self.converge_time = self.end_time  - self.start_time
        return self.uniqueMFEs,self.average_MFEs,self.minimumMFEs,self.final_solutions,self.all_solutions,self.converge_time

def main(gene20ntInfo,insulatorSeq,repeatSeq,maxIterations,numSolutions,cooling_rate):
    sa_instance = SA(gene20ntInfo,insulatorSeq,repeatSeq,numSolutions,maxIterations,cooling_rate)
    for _ in range(int(maxIterations-1)):
        MFEs,Average_MFEs,minimumMFEs,final_solutions,all_solutions,converge_time = sa_instance.iter_once()
    return MFEs,Average_MFEs,minimumMFEs,final_solutions,all_solutions,converge_time
        
