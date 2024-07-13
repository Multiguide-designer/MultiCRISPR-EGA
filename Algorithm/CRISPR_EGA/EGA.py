from functools import partial 
import random
import subprocess
import multiprocessing
import copy
import numpy as np
import time



class Individual:
    def __init__(individual,self):
        individual.index = [random.randint(1, max_index) for max_index in self.GeneMaxIndices]
    
    def combineSeq(individual,self):
        sequence_parts = [self.gene20ntInfo[gene][index - 1] for gene, index in zip(list(self.gene20ntInfo.keys()), individual.index)]
        individual.seq = self.RIBOJ + self.REPEAT + self.REPEAT.join(sequence_parts) + self.REPEAT


def run_RNAfold(sequence):
    cmd = ["RNAfold", "--noPS"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def evaluate_individual(individual):
    output = run_RNAfold(individual.seq)
    if output:
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
    else:
        energy = 0
    return energy


# Define tournament selection function
def tournament_selection(population, tournamentSize):
    selected_parents = []
    while len(selected_parents) < 1 * len(population):
        participants = random.sample(population, tournamentSize)
        winner = min(participants, key=lambda individual: individual.fitness)
        selected_parents.append(winner)
    return selected_parents


# Define roulette selection function
def roulette_selection(population):
    fitness_values = [ind.fitness for ind in population]
    total_fitness = sum(fitness_values)
    selection_probabilities = [fitness / total_fitness for fitness in fitness_values]
    selected_parents = np.random.choice(population, size=len(population), p=selection_probabilities).tolist()
    return selected_parents


def cxTwoPoint(individual1, individual2):
    size = min(len(individual1.index), len(individual2.index))
    cxpoint1 = random.randint(0, size)
    cxpoint2 = random.randint(0, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else:
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    individual1.index[cxpoint1:cxpoint2], individual2.index[cxpoint1:cxpoint2] = individual2.index[
                                                                               cxpoint1:cxpoint2], individual1.index[
                                                                                                   cxpoint1:cxpoint2]
    return individual1, individual2


def mutUniformInt(individual, mutation_probability, GENE_MAX_INDICES):
    for i in range(len(individual.index)):
        if random.random() < mutation_probability:  # Control mutation probability
            individual.index[i] = random.randint(1, GENE_MAX_INDICES[i])
    return individual


def get_uniqueMFE(population):
    unique_MFEs = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_MFEs.append(individual.fitness)
    return unique_MFEs
    

def remove_duplicate_individual(population):
    unique_individuals = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_individuals.append(individual)
    return unique_individuals

class Ga:
    def __init__(self,gene20ntInfo,insulatorSeq,repeatSeq,populationSize,maxIterations,crossover,mutation,tournamentSize):
        # EGA parameters
        self.populationSize = populationSize
        self.crossover = crossover
        self.mutation = mutation
        self.tournamentSize = tournamentSize
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
        self.final_individuals = []
        self.all_individuals = []

        self.converge_time = 0
        self.start_time = time.time()
        self.end_time = self.start_time
        
        # Create a pool of worker processes
        self.pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        self.population = [Individual(self) for _ in range(self.populationSize)]
        for individual in self.population:
            individual.combineSeq(self)

        # Parallelize the calculation of MFE for each individual in the population
        self.energies = self.pool.map(partial(evaluate_individual),self.population)
        
        # Assign fitness values to individuals
        for individual, energy in zip(self.population, self.energies):
            individual.fitness = energy  
        
        self.uniqueMFEs.append(get_uniqueMFE(self.population))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.best_parent_individual = min(self.population, key=lambda individual: individual.fitness)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.population)},best individual index:{self.best_parent_individual.index},best MFE:{self.best_parent_individual.fitness},average MFE:{self.average_MFEs[-1]}")
        

    def get_first_minMFE(self):
        self.best_parent_individual = min(self.population, key=lambda individual: individual.fitness)
        return  self.best_parent_individual.fitness
        
        
    def iter_once(self):
        
        # Select a parent (use tournament selection)
        self.parents = tournament_selection(self.population, self.tournamentSize)

        # Select parent (use roulette wheel selection)
        # parents = roulette_selection(population)

        # Select the next generation
        self.offspring = []
        for _ in range(int(self.populationSize/2)):
            self.parent1, self.parent2 = random.sample(self.parents, 2)
            if random.random() < self.crossover:
                self.child1, self.child2 = cxTwoPoint(copy.deepcopy(self.parent1), copy.deepcopy(self.parent2))
            else:
                self.child1, self.child2 = copy.deepcopy(self.parent1), copy.deepcopy(self.parent2)
            self.child1 = mutUniformInt(self.child1, self.mutation, self.GeneMaxIndices)
            self.child2 = mutUniformInt(self.child2, self.mutation, self.GeneMaxIndices)
            self.offspring.extend([self.child1, self.child2])

        for individual in self.offspring:
            individual.combineSeq(self)
            
        self.energies = self.pool.map(partial(evaluate_individual),self.offspring)
        
        # Assign fitness values to individuals
        for individual, energy in zip(self.offspring, self.energies):
            individual.fitness = energy  

        self.best_offspring_individual = min(self.offspring, key=lambda ind: ind.fitness)
        self.best_offspring_index = self.offspring.index(self.best_offspring_individual)

        # If the lowest value of the parent is less than the lowest value of the offspring,
        #  replace the lowest individual in the offspring with the lowest individual in the parent.
        if self.best_parent_individual.fitness < self.best_offspring_individual.fitness:
            self.offspring[self.best_offspring_index] = self.best_parent_individual
            self.energies[self.best_offspring_index] = self.best_parent_individual.fitness
            
        # Replace the old generation with the offspring (using Elitist GA)
        self.population[:] = self.offspring

        self.uniqueMFEs.append(get_uniqueMFE(self.population))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        if self.minimumMFEs[-1] != min(self.uniqueMFEs[-1]):
            self.end_time = time.time() # Update convergence end time
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.generation = self.generation + 1
        
        if self.generation >= self.maxIterations-10:
            self.all_individuals.extend(sorted(self.population, key=lambda ind: ind.fitness))
            self.all_individuals = remove_duplicate_individual(self.all_individuals)
        
        self.best_parent_individual = min(self.population, key=lambda individual: individual.fitness)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.population)},best individual index:{self.best_parent_individual.index},best MFE:{self.best_parent_individual.fitness},average MFE:{self.average_MFEs[-1]}")
        
        if self.generation == self.maxIterations:
            self.final_individuals = remove_duplicate_individual(self.population)
            self.converge_time = self.end_time  - self.start_time
        return self.uniqueMFEs,self.average_MFEs,self.minimumMFEs,self.final_individuals,self.all_individuals,self.converge_time
    

def main(gene20ntInfo,insulatorSeq,repeatSeq,maxIterations,populationSize,crossover,mutation,tournamentSize):
    ga_instance = Ga(gene20ntInfo,insulatorSeq,repeatSeq,populationSize,maxIterations,crossover,mutation,tournamentSize)
    for _ in range(int(maxIterations-1)):
        MFEs,Average_MFEs,minimumMFEs,final_individuals,all_individuals,converge_time = ga_instance.iter_once()
    return MFEs,Average_MFEs,minimumMFEs,final_individuals,all_individuals,converge_time

