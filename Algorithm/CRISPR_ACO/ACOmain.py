from functools import partial  
import random
import subprocess
import multiprocessing
import time
import numpy as np
import copy



class Ant:
    def __init__(ant,self):
        ant.index = [random.randint(1, max_index) for max_index in self.GeneMaxIndices]
    
    def combineSeq(ant,self):
        sequence_parts = [self.gene20ntInfo[gene][index - 1] for gene, index in zip(list(self.gene20ntInfo.keys()), ant.index)]
        ant.seq = self.RIBOJ + self.REPEAT + self.REPEAT.join(sequence_parts) + self.REPEAT


def run_RNAfold(sequence):
    cmd = ["RNAfold", "--noPS"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def evaluate_ant(ant):
    output = run_RNAfold(ant.seq)
    if output:
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
    else:
        energy = 0
    return energy


'''Get heuristic information'''
def get_heuristic(gene20ntInfo,insulatorSeq,repeatSeq,):
    
    heuristic =  {(gene, spacer): 0 for gene, spacers in gene20ntInfo.items() for spacer in spacers}  
    index = 0
    
    for gene, spacers in gene20ntInfo.items():
        energies = []
        index = index + 1
        for spacer in spacers:
            if index > 1:
                sequence = repeatSeq + spacer + repeatSeq
            else:
                sequence = insulatorSeq + spacer + repeatSeq
            
            output = run_RNAfold(sequence)
            if output:
                energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
            else:
                energy = 0
                
            energies.append(energy)
        max_energy = np.max(energies)
        min_energy = np.min(energies)
        # Use Min-Max normalized free energy as heuristic information. Note that the smaller the free energy, the better.
        normalized_energies = [(max_energy - energy) / (max_energy - min_energy + 1) for energy in energies] 
        for spacer, normalized_energy in zip(spacers, normalized_energies):
            heuristic[(gene, spacer)] = normalized_energy
        
    return heuristic
            

'''Update total pheromones'''  
def update_global_pheromones(self):

    energies = [ant.energy for ant in self.antPopulation]
    
    mean_energy = np.mean(energies)
    max_energy = np.max(energies)
    min_energy = np.min(energies)

    # Pheromone volatilization
    for key in self.pheromones.keys():
        self.pheromones[key] *= (1 - self.Global_evaporationRate)

    # Update pheromones
    # Global increase, only ants smaller than the average are allowed to release pheromone
    for ant in self.antPopulation:
        if ant.energy < mean_energy:
            for gene, index in zip(list(self.gene20ntInfo.keys()), ant.index):
                spacer = self.gene20ntInfo[gene][index - 1]
                num_spacers = len(self.gene20ntInfo[gene]) 
                self.pheromones[(gene, spacer)] += ((max_energy - ant.energy) / (max_energy - min_energy+1))/num_spacers

    # (local increase)
    best_ant = min(self.antPopulation, key=lambda ant: ant.energy)
    for gene, index in zip(list(self.gene20ntInfo.keys()), best_ant.index):
        spacer = self.gene20ntInfo[gene][index - 1]
        num_spacers = len(self.gene20ntInfo[gene])  
        self.pheromones[(gene, spacer)] +=  (((max_energy - best_ant.energy ) / (max_energy - min_energy+1)))/num_spacers

    return self.pheromones

'''Update the pheromone of a single path'''  


def update_single_pheromones(self,gene,next_spacer_index):

    spacer = self.gene20ntInfo[gene][next_spacer_index]
    self.pheromones[(gene, spacer)]  =  self.pheromones[(gene, spacer)]*(1-self.Single_evaporationRate)+self.Single_evaporationRate * self.pheromones_ori[(gene, spacer)]
    

'''Roulette selection function'''
def roulette_selection(self,gene,spacers,q0):
    if random.random() <= q0:
        # Exploitation strategy: directly select the next node with the greatest appeal
        max_attractiveness_index = max(range(len(spacers)),
                                       key=lambda i: (self.pheromones[(gene, spacers[i])] ** self.alpha) *
                                                     (self.heuristic[(gene, spacers[i])] ** self.beta))
        return max_attractiveness_index
    
    else:

        attractiveness_values = []

        for spacer in spacers:
            pheromone_level = self.pheromones[(gene, spacer)] ** self.alpha
            heuristic_info = self.heuristic[(gene, spacer)] ** self.beta
            attractiveness = pheromone_level * heuristic_info
            attractiveness_values.append(attractiveness)
        
        total_attractiveness = sum(
            (self.pheromones[(gene, spacer)] ** self.alpha) * (self.heuristic[(gene, spacer)] ** self.beta)
        for spacer in spacers
        )
        random_number = random.uniform(0, total_attractiveness)
        accumulated_attractiveness = 0

        for next_spacer_index,spacer in enumerate(spacers):
            pheromone_level = self.pheromones[(gene, spacer)] ** self.alpha
            heuristic_info = self.heuristic[(gene, spacer)]  ** self.beta
            attractiveness = pheromone_level * heuristic_info
            accumulated_attractiveness += attractiveness

            if accumulated_attractiveness > random_number:
                break
        
        return next_spacer_index


def update_index(self,gene,spacers,ant,pos):
    
    next_spacer_index = roulette_selection(self,gene,spacers,self.q0)
    update_single_pheromones(self,gene,next_spacer_index)
    
    ant.index[pos] = next_spacer_index+1
   
    
def construct_ant_path(self):

    pos = -1
    for gene, spacers in self.gene20ntInfo.items():
        pos = pos + 1
        for ant in self.newAntPopulation:
            update_index(self,gene,spacers,ant,pos) 


def get_uniqueMFE(population):

    unique_MFEs = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_MFEs.append(individual.energy)
    return unique_MFEs


def remove_duplicate_ant(population):
    unique_ants = []
    seen_indices = set()
    for ant in population:
        index_tuple = tuple(ant.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_ants.append(ant)
    return unique_ants



class ACO:
    def __init__(self,gene20ntInfo,insulatorSeq,repeatSeq,numAnts,maxIterations,alpha, beta, Global_evaporationRate,Single_evaporationRate,q0):
        # ACO parameters
        self.numAnts = numAnts
        self.alpha = alpha
        self.beta = beta
        self.Global_evaporationRate = Global_evaporationRate
        self.Single_evaporationRate = Single_evaporationRate
        self.q0 = q0
        
        # Initialize pheromone
        pheromones = {(gene, spacer): 1/len(spacers) for gene, spacers in gene20ntInfo.items() for spacer in spacers} 
        pheromones_ori = copy.copy(pheromones)

        # Heuristic information
        heuristic = get_heuristic(gene20ntInfo,insulatorSeq,repeatSeq,)
        
        self.pheromones_ori = pheromones_ori
        self.pheromones = pheromones
        self.heuristic = heuristic
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
        self.final_ants = []
        self.all_ants = []

        self.converge_time = 0
        self.start_time = time.time()
        self.end_time = self.start_time
        
        # Create a pool of worker processes
        self.pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        self.antPopulation = [Ant(self) for _ in range(self.numAnts)]
        for ant in self.antPopulation:
            ant.combineSeq(self)

        # Parallelize the calculation of MFE for each individual in the population
        self.energies = self.pool.map(partial(evaluate_ant),self.antPopulation)
    
        for ant, energy in zip(self.antPopulation, self.energies):
            ant.energy = energy  

        self.uniqueMFEs.append(get_uniqueMFE(self.antPopulation))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.bestAnt = min(self.antPopulation, key=lambda ant: ant.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.antPopulation)},best ant index:{self.bestAnt.index},best MFE:{self.bestAnt.energy},average MFE:{self.average_MFEs[-1]}")
        
        
    def get_first_minMFE(self):
        self.bestAnt = min(self.antPopulation, key=lambda ant: ant.energy)
        return  self.bestAnt.energy
        

    def iter_once(self):

        # Update pheromones
        self.pheromones = update_global_pheromones(self)
        
        self.newAntPopulation = [Ant(self) for _ in range(self.numAnts)]
        construct_ant_path(self)
        for newAnt in self.newAntPopulation:
            newAnt.combineSeq(self)
        
        self.energies = self.pool.map(partial(evaluate_ant),self.newAntPopulation)
        
        for ant, energy in zip(self.newAntPopulation, self.energies):
            ant.energy = energy  

        self.best_newAnt = min(self.newAntPopulation, key=lambda ant: ant.energy)
        self.best_newAnt_index = self.newAntPopulation.index(self.best_newAnt)

        # If the lowest value of the original ant colony is less than the lowest value of the new ant colony, 
        # replace the lowest individual in the new ant colony with the lowest individual in the original ant colony.
        if self.bestAnt.energy < self.best_newAnt.energy:
            self.newAntPopulation[self.best_newAnt_index] = self.bestAnt
            self.energies[self.best_newAnt_index] = self.bestAnt.energy
    
        # Replace the old generation with the offspring 
        self.antPopulation[:] = self.newAntPopulation

        self.uniqueMFEs.append(get_uniqueMFE(self.antPopulation))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        if self.minimumMFEs[-1] != min(self.uniqueMFEs[-1]):
            self.end_time = time.time() # Update convergence end time
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        if self.generation >= self.maxIterations-10:
            self.all_ants.extend(sorted(self.antPopulation, key=lambda ind: ind.energy))
            self.all_ants = remove_duplicate_ant(self.all_ants)
        
        self.generation = self.generation + 1
        self.bestAnt = min(self.antPopulation, key=lambda ant: ant.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.antPopulation)},best ant index:{self.bestAnt.index},best MFE:{self.bestAnt.energy},average MFE:{self.average_MFEs[-1]}")
        
        if self.generation == self.maxIterations:
            self.final_ants = remove_duplicate_ant(self.antPopulation)
            self.converge_time = self.end_time  - self.start_time
        return self.uniqueMFEs,self.average_MFEs,self.minimumMFEs,self.final_ants,self.all_ants,self.converge_time


def main(gene20ntInfo,insulatorSeq,repeatSeq,maxIterations,numAnts,alpha, beta, Global_evaporationRate,Single_evaporationRate,q0):
    aco_instance = ACO(gene20ntInfo,insulatorSeq,repeatSeq,numAnts,maxIterations,alpha, beta, Global_evaporationRate,Single_evaporationRate,q0)
    for _ in range(int(maxIterations-1)):
        MFEs,Average_MFEs,minimumMFEs,final_ants,all_ants,converge_time = aco_instance.iter_once()
    return MFEs,Average_MFEs,minimumMFEs,final_ants,all_ants,converge_time
        
