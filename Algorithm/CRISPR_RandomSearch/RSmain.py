from functools import partial  
import random
import subprocess
import multiprocessing
import time
import numpy as np
import copy



class RandomSearchAgent:
    def __init__(Agent,self):
        Agent.index = [random.randint(1, max_index) for max_index in self.GeneMaxIndices]
    
    def combineSeq(Agent,self):
        sequence_parts = [self.gene20ntInfo[gene][index - 1] for gene, index in zip(list(self.gene20ntInfo.keys()), Agent.index)]
        Agent.seq = self.RIBOJ + self.REPEAT + self.REPEAT.join(sequence_parts) + self.REPEAT


def run_RNAfold(sequence):
    cmd = ["RNAfold", "--noPS"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def evaluate_Agent(Agent):
    output = run_RNAfold(Agent.seq)
    if output:
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
    else:
        energy = 0
    return energy


def get_uniqueMFE(population):
    unique_MFEs = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.index)  
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_MFEs.append(individual.energy)
    return unique_MFEs


def remove_duplicate_Agent(population):
    unique_Agents = []
    seen_indices = set()
    for Agent in population:
        index_tuple = tuple(Agent.index) 
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_Agents.append(Agent)
    return unique_Agents



class RandomSearch:
    def __init__(self,gene20ntInfo,insulatorSeq,repeatSeq,numAgents,maxIterations):
        # RS parameters
        self.numAgents = numAgents
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
        self.final_Agents = []
        self.all_Agents  =[]

        self.converge_time = 0
        self.start_time = time.time()
        self.end_time = self.start_time
        
        # Create a pool of worker processes
        self.pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        self.Agents = [RandomSearchAgent(self) for _ in range(self.numAgents)]
        for Agent in self.Agents:
            Agent.combineSeq(self)

        # Parallelize the calculation of MFE for each individual in the population
        self.energies = self.pool.map(partial(evaluate_Agent),self.Agents)
    
        for Agent, energy in zip(self.Agents, self.energies):
            Agent.energy = energy  

        self.uniqueMFEs.append(get_uniqueMFE(self.Agents))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.bestAgent = min(self.Agents, key=lambda Agent: Agent.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size::{len(self.Agents)},Best agent index:{self.bestAgent.index},best MFE:{self.bestAgent.energy},average MFE:{self.average_MFEs[-1]}")

        
    def get_first_minMFE(self):
        self.bestAgent = min(self.Agents, key=lambda Agent: Agent.energy)
        return  self.bestAgent.energy
        
        
    def iter_once(self):

        self.newAgents = [RandomSearchAgent(self) for _ in range(self.numAgents)]

        for newAgent in self.newAgents:
            newAgent.combineSeq(self)
        
        self.energies = self.pool.map(partial(evaluate_Agent),self.newAgents)
        
        for Agent, energy in zip(self.newAgents, self.energies):
            Agent.energy = energy  

        combined_agents = self.Agents + self.newAgents
        sorted_agents = sorted(combined_agents, key=lambda agent: agent.energy)
        self.newAgents = sorted_agents[:self.numAgents]

        self.best_newAgent = min(self.newAgents, key=lambda Agent: Agent.energy)
        self.best_newAgent_index = self.newAgents.index(self.best_newAgent)

        # If the lowest value of the original ant colony is less than the lowest value of the new ant colony, 
        # replace the lowest individual in the new ant colony with the lowest individual in the original ant colony
        if self.bestAgent.energy < self.best_newAgent.energy:
            self.newAgents[self.best_newAgent_index] = self.bestAgent
            self.energies[self.best_newAgent_index] = self.bestAgent.energy
            
        # Replace the old generation with the offspring 
        self.Agents[:] = self.newAgents

        self.uniqueMFEs.append(get_uniqueMFE(self.Agents))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        if self.minimumMFEs[-1] != min(self.uniqueMFEs[-1]):
            self.end_time = time.time() # Update convergence end time
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        if self.generation >= self.maxIterations-10:
            self.all_Agents.extend(sorted(self.Agents, key=lambda ind: ind.energy))
            self.all_Agents = remove_duplicate_Agent(self.all_Agents)
        
        self.generation = self.generation + 1
        self.bestAgent = min(self.Agents, key=lambda Agent: Agent.energy)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size:{len(self.Agents)},best agent index:{self.bestAgent.index},best MFE:{self.bestAgent.energy},average MFE:{self.average_MFEs[-1]}")

        if self.generation == self.maxIterations:
            self.final_Agents = remove_duplicate_Agent(self.Agents)
            self.converge_time = self.end_time  - self.start_time
        return self.uniqueMFEs,self.average_MFEs,self.minimumMFEs,self.final_Agents,self.all_Agents,self.converge_time


def main(gene20ntInfo,insulatorSeq,repeatSeq,maxIterations,numAgents):
    rs_instance = RandomSearch(gene20ntInfo,insulatorSeq,repeatSeq,numAgents,maxIterations)
    for _ in range(int(maxIterations-1)):
        MFEs,Average_MFEs,minimumMFEs,final_Agents,all_Agents,converge_time = rs_instance.iter_once()
    return MFEs,Average_MFEs,minimumMFEs,final_Agents,all_Agents,converge_time
        
     