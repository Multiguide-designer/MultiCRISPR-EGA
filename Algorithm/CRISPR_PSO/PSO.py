from functools import partial  
import random
import subprocess
import multiprocessing
import copy
import numpy as np
import time



class Particle:
    def __init__(particle,self):
        particle.position = [random.randint(1, max_position) for max_position in self.GeneMaxIndices]
        particle.velocity = np.zeros(len(self.GeneMaxIndices))
        particle.best_position = particle.position[:]
        particle.best_fitness = float('inf')

    def combineSeq(particle,self):
        sequence_parts = [self.gene20ntInfo[gene][position - 1] for gene, position in zip(list(self.gene20ntInfo.keys()),  particle.position)]
        particle.seq = self.RIBOJ + self.REPEAT + self.REPEAT.join(sequence_parts) + self.REPEAT


def run_RNAfold(sequence):
    cmd = ["RNAfold", "--noPS"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def evaluate_particle(particle):
    output = run_RNAfold(particle.seq)
    if output:
        energy = float(output.split(" ")[-1].replace("(", "").replace(")", "").strip())
    else:
        energy = 0
    return energy


def update_velocity(self):
    for particle in self.new_particles:
        for i in range(len(particle.velocity)):
            r1 = random.random()
            r2 = random.random()
            cognitive = self.cognitive_parameter * r1 * (particle.best_position[i] - particle.position[i])
            social = self.social_parameter * r2 * (self.global_best_position[i] - particle.position[i])
            particle.velocity[i] = self.inertia_weight * particle.velocity[i] + cognitive + social


def update_position(self):
    for particle in self.new_particles:
        for i in range(len(particle.position)):
            particle.position[i] = max(1, min(round(particle.position[i] + particle.velocity[i]), self.GeneMaxIndices[i]))

    
def update_particle(self):
    for particle in self.particles:
        if particle.fitness < particle.best_fitness:
            particle.best_position = particle.position[:]
            particle.best_fitness = particle.fitness

        if particle.fitness < self.global_best_fitness:
            self.global_best_position = particle.position[:]
            self.global_best_fitness = particle.fitness
            
    self.new_particles = copy.deepcopy(self.particles)
    update_velocity(self)
    update_position(self)
    
    return self.new_particles


def get_uniqueMFE(population):
    unique_MFEs = []
    seen_indices = set()
    for individual in population:
        index_tuple = tuple(individual.position) 
        if index_tuple not in seen_indices:
            seen_indices.add(index_tuple)
            unique_MFEs.append(individual.fitness)
    return unique_MFEs


def remove_duplicate_particle(particles):
    unique_particles = []
    seen_indices = set()
    for particle in particles:
        position_tuple = tuple(particle.position)  
        if position_tuple not in seen_indices:
            seen_indices.add(position_tuple)
            unique_particles.append(particle)
    return unique_particles
    

class PSO:
    def __init__(self, gene20ntInfo,insulatorSeq,repeatSeq,particleSize,maxIterations,inertia_weight=0.9,cognitive_parameter=1.5,social_parameter=2.0):
        # PSO parameters
        self.particleSize = particleSize
        self.inertia_weight = inertia_weight  
        self.cognitive_parameter = cognitive_parameter 
        self.social_parameter = social_parameter
        self.generation = 1
        self.maxIterations = maxIterations 
        
        # Gene sequence parameters
        self.gene20ntInfo = gene20ntInfo
        self.GeneMaxIndices = [len(spacers) for spacers in gene20ntInfo.values()]
        
        self.RIBOJ = insulatorSeq
        self.REPEAT = repeatSeq

        random.seed(42)
        self.num_dimensions = len(self.GeneMaxIndices)
        
        # Collect the fitness value of each generation, 
        # as well as the average and minimum fitness value of each generation group, 
        # and the unique individuals of the last generation
        self.uniqueMFEs = []
        self.average_MFEs = []
        self.minimumMFEs = []
        self.final_particles = []
        self.all_particles = []
        
        self.converge_time = 0
        self.start_time = time.time()
        self.end_time = self.start_time
        
        # Create a pool of worker processes
        self.pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        self.particles = [Particle(self) for _ in range(self.particleSize)]
        for particle in self.particles:
            particle.combineSeq(self)

        # Parallelize the calculation of MFE for each individual in the population
        self.energies = self.pool.map(partial(evaluate_particle),self.particles)
        
        # Assign fitness values to particles
        for particle, energy in zip(self.particles, self.energies):
            particle.fitness = energy  
        
        self.uniqueMFEs.append(get_uniqueMFE(self.particles))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0])))  
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
        
        self.best_particle = min(self.particles, key=lambda particle: particle.fitness)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size::{len(self.particles)},best individual index:{self.best_particle.position},best MFE:{self.best_particle.fitness},Average MFE:{self.average_MFEs[-1]}")
        
    def iter_once(self):
       
        self.global_best_position = self.best_particle.position
        self.global_best_fitness = self.best_particle.fitness
        
        # Update particle swarm
        self.new_particles = update_particle(self)
        for particle in self.new_particles:
            particle.combineSeq(self)
        
        # Use adaptive PSO
        self.inertia_weight *= 0.99  
        if random.random() < 0.2:  # update weight
            self.inertia_weight = max(0.4, min(1.0, self.inertia_weight))
            
        self.energies = self.pool.map(partial(evaluate_particle),self.new_particles)
        
        # Assign fitness values to particles
        for particle, energy in zip(self.new_particles, self.energies):
            particle.fitness = energy  

        self.best_newParticle = min(self.new_particles, key=lambda ind: ind.fitness)
        self.best_newParticle_position = self.new_particles.index(self.best_newParticle)

        # If the lowest value of the parent is less than the lowest value of the offspring, r
        # eplace the lowest individual in the offspring with the lowest individual in the parent.
        if self.best_particle.fitness < self.best_newParticle.fitness:
            self.new_particles[self.best_newParticle_position] = self.best_particle
            self.energies[self.best_newParticle_position] = self.best_particle.fitness
            
        # Replace the old generation with the new_particles
        self.particles[:] = self.new_particles

        self.uniqueMFEs.append(get_uniqueMFE(self.particles))
        self.average_MFEs.append(np.mean(np.array([energy for energy in self.uniqueMFEs[-1] if energy != 0]))) 
        if self.minimumMFEs[-1] != min(self.uniqueMFEs[-1]):
            self.end_time = time.time() # Update convergence end time
        self.minimumMFEs.append(min(self.uniqueMFEs[-1]))
           
        if self.generation >= self.maxIterations-10:
            self.all_particles.extend(sorted(self.particles, key=lambda ind: ind.fitness))
            self.all_particles = remove_duplicate_particle(self.all_particles)
           
        self.generation = self.generation + 1
        self.best_particle = min(self.particles, key=lambda particle: particle.fitness)
        print(
            f"Generation {self.generation}/{self.maxIterations},population size::{len(self.particles)},best individual index:{self.best_particle.position},best MFE:{self.best_particle.fitness},Average MFE:{self.average_MFEs[-1]}")
        
        if self.generation == self.maxIterations:
            self.final_particles = remove_duplicate_particle(self.particles)
            self.converge_time = self.end_time  - self.start_time
        return self.uniqueMFEs,self.average_MFEs, self.minimumMFEs,self.final_particles,self.all_particles,self.converge_time
    
        
def main(gene20ntInfo,insulatorSeq,repeatSeq,maxIterations, particleSize,inertia_weight,cognitive_parameter,social_parameter):
    pso_instance = PSO(gene20ntInfo,insulatorSeq,repeatSeq,particleSize,maxIterations,inertia_weight,cognitive_parameter,social_parameter)
    for _ in range(int(maxIterations-1)):
        MFEs,Average_MFEs,minimumMFEs,final_particles,all_particles,converge_time = pso_instance.iter_once()
    return MFEs,Average_MFEs,minimumMFEs,final_particles,all_particles,converge_time

