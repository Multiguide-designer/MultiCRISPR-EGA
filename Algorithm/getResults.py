import sys
import os

current_directory = os.path.dirname(__file__)
parent_directory = os.path.dirname(current_directory)

from CRISPR_EGA import EGA
from CRISPR_ACO import ACOmain
from CRISPR_PSO import PSO
from CRISPR_SA import SAmain
from CRISPR_RandomSearch import RSmain
import time,os
import csv
from datetime import datetime


if __name__ == '__main__':
    maxIterations = 100
    populationSize = 1000
    numAnts = 1000
    particleSize = 1000
    numSolutions = 1000

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    folder_name = f"Result_{timestamp}"
    folder_path = os.path.join(current_directory, 'Result_all',folder_name)
    os.makedirs(folder_path)

    # Extract spacer information
    geneSpacerInfo = {}
    gene_spacer_file = os.path.join(parent_directory, 'Gene_Info','Gene_Info8','spacer.csv')
    with open(gene_spacer_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row['geneName']
            spacer_sequence = row['spacer']

            if gene_name in geneSpacerInfo:
                geneSpacerInfo[gene_name].append(spacer_sequence)
            else:
                geneSpacerInfo[gene_name] = [spacer_sequence]

    #prompter and direct repeat seq
    insulatorSeq = 'AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA'.replace('T','U')
    repeatSeq = 'GUUUGAGAGUUGUGUAAUUUAAGAUGGAUCUCAAAC' 
    
    #EGA
    print('Running the Elite Genetic Algorithm...')
    crossover = 0.8
    mutation = 0.05
    tournamentSize = 3
    start_time = time.time()
    EGA_MFEs,EGA_Average_MFEs,EGA_minimumMFEs,EGA_final_solutions,EGA_all_solutions,EGA_converge_time= EGA.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations,populationSize,crossover,mutation,tournamentSize)
    end_time = time.time()
    EGA_time = end_time - start_time
    print('\n')

    #ACO
    print('Running the Ant Colony Algorithm...')
    alpha = 3.0 # Importance of pheromones
    beta = 2.5 # Importance of Heuristic information
    Global_evaporationRate = 0.3
    Single_evaporationRate = 0.5
    q0 = 0.1 # Exploitation and exploration-biased control ratios
    start_time = time.time()
    ACO_MFEs,ACO_Average_MFEs,ACO_minimumMFEs,ACO_final_solutions,ACO_all_solutions,ACO_converge_time = ACOmain.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations, numAnts, alpha, beta, Global_evaporationRate, Single_evaporationRate, q0)
    end_time = time.time()
    ACO_time = end_time - start_time
    print('\n')


    #PSO
    print('Running the Particle Swarm Algorithm...')
    inertia_weight = 1.8 
    cognitive_parameter = 1.6
    social_parameter = 1.6
    start_time = time.time()
    PSO_MFEs,PSO_Average_MFEs,PSO_minimumMFEs,PSO_final_solutions,PSO_all_solutions,PSO_converge_time= PSO.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations, particleSize,inertia_weight,cognitive_parameter,social_parameter)
    end_time = time.time()
    PSO_time = end_time - start_time
    print('\n')

    #SA
    print('Running the Simulated Annealing...')
    cooling_rate = 0.3
    start_time = time.time()
    SA_MFEs,SA_Average_MFEs,SA_minimumMFEs,SA_final_solutions,SA_all_solutions,SA_converge_time = SAmain.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations,numSolutions,cooling_rate)
    end_time = time.time()
    SA_time = end_time - start_time
    print('\n')
    
    # RS
    print('Running the Random Search...')
    start_time = time.time()
    RS_MFEs,RS_Average_MFEs,RS_minimumMFEs,RS_final_solutions,RS_all_solutions,RS_converge_time = RSmain.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations,numSolutions)
    end_time = time.time()
    RS_time = end_time - start_time
    print('\n')
    
    # EGA results saved as CSV file
    EGA_results = {'MFEs': EGA_MFEs, 'Final Average MFE': EGA_Average_MFEs, 'Minimum MFEs': EGA_minimumMFEs}
    EGA_filename = f"EGA_MFE_result.csv"
    EGA_filepath = os.path.join(folder_path, EGA_filename)
    with open(EGA_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['MFEs'] + [''] * (populationSize-1) + ['Average MFEs', 'Minimum MFEs'])
        for i in range(maxIterations):
            new_row = EGA_MFEs[i] 
            new_row.append(EGA_Average_MFEs[i]) 
            new_row.append(EGA_minimumMFEs[i])  
            writer.writerow(new_row)
        writer.writerow(['Time taken:'+str(EGA_time) + 's'])
        writer.writerow(['Obtain the lowest MFE sequence time:'+str(EGA_converge_time) + 's'])
    print("EGA MFE information has been saved to",EGA_filepath)
    EGA_available_solutions_final = os.path.join(folder_path, 'EGA_Available_Sol_final.csv')
    with open(EGA_available_solutions_final, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in EGA_final_solutions:
            writer.writerow([solution.seq, solution.fitness])
    print("The last iteration result of the EGA has been saved to",EGA_available_solutions_final)
    EGA_available_solutions_all = os.path.join(folder_path, 'EGA_Available_Sol_all.csv')
    with open(EGA_available_solutions_all, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in EGA_all_solutions:
            writer.writerow([solution.seq, solution.fitness])
    print("The sequence of the last ten iterations found by the EGA has been saved to",EGA_available_solutions_all)
    print('\n')
            
    # ACO results saved as CSV file
    aco_results = {'MFEs': ACO_MFEs, 'Final Average MFE': ACO_Average_MFEs, 'Minimum MFEs': ACO_minimumMFEs}
    aco_filename = f"ACO_MFE_result.csv"
    aco_filepath = os.path.join(folder_path, aco_filename)
    with open(aco_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['MFEs'] + [''] * (populationSize-1) + ['Final Average MFE', 'Minimum MFEs'])
        for i in range(maxIterations):
            new_row = ACO_MFEs[i] 
            new_row.append(ACO_Average_MFEs[i]) 
            new_row.append(ACO_minimumMFEs[i])  
            writer.writerow(new_row)
        writer.writerow(['Time taken:'+str(ACO_time) + 's'])
        writer.writerow(['Obtain the lowest MFE sequence time:'+str(ACO_converge_time) + 's'])
    print("ACO MFE information has been saved to",aco_filepath)
    aco_available_solutions_final = os.path.join(folder_path, 'ACO_Available_Sol_final.csv')
    with open(aco_available_solutions_final, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in ACO_final_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The last iteration result of the ACO has been saved to",aco_available_solutions_final)
    aco_available_solutions_all = os.path.join(folder_path, 'ACO_Available_Sol_all.csv')
    with open(aco_available_solutions_all, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in ACO_all_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The sequence of the last ten iterations found by the ACO has been saved to",aco_available_solutions_all)
    print('\n')

    # PSO results saved as CSV file
    pso_results = {'MFEs': PSO_MFEs, 'Final Average MFE': PSO_Average_MFEs, 'Minimum MFEs': PSO_minimumMFEs}
    pso_filename = f"PSO_MFE_result.csv"
    pso_filepath = os.path.join(folder_path, pso_filename)
    with open(pso_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['MFEs'] + [''] * (populationSize-1) + ['Final Average MFE', 'Minimum MFEs'])
        for i in range(maxIterations):
            new_row = PSO_MFEs[i] 
            new_row.append(PSO_Average_MFEs[i]) 
            new_row.append(PSO_minimumMFEs[i])  
            writer.writerow(new_row)
        writer.writerow(['Time taken:'+str(PSO_time) + 's'])
        writer.writerow(['Obtain the lowest MFE sequence time:'+str(PSO_converge_time) + 's'])
    print("PSO MFE information has been saved to", pso_filepath)
    pso_available_solutions_final = os.path.join(folder_path, 'PSO_Available_Sol_final.csv')
    with open(pso_available_solutions_final, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in PSO_final_solutions:
            writer.writerow([solution.seq, solution.fitness])
    print("The last iteration result of the PSO has been saved to",pso_available_solutions_final)
    pso_available_solutions_all = os.path.join(folder_path, 'PSO_Available_Sol_all.csv')
    with open(pso_available_solutions_all, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in PSO_all_solutions:
            writer.writerow([solution.seq, solution.fitness])
    print("The sequence of the last ten iterations found by the PSO has been saved to",pso_available_solutions_all)
    print('\n')

    # SA results saved as CSV file
    sa_results = {'MFEs': SA_MFEs, 'Final Average MFE': SA_Average_MFEs, 'Minimum MFEs': SA_minimumMFEs}
    sa_filename = f"SA_MFE_result.csv"
    sa_filepath = os.path.join(folder_path, sa_filename)
    with open(sa_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['MFEs'] + [''] * (populationSize-1) + ['Final Average MFE', 'Minimum MFEs'])
        for i in range(maxIterations):
            new_row = SA_MFEs[i] 
            new_row.append(SA_Average_MFEs[i]) 
            new_row.append(SA_minimumMFEs[i])  
            writer.writerow(new_row)
        writer.writerow(['Time taken:'+str(SA_time) + 's'])
        writer.writerow(['Obtain the lowest MFE sequence time:'+str(SA_converge_time) + 's'])
    print("SA MFE information has been saved to", sa_filepath)
    sa_available_solutions_final = os.path.join(folder_path, 'SA_Available_Sol_final.csv')
    with open(sa_available_solutions_final, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in SA_final_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The last iteration result of the SA has been saved to",sa_available_solutions_final)
    sa_available_solutions_all = os.path.join(folder_path, 'SA_Available_Sol_all.csv')
    with open(sa_available_solutions_all, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in SA_all_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The sequence of the last ten iterations found by the SA has been saved to",sa_available_solutions_all)
    print('\n')

    # RS results saved as CSV file
    rs_results = {'MFEs': RS_MFEs, 'Final Average MFE': RS_Average_MFEs, 'Minimum MFEs': RS_minimumMFEs}
    rs_filename = f"RS_MFE_result.csv"
    rs_filepath = os.path.join(folder_path, rs_filename)
    with open(rs_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['MFEs'] + [''] * (populationSize-1) + ['Final Average MFE', 'Minimum MFEs'])
        for i in range(maxIterations):
            new_row = RS_MFEs[i] 
            new_row.append(RS_Average_MFEs[i]) 
            new_row.append(RS_minimumMFEs[i])  
            writer.writerow(new_row)
        writer.writerow(['Time taken:'+str(RS_time) + 's'])
        writer.writerow(['Obtain the lowest MFE sequence time:'+str(RS_converge_time) + 's'])
    print("RS MFE information has been saved to", rs_filepath)
    rs_available_solutions_final = os.path.join(folder_path, 'RS_Available_Sol_final.csv')
    with open(rs_available_solutions_final, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in RS_final_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The last iteration result of the RS has been saved to",rs_available_solutions_final)
    rs_available_solutions_all = os.path.join(folder_path, 'RS_Available_Sol_all.csv')
    with open(rs_available_solutions_all, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for solution in RS_all_solutions:
            writer.writerow([solution.seq, solution.energy])
    print("The sequence of the last ten iterations found by the RS has been saved to",rs_available_solutions_all)
    print('\n')
