from hyperopt import hp, fmin, tpe, Trials
import EGA
import os,csv


current_directory = os.path.dirname(__file__)
parent_directory = os.path.dirname(current_directory)
grandparent_directory = os.path.dirname(parent_directory)
# Extract spacer information
geneSpacerInfo = {}
gene_spacer_file = os.path.join(grandparent_directory, 'Gene_Info','Gene_Info8','spacer.csv')
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

space = {
    'crossover': hp.choice('alpha', [0.5, 0.6, 0.7,0.8,0.9]),  
    'mutation': hp.choice('beta', [0.05, 0.1,0.15,0.2,0.25,0.3]),   
    'tournamentSize': hp.choice('tournamentSize',[2,3,4])  
}

def objective(params):
    crossover = params['crossover']
    mutation = params['mutation']
    tournamentSize= params['tournamentSize']
    
    maxIterations = 100  
    populationSize = 1000 
    EGA_MFEs,EGA_Average_MFEs,EGA_minimumMFEs,EGA_final_solutions,EGA_all_solutions,EGA_converge_time= EGA.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations,populationSize,crossover,mutation,tournamentSize)
    minimum_energy = EGA_minimumMFEs[-1]
    final_average_MFE = EGA_Average_MFEs[-1]
    target_metric = 0.1*final_average_MFE + 0.9*minimum_energy
    return {'loss': target_metric,'params': params,'status': 'ok'} 

if __name__ == '__main__':
    trials = Trials()

    best = fmin(objective, space, algo=tpe.suggest, max_evals=50, trials=trials)

    print("Optimal hyperparameter configuration:", best)

    # Get the optimal N parameter combinations
    N = 20 
    best_trials = sorted(trials.results, key=lambda x: x['loss'])[:N]
    print("Optimal {} parameter combinations:".format(N))
    for trial in best_trials:
        print(trial['params'])
        
    with open("recordEGA.txt", "w") as f:
        f.write("Optimal hyperparameter configuration: {}\n".format(best))
        f.write("The optimal 20 parameter combinations:\n")
        for trial in best_trials:
            f.write(str(trial['params']) + "\n")
