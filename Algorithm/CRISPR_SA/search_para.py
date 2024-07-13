from hyperopt import hp, fmin, tpe, Trials
import SAmain
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
    'cooling_rate': hp.choice('cooling_rate', [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9]),  
}

def objective(params):
    cooling_rate = params['cooling_rate']

    maxIterations = 100 
    numSolutions = 1000

    SA_MFEs,SA_Average_MFEs,SA_minimumMFEs,SA_final_solutions,SA_all_solutions,SA_converge_time = SAmain.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations,numSolutions,cooling_rate)
    minimum_energy = SA_minimumMFEs[-1]
    final_average_MFE = SA_Average_MFEs[-1]
    target_metric = 0.1*final_average_MFE + 0.9*minimum_energy
    return {'loss': target_metric,'params': params,'status': 'ok'} 

if __name__ == '__main__':
    trials = Trials()

    best = fmin(objective, space, algo=tpe.suggest, max_evals=100, trials=trials)

    print("Optimal hyperparameter configuration:", best)

    # Get the optimal N parameter combinations
    N = 20 
    best_trials = sorted(trials.results, key=lambda x: x['loss'])[:N]
    print("Optimal {} parameter combinations:".format(N))
    for trial in best_trials:
        print(trial['params'])
        
    with open("recordSA.txt", "w") as f:
        f.write("Optimal hyperparameter configuration: {}\n".format(best))
        f.write("The optimal 20 parameter combinations:\n")
        for trial in best_trials:
            f.write(str(trial['params']) + "\n")
        