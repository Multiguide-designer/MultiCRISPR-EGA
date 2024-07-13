from hyperopt import hp, fmin, tpe, Trials
import PSO
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
    'inertia_weight': hp.choice('inertia_weight', [1.2, 1.5, 1.8,2.1]),  
    'cognitive_parameter': hp.choice('cognitive_parameter', [1.2, 1.6, 2.0,2.5]), 
    'social_parameter': hp.choice('social_parameter', [1.2, 1.6, 2.0,2.5]),  
   }


def objective(params):
    inertia_weight = params['inertia_weight']
    cognitive_parameter = params['cognitive_parameter']
    social_parameter= params['social_parameter']
    
    maxIterations = 100  
    particleSize = 1000

    PSO_MFEs,PSO_Average_MFEs,PSO_minimumMFEs,PSO_final_solutions,PSO_all_solutions,PSO_converge_time= PSO.main(geneSpacerInfo,insulatorSeq,repeatSeq,maxIterations, particleSize,inertia_weight,cognitive_parameter,social_parameter)
    minimum_energy = PSO_minimumMFEs[-1]
    final_average_MFE = PSO_Average_MFEs[-1]
    target_metric = 0.1*final_average_MFE + 0.9*minimum_energy
    return {'loss': target_metric,'params': params,'status': 'ok'} 

if __name__ == '__main__':
    trials = Trials()

    best = fmin(objective, space, algo=tpe.suggest, max_evals=20, trials=trials)

    print("Optimal hyperparameter configuration:", best)

    # Get the optimal N parameter combinations
    N = 20 
    best_trials = sorted(trials.results, key=lambda x: x['loss'])[:N]
    print("Optimal {} parameter combinations:".format(N))
    for trial in best_trials:
        print(trial['params'])
        
    with open("recordPSO.txt", "w") as f:
        f.write("Optimal hyperparameter configuration: {}\n".format(best))
        f.write("The optimal 20 parameter combinations:\n")
        for trial in best_trials:
            f.write(str(trial['params']) + "\n")
