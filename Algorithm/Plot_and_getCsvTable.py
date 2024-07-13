import csv
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import math
import statsmodels.api as sm

def read_MFE_csv(filename):
    MFEs = []
    Average_MFEs = []
    minimumMFEs = []
    with open(filename, 'r', newline='') as file:
        reader = csv.reader(file)
        next(reader) 
        for row in reader:
            if len(row) != 1:
                MFEs.append(row[:-2])  
                Average_MFEs.append(row[-2])  
                minimumMFEs.append(row[-1]) 
            else:
                time_str = row[0].split(':')[-1].strip()  
                algorithm_time = float(time_str[:-1]) 
                next_row = next(reader, None)
                if next_row:
                    time_str = next_row[0].split(':')[-1].strip() 
                    converge_time = float(time_str[:-1]) 
                
    MFEs = [[float(num) for num in sublist] for sublist in MFEs]
    Average_MFEs = [float(num) for num in Average_MFEs]
    minimumMFEs = [float(num) for num in minimumMFEs]
    return MFEs, Average_MFEs, minimumMFEs,algorithm_time,converge_time 
    
def read_available_sols_csv(filename):
    available_sols = []
    available_MFEs = []
    with open(filename, 'r', newline='') as file:
        reader = csv.reader(file)
        next(reader) 
        for row in reader:
            available_sols.append(row[0])  
            available_MFEs.append(row[1]) 
    available_MFEs = [float(num) for num in available_MFEs]   
    return available_sols, available_MFEs
    
def process_data_all_iteration():
    MFEs = [ega_MFEs, aco_MFEs, pso_MFEs, sa_MFEs,rs_MFEs]
    std_values = [[np.std(lst) for lst in mfe_lst] for mfe_lst in MFEs]  #np.std function calculates the standard error by default. By default, ddof=0
    percentile_values = [[np.percentile(lst, [25, 75]) for lst in mfe_lst] for mfe_lst in MFEs]
    lower_bounds = [[percentile[0] for percentile in percentiles] for percentiles in percentile_values]
    upper_bounds = [[percentile[1] for percentile in percentiles] for percentiles in percentile_values]
    average_MFEs = [[EGA_Average_MFE,aco_Average_MFE,pso_Average_MFE,sa_Average_MFE,rs_Average_MFE] for EGA_Average_MFE,aco_Average_MFE,pso_Average_MFE,sa_Average_MFE,rs_Average_MFE in zip(ega_Average_MFEs, aco_Average_MFEs, pso_Average_MFEs, sa_Average_MFEs,rs_Average_MFEs)]
    minimum_MFEs = [ega_minimum_MFEs, aco_minimum_MFEs, pso_minimum_MFEs, sa_minimum_MFEs,rs_minimum_MFEs]

    # Convert the data to dataframe long format and add new variables "Iteration" and "Standard Deviation"
    df = pd.DataFrame(average_MFEs , columns=algorithm_name)
    df_long = df.melt(var_name="Algorithm", value_name="Average MFE")
    df_long.insert(0, "Iteration", df_long.groupby("Algorithm").cumcount() + 1)
    df_long["Standard Error"] = np.concatenate(std_values)
    df_long["Lower Bound"] = np.concatenate(lower_bounds)
    df_long["Upper Bound"] = np.concatenate(upper_bounds)
    df_long["Minimum MFE"] = np.concatenate(minimum_MFEs)
    print(df_long)
    return df_long

    
def process_data_final_MFE_forViolin():
    determine_data = []
    
    for comparision_algorithm in algorithm_name:
        if comparision_algorithm != 'EGA':
            determine_data.extend(['EGA'] * (len(ega_MFEs[-1])))
            determine_data.extend(['Other Method'] * len(eval(f"{comparision_algorithm.lower()}_MFEs[-1]")))
    
    algorithm_name_data = []
    for comparision_algorithm in algorithm_name:
        if comparision_algorithm != 'EGA':
            algorithm_name_data.extend([comparision_algorithm] * (len(ega_MFEs[-1]) + len(eval(f"{comparision_algorithm.lower()}_MFEs[-1]"))))
    
    final_MFE_data = []
    for comparision_algorithm in algorithm_name:
        if comparision_algorithm != 'EGA':
            final_MFE_data.extend(ega_MFEs[-1])
            final_MFE_data.extend(eval(f"{comparision_algorithm.lower()}_MFEs[-1]"))
            
    data = {
        'Is EGA': determine_data,
        'Comparison algorithm': algorithm_name_data,
        'Final MFE': final_MFE_data
    }

    df = pd.DataFrame(data)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    return df
    
    
def line_graph_average(df_allIteration,figure_path):
    plt.figure(figsize=(8, 6), dpi=300)
    
    palette = sns.color_palette("muted")
    palette[0], palette[3] = palette[3], palette[0]
    markers = {'EGA': '^', 'ACO': '*', 'PSO': 'D', 'SA': 'p', 'RS': 'o'}
    
   # Loop to draw polylines, error bands and special markers for each algorithm
    for idx, algorithm in enumerate(algorithm_name):
        group = df_allIteration[df_allIteration['Algorithm'] == algorithm]
        plt.plot(group['Iteration'], group['Average MFE'], label=algorithm, linestyle='-', linewidth=2, 
                 marker=markers[algorithm], markevery=(19,20), markersize=10,markeredgewidth=1,
                 markeredgecolor='white',color=palette[idx])
        
        # Use LOWESS (Locally Weighted Scatterplot Smoothing)
        lowess = sm.nonparametric.lowess

        # Smooth the lower and upper bounds
        smoothed_lower = lowess(group['Lower Bound'], group['Iteration'], frac=0.2)
        smoothed_upper = lowess(group['Upper Bound'], group['Iteration'], frac=0.2)
        plt.fill_between(group['Iteration'], smoothed_lower[:, 1], smoothed_upper[:, 1], color=palette[idx], alpha=0.2)
           
    plt.xlim(left=1)
   
    plt.xlabel("Iteration",fontsize=16)
    plt.ylabel("Average MFE",fontsize=16)

    plt.title('Average Minimum Free Energy(MFE) Comparison',fontsize=16)
    plt.legend(fontsize=16, loc='upper right')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.tight_layout() 

    os.makedirs(figure_path, exist_ok=True) 
    plt.savefig(os.path.join(figure_path, "Average_MFE_line_graph.png"))
    
    plt.close()
    print("Chart saved to", os.path.join(figure_path, "Average_MFE_line_graph.png"))
    return True


def line_graph_minimum(df_allIteration,figure_path):
    plt.figure(figsize=(8, 6), dpi=300)
    palette = sns.color_palette("muted")
    palette[0], palette[3] = palette[3], palette[0]
    markers = {'EGA': '^', 'ACO': '*', 'PSO': 'D', 'SA': 'p', 'RS': 'o'}
    
    for idx, algorithm in enumerate(algorithm_name):
        group = df_allIteration[df_allIteration['Algorithm'] == algorithm]
        plt.plot(group['Iteration'], group['Minimum MFE'], label=algorithm, linestyle='-', linewidth=2, 
                 marker=markers[algorithm], markevery=(19,20), markersize=10,markeredgewidth=1,
                 markeredgecolor='white',color=palette[idx])
    
    plt.xlim(left=1)
            
    plt.title('Minimum Minimum Free Energy(MFE) Comparison',fontsize=16)
    plt.xlabel('Iteration',fontsize=16)
    plt.ylabel('Minimum MFE',fontsize=16)
    plt.legend(fontsize=16, loc='upper right')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.tight_layout() 

    os.makedirs(figure_path, exist_ok=True) 
    plt.savefig(os.path.join(figure_path, "Minimum_MFE_line_graph.png"))
    
    plt.close()
    print("Chart saved to", os.path.join(figure_path, "Minimum_MFE_line_graph.png"))
    return True
    
    
def violin_graph(df_final_MFE,figure_path):
    plt.figure(figsize=(8, 6), dpi=300)

    sns.violinplot(data=df_final_MFE, x="Comparison algorithm", y="Final MFE", hue="Is EGA",
               split=True, inner="quart", fill=False,
               palette={"EGA": "g", "Other Method": ".35"},
               density_norm='area',cut=0)
    
    plt.title('Final Iteration Minimum Free Energy Comparison',fontsize=16)
    plt.xlabel('Algotithm',fontsize=16)
    plt.ylabel('MFE',fontsize=16)
    plt.legend(fontsize=16)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.tight_layout() #是 Matplotlib 中的一个函数，用于自动调整图表布局，以确保图表元素之间的间距合适，避免它们重叠或太靠近边界。

    os.makedirs(figure_path, exist_ok=True) 
    plt.savefig(os.path.join(figure_path, "Final_Iteration_MFE_violin_graph.png"))
    
    plt.close()
    print("Chart saved to", os.path.join(figure_path, "Final_Iteration_MFE_violin_graph.png"))
    
    

def convergence_time_barplot(convergence_times, algorithm_names, figure_path):
    plt.figure(figsize=(8, 6), dpi=300)
    
    palette = sns.color_palette("muted")
    palette[0], palette[3] = palette[3], palette[0]
    ax = sns.barplot(x=algorithm_names, y=convergence_times, palette=palette)
    plt.title("Convergence Time Comparison",fontsize=16)
    plt.xlabel("Algorithm",fontsize=16)
    plt.ylabel("Convergence Time (seconds)",fontsize=16)
    plt.xticks(rotation=45, ha='right')  
    
    # Add a label to the top of each bar
    for i, v in enumerate(convergence_times):
        ax.text(i, v + 0.05, str(round(v, 2)), ha='center', va='bottom')
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
        
    plt.tight_layout()  
    
    os.makedirs(figure_path, exist_ok=True)
    plt.savefig(os.path.join(figure_path, "Convergence_Time_Barplot.png"))

    plt.close()
    print("Chart saved to", os.path.join(figure_path, "Convergence_Time_Barplot.png"))


def get_infoTable_and_validSeq(validSeq_path,comparisonTable_path):
    initialSeqNum = []
    finalSeqNum = []
    minimum_MFE = []
    Average_MFE = []
    Time = []
    for algorithm in algorithm_name:
        initialSeqNum.append(len(eval(f"{algorithm.lower()}_MFEs[0]")))
        finalSeqNum.append(len(eval(f"{algorithm.lower()}_MFEs[-1]")))
        minimum_MFE.append(eval(f"{algorithm.lower()}_minimum_MFEs[-1]"))
        Average_MFE.append(eval(f"{algorithm.lower()}_Average_MFEs[-1]"))
        Time.append(eval(f"{algorithm.lower()}_time"))
        
    # Calculate the number and proportion of effective sequences
    Effective_energy = rs_minimum_MFEs[-1] #The randomly generated minimum MFE is set as the effective energy, and sequences not higher than this energy are valid sequences.
    all_effective_sequences = set()
    algorithm_effective_counts = {}
    
    for algorithm in algorithm_name:
        available_sols = eval(f'{algorithm.lower()}_available_all_sols')
        available_MFEs = eval(f'{algorithm.lower()}_available_all_MFEs')
        effective_sequences = set()
        for i, mfe in enumerate(available_MFEs):
            if mfe <= Effective_energy:
                effective_sequences.add((available_sols[i],mfe,algorithm))
                all_effective_sequences.add((available_sols[i],mfe,algorithm))    
        algorithm_effective_counts[algorithm] = len(effective_sequences)
    
    total_effective_count = len(all_effective_sequences)
        
    all_effective_sequences = list(all_effective_sequences)        
    all_effective_sequences.sort(key=lambda x: x[1]) 
    with open(validSeq_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow( ['Sequence','MFE','Method'])
        for sequence_info in all_effective_sequences:
            writer.writerow(sequence_info)
    print("Valid sequence saved to",validSeq_path)
    
    # Calculate the proportion of the top 10%, top 30% and top 50% effective sequences
    top_10_percent_index = math.ceil(total_effective_count * 0.1)
    top_30_percent_index = math.ceil(total_effective_count * 0.3)
    top_50_percent_index = math.ceil(total_effective_count * 0.5)
    
    top_10_percent_mfe = all_effective_sequences[top_10_percent_index][1]
    top_30_percent_mfe = all_effective_sequences[top_30_percent_index][1]
    top_50_percent_mfe = all_effective_sequences[top_50_percent_index][1]
    
    top_10_percent_sequences = [seq for seq in all_effective_sequences if seq[1] <= top_10_percent_mfe]
    top_30_percent_sequences = [seq for seq in all_effective_sequences if seq[1] <= top_30_percent_mfe]
    top_50_percent_sequences = [seq for seq in all_effective_sequences if seq[1] <= top_50_percent_mfe]
    
    top_10_percent_counts = {algorithm: sum(1 for seq in top_10_percent_sequences if seq[2] == algorithm) for algorithm in algorithm_name}
    top_30_percent_counts = {algorithm: sum(1 for seq in top_30_percent_sequences if seq[2] == algorithm) for algorithm in algorithm_name}
    top_50_percent_counts = {algorithm: sum(1 for seq in top_50_percent_sequences if seq[2] == algorithm) for algorithm in algorithm_name}
    
    with open(comparisonTable_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow( ['']  + algorithm_name)
        writer.writerow(['Minimum MFE '] + minimum_MFE)
        writer.writerow(['Average MFE'] + Average_MFE)
        writer.writerow(['Effective sequence count'] + list(algorithm_effective_counts.values()))
        if len(top_10_percent_sequences)!=0:
            writer.writerow(['Top 10% effective sequence percentage'] + [top_10_percent_counts[alg] / len(top_10_percent_sequences) for alg in algorithm_name])
        else:
            writer.writerow(['Top 10% effective sequence percentage'] + [0 for alg in algorithm_name])
        if len(top_30_percent_sequences)!=0:
            writer.writerow(['Top 30% effective sequence percentage'] + [top_30_percent_counts[alg] / len(top_30_percent_sequences) for alg in algorithm_name])
        else:
            writer.writerow(['Top 30% effective sequence percentage'] + [0 for alg in algorithm_name])
        if len(top_10_percent_sequences)!=0:
            writer.writerow(['Top 50% effective sequence percentage'] + [top_50_percent_counts[alg] / len(top_50_percent_sequences) for alg in algorithm_name])
        else:
            writer.writerow(['Top 50% effective sequence percentage'] + [0 for alg in algorithm_name])
        writer.writerow(['Time taken'] + Time)
    print("Algorithm comparison table information has been saved to",comparisonTable_path)
    return 
    

    
current_directory = os.path.dirname(__file__)
folder_path = os.path.join(current_directory,'Result_all')
# Use the modification time of the files to sort and get the time of the last file that was recently modified.
files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]
files.sort(key=lambda x: os.path.getmtime(x))
if files:
    last_file_by_time = files[-1]
    print("Path to the last modified file:", last_file_by_time)
    last_modified_time = os.path.basename(last_file_by_time)[len('Result_'):]
else:
    print("No files found in the directory.")
    
folder_path = os.path.join(folder_path,last_file_by_time)


# Read CSV file
ega_MFEfilename = "EGA_MFE_result.csv"
ega_MFEfilepath = os.path.join(folder_path, ega_MFEfilename)
ega_MFEs, ega_Average_MFEs, ega_minimum_MFEs,ega_time,ega_converge_time = read_MFE_csv(ega_MFEfilepath)
ega_final_Solfilename = "EGA_Available_Sol_final.csv"
ega_final_Solfilepath = os.path.join(folder_path, ega_final_Solfilename)
ega_available_final_sols, ega_available_final_MFEs = read_available_sols_csv(ega_final_Solfilepath)
ega_all_Solfilename = "EGA_Available_Sol_all.csv"
ega_all_Solfilepath = os.path.join(folder_path, ega_all_Solfilename)
ega_available_all_sols, ega_available_all_MFEs = read_available_sols_csv(ega_all_Solfilepath)

aco_MFEfilename = "ACO_MFE_result.csv"
aco_MFEfilepath = os.path.join(folder_path, aco_MFEfilename)
aco_MFEs, aco_Average_MFEs, aco_minimum_MFEs,aco_time,aco_converge_time = read_MFE_csv(aco_MFEfilepath)
aco_final_Solfilename = "ACO_Available_Sol_final.csv"
aco_final_Solfilepath = os.path.join(folder_path, aco_final_Solfilename)
aco_available_final_sols, aco_available_final_MFEs = read_available_sols_csv(aco_final_Solfilepath)
aco_all_Solfilename = "ACO_Available_Sol_all.csv"
aco_all_Solfilepath = os.path.join(folder_path, aco_all_Solfilename)
aco_available_all_sols, aco_available_all_MFEs = read_available_sols_csv(aco_all_Solfilepath)

pso_MFEfilename = "PSO_MFE_result.csv"
pso_MFEfilepath = os.path.join(folder_path, pso_MFEfilename)
pso_MFEs, pso_Average_MFEs, pso_minimum_MFEs,pso_time,pso_converge_time = read_MFE_csv(pso_MFEfilepath)
pso_final_Solfilename = "PSO_Available_Sol_final.csv"
pso_final_Solfilepath = os.path.join(folder_path, pso_final_Solfilename)
pso_available_final_sols, pso_available_final_MFEs = read_available_sols_csv(pso_final_Solfilepath)
pso_all_Solfilename = "PSO_Available_Sol_all.csv"
pso_all_Solfilepath = os.path.join(folder_path, pso_all_Solfilename)
pso_available_all_sols, pso_available_all_MFEs = read_available_sols_csv(pso_all_Solfilepath)

sa_MFEfilename = "SA_MFE_result.csv"
sa_MFEfilepath = os.path.join(folder_path, sa_MFEfilename)
sa_MFEs, sa_Average_MFEs, sa_minimum_MFEs,sa_time,sa_converge_time = read_MFE_csv(sa_MFEfilepath)
sa_final_Solfilename = "SA_Available_Sol_final.csv"
sa_final_Solfilepath = os.path.join(folder_path, sa_final_Solfilename)
sa_available_final_sols, sa_available_final_MFEs = read_available_sols_csv(sa_final_Solfilepath)
sa_all_Solfilename = "SA_Available_Sol_all.csv"
sa_all_Solfilepath = os.path.join(folder_path, sa_all_Solfilename)
sa_available_all_sols, sa_available_all_MFEs = read_available_sols_csv(sa_all_Solfilepath)

rs_MFEfilename = "RS_MFE_result.csv"
rs_MFEfilepath = os.path.join(folder_path, rs_MFEfilename)
rs_MFEs, rs_Average_MFEs, rs_minimum_MFEs,rs_time,rs_converge_time = read_MFE_csv(rs_MFEfilepath)
rs_final_Solfilename = "RS_Available_Sol_final.csv"
rs_final_Solfilepath = os.path.join(folder_path, rs_final_Solfilename)
rs_available_final_sols, rs_available_final_MFEs = read_available_sols_csv(rs_final_Solfilepath)
rs_all_Solfilename = "RS_Available_Sol_all.csv"
rs_all_Solfilepath = os.path.join(folder_path, rs_all_Solfilename)
rs_available_all_sols, rs_available_all_MFEs = read_available_sols_csv(rs_all_Solfilepath)

algorithm_name = ['EGA', 'ACO', 'PSO', 'SA','RS']

# Data preprocessing
df_allIteration = process_data_all_iteration()
df_final_MFE_violin = process_data_final_MFE_forViolin()
convergence_times = [ega_converge_time, aco_converge_time, pso_converge_time, sa_converge_time, rs_converge_time]

# Get figures
folder_path_fig = os.path.join(current_directory,'Result_fig_and_table')
folder_path_fig = os.path.join(folder_path_fig,last_modified_time)
line_graph_average(df_allIteration,folder_path_fig)
line_graph_minimum(df_allIteration,folder_path_fig)
violin_graph(df_final_MFE_violin,folder_path_fig)
convergence_time_barplot(convergence_times, algorithm_name, folder_path_fig)

# Get the method comparison table and the valid sequences of all methods
folder_path_table = os.path.join(current_directory, 'Result_fig_and_table',last_modified_time)
comparisonTable_filename = "infoTable.csv"
comparisonTable_path = os.path.join(folder_path_table,comparisonTable_filename)
folder_path_seq = os.path.join(current_directory, 'Result_fig_and_table',last_modified_time )
validSeq_filename = "validSeq.csv"
validSeq_path = os.path.join(folder_path_seq,validSeq_filename)
get_infoTable_and_validSeq(validSeq_path,comparisonTable_path)


    
    






