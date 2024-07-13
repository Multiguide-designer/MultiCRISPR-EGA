import csv
import subprocess
import os
import matplotlib.pyplot as plt

def obtainFile(all_individuals,filePath):
    all_individuals_sorted = sorted(all_individuals, key=lambda ind: ind.fitness)
    with open(filePath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence', 'Energy'])
        for individual in all_individuals_sorted:
            writer.writerow([individual.seq,individual.fitness])
    print("迭代结果已保存至",filePath)


def line_graph_average(average_MFEs, folder_path):
    plt.figure(figsize=(8, 6), dpi=300)
    
    # 循环绘制每个算法的折线和误差带
    plt.plot(range(len(average_MFEs)), average_MFEs, color='purple', linestyle='-')

    plt.xlabel("Iteration")
    plt.ylabel("MFE")

    # 添加标题和标签
    plt.title('Average Minimum Free Energy (MFE) Comparison')

    # 显示图表
    plt.tight_layout()

    # 保存图表
    os.makedirs(folder_path, exist_ok=True)
    plt.savefig(os.path.join(folder_path, "Average_MFE_line_graph.png"))

    plt.close()
    print("图表已保存至", os.path.join(folder_path, "Average_MFE_line_graph.png"))
    return True


def line_graph_minimum(minimumMFEs, folder_path):
    plt.figure(figsize=(8, 6), dpi=300)
    plt.plot(range(len(minimumMFEs)), minimumMFEs, color='purple', linestyle='-')

    # 添加标题和标签
    plt.title('Minimum Minimum Free Energy (MFE) Comparison')
    plt.xlabel('Iteration')
    plt.ylabel('Minimum MFE')

    # 显示图表
    plt.tight_layout()

    # 保存图表
    os.makedirs(folder_path, exist_ok=True)
    plt.savefig(os.path.join(folder_path, "Minimum_MFE_line_graph.png"))

    plt.close()
    print("图表已保存至", os.path.join(folder_path, "Minimum_MFE_line_graph.png"))
    return True


# 使用RNAplot可视化RNA二级结构
def run_RNAfold(sequence):
    cmd = ["RNAfold"]
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=sequence)
    return stdout


def visualize_rna(sequence, output_directory):
    output = run_RNAfold(sequence)
    with open('temp_output.txt', "w") as f:
        f.write(output)  # output 是表示RNA二级结构的文本数据,分别是（序列 \n 二级结构 （自由能））
    cmd = ["RNAplot", "-i", "temp_output.txt", "-o", "ps"]
    subprocess.run(cmd)
    os.remove('temp_output.txt')
    
    # 检查.ps文件是否已存在，如果存在则重命名
    output_ps_file = "sequence_1.ps"
    index = 1
    while os.path.exists(os.path.join(output_directory, output_ps_file)):
        index += 1
        output_ps_file = f"sequence_{index}.ps"    
    os.rename("rna.ps", os.path.join(output_directory, output_ps_file))

    # 将.ps文件转换为.pdf文件
    output_pdf_file = f"sequence_{index}.pdf"
    cmd = ["ps2pdf", os.path.join(output_directory, output_ps_file),
           os.path.join(output_directory, output_pdf_file)]
    subprocess.run(cmd)

    print("pdf file has been generated in ",output_directory)
    
