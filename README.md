# MultiCRISPR-EGA
The supplementary folder contains the supplemental materials for the paper, which include detailed descriptions of how to run the library. 

Datasets of varying sizes for Escherichia coli are available in the '**Gene_Info**' folder, and the generated multiplexed gRNA sequences can be found in the '**Result_xxx**' folders.


# Preparation
run:

```
conda create --name multicrisprEGA python==3.10

conda activate multicrisprEGA

conda install --file requirements.txt
```

# Additional Software Installation
Install ViennaRNA and Ghostscript

1. Download the ViennaRNA from [ViennaRNA official site](https://www.tbi.univie.ac.at/RNA/) and Ghostscript from [Ghostscript download page](https://www.ghostscript.com/download/gsdnld.html).

2. Follow the installation instructions provided in the downloaded package.

3. Add their bin directory to your system's environment variables.

We also have provided a pre-downloaded ViennaRNA package in the '**Additional Software and Package**' folder

# Run the Application
This application can be run in two different modes depending on your preference for using a graphical user interface (GUI) or running it via a script.

## Run via GUI
If you prefer to use the GUI for a more interactive experience, simply run the following command in your terminal:


```
python GUI.py
```
This will launch the graphical user interface, allowing you to input parameters and execute tasks through a user-friendly interface.

## Run via CLI
If you prefer to run the application via comand-line for automation or integration into workflows, you can execute the script version by running:

``` 
python MultiCRIPSR_EGA_CLI.py --insulatorSeq 'AGCUGUCACCGGAUGTGCUUUCCGGUCUGAUGAGUCCGUGAGGACGAAACAGCCUCUACAAAUAAUUUUGUUUAA' --repeatSeq 'GUGUCAUAGCCCAGCUUGGCGGGCGAAGGCCAAGAC' --geneSpacerFile .\Gene_Info\Gene_Info8\spacer.csv --populationSize 100 --maxIterations 10 --crossover 0.8 --mutation 0.05 --tournamentSize 3
```

## Dataset Preparation Guide
To extract genes from a specific species and generate CRISPR target sequences, follow these steps:

### 1. Extract Gene Data from NCBI
First, you can use the  **Extract_form_NCBI.ipynb**  notebook to retrieve gene information from the NCBI database for a specific species.

### 2. Generate CRISPR Target Sequences
Once you have the your_dataset_file.csv file, use the following command to generate CRISPR target sequences (spacers):
```
cd CRISPR_library

python spacer_extract.py --geneInfoFile data/your_dataset_file.csv --PAM [PAM sequences] --pamPosition [before|after] --targetLength [target length] --outputFile [desired output folder]
```
