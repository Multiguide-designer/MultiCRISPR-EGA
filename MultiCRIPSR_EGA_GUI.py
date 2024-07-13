import os,csv
import webview
import CRISPR_library.EGA as ga_class
import CRISPR_library.spacer_extract as extract
import CRISPR_library.outputResult as outputResult
import datetime


config = {
            'proteinType':'Cas9',
            'geneNumber': 4,
            'pamSeq': ['AGG', 'TGG'],
            'insulatorSeq': '',
            'spacerSeq': '',
            'geneList': [],
            'spacerList': [],
            'populationSize': 0,
            'termination': 0,
            'crossover': 0.0,
            'mutation': 0.0,
        }

resultPath = os.path.join(os.path.dirname(__file__),'Result_GUI')
if not os.path.exists(resultPath):
    os.makedirs(resultPath)


class Api():
    def __init__(self) -> None:
        self.resultPath = resultPath
        current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.folder_path = os.path.join(self.resultPath, current_time)
        pass
    
    def setInitParam(self,geneNumber,pamSeq,insulatorSeq,spacerSeq,proteinType,targetLength,pamPosition):
        config['geneNumber'] = geneNumber
        config['pamSeq'] = pamSeq.split()
        config['insulatorSeq'] = insulatorSeq.replace('T', 'U')
        config['spacerSeq'] = spacerSeq
        config['targetLength'] = targetLength
        config['pamPosition'] = pamPosition
        config['proteinType'] = proteinType
        pass
    
    
    def setGeneInfo(self, geneList):
        # geneList: [
        # {
        #   'geneName': string,
        #   'geneSeq': string,
        #   '20ntonNegative': bool,
        # },
        # ....
        # ]
        config['geneList'] = geneList
        pass


    def getSpacer(self):
        
        spacerList = extract.generate_spacerlist(config['geneList'],config['pamSeq'],config['pamPosition'],config['targetLength'])
        # [
        #     {
        #         'geneName': 'gene1',
        #         'sequence': 'ATCGATCGATCGATCGATCG',
        #         'position' : '+'
        #     },
        #     {
        #         'geneName': 'gene1',
        #         'sequence': 'ATCGATCGATCGATCGATCG',
        #         'position' : '-'
        #     },
        #     {
        #         'geneName': 'gene1',
        #         'sequence': 'ATCGATCGATCGATCGATCG',
        #         'position' : '+'
        #     },
        #     {
        #         'geneName': 'gene1',
        #         'sequence': 'ATCGATCGATCGATCGATCG',
        #         'position' : '+'
        #     },
        #     {
        #         'geneName': 'gene2',
        #         'sequence': 'ATCGATCGATCGATCGATCG',
        #         'position' : '-'
        #     },
        #      ...
        # ]
        return spacerList
    
    
    def updateSpacer(self, spacerList):
        config['spacerList'] = spacerList
         #[
            # {
            #     'geneName': 'gene1',
            #     'sequence': 'ATCGATCGATCGATCGATCG'
            
            # },
            # ...
        #]
        
        spacerDist = {}
        for geneInfo in config['spacerList']:    
            spacer = geneInfo['spacer']
            geneName = geneInfo['geneName']
            if geneInfo['geneName'] in spacerDist :
                spacerDist[geneName].append(spacer)
            else:
                spacerDist[geneName] = [spacer]
        config['spacerDist']  = spacerDist 
        
        # {
            #'gene1': ['CGCTATGGTTATGCGTAAGC','GTCGCGGCAATATAATGAGA',...],
            #'gene2': ... 
        # }
        
        numSpacerList = [len(seq_list) for seq_list in spacerDist.values()]
        config['numSpacerList']  = numSpacerList
        
        self.get20ntFlag = False
        pass
    
    
    def getSpacerFile(self,spacerList):
        index = 1
        output_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'crispr_targets')
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'crispr_targets', 'spacer.csv'), "w+",
              newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Gene Name", "spacer","Position(w.r.t given sequence)", 
                         "Position", "Index"])
            
            geneName = spacerList[0]['geneName']
            for geneInfo in spacerList:
                if  geneName== geneInfo['geneName']:
                    writer.writerow([geneInfo['geneName'], geneInfo['spacer'], geneInfo['position'], index])
                    index += 1
                else:
                    index = 1
                    geneName = geneInfo['geneName']
                    writer.writerow([geneInfo['geneName'], geneInfo['spacer'], geneInfo['position'], index])
                    index += 1
        
        print("The extracted spacer file has been generated in crispr_targets .")
        pass
            
               
    def setGaParam(self, populationSize, termination, crossover, mutation,tournamentSize):
        if  self.get20ntFlag is False:
            self.getSpacerFile(config["spacerList"])
            self.get20ntFlag = True
        
        config['populationSize'] = populationSize
        config['termination'] = termination
        config['crossover'] = crossover
        config['mutation'] = mutation
        config['tournamentSize'] = tournamentSize
    
        self.ga_instance = None
        pass


    def runGaIter(self):
        if self.ga_instance is None:
            self.ga_instance = ga_class.Ga(config)
            return self.ga_instance.get_first_minMFE()

        minimumMFEs = self.ga_instance.iter_once()
        return minimumMFEs[-1]


    def getResult(self):
        return self.ga_instance.sort()
# [
#   {
#     "freeEnergy": -276,
#     "geneList": [
#       {
#         "type": "insulator",
#         "sequence": "AGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#       {
#         "type": "gene1",
#         "sequence": "CCAGCACATTCCACAGGGTA"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#       {
#         "type": "gene2",
#         "sequence": "CAGCTAAAATTAGTCGCTTT"
#       },
#       {
#         "type": "spacer",
#         "sequence": "GTGTCATAGCCCAGCTTGGCGGGCGAAGGCCAAGAC"
#       },
#     ]
#   },
#  {...


    def exportToFile(self):
        
        os.makedirs(self.folder_path, exist_ok=True)
        SeqFilePath = os.path.join(self.folder_path, 'AllSequence_Iteration{}.csv'.format(self.ga_instance.generation))
        outputResult.obtainFile(self.ga_instance.all_individuals, SeqFilePath)
        return SeqFilePath

    
    def getRnafig(self, seq):
        RnaFigFilePath = os.path.join(self.folder_path,'RNA_visualize')
        os.makedirs(RnaFigFilePath, exist_ok=True)
        outputResult.visualize_rna(seq,RnaFigFilePath)
        return None


if __name__ == '__main__':
    api = Api()
    webview.create_window('Bioview', 'dist/index.html', js_api=api, min_size=(800, 650))
    webview.start()
    # webview.start(debug=True)