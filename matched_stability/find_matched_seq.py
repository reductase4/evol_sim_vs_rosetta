
import shutil
import os
import fnmatch
import math

PROJECT_DIR = os.getcwd()
SCORES_DIR = PROJECT_DIR + "/designed_scores"
SEQ_DIR = PROJECT_DIR + "/designed_sequences_fasta"
RESULT_SCORE_DIR = PROJECT_DIR + "/results_scores"
RESULT_FASTA_DIR = PROJECT_DIR + "/result_fasta"
RESULT_CSV_DIR = PROJECT_DIR + "/result_csv"

SCORE_PREFIX = "results_"
DESIGNED_SCORE_POSTFIX = "_rosetta.csv"

def readScore(path, skip_title=False, delimiter=" "):
    with open(path) as f:
        if skip_title:
            f.readline()
        content = f.readlines()
        scores = {}
        for line in content:
            strs = [x.strip() for x in line.split(delimiter)]
            scores[strs[0]] = strs[1]
        return scores
        
def initDir(path):
    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path) 

initDir(RESULT_SCORE_DIR)
initDir(RESULT_FASTA_DIR)
initDir(RESULT_CSV_DIR)

for root, directories, files in os.walk(SCORES_DIR):
    for one_file in fnmatch.filter(files, "*scores.csv"):
        name = one_file[0:6]
        
        evolved_scores = readScore(path = SCORES_DIR + "/" + one_file, skip_title = True, delimiter = ",")
        sorted_evolved_scores = sorted(evolved_scores.items(), key=lambda x: x[1])

        max_evolved_score = float(sorted_evolved_scores[0][1])
        min_evolved_score = float(sorted_evolved_scores[len(sorted_evolved_scores) - 1][1])

        designed_scores = readScore(SCORES_DIR + "/" + name + "_score_rosetta.csv")
        sorted_designed_scores = sorted(designed_scores.items(), key=lambda x: x[1])
        
        index_match = {}
        for designed_iter in sorted_designed_scores:
            designed_index = designed_iter[0]
            designed_score = float(designed_iter[1])
            
            if designed_score >= min_evolved_score and designed_score <= max_evolved_score:
                closest_score = float('nan')
                closest_index = float('nan')
                for evolved_iter in sorted_evolved_scores:
                    evolved_index = evolved_iter[0]
                    evolved_score = float(evolved_iter[1])
                    
                    if evolved_index not in index_match:
                        if math.isnan(closest_score) or abs(closest_score - designed_score) > abs(evolved_score - designed_score):
                            closest_score = evolved_score
                            closest_index = evolved_index

                index_match[closest_index] = designed_index
                
        for evolved_index in index_match:
            with open(RESULT_SCORE_DIR + "/" + name + "_evolved_scores.csv", "a") as f:
                f.write(evolved_index + "," + str(evolved_scores[evolved_index]) + "\r\n")
            
            with open(RESULT_SCORE_DIR + "/" + name + "_rosetta_scores.csv", "a") as f:
                    f.write(index_match[evolved_index] + "," + str(designed_scores[index_match[evolved_index]]) + "\r\n")
            
        with open(SEQ_DIR + "/" + name + "_evolved_seqs.fasta") as f:
            content = f.readlines()
            for i in range(len(content)):
                name_and_idx = content[i][1:].split(".")[0]
                if name_and_idx in index_match:
                    with open(RESULT_FASTA_DIR + "/" + name + "_evolved_seqs.fasta", "a") as f1:
                        f1.write(content[i])
                        f1.write(content[i + 1])
                    
                    with open(RESULT_CSV_DIR + "/" + SCORE_PREFIX + name + "_evolved.csv", "a") as f1:
                        f1.write(content[i + 1])
                i += 2
                
        with open(SEQ_DIR + "/" + name + "_designed_seqs.fasta") as f:
            content = f.readlines()
            for i in range(len(content)):
                name_and_idx = content[i][1:].split(".")[0]
                if name_and_idx in index_match.values():
                    with open(RESULT_FASTA_DIR + "/" + name + "_designed_seqs.fasta", "a") as f1:
                        f1.write(content[i])
                        f1.write(content[i + 1])
                    
                    with open(RESULT_CSV_DIR + "/" + SCORE_PREFIX + name + "_rosetta.csv", "a") as f1:
                        f1.write(content[i + 1])
                i += 2
