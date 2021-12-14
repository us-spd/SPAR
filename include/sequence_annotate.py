
#from Bio.SubsMat import MatrixInfo
from Bio.Align import substitution_matrices
import json
import math
import numpy as np
import os
import re
import shutil
from sklearn.cluster import SpectralClustering
#import statistics as stats
import subprocess
import sys


if __name__ == "__main__":
    from misc import alignment, fasta, message, nested, string
    from settings import bash_path, translational_frameshift_di
else:
    from .misc import alignment, fasta, message, nested, string
    from .settings import bash_path, translational_frameshift_di
    
################################################################################
        
spar_path = re.sub("/include$", "", os.path.dirname(os.path.realpath(__file__)))
valid_organism_li = [d for d in os.listdir(spar_path+"/required/") if os.path.isdir(spar_path+"/required/"+d) and d != "blastdb"]
if bash_path != None:
    os.environ["PATH"] = bash_path

################################################################################   
################################################################################

class sequence_annotate:
    
    def __init__():
        '''The sequence_annotate class is initialized'''
    
    class prepare:
        
        def __init__(self):
            '''The prepare subclass is initialized'''
            
            # check required bash utilities
            required_program_li = ["hmmbuild", "mafft"]
            missing_programs = string.readable_list([x for x in required_program_li if shutil.which(x) == None])
            if len(missing_programs) > 0:
                raise Exception("The following dependencies could not be found: "+missing_programs+". Please install and re-run.")
                
        
        def cluster_expression(self, sequence_li, cropped_msa_li):
            '''Iterative clustering of residues to produce selective (minimize incorrect matches) and specific (maximize correct matches) regular expressions'''
            # https://stats.stackexchange.com/questions/123060/clustering-a-long-list-of-strings-words-into-similarity-groups
            # https://www.programcreek.com/python/example/85778/sklearn.cluster.AffinityPropagation
            # stepwise process of increasing selectivity whle decreasing specificity
            # input sequences are grouped by similarity; residues may be removed
            
            aa_li = [x for x in "ARNDCEQGHILKMNFPSTWYV"]
            ambiguous_aa_di = {"B":["N", "D"], "Z":["E", "Q"], "J":["I", "L"], "X":aa_li}
            #blosum_di = MatrixInfo.blosum62
            b_matrix = substitution_matrices.load("BLOSUM62")
            blosum_di = {(b_matrix.alphabet[i1], b_matrix.alphabet[i2]): b_matrix[i1, i2]  for i1, x in enumerate(b_matrix) for i2, y in enumerate(x)}
            
            read_fasta_li = [x.replace("-", "") for x in sequence_li]    
            for crop_index, crop in enumerate(cropped_msa_li):
                crop = crop.replace("-", "") # short regex unlikely to pass duplicate check
                if len([x for x in crop if x not in ambiguous_aa_di]) >= min(10, len(cropped_msa_li[0])): # expression must be at least 10 non-ambiguous residues in length
                    for read_index, read in enumerate(read_fasta_li):
                        if read.count(crop) > 1: # target sequence likely contains repeat region
                            read_fasta_li[read_index] = "removed"
                else:
                    cropped_msa_li[crop_index] = "removed"
            read_fasta_li = [x for x in read_fasta_li if x != "removed"]
            cropped_msa_li = [x for x in cropped_msa_li if x != "removed"]
            output = None
                        
            # remove low occurrence residues prior to clustering
            # rare residues can be a result of mutation or poor alignment
            remove_threshold, remove_count = 0, 0
            word_li = list(set(cropped_msa_li))
            copy_word_li = word_li[:]
            if len(word_li) == 1:
                output = word_li[0], 0
            while len(word_li) > 1:
                if output != None and output[1] <= remove_count:
                    return(output)
                distance_similarity = [[blosum_di[w1, w2] if (w1, w2) in blosum_di else \
                                        [min([blosum_di[x] for x in blosum_di]) if w1 == w2 else max([blosum_di[x] for x in blosum_di])][0] \
                                        for w1 in word_li] for w2 in word_li]
                
                max_n = min(10, len(word_li)-1)
                for n_cluster in range(1, max_n+1):
                    # random error may occur during fit
                    # grouping is random; difficult to optomize
                    try:
                        affprop = SpectralClustering(n_clusters=n_cluster, affinity="precomputed").fit_predict(distance_similarity)
                        segment_li2, cluster_li2 = [], []
                        
                        for cluster_index in list(set(affprop)):
                            cluster = [word_li[x] for x, y in enumerate(affprop) if y == cluster_index]
                            segment_li = ["".join(set([x[index] for x in cluster])) for index in range(len(cluster[0]))]
                            segment_li2.append(segment_li)
                            cluster_li2.append(cluster)
                            
                        for index in range(len(segment_li2)):
                            product_li = []
                            for column in segment_li2[index]:
                                for character in ambiguous_aa_di:
                                    column.replace(character, "".join(ambiguous_aa_di[character]))
                                if "-" in column:
                                    product_li.append(1)
                                else:
                                    product_li.append(len(column)/20)
                            if np.prod(product_li) > 10**-11:
                                segment_li2[index] = []
                            else:
                                cluster_li2[index] = []
                            # replace beginning positions containing "-" with "." to maintain start position
                            segment_li2[index] = ["["+x.replace("-", "")+"]?" if "-" in x else "["+x+"]" for x in segment_li2[index]]
                            segment_li2[index] = [re.sub("["+"".join([y for y in ambiguous_aa_di])+"]", "", x) for x in segment_li2[index]] #remove ambiguous characters
                            segment_li2[index] = [x.replace("[]?", "").replace("[]", ".") for x in segment_li2[index]] #remove "[]" (former ambiguous character) and "[]?" (former gap)
                            segment_li2[index] = [re.sub(r"\[([A-Z\.]{1})\]", r"\1", x) for x in segment_li2[index]] #remove brackets around single characters
                        
                        miss_count = 0
                        regex = "|".join(["".join(x) for x in segment_li2 if x != []])
                        for sequence in read_fasta_li:
                            # detect misalignment in protein MSA
                            if len(re.findall(regex, sequence)) > 1:
                                miss_count = len(read_fasta_li) # prevent output if multiple matches detected (reiterate)
                                break
                            elif len(re.findall(regex, sequence)) == 0:
                                miss_count += 1
                        
                        if miss_count == 0:
                            return(regex, miss_count)
                        elif miss_count < len(read_fasta_li) and (output == None or miss_count < output[1]):
                            output = regex, miss_count
                    # multiple consecutive failures unlikely
                    except:
                        #print(n_cluster, word_li, file=sys.stderr)
                        pass
                            
                # rebuilds cropped sequences to reflect highest character (including gaps/ambiguities) homology
                while word_li == copy_word_li:
                    cropped_column_li2 = [[cropped_msa_li[x][y] for x in range(len(cropped_msa_li))] for y in range(len(cropped_msa_li[0]))]
                    for column_index in range(len(cropped_column_li2)):
                        column_li = cropped_column_li2[column_index]
                        max_residue = max(set(column_li), key = column_li.count)
                        replace_li = []
                        for residue in set(column_li):
                            if column_li.count(residue) == remove_threshold:
                                replace_li.append(residue)
                        for row_index in range(len(cropped_column_li2[column_index])):
                            if cropped_column_li2[column_index][row_index] in replace_li:
                                split_row_li = [x for x in cropped_msa_li[row_index]]
                                split_row_li[column_index] = max_residue
                                cropped_msa_li[row_index] = "".join(split_row_li)
                    remove_threshold += 1
                    word_li = list(set(cropped_msa_li))
                copy_word_li = word_li[:]
                remove_count = len([1 for x in read_fasta_li if len(re.findall("|".join([y.replace("-", "") for y in word_li]), x)) == 0])
                        
            return(output)
                
                    
        def calculate_shannon(self, input_fasta):
            '''Calculates Shannon diversity index for each column in MSA fasta to gauge homology/conservation'''
            # https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/shannon.htm
            
            read_fasta_li = fasta.read(input_fasta)[1]
            shannon_li = []
            for column_index in range(len(read_fasta_li[0])):
                column = [read_fasta_li[x][column_index] for x in range(len(read_fasta_li))]
                gap_count = [read_fasta_li[x][column_index] for x in range(len(read_fasta_li))].count("-")
                column = [x for x in column if x != "-"] # remove gaps from target column
                proportion_li = [column.count(x)/len(column) for x in set(column)]
                if len(proportion_li) > 1:
                    shannon_diversity = -sum([x*math.log(x) for x in proportion_li])
                    shannon_equitability = shannon_diversity/(math.log(len(proportion_li)))
                else:
                    shannon_equitability = 0
                shannon_li.append((1-shannon_equitability)*len(column)/(len(column)+gap_count))
            return(shannon_li)
            
            
        def regex_select(self, protein_directory, cropped_regex_json):
            # from descending order of highest homology/conservation, regions are converted into regular expressions
            # regulate maximum distance/overlap between regions
            
            fasta_li = [f for f in os.listdir(protein_directory+"/msa") if os.path.isfile(os.path.join(protein_directory+"/msa", f)) and ".fa" in f]    
            for fasta_file in fasta_li:
                annotation = fasta_file.split(".")[0]
                msa_save_path = protein_directory.replace("protein", "")+"/hmm_profiles/msa_save/"+fasta_file
                if annotation not in cropped_regex_json:
                    cropped_regex_json[annotation] = []
                    protein_fasta = protein_directory+"/msa/"+fasta_file
                    protein_fasta_li2 = fasta.read(protein_fasta)
                    shannon_li = self.calculate_shannon(protein_fasta)
                    shannon_window_li = [sum(shannon_li[x:x+15]) for x in range(1, len(shannon_li)-16)] # allow near complete overlap with terminal ends
                    shannon_window_li2 = [[x+1, y] for x, y in enumerate(shannon_window_li)] # index, value
                    shannon_window_li2.sort(reverse = True, key = lambda x: x[1]) # sort by highest window sum
                    # identify gaps removed from msa_save
                    save_fasta_li2 = fasta.read(msa_save_path)
                    modify_li = fasta.remove_column(save_fasta_li2[1], 0.5)[1] # gap proportion should be the same as in check_msa()
                    
                    index_li = [0, max(0, len(shannon_li)-15)]
                    while len(shannon_window_li2) > 0:
                        window_index = shannon_window_li2[0][0]
                        min_less_int = min([x for x in index_li if x < window_index], key=lambda x:abs(x-window_index))
                        min_more_int = min([x for x in index_li if x > window_index], key=lambda x:abs(x-window_index))
                        overlap = max(0, min_less_int+15-window_index) + max(0, window_index+15-min_more_int)
                        # prior windows are assumed to be 15 long regardless of regex length; match redundancy possible but unlikely
                        if overlap <= 10: # maximum allowed overlap: 10 (2/3)
                            index_li.append(window_index)
                        shannon_window_li2 = shannon_window_li2[1:]
                    
                    for index, index_value in enumerate(index_li):
                        cropped_msa_li = [x[index_value:index_value+15] for x in protein_fasta_li2[1]]
                        cluster = self.cluster_expression(protein_fasta_li2[1], cropped_msa_li)
                        if cluster != None: # failure when region is highly diverse, resulting in complex, non-specific regex
                            # regex_select() requires indices of residues relative to save_msa
                            frame_mod = re.findall("#frameshift_modification=\((.*)\)", protein_fasta_li2[0][0])
                            if len(frame_mod) > 0 and len([x for x in frame_mod[0].split(",") if x.strip().isdigit()]) == 3:
                                adjust_li = [int(num.strip()) for num in frame_mod[0].split(",")]
                                adjust_li = adjust_li[:-1] # remove value used in modify_hmmalign()
                                save_fasta_li2[1] = [x[:adjust_li[1]]+x[adjust_li[1]-adjust_li[0]:] for x in save_fasta_li2[1]]
                                modify_li = modify_li[:adjust_li[1]]+modify_li[adjust_li[1]-adjust_li[0]:]
                            # cluster may contain gaps (should be forbidden; fails to match)
                            protein_match_regex_li = [re.findall(cluster[0], x.replace("-", "")) for x in cropped_msa_li] # only consider sequences that match regex
                            combine_li = [len(x[:index_value].replace("-", ""))*3 if len(protein_match_regex_li[i]) > 0 else "" for i, x in enumerate(protein_fasta_li2[1])]
                            # assumes that protein_fasta_li2 and save_fasta_li2 have the same header order
                            # correction for gaps in msa_save nucleotide alignment
                            for combine_index, combine_value in enumerate(combine_li):
                                if combine_value != "":
                                    save_value = combine_value
                                    # correction for gaps that were not replaced by fasta.remove_column()
                                    while len(save_fasta_li2[1][combine_index][:save_value].replace("-", "")) < combine_value and \
                                    len(save_fasta_li2[1][combine_index][:save_value].replace("-", "")) < len(save_fasta_li2[1][combine_index].replace("-", "")):
                                        save_value += 1
                                    combine_li[combine_index] = save_value
                                    
                            # correction for gaps removed by fasta.remove_column()
                            combine_li = [x - modify_li[x] for x in combine_li if x != ""]
                            if len(combine_li) > 0:
                                aa_position = round(index_value/(len(shannon_li)-15), 6)
                                if index == 0:
                                    aa_position = 0
                                elif index == 1:
                                    aa_position = 1                        
                                cropped_regex_json[annotation].append({
                                    "regex": cluster[0],
                                    "aa_position": aa_position,
                                    "nucleotide_position": max(combine_li)
                                })
    
                    with open(protein_directory+'/protein.json', 'w') as outfile:
                        for cds in cropped_regex_json:
                            cropped_regex_json[cds].sort(key=lambda x: x["aa_position"])
                        json.dump(cropped_regex_json, outfile)
                            
            
        def check_msa(self, organism, dir_path=spar_path):
            '''Checks msa_save for consistent logic and constructs msa_build file (reference for hmm profile construction)'''
            #requires manual modification for translational frameshifts, limited to -1 frameshifts (not user friendly)
            
            input_nucleotide_msa_directory = dir_path + "/required/" + organism + "/hmm_profiles/msa_save"
            output_nucleotide_msa_directory = dir_path + "/required/" + organism + "/hmm_profiles/msa_build"
            output_protein_msa_directory = dir_path + "/required/" + organism + "/protein/msa"
            target_files = [f for f in os.listdir(input_nucleotide_msa_directory) if os.path.isfile(os.path.join(input_nucleotide_msa_directory, f))]
            for file in target_files:
                if ".fa" in file:
                    annotation_str = file.split(".")[0]
                    read_fasta_li2 = fasta.read(input_nucleotide_msa_directory + "/" + file) # nucleotide
                    read_fasta_li2[1] = [x.upper() for x in read_fasta_li2[1]]
                    
                    # location of modified residue/s for frameshift correction
                    frameshift_li = [] # slide distance/direction, msa_save index, msa_build index (gaps removed)
                    consensus_str = "".join([max([y[x] for y in read_fasta_li2[1]], key=[y[x] for y in read_fasta_li2[1]].count) for x in range(len(read_fasta_li2[1][0]))])
                    # find translation frameshift location
                    if organism in translational_frameshift_di and annotation_str in translational_frameshift_di[organism]["fs_annotation"]:
                        fs_di = translational_frameshift_di[organism]
                        fs_regex, fs_index, fs_type = fs_di["fs_regex"], fs_di["fs_index"], fs_di["fs_type"]
                        fs_regex = "-*".join([x for x in fs_regex])
                        fs_match_regex = re.findall(fs_regex, consensus_str)
                        # check that regex matches and premature stop codons exist (frame change necessary)
                        if len(fs_match_regex) > 0 and min([string.translate_nucleotides(read_fasta_li2[1][x], 0)[:-1].count("x") for x in range(len(read_fasta_li2[1]))]) > 0:
                            add_adjust = len(re.findall("^"+"-*".join(["[A-Z]" for x in range(fs_index)]), fs_match_regex[0])[0])
                            fs_regex_adjusted_index = consensus_str.index(fs_match_regex[0]) + add_adjust
                            frameshift_li = [fs_type, fs_regex_adjusted_index]
                            # track frequency of nucleotides adjacent to translational frameshift location
                            adjacent_nucleotide_di = {x:0 for x in "ATGC"} #"AUGC"}
                            for x in range(len(read_fasta_li2[1])):
                                for nucleotide in read_fasta_li2[1][x][fs_regex_adjusted_index-1:fs_regex_adjusted_index+1]:
                                    if nucleotide in adjacent_nucleotide_di:
                                        adjacent_nucleotide_di[nucleotide] += 1
                            # compensating for translational frameshift by inserting least likely nucleotide match as lower()
                            for x in range(len(read_fasta_li2[1])):
                                copy_adjacent_nucleotide_di = {x:adjacent_nucleotide_di[x] for x in adjacent_nucleotide_di}
                                for nucleotide in read_fasta_li2[1][x][fs_regex_adjusted_index-1:fs_regex_adjusted_index+1]:
                                    if nucleotide in copy_adjacent_nucleotide_di:
                                        del copy_adjacent_nucleotide_di[nucleotide]
                                read_fasta_li2[1][x] = read_fasta_li2[1][x][:fs_regex_adjusted_index] + min(copy_adjacent_nucleotide_di, key=copy_adjacent_nucleotide_di.get).lower() + read_fasta_li2[1][x][fs_regex_adjusted_index:]
                    # remove sequences from MSA that cannot be translated (premature stop codon/s)
                    protein_li = [string.translate_nucleotides(x, 0) for x in read_fasta_li2[1]]
                    start_li, end_li = [x[0] for x in protein_li], [x[-1] for x in protein_li]
                    start_aa, end_aa = max(start_li, key=start_li.count), max(end_li, key=end_li.count)
                    for index, protein_str in enumerate(protein_li):
                        if ("x" in protein_str[:-1] or                      # premature stop codons
                            (start_aa != "M" and protein_str[0] == "M") or  # consistent start
                            (end_aa != "x" and protein_str[-1] == "x")):    # consistent end
                            protein_li[index] = ""
                            for l in range(2):
                                read_fasta_li2[l][index] = ""
                    if protein_li.count("") > 0:
                        print(protein_li.count(""), "sequences eliminated from", file, file=sys.stderr)
                        protein_li = [x for x in protein_li if x != ""]
                        for l in range(2):
                            read_fasta_li2[l] = [x for x in read_fasta_li2[l] if x != ""]
                            
                    # remove columns by gap proportion (avoid problems with hmmalign)
                    remove_gap_out_li = fasta.remove_column(read_fasta_li2[1], 0.5) # translate indices between save_msa and build_msa
                    read_fasta_li2[1] = remove_gap_out_li[0]
                    if len(frameshift_li) > 0:
                        frameshift_li.append(frameshift_li[1] - remove_gap_out_li[1][frameshift_li[1]])
                        # add location of modified residue/s for frameshift correction
                        # modify_hmmalign() requires indices of modified residues relative to build_msa
                        # regex_select() requires indices of residues relative to save_msa
                        for i, head in enumerate(read_fasta_li2[0]):
                            read_fasta_li2[0][i] += "#frameshift_modification=("+",".join([str(x) for x in frameshift_li])+")"
                    # write protein fasta and align
                    fasta.write(output_protein_msa_directory+"/"+file, [read_fasta_li2[0], protein_li])
                    alignment.mafft.align(output_protein_msa_directory+"/"+file, output_protein_msa_directory+"/"+file)
                    # write msa_build; do not re-align after removing sequences/columns (avoid disrupting frameshift index)
                    fasta.write(output_nucleotide_msa_directory + "/" + file, read_fasta_li2)
            
            
        def build_reference_gff(self, organism, dir_path=spar_path):
            '''Builds reference template gff3 file (if one does not already exist)'''
            # user may make modifications to the template file or provide their own
            
            organism_path = dir_path+"/required/"+organism
            if not os.path.isdir(organism_path+"/gff3"):
                subprocess.call("mkdir -p "+organism_path+"/gff3", shell=True)
            if not os.path.isfile(organism_path+"/gff3/template.gff3"):
                print("Required sample template could not be found for "+organism+". Generating..." )
                gff3_di = {}
                for protein_file in os.listdir(organism_path+"/protein/msa"):
                    if ".fa" in protein_file:
                        protein = re.sub("\.fa.*$", "", protein_file)
                        protein_li = fasta.read(organism_path+"/protein/msa/"+protein_file)[1]
                        # alignment unnecessary; variable length may impede appropriate treatment
                        protein_li = [x.replace("-", "") for x in protein_li]
                        # all gene and CDS are required to have start/stop codons
                        if len([x for x in protein_li if x[0] == "M" and x[-1] == "x"]) > len(protein_li)/2:
                            if "gene/CDS" not in gff3_di:
                                gff3_di["gene/CDS"] = []
                            gff3_di["gene/CDS"].append(protein)
                        else:
                            if "mat_peptide" not in gff3_di:
                                gff3_di["mat_peptide"] = []
                            gff3_di["mat_peptide"].append(protein)      
                # find parent of all mat_peptide
                parent_di = {}
                if "gene/CDS" in gff3_di and "mat_peptide" in gff3_di:
                    # compare unaligned sequences in required reference files
                    for mat_peptide in gff3_di["mat_peptide"]:
                        sub_di = {re.sub("\.fa.*$", "", f):f for f in os.listdir(organism_path+"/hmm_profiles/msa_save")}
                        if mat_peptide in sub_di:
                            parent_di[mat_peptide] = []
                            mat_peptide_li = list(set([x.replace("-", "") for x in fasta.read(organism_path+"/hmm_profiles/msa_save/"+sub_di[mat_peptide])[1]]))
                            for gene_cds in gff3_di["gene/CDS"]:
                                if gene_cds in sub_di:
                                    match_count, cut_count_li = 0, []
                                    gene_cds_li = [x.replace("-", "") for x in fasta.read(organism_path+"/hmm_profiles/msa_save/"+sub_di[gene_cds])[1]]
                                    for mat_peptide_str in mat_peptide_li:
                                        for gene_cds_str in gene_cds_li:
                                            find_mat_peptide = gene_cds_str.find(mat_peptide_str)
                                            if find_mat_peptide != -1:
                                                match_count += 1
                                                # number of cuts necessary to produce mat_peptide
                                                cut_count_li.append(0)
                                                if find_mat_peptide > 0:
                                                    cut_count_li[-1] += 1
                                                if len(gene_cds_str) - find_mat_peptide - len(mat_peptide_str) > 3: # ignores stop codon
                                                    cut_count_li[-1] += 1
                                                break
                                    if match_count > len(mat_peptide_li)/2:
                                        parent_di[mat_peptide].append([gene_cds, len(gene_cds_li[0]), round(sum(cut_count_li)/len(cut_count_li)), match_count/len(mat_peptide_li)])
                gff3_li2 = []
                if "gene/CDS" not in gff3_di:
                    gff3_di["gene/CDS"] = ["polyprotein"]
                for x in range(2):
                    for gene_cds in gff3_di["gene/CDS"]:
                        if x == 0: 
                            gff3_li2.append(["sample", ".", "gene", "start", "end", ".", "+", ".", ";".join(["ID="+gene_cds+"-gene", "Name="+gene_cds])])
                        else: 
                            gff3_li2.append(["sample", ".", "CDS", "start", "end", ".", "+", ".", ";".join(["ID="+gene_cds+"-cds", "Name="+gene_cds+" CDS", "Parent="+gene_cds+"-gene"])])
                if "mat_peptide" in gff3_di:
                    for mat_peptide in gff3_di["mat_peptide"]:
                        mat_peptide_value_li = ["ID="+mat_peptide, "Name="+mat_peptide, "Parent="]
                        if mat_peptide in parent_di and len(parent_di[mat_peptide]) > 0:
                            # remove parents that have more than the minimum cut count (fewer cuts, more likely parent)
                            min_cut_count = min([x[2] for x in parent_di[mat_peptide]])
                            parent_di[mat_peptide] = [x for x in parent_di[mat_peptide] if x[2] == min_cut_count]
                            # add longest parent only
#                           mat_peptide_value_li[-1] += max(parent_di[mat_peptide], key = lambda x: x[1])[0]+"-gene"
                            # add all parents (sorted by length) separated by comma
                            mat_peptide_value_li[-1] += ",".join([x[0]+"-gene" for x in sorted(parent_di[mat_peptide], key = lambda x: x[1])])
                        elif len(parent_di) == 0:
                            mat_peptide_value_li[-1] += "polyprotein-gene"
                        else:
                            print("Missing parent for", mat_peptide, file=sys.stderr)
                            mat_peptide_value_li[-1] += "MISSING"
                        gff3_li2.append(["sample", ".", "mat_peptide", "start", "end", ".", "+", ".", ";".join(mat_peptide_value_li)])
                #write file
                write_file = open(organism_path+"/gff3/template.gff3", "w")
                write_file.write("##gff-version 3\n")
                write_file.write("# organsim "+organism+"\n")
                write_file.write("\t".join(["sample", ".", "region", "start", "end", ".", "+", ".", ";".join(["ID="+organism, "Name="+organism])]) + "\n")
                for line in gff3_li2:
                    write_file.write("\t".join(line) + "\n")
                write_file.close()
                
        
        def main(self, organism):
            '''Checks that all necessary dependencies exist; prompts user to install missing programs or build missing reference files'''
            
            if organism in valid_organism_li:
                msa_save_li, msa_build_li, hmm_build_li, protein_msa_li = [], [], [], []
                if os.path.isdir(spar_path+"/required/"+organism+"/hmm_profiles/msa_save"):
                    msa_save_li =    [re.sub("\.fa.*", "", f) for f in os.listdir(spar_path+"/required/"+organism+"/hmm_profiles/msa_save") \
                                      if os.path.isfile(spar_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f) and f[0] != "."]
                if os.path.isdir(spar_path+"/required/"+organism+"/hmm_profiles/msa_build"):
                    msa_build_li =   [re.sub("\.fa.*", "", f) for f in os.listdir(spar_path+"/required/"+organism+"/hmm_profiles/msa_build") \
                                      if os.path.isfile(spar_path+"/required/"+organism+"/hmm_profiles/msa_build/"+f)]
                if os.path.isdir(spar_path+"/required/"+organism+"/hmm_profiles/hmm_build"):
                    hmm_build_li =   [re.sub("\.hmm.*", "", f) for f in os.listdir(spar_path+"/required/"+organism+"/hmm_profiles/hmm_build") \
                                      if os.path.isfile(spar_path+"/required/"+organism+"/hmm_profiles/hmm_build/"+f)]
                if os.path.isdir(spar_path+"/required/"+organism+"/protein/msa"):
                    protein_msa_li = [re.sub("\.fa.*", "", f) for f in os.listdir(spar_path+"/required/"+organism+"/protein/msa") \
                                      if os.path.isfile(spar_path+"/required/"+organism+"/protein/msa/"+f)]            
                try:
                    with open(spar_path+"/required/"+organism+"/protein/protein.json") as json_file:  
                        cropped_regex_json = json.load(json_file)
                except:
                    cropped_regex_json = {}        
                if (len([x for x in msa_save_li if x not in msa_build_li]) != 0 or \
                    len([x for x in msa_save_li if x not in hmm_build_li]) != 0 or \
                    len([x for x in msa_save_li if x not in protein_msa_li]) != 0 or \
                    len([x for x in msa_save_li if x not in cropped_regex_json]) != 0) or \
                    not os.path.isfile(spar_path+"/required/"+organism+"/gff3/template.gff3"):
                    if message.proceed(input_value="Additional setup required for "+organism+". Continue? (y/n) ", \
                                       response_y="You have selected to continue. Please wait. Notice: this is a slow process.", \
                                       response_n="You have chosen to skip setup for "+organism+"."):
                        if len([x for x in msa_save_li if x not in msa_build_li]) != 0:
                            self.check_msa(organism, spar_path)
                        if len([x for x in msa_save_li if x not in hmm_build_li]) != 0:
                            #use hmmbuild to create hmm profiles (requires bash)
                            target_dir = spar_path + "/required/" + organism + "/hmm_profiles/msa_build"
                            build_dir = spar_path + "/required/" + organism + "/hmm_profiles/hmm_build"
                            subprocess.call("mkdir -p "+build_dir, shell=True)
                            target_file_li = [f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f)) and ".fa" in f]
                            for target_file in target_file_li:
                                msa_input = target_dir+"/"+target_file
                                hmm_profile = build_dir+"/"+re.findall("^[^.]*", target_file)[0]+".hmm"
                                alignment.hmm.build(msa_input, hmm_profile)
                        if len([x for x in msa_save_li if x not in cropped_regex_json]) != 0:
                            self.regex_select(spar_path+"/required/"+organism+"/protein", cropped_regex_json)
                        if not os.path.isfile(spar_path+"/required/"+organism+"/gff3/template.gff3"):
                            self.build_reference_gff(organism, spar_path)
                        print("Setup for "+organism+" finished successfully.", file=sys.stderr)
                        return(True)
                else:
                    return(True)
                
    ################################################################################   
    ################################################################################
    
    class process:
        
        # classify all functions that call global variable/s
        def __init__(self, organism, temp_dir=""):
            '''The process subclass is initialized'''
            # value appended to file names, avoids overwriting files when running in parallel
            
            # check required bash utilities
            required_program_li = ["hmmalign"]
            missing_programs = string.readable_list([x for x in required_program_li if shutil.which(x) == None])
            if len(missing_programs) > 0:
                raise Exception("The following dependencies could not be found: "+missing_programs+". Please install and re-run.")
            
            self.required_path = spar_path + "/required"
            self.temporary_path = spar_path + "/temporary"
            if os.path.isdir(temp_dir):
                self.temporary_path = temp_dir
            if organism in valid_organism_li:
                self.organism = organism
            else:
                raise Exception("Invalid input supplied for organism field:", organism)
                
            self.total_cds = [re.findall("^[^.]*", f)[0] for f in os.listdir(self.required_path+"/"+self.organism+"/hmm_profiles/msa_save") if ".fa" in f]
            with open(self.required_path+"/"+self.organism+"/protein/protein.json") as json_file:  
                cropped_regex_json = json.load(json_file)
            for cds in cropped_regex_json:
                cropped_regex_json[cds].sort(key=lambda x: x["aa_position"])
            self.match_residues_di = cropped_regex_json
            self.msa_len_li, self.translational_frameshift_di = [], {}
            
            for cds in self.total_cds:
                self.msa_len_li.append(len(fasta.read(self.required_path+"/"+self.organism+"/protein/msa/"+cds+".fasta")[1][0]))
                
                translational_frameshift_li = []
                msa_build = self.required_path+"/"+self.organism+"/hmm_profiles/msa_build/"+cds+".fasta"
                open_file = open(msa_build, "r")
                head_line = open_file.readline()
                open_file.close()
                # modify_hmmalign() requires indices of modified residues relative to build_msa
                frame_mod = re.findall("#frameshift_modification=\((.*)\)", head_line)
                if len(frame_mod) > 0:
                    translational_frameshift_li = [int(num.strip()) for num in frame_mod[0].split(",")]
                    translational_frameshift_li = [-translational_frameshift_li[0], translational_frameshift_li[-1]] #remove value used in regex_select()
                    if translational_frameshift_li[0] < 0:
                        translational_frameshift_li = ["[A-Z][a-z]{"+str(translational_frameshift_li[0])+"}$"] + translational_frameshift_li
                    else:
                        translational_frameshift_li = ["[A-Z]\-{"+str(translational_frameshift_li[0])+"}$"] + translational_frameshift_li
                self.translational_frameshift_di[cds] = translational_frameshift_li
                
    
        # https://stackoverflow.com/questions/3519565/find-the-indexes-of-all-regex-matches
    
        def match_regex(self, dna_str, choose_match_li):
            '''Match regular expressions in protein.json to determine position of features'''
            # Note - regex can return more than one positive match
            
            total_match_di = {}
            # merge identical regex, which is possible from overlapping features
            regex_di3 = {}
            for cds in choose_match_li:
                if cds in self.match_residues_di:
                    for cds_di in self.match_residues_di[cds]:
                        if cds_di["regex"] not in regex_di3:
                            regex_di3[cds_di["regex"]] = {}
                        regex_di3[cds_di["regex"]][cds] = dict(cds_di)
                        del regex_di3[cds_di["regex"]][cds]["regex"]
            # match regex against translated (three frames) input nucleotide sequence     
            protein_li = [string.translate_nucleotides(dna_str, i) for i in range(3)]
            for regex in regex_di3:
                for reading_frame, protein in enumerate(protein_li):
                    match_li2 = [(m.start(), m.group()) for m in re.finditer(regex, protein)]
                    for cds in regex_di3[regex]:
                        if cds not in total_match_di:
                            total_match_di[cds] = []
                        for match in match_li2:
                            total_match_di[cds].append([reading_frame, match[0], regex_di3[regex][cds]["aa_position"], regex_di3[regex][cds]["nucleotide_position"], match[1]])
            
            # Attempt to force fit for start/stop if not already matched
            for regex in regex_di3:
                for cds in regex_di3[regex]:
                    aa_position = regex_di3[regex][cds]["aa_position"]
                    if aa_position == self.match_residues_di[cds][0]["aa_position"] or aa_position == self.match_residues_di[cds][-1]["aa_position"]:
                        if cds in total_match_di and aa_position not in [x[2] for x in total_match_di[cds]]: #minimum regex match count requirement, 1
                            for reading_frame, protein in enumerate(protein_li):
                                modify_regex_li = [re.findall("\[[A-Za-z]+\]|[A-Za-z]\??", x) for x in regex.split("|")]
                                for split_index in range(15): # maximum window size is 15; regular expression may be shorter and vary in length
                                    mod_regex = "|".join(["".join(x[:split_index]+["."]+x[split_index+1:]) for x in modify_regex_li])
                                    match_li2 = [(m.start(), m.group()) for m in re.finditer(mod_regex, protein)]
                                    for cds in regex_di3[regex]:
                                        for match in match_li2:
                                            add = [reading_frame, match[0], regex_di3[regex][cds]["aa_position"], regex_di3[regex][cds]["nucleotide_position"], match[1]]
                                            if add not in total_match_di[cds]:
                                                total_match_di[cds].append(add)
            # reading frame, translated aa position (as index), position relative to aa conseq, position relative to hmm profile, match string
            return(total_match_di)
    
        
        def modify_hmmalign(self, segment_key, dna_str, separate_match_li2, start, stop):
            '''Prepares and processes sequence for hmmalign, interprets result'''
            
            def end_overlap(s1, s2):
                # match end of s1 to beginning of s2
                for i in range(1, len(s2)+1):
                    if s1[-i:] == s2[:i]:
                        return(i)
                    elif i == len(s1):
                        break
                return(0)
                
            dna_str = dna_str.upper().replace("U", "T").replace("-", "")
            # check for reading frame changes in matched region/s
            match_li2, match_li3 = [separate_match_li2[0]], []
            for i, match in enumerate(separate_match_li2[1:]):
                if match[0] == match_li2[0][0]:
                    match_li2.append(match)
                else:
                    # check that no overlap occurs when reading frame changes
                    front_index = match_li2[-1][1] + len(match_li2[-1][-1]) - end_overlap(match_li2[-1][-1], match[-1])
                    back_index = None
                    for x in separate_match_li2[i+1:]:
                        if x[0] == match_li2[-1][0]:
                            back_index = x[1] - len(x[-1]) + end_overlap(match[-1], x[-1])
                    if match[1] >= front_index and (back_index == None or match[1] <= front_index):
                        match_li3.append(match_li2)
                        match_li2 = [match]
            match_li3.append(match_li2)
            translate_str = "".join([string.translate_nucleotides(dna_str[x[0][1]*3+x[0][0]:(x[-1][1]+len(x[-1][4]))*3+x[-1][0]], 0) for x in match_li3])
                            
            translational_frameshift_li = self.translational_frameshift_di[segment_key][:]
            hmm_build = self.required_path+"/"+self.organism+"/hmm_profiles/hmm_build/"+segment_key+".hmm"
            hmmalign_input = self.temporary_path+"/hmmalign_input.fasta"
            
            # tracks the corrected homologous sequence length relative to cropped alignments in else statement
            # necessary for finding location of translational frameshift/s; unable to reliably determine reading frame
            crop_li3 = [] #crop_dna_li, crop_msa_build_li
            trim_sequence_str, end_li = dna_str[start:stop], ["", ""]
            if translate_str[:-1].count("x") != 0 or separate_match_li2[0][2] != 0 or separate_match_li2[-1][2] != 1:
                read_fasta_str = ">Sequence\n"+dna_str[start:stop]
                function_di = {
                    alignment.hmm.align: [hmmalign_input, hmm_build, 3] # force align up to 3 unmatched terminal residues
                }
                trim_sequence_str, end_li = nested.temporary(hmmalign_input, read_fasta_str, function_di)[0]
                
            elif len(match_li3) > 1:
                in_sequence_li = []
                for match_index, match_li2 in enumerate(match_li3[:-1]):
                    # crop input and crop/rebuild temporary hmm profile
                    crop_li2 = [match_li2[-1], match_li3[match_index+1][0]]
                    # crop_dna_li and crop_msa_build_li should be approximately equal lengths (treatments for when either is longer)
                    crop_dna_li = [crop_li2[0][1]*3+crop_li2[0][0], (crop_li2[1][1]+len(crop_li2[1][4]))*3+crop_li2[1][0]]
                    dna_segment_str = dna_str[crop_dna_li[0]:crop_dna_li[1]]
                    crop_msa_build_li = [crop_li2[0][3], crop_li2[1][3]+len(crop_li2[1][-1])*3]
                    crop_li3.append([crop_dna_li, crop_msa_build_li])
                                    
                    read_fasta_li2 = fasta.read(self.required_path+"/"+self.organism+"/hmm_profiles/msa_build/"+segment_key+".fasta")
                    read_fasta_li2[1] = [x[crop_msa_build_li[0]:crop_msa_build_li[1]] for x in read_fasta_li2[1]]
                    read_fasta_str = fasta.list2_to_string(read_fasta_li2)
                    msa_build = self.temporary_path+"/crop.fasta"
                    hmm_build = self.temporary_path+"/crop.hmm"
                    function_di = {
                        alignment.hmm.build: [msa_build, hmm_build], 
                        nested.temporary: [hmmalign_input, ">Sequence\n"+dna_segment_str, {
                            # force align up to all unmatched terminal residues
                            alignment.hmm.align: [hmmalign_input, hmm_build, crop_msa_build_li[1]-crop_msa_build_li[0]]    
                        }]
                    }
                    hmm_out = nested.temporary(msa_build, read_fasta_str, function_di)[-1][0]
                    in_sequence_str = hmm_out[0].join(hmm_out[1]) # combine trim_sequence_str and end_li
                
                    # ignore/remove terminal gaps and update crop_li3 (crop_msa_build_li section)
                    for i, pattern in enumerate(["^\-*", "\-*$"]):
                        	terminal_match, scale_end_li = re.findall(pattern, in_sequence_str)[0], [1, -1]
                        	if len(terminal_match) > 0: #hmm_build too long
                        		crop_li3[-1][1][i] += scale_end_li[i] * len(terminal_match)
                        		in_sequence_str = re.sub(pattern, "", in_sequence_str)
                    # make regex matched residues uppercase
                    for i, pattern in enumerate([["^(?:[^a-zA-Z]*[a-zA-Z]){"+str(len(crop_li2[0][-1])*3)+"}", "^(?:[\-]*[a-z])*"],
                                                 ["(?:[^a-zA-Z]*[a-zA-Z]){"+str(len(crop_li2[1][-1])*3)+"}$", "(?:[\-]*[a-z])*$"]]):
                        # find characters matched by match_li3 regex
                        terminal_match_li = re.findall(pattern[0], in_sequence_str)
                        if len(terminal_match_li) > 0:
                            	terminal_match, scale_end_li = re.findall(pattern[1], terminal_match_li[0])[0], [1, -1]
                            	lower_count = len([x for x in terminal_match if x.islower()])
                            	if lower_count > 0: # dna_segment too long
                            		in_sequence_str = re.sub("^"*abs(i-1)+terminal_match+"$"*i, terminal_match.upper(), in_sequence_str)
                    # compares total potential frame changes in aligned area to known reading frame of matched regex
                    sub_gap = (in_sequence_str.count("-") - len([x for x in in_sequence_str if x.islower()]))%3
                    frame_change = crop_li2[0][0] - crop_li2[1][0]
                    if sub_gap != frame_change and sub_gap != 3-frame_change and sub_gap != 3+frame_change:
                        in_sequence_str = in_sequence_str.replace("-", "").upper()
                        sub_gap = frame_change
                    # adds alignment correction (end treatments removed aligned information)
                    if in_sequence_str.replace("-", "").upper() == in_sequence_str: # no misaligned nucleotides for frameshift
                        sub_gap = crop_li2[0][0] - crop_li2[1][0]
                        crop_location = len(crop_li2[0][4])*3
                        if sub_gap == -1 or sub_gap == 2: # add gap
                            in_sequence_str = in_sequence_str[:crop_location] + in_sequence_str[crop_location].lower() + in_sequence_str[crop_location+1:]
                        elif sub_gap == 1 or sub_gap == -2: # sub lower
                            in_sequence_str = in_sequence_str[:crop_location] + "-" + in_sequence_str[crop_location:]
                    in_sequence_li.append(in_sequence_str)
                for index, in_sequence_str in enumerate(in_sequence_li[:-1]):
                    # overlap correction
                    overlap_int = crop_li3[index][0][1] - crop_li3[index+1][0][0]
                    if overlap_int > 0:
                        if in_sequence_str[-overlap_int:].isupper(): # end uppercase
                            in_sequence_li[index] = in_sequence_str[:-overlap_int]
                            crop_li3[index][0][1] -= overlap_int
                        elif in_sequence_li[index+1][:overlap_int].isupper(): # beginning uppercase
                            in_sequence_li[index+1] = in_sequence_li[index+1][overlap_int:]
                            crop_li3[index+1][0][0] += overlap_int
                        else:
                            print("overlap correction error: evidence of frameshift found in both sequences", file=sys.stderr) #?
                # merge aligned and unaligned sections of dna_str
                flat_crop_dna_li = [start] + [y for x in crop_li3 for y in x[0]] + [stop]
                out_sequence_li = [dna_str[flat_crop_dna_li[x]:flat_crop_dna_li[x+1]] for x in range(0, len(flat_crop_dna_li), 2)]
                trim_sequence_str = out_sequence_li[0]
                for index in range(len(in_sequence_li)):
                    trim_sequence_str += in_sequence_li[index]
                    trim_sequence_str += out_sequence_li[index+1]
                # extract end_li from trim_sequence_str
                end_li = [re.findall("^[-a-z]*", trim_sequence_str)[0], re.findall("[-a-z]*$", trim_sequence_str)[0]]
                trim_sequence_str = trim_sequence_str[len(end_li[0]):len(trim_sequence_str)-len(end_li[1])]
                
            # fragmentation of the modified input nucleotide sequence
            # lowercase characters and gaps ("-") in the hmm alignemnt are indicative of low homology regions where reading frame shifts are most likely
            modified_nucleotide_li = re.findall("[A-Z]+[a-z\-]*", trim_sequence_str)
            # adjust start value
            start += len("".join(re.findall("[a-zA-Z]", end_li[0])))
            # determine reading frame based on end_li[0] (front)
            match_li2 = [x for x in separate_match_li2 if x[1]*3+x[0] > start]
            if len(match_li2) == 0:
                match_li2 = [separate_match_li2[-1]]
            
            # determine reading frame based on HMMER alignment
            reading_frame = (3 - len([x for x in end_li[0] if not x.islower()]) % 3) % 3
            """
            # determine reading frame based on first regex match
            reading_frame = match_li2[0][0]
            start_to_match = trim_sequence_str[:match_li2[0][1]*3+reading_frame - start]
            reading_frame = (reading_frame + len([x for x in start_to_match if not x.islower()])) % 3
            """
            frameshift_li2 = [] # start, stop, type (translational=0, error=1)
            read_len_correction = 0
            for modified_nucleotide_index, modified_nucleotide_str in enumerate(modified_nucleotide_li[:-1]):
                # homolog frame must be preserved; lower cannot be part of substitution with upper or gap
                homolog_len = len(re.sub("[a-z]", "", "".join(modified_nucleotide_li[:modified_nucleotide_index+1])))
                read_len = len("".join(modified_nucleotide_li[:modified_nucleotide_index+1]).replace("-", ""))
                corrected_read_len = read_len + read_len_correction 
                # correct_homolog_len references crop_msa_build_li in crop_li3; parts of modified_nucleotide_li may not have been hmmalign-ed
                pseudo_align_len = corrected_read_len
                for crop_index, crop_li2 in enumerate(crop_li3[::-1]): #[crop_dna_li, crop_msa_build_li]
                    if crop_li2[0][0]-start < pseudo_align_len:
                        pseudo_align_len -= crop_li2[0][0]-start - crop_li2[1][0]
                        break
                # unique condition for translational frameshits
                if len(translational_frameshift_li) == 3 \
                and len(re.findall(translational_frameshift_li[0], modified_nucleotide_str)) > 0 \
                and abs(translational_frameshift_li[2] - pseudo_align_len) <= 1:  # aligned length, target location, indel length correction
                    location_difference = translational_frameshift_li[2] - pseudo_align_len
                    slip_len = translational_frameshift_li[1]
                    if slip_len > 0: # check slip direction (alternative direction not currently needed)
                        modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_str[:-slip_len] #crop out gap/s
                        # corrective measure if hmmalignment slightly misplaces site of translational frameshift
                        if location_difference > 0: # move frameshift upstream
                            modified_nucleotide_li[modified_nucleotide_index] += modified_nucleotide_li[modified_nucleotide_index+1][:location_difference]
                            modified_nucleotide_li[modified_nucleotide_index+1] = modified_nucleotide_li[modified_nucleotide_index+1][location_difference:]
                        elif location_difference < 0: # move frameshift downstream
                            modified_nucleotide_li[modified_nucleotide_index+1] = modified_nucleotide_li[modified_nucleotide_index][location_difference:] + modified_nucleotide_li[modified_nucleotide_index+1]
                            modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_li[modified_nucleotide_index][:location_difference]
                        corrected_read_len += location_difference
                        modified_nucleotide_li[modified_nucleotide_index] += modified_nucleotide_li[modified_nucleotide_index][-slip_len:]
                        frameshift_li2.append([corrected_read_len, corrected_read_len-slip_len, 0])
                        read_len_correction -= slip_len
                else: # must read next segment for presence of stop codons
                    # reading frame determination
                    sub_gap = (homolog_len-read_len)%3
                    if sub_gap != 0:
                        end_gap_str = re.findall("[a-z\-]*$", modified_nucleotide_str)[0]
                        if len(end_gap_str) > 0:
                            modified_nucleotide_str = modified_nucleotide_str[:-len(end_gap_str)]
                        # add gap/s to preserve homolog_len; track changes in read_len with slip_len variable
                        if sub_gap == 2: # or sub_gap == -1: # skip nucleotide
                            modified_nucleotide_str = modified_nucleotide_str[:-1]+"-" + end_gap_str
                            frameshift_li = [corrected_read_len-1, corrected_read_len, 1]
                            slip_len = 1
                        elif sub_gap == 1: # or sub_gap == -2: #re-read nucleotide
                            modified_nucleotide_str += modified_nucleotide_str[-1].lower() + end_gap_str
                            frameshift_li = [corrected_read_len, corrected_read_len-1, 1]
                            slip_len = -1
                        read_str = "".join(modified_nucleotide_li[:modified_nucleotide_index+2])
                        sub_str = "".join(modified_nucleotide_li[:modified_nucleotide_index])+modified_nucleotide_str+modified_nucleotide_li[modified_nucleotide_index+1]
                        translate_read_str = string.translate_nucleotides(read_str.replace("-", ""), reading_frame)
                        translate_sub_str = string.translate_nucleotides(sub_str.replace("-", ""), reading_frame)
                        if translate_sub_str[:-1].count("x") < translate_read_str[:-1].count("x"):
                            modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_str
                            frameshift_li2.append(frameshift_li)
                            read_len_correction += slip_len
            modified_nucleotide_str = "".join(modified_nucleotide_li).replace("-", "").upper()
            # correct frameshift locations based on adjusted start
            frameshift_li2 = [[x[0]+start, x[1]+start, x[2]] for x in frameshift_li2]
            nucleotide_out_li = end_li
            return(modified_nucleotide_str, frameshift_li2, nucleotide_out_li, reading_frame)
    
        
        def main(self, dna_str, choose_match_li):
            '''Processes input sequence based on provided list of possible feature/s'''
            fix_dna_str = dna_str.upper().replace("U", "T").replace("-", "")
            total_annotate_li2 = []
            total_match_di = self.match_regex(fix_dna_str, choose_match_li)
            for key, match_li2 in total_match_di.items():
                match_li2.sort(key=lambda x: x[1])
                if len(match_li2) > 1: # check if enough conserved matches, minimum 2
                    # separate values that fall out of order (possible duplicate sequences or false matches)
                    separate_match_li3 = [[match_li2[0]]]
                    for value in match_li2:
                        if value not in separate_match_li3[-1]:
                            if value[2] > separate_match_li3[-1][-1][2]:
                                separate_match_li3[-1].append(value)            
                            else:
                                separate_match_li3.append([value])
                    
                    location_li3 = []
                    msa_len = self.msa_len_li[self.total_cds.index(key)]
                    for separate_match_index, separate_match_li2 in enumerate(separate_match_li3):
                        if len(separate_match_li2) > 1: # recheck count requirement post separation, minimum 2
                            match_start = separate_match_li2[0][0]+(separate_match_li2[0][1]*3) # reading frame + amino acid index
                            match_stop = separate_match_li2[-1][0]+(separate_match_li2[-1][1]+len(separate_match_li2[-1][-1]))*3 # reading frame + amino acid index + adjustment to end of regex pattern
                            # estimated start assuming feature is full length; helps match discontinuity (frame change/s)
                            estimate_start = match_start - msa_len*3*(0.10 + separate_match_li2[0][2])
                            estimate_stop = match_stop + msa_len*3*(0.10 + 1-separate_match_li2[-1][2])
                            location_li3.append([[match_start, match_stop], [estimate_start, estimate_stop]])
                        else:
                            separate_match_li3[separate_match_index] = []
                    location_li3 = [[[int(z) for z in y] for y in x] for x in location_li3] # estimated values may be decimals
                    separate_match_li3 = [x for x in separate_match_li3 if x != []]
                    
                    # merge overlapping locations
                    len_location_li3 = None
                    while len_location_li3 == None or len_location_li3 != len(location_li3):
                        for loc1_index, loc1 in enumerate(location_li3[:-1]):
                            loc2 = location_li3[loc1_index+1]
                            if loc1 != [] and loc2 != []:
                                if loc1[0][0] >= loc2[1][0] and loc1[0][1] <= loc2[1][1] or loc2[0][0] >= loc1[1][0] and loc2[0][1] <= loc1[1][1]:
                                    location_li3[loc1_index] = [[min(loc1[0][0], loc2[0][0]), max(loc1[0][1], loc2[0][1])], \
                                                                [min(loc1[1][0], loc2[1][0]), max(loc1[1][1], loc2[1][1])]]
                                    separate_match_li3[loc1_index] = separate_match_li3[loc1_index]+separate_match_li3[loc1_index+1]
                                    location_li3[loc1_index+1], separate_match_li3[loc1_index+1] = [], [] 
                                    break
                        len_location_li3 = len(location_li3)
                        location_li3 = [x for x in location_li3 if x != []]
                        separate_match_li3 = [x for x in separate_match_li3 if x != []]
                    for loc_index, separate_match_li2 in enumerate(separate_match_li3):
                        total_annotate_li = ["" for x in range(7)]
                        total_annotate_li[0] = key+" full"
                        start, stop = location_li3[loc_index][0][0], location_li3[loc_index][0][1]
                        if separate_match_li2[0][2] != 0:
                            # estimated start must be in same reading frame as separate_match_li2[0] value
                            start = location_li3[loc_index][1][0]
                            while start%3 != separate_match_li2[0][0]:
                                start -= 1
                            if start < 0:
                                start = 0
                        if separate_match_li2[-1][2] != 1:
                            stop = location_li3[loc_index][1][1]
                            if stop > len(fix_dna_str):
                                stop = len(fix_dna_str)
                        modified_nucleotide_str, frameshift_li2, nucleotide_out_li, reading_frame = self.modify_hmmalign(key, fix_dna_str, separate_match_li2, start, stop)
                        if len(modified_nucleotide_str) > 45: # output must be at least as long as single regex match (15 aa); shorter indicates poor match
                            start += len(nucleotide_out_li[0].replace("-", ""))
                            stop -= len(nucleotide_out_li[1].replace("-", ""))
                            translate_str = string.translate_nucleotides(modified_nucleotide_str, reading_frame)
                            align_len, crop_len = len(re.sub("[a-z]", "", modified_nucleotide_str)), len(re.sub("[a-z]", "", "".join(nucleotide_out_li)))
                            if crop_len > 0:
                                if crop_len/(align_len+crop_len) < 0.05:
                                    total_annotate_li[0] = key+" near complete"
                                else:
                                    total_annotate_li[0] = key+" partial"                           
                            for i in range(2):
                                total_annotate_li[i+3] = "; ".join([",".join([str(x[0]), str(x[1]+1)]) for x in frameshift_li2 if x[-1] == i])
                            
                            # scale output values if input sequence contains gaps
                            if dna_str.count("-") > 0:
                                location_li = [start] + [y for x in frameshift_li2 for y in x[:2]] + [stop]
                                for index, location in enumerate(location_li):
                                    while len(dna_str[:location].replace("-", "")) != location_li[index]:
                                        location += 1
                                    location_li[index] = location
                                start, stop = location_li[0], location_li[-1]
                                frameshift_li2 = [location_li[1:-1][x:x+2]+[frameshift_li2[i][-1]] for i, x in enumerate([y for y in range(0, len(location_li[1:-1]), 2)])]
                            
                            total_annotate_li[1] = str(start+1)+","+str(stop) # add one to start to account for genbank format
                            total_annotate_li[2] = str(reading_frame)
                            for i in range(2):
                                total_annotate_li[i+3] = "; ".join([",".join([str(x[0]), str(x[1]+1)]) for x in frameshift_li2 if x[-1] == i])
                            if translate_str[-1] == "x":
                                translate_str = translate_str[:-1]
                            total_annotate_li[-2:] = str(len(translate_str)), translate_str
                            total_annotate_li2.append(total_annotate_li)
    
            # 0:translation type, 1:start/stop range, 2:reading frame, 3:location of translational frameshift, 4:homology adjustment (error correction), 5:sequence length, 6:AA sequence
            return(total_annotate_li2)
            
    