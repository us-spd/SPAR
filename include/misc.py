
from datetime import datetime
import os
import pandas as pd
import random
import re
from string import ascii_uppercase, ascii_lowercase, digits
import subprocess
import sys
import time
import unicodedata
from functools import lru_cache
from operator import itemgetter

################################################################################

#Stop codons are currently represented by 'X' (ambig residue) to prevent MAFFT auto deletion
dna_convert_aa_di = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'x', 'TAG':'x',
    'TGC':'C', 'TGT':'C', 'TGA':'x', 'TGG':'W'}

ambiguous_dna_di = {
    'Y':['C', 'T'],         'R':['A', 'G'],         'W':['A', 'T'],
    'S':['G', 'C'],         'K':['T', 'G'],         'M':['C', 'A'],
    'D':['A', 'G', 'T'],    'V':['A', 'C', 'G'],    'H':['A', 'C', 'T'],
    'B':['C', 'G', 'T'],    'N':['A', 'C', 'T', 'G']
}

aa_convert_codon_di =  {
    'A':['[GRSK][CYSM].'],
    'B':['[ARWM][ARWM][CTYWKSM]', '[GRSK][ARWM][TCYWKSM]'],
    'C':['[TYWK][GRSK][TCYWKSM]'],
    'D':['[GRSK][ARWM][TCYWKSM]'],
    'E':['[GRSK][ARWM][AGRSKWM]'],
    'F':['[TYWK][TYWK][CTYWKSM]'],
    'G':['[GRSK][GRSK].'],
    'H':['[CYSM][ARWM][TCYWKSM]'],
    'I':['[ARWM][TYWK][^G]'],
    'J':['[ARWM][TYWK][^G]', '[CYSM][TYWK].', '[TYWK][TYWK][AGRSKWM]'],
    'K':['[ARWM][ARWM][AGRSKWM]'],
    'L':['[CYSM][TYWK].', '[TYWK][TYWK][AGRSKWM]'],
    'M':['[ARWM][TYWK][GRSK]'],
    'N':['[ARWM][ARWM][CTYWKSM]'],
    'O':['[TYWK][ARWM][GRSK]'],
    'P':['[CYSM][CYSM].'],
    'Q':['[CYSM][ARWM][AGRSKWM]'],
    'R':['[CYSM][GRSK].', '[ARWM][GRSK][GARSKWM]'],
    'S':['[TYWK][CYSM].', '[ARWM][GRSK][CTYWKSM]'],
    'T':['[ARWM][CYSM].'],
    'U':['[TYWK][GRSK][ARWM]'],
    'V':['[GRSK][TYWK].'],
    'W':['[TYWK][GRSK][GRSK]'],
    'X':['...'],
    'Y':['[TYWK][ARWM][CTYWKSM]'],
    'Z':['[CYSM][ARWM][AGRSKWM]','[GRSK][ARWM][AGRSKWM]'],
    '_':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]'],
    '*':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]'],
    'x':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]']}

################################################################################
################################################################################

class fasta:
    
    def __init__():
        '''The fasta class is initialized'''
        
    def string_to_list2(fasta_str):
        fasta_match = [list(x) for x in re.findall("(>.*\n)([^>]*)", fasta_str)]
        for index, match in enumerate(fasta_match):
            if match[1].replace("\n", "").islower(): #convert all characters to uppercase only if all characters are lowercase
                fasta_match[index][1] = match[1].upper()
        #convert RNA to DNA
        read_fasta_li2 = [[match[0].replace("\n", "").replace(">", "") for match in fasta_match], [match[1].replace("\n", "").replace("U", "T") for match in fasta_match]]
        return read_fasta_li2
    
        
    def list2_to_string(fasta_li2):
        return "".join([">" + fasta_li2[0][x].replace(" ", "_") + "\n" + fasta_li2[1][x] + "\n" for x in range(len(fasta_li2[0]))])
    
    
    def read(input_filepath):
        '''Reads input fasta file and returns read_fasta_li2 structure [[header1, header2, ...], [sequence1, sequence2, ...]]'''
        #supports mixed case (upper/lower); can be used to indicate residues of interest
        #example: "x" is used as (unofficial) notation for stop codons in aminio acid sequences while "X" is for ambiguous residues
        #altered notation prevents errors when using mafft, which does not support standard notation of stop codons ("*")
        read_fasta_li2 = [[], []]
        if os.path.isfile(input_filepath):
            read_file = open(input_filepath, "r")
            read_file_str = "".join([line for line in read_file])
            read_file.close()
            read_fasta_li2 = fasta.string_to_list2(read_file_str)
        else:
            print("file does not exist: "+input_filepath, file=sys.stderr)
        return(read_fasta_li2)
        
    
    def write(output_filepath, read_fasta_li2):
        '''Writes read_fasta_li2 structure [[header1, header2, ...], [sequence1, sequence2, ...]] to provided output filepath''' 
        #create missing directories if path does not exist (requires bash)
        subprocess.call(["mkdir", "-p", re.sub("[^\/]*$", "", output_filepath)])
        write_file = open(output_filepath, "w")
        write_file.write(fasta.list2_to_string(read_fasta_li2))
        write_file.close()
        

    def remove_column(sequence_li, gap_proportion):
        '''removes gaps from multiple sequence alignment file based on specified gap proportion'''
        #remove columns based on allowed gap proportion; track alterations
        #mafft alignments may produce columns that contain all gaps
        #hmm profile alignments do not track excessively gapped (>50%) regions
        column_index = 0 #index conversion from input to output: index of input minus indexed value of modify_li equals the index of output
        modify_li = [0]
        #avoid overwriting input list
        read_fasta_sequences = sequence_li[:]
        while column_index < len(read_fasta_sequences[0]):
            column_li = [read_fasta_sequences[x][column_index] for x in range(len(read_fasta_sequences))]
            if column_li.count("-") >= len(column_li)*gap_proportion:
#               print("column removed", {x:column_li.count(x) for x in set(column_li)}, file=sys.stderr)
                for row_index in range(len(read_fasta_sequences)):
                    read_fasta_sequences[row_index] = read_fasta_sequences[row_index][:column_index] + read_fasta_sequences[row_index][column_index+1:]
                modify_li[-1] += 1
            else:
                column_index += 1
            modify_li.append(modify_li[-1])
        return(read_fasta_sequences, modify_li[:-1])
    
################################################################################
        
class excel:
    
    def __init__():
        '''The excel class is initialized'''

    def read(input_file):
        '''Reads input excel file and returns pandas dataframe'''
        return pd.read_excel(input_file, dtype=object, keep_default_na=False).astype(str)
                
        
    def write(output_file, dataframe):
        '''Writes pandas dataframe to provided output filepath''' 
        writer = pd.ExcelWriter(output_file, engine="xlsxwriter", options={"strings_to_numbers": False})
        dataframe.fillna("").to_excel(writer, index=False)
        writer.save()
        
################################################################################
    
class multiple:
    
    def __init__():
        '''The multiple class is initialized'''

    def merge(excel_file, fasta_file, fasta_head):
        '''Combines input excel file and fasta file into pandas dataframe'''
        if os.path.exists(excel_file) and os.path.exists(fasta_file):
            read_fasta_li2 = fasta.read(fasta_file)
            transpose_read_fasta_li2 = [[row[i] for row in read_fasta_li2] for i in range(len(read_fasta_li2[0]))]
            sequence_df = pd.DataFrame(transpose_read_fasta_li2, columns=[fasta_head, "sequence"])
            metadata_df = excel.read(excel_file)
            merge_df = pd.merge(metadata_df, sequence_df, on=fasta_head)
        else:
            merge_df = None
        return(merge_df.fillna(""))
        
        
    def split(excel_file, fasta_file, merge_df, fasta_head):
        '''Separates pandas dataframe into excel file and fasta file'''
        copy_df = merge_df.copy().fillna("")
        read_fasta_li2 = [list(copy_df[fasta_head]), list(copy_df["sequence"])]
        fasta.write(fasta_file, read_fasta_li2)
        copy_df.drop(["sequence"], axis=1, inplace=True)
        excel.write(excel_file, copy_df)

################################################################################
        
class gff3:
    
    def __init__():
        '''The gff3 class is initialized'''
    
    # https://m.ensembl.org/info/website/upload/gff3.html
    # http://gmod.org/wiki/GFF3#GFF3_Format
    # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    # https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
        
    def write(read_fasta_li2, annotate_li4, required_path):
        '''Returns gff3 format string structure from annotated input and SPAR required path'''
######### does not distinguish duplicate features (should recieve unique ids)
        
        # accepts annotated output (non-standard format) and uses reference template to produce gff3 file (standard format)
        organism_li = list(set([x[0] for x in annotate_li4 if len(x) == 2]))
        # http://www.sequenceontology.org/so_wiki/index.php/GFF3_character_encoding
        # SO calls for limited percent encoding and escaping specific reserved characters
        copy_read_fasta_li2 = [x[:] for x in read_fasta_li2]
        copy_read_fasta_li2[0] = [string.escape(x, [";", "=", "&", ","]) for x in copy_read_fasta_li2[0]]
        copy_read_fasta_li2[0] = [string.percent_encode(x, ["\t", "\n", "%"]) for x in copy_read_fasta_li2[0]]
        gff_di = {}
        for organism in organism_li:
            gff_path = required_path+organism+"/gff3/template.gff3"
            if os.path.isfile(gff_path):
                read_gff = open(gff_path, "r")
                read_gff_li2 = [line.strip().split("\t") for line in read_gff]
                read_gff.close()
                gff_di[organism] = read_gff_li2
            else:
                print("GFF3 template file could not be found for "+organism, file=sys.stderr)        
        gff_output_str = "##gff-version 3\n"
        for fasta_index, annotate_li3 in enumerate(annotate_li4):
            if len(annotate_li3) == 2 and annotate_li3[0] in gff_di:
                organism, annotate_li2 = annotate_li3
                # add organism header/region
                gff_output_str += "# organsim "+organism+"\n"
                annotate_li2.append([organism, str(1)+","+str(len(read_fasta_li2[1][fasta_index])), "0"]+[""]*4)
                feature_li = [re.sub(" full| near complete| partial", "", x[0]) for x in annotate_li2]
                gff_li3 = [[x] for x in gff_di[organism]]
                for gff_index, gff_li2 in enumerate(gff_li3):
                    gff_template_li = gff_li2[0][:]  #first list always template (maintain order)
                    #0:translation type, 1:start/stop range, 2:reading frame, 3:location of translational frameshift, 4:homology adjustment (error correction), 5:sequence length, 6:AA sequence
                    for feature_index, annotate_li in enumerate(annotate_li2):
                        feature = feature_li[feature_index]
                        note_str = re.sub(feature, "", annotate_li[0]).strip() # completeness
                        find_parent = re.findall("parent\=([^\;]*)|$", gff_template_li[-1], flags=re.I)[0]
                        if len(re.findall("[^a-zA-Z0-9]"+feature+"[^a-zA-Z0-9]", gff_template_li[-1].replace(find_parent, ""))) > 0: # regex should be specific enough for single character feature names
                            # apply feature annotations to gff file line
                            gff_template_li[0] = copy_read_fasta_li2[0][fasta_index]
                            gff_template_li[1] = "Swine Pathogen Analysis Resource (SPAR)"
                            # edit/modify note
                            if len(annotate_li[3]) > 0:
                                translational_frameshift_li = [int(x.strip()) for x in annotate_li[3].split(",")] #assumes maximum of one translational frameshift per genomic feature
                                note_str += ", "+"{:+}".format(translational_frameshift_li[0]-translational_frameshift_li[1]-1)+" translational frameshift at "+str(translational_frameshift_li[0])
                            if len(annotate_li[4]) > 0:
                                note_str += ", adjustments made to account for premature stop codon/s"
                            note_index = gff_template_li[-1].upper().find("NOTE=")
                            if note_index >= 0:
                                gff_template_li[-1] = gff_template_li[-1][:note_index]+note_str+", "+gff_template_li[-1][note_index:]
                            elif len(note_str) > 0:
                                gff_template_li[-1] += ";Note="+note_str
                                
                            #translational frameshift and homology adjustment treatment
                            frameshift_str = annotate_li[3] + "; " + annotate_li[4]
                            frameshift_li2 = [[y.strip() for y in x.split(",")] for x in frameshift_str.split(";") if x.strip() != ""]
                            frameshift_li2 = sorted(frameshift_li2, key=lambda x: int(x[0]))
                            
                            start, stop = tuple([x.strip() for x in annotate_li[1].split(",")])
                            frameshift_li = [start]+[x for y in frameshift_li2 for x in y]+[stop]
                            frameshift_li2 = [[frameshift_li[x], frameshift_li[x+1]] for x in range(0, len(frameshift_li), 2)]
                            index_gff_li2 = [gff_template_li[:3]+x+gff_template_li[5:] for x in frameshift_li2]
                            if annotate_li[2] != "0":
                                index_gff_li2[0][-2] = str(annotate_li[2])
                            gff_li3[gff_index] += index_gff_li2
                           
                            # missing parent treatment (inference)
                            if len(find_parent) > 0:
                                parent_li = find_parent.split(",")
                                for parent in parent_li:
                                    #check that parent does not exist in annotate_li2
                                    if len([x for x in feature_li if len(re.findall("(?:^|[^a-zA-Z0-9])"+x+"(?:[^a-zA-Z0-9]|$)", parent)) > 0]) == 0:
                                        #find appropriate line in gff3 for gene
                                        for index, li2 in enumerate(gff_li3):
                                            template_li = li2[0][:]
                                            if parent in re.sub("parent\=([^\;]*)", "", template_li[-1], flags=re.I):
                                                if template_li[0] != copy_read_fasta_li2[0][fasta_index]:
                                                    template_li[0] = copy_read_fasta_li2[0][fasta_index]
#                                                   template_li[1] = "United States Swine Pathogen Database" #not supported
                                                    #edit/modify note
                                                    note_index = template_li[-1].upper().find("NOTE=")
                                                    if note_index >= 0:
                                                        template_li[-1] = template_li[-1][:note_index]+"Infered, "+template_li[-1][note_index:]
                                                    else:
                                                        template_li[-1] += ";Note=Infered"
                                                
                                                index_gff_li2 = [template_li[:3]+x+template_li[5:] for x in frameshift_li2]
                                                if annotate_li[2] != 0:
                                                    index_gff_li2[0][-2] = annotate_li[2]
                                                gff_li3[index] += index_gff_li2
                                                #order and merge values
                                                gff_li3[index] = [gff_li3[index][0]] + sorted(gff_li3[index][1:], key=lambda x: int(x[3]))
                                                i = 1 #first list is a reference that should not be changed/deleted
                                                while i < len(gff_li3[index])-1:
                                                    if int(gff_li3[index][i+1][3]) - int(gff_li3[index][i][4]) == 1:
                                                        gff_li3[index][i][4] = gff_li3[index][i+1][4]
                                                        gff_li3[index] = gff_li3[index][:i+1]+gff_li3[index][i+2:]
                                                    else:
                                                        i += 1
                # rearrange order; parents first (assumes parents cannot also be children)
                parent_li2 = [x for x in gff_li3 if "parent=" not in x[0][-1].lower()]
                child_li2 = [x for x in gff_li3 if "parent=" in x[0][-1].lower()]
                gff_li3 = parent_li2+child_li2
                for gff_li2 in gff_li3:
                    if len(gff_li2) > 1:
                        # add accession to id/parent
                        for i, gff_li in enumerate(gff_li2[1:]):
                            accession = gff_li[0].split(" ")[0]
                            gff_li2[i+1][0] = accession
                            value_match_li = re.findall("id=.*?;|parent=.*?;", gff_li[-1], re.IGNORECASE)
                            for value_match in value_match_li:
                                split_eq = value_match.split("=")
                                split_eq[1] = ",".join([accession+":"+x for x in split_eq[1].split(",")])
                                gff_li2[i+1][-1] = re.sub(value_match, "=".join(split_eq), gff_li[-1])
                        gff_output_str += "\n".join(["\t".join(x) for x in gff_li2[1:]])+"\n"
        return(gff_output_str)
 
################################################################################
        
class string:

    def __init__():
        '''The string_mod class is initialized'''
        
    def escape(input_str, char_li=[]):
        input_li = [x for x in input_str]
        for i, char in enumerate(input_li):
            if char in char_li:
                input_li[i] = "\\"+char
        return("".join(input_li))
    
    
    def longest_common_substring(x: str, y: str) -> (int, int, int):
        # https://www.geeksforgeeks.org/longest-common-substring-dp-29/
    
        # function to find the longest common substring
    
        # Memorizing with maximum size of the memory as 1
        @lru_cache(maxsize=1) 
    
        # function to find the longest common prefix
        def longest_common_prefix(i: int, j: int) -> int:
    
            if 0 <= i < len(x) and 0 <= j < len(y) and x[i] == y[j]:
                return 1 + longest_common_prefix(i + 1, j + 1)
            else:
                return 0
    
        # digonally computing the subproplems
        # to decrease memory dependency
        def digonal_computation():
    
            # upper right triangle of the 2D array
            for k in range(len(x)):       
                yield from ((longest_common_prefix(i, j), i, j)
                            for i, j in zip(range(k, -1, -1),
                                        range(len(y) - 1, -1, -1)))
    
            # lower left triangle of the 2D array
            for k in range(len(y)):       
                yield from ((longest_common_prefix(i, j), i, j)
                            for i, j in zip(range(k, -1, -1),
                                        range(len(x) - 1, -1, -1)))
    
        # returning the maximum of all the subproblems
        length, i, j = max(digonal_computation(), key=itemgetter(0), default=(0, 0, 0))
        return x[i: i + length]
            
            
    def percent_encode(input_str, char_li=[]):
        # https://en.wikipedia.org/wiki/Percent-encoding
        
        # reserved characters and line spacing
        percent_encode_di = {"!": "%21", "#": "%23", "$": "%24", "%": "%25", "&": "%26", "'": "%27", "(": "%28", ")": "%29", "*": "%2A", "+": "%2B", ",": "%2C", "/": "%2F", ":": "%3A", ";": "%3B", "=": "%3D", "?": "%3F", "@": "%40", "[": "%5B", "]": "%5D",
                             "\n": "%0A", "\t": "%09", " ": "%20"}
        input_li = [x for x in input_str]
        for i, char in enumerate(input_li):
            if char in char_li and char in percent_encode_di:
                input_li[i] = percent_encode_di[char]
        return("".join(input_li))

        
    def remove_accent(input_str, force=False):
        '''remove accents from input string; fails for non-ascii characters'''
        # https://stackoverflow.com/questions/517923/what-is-the-best-way-to-remove-accents-normalize-in-a-python-unicode-string
        
        def method1():
            '''fails to identify strokes'''
            nfkd_form = unicodedata.normalize("NFKD", input_str)
            return("".join([c for c in nfkd_form if not unicodedata.combining(c)]))
        def method2():
            '''able to identify strokes'''
            # break sybols: 
            def sub1(char):
                try:
                    desc = unicodedata.name(char)
                    cutoff = desc.find(" WITH ")
                    if cutoff != -1:
                        desc = desc[:cutoff]
                        char = unicodedata.lookup(desc)
                except:
                    pass
                return(char)
            return("".join([sub1(x) for x in input_str]))
        translate_di = {"Ð": "Dj", "ð": "o", "Þ": "Th", "þ": "th", "æ": "ae"}
        input_str = "".join(translate_di[x] if x in translate_di else x for x in input_str)
        best = input_str.encode("ascii", "ignore").decode("utf-8").strip()
        for method in [method1(), method2()]:
            remove = method.encode("ascii", "ignore").decode("utf-8").strip()
            if len(remove) > len(best):
                best = remove
            if len(method) == len(remove):
                return(remove)
        if force:
            return(best)
        
    
    def readable_list(string_li):
        '''Punctuates list following oxford comma format'''
        #https://stackoverflow.com/questions/53981845/grammatically-correct-human-readable-string-from-list-with-oxford-comma
        string_li = [x.strip() for x in string_li]
        if len(string_li) < 3:
            return " and ".join(map(str, string_li))
        *a, b = string_li
        return f"{', '.join(map(str, a))}, and {b}"

    
    def reverse_compliment(input_dna, force=True):
        reverse_compliment_di = {"A": "T", "T": "A", "G": "C", "C": "G",
                                 "a": "t", "t": "a", "g": "c", "c": "g", 
                                 "Y": "R", "R": "Y", "K": "M", "M": "K",
                                 "y": "r", "r": "y", "k": "m", "m": "k",
                                 "D": "H", "H": "D", "V": "B", "B": "V",
                                 "d": "h", "h": "d", "v": "b", "b": "v",
                                 "W": "W", "S": "S", "N": "N", "-": "-",
                                 "w": "w", "s": "s", "n": "n"}
        input_dna = input_dna.replace("U", "T").replace("u", "t")
        invalid_li = [x for x in input_dna if x not in reverse_compliment_di]
        if len(invalid_li) > 0:
            print("input contains invalid characters")
            if force:
                input_dna = [x if x in reverse_compliment_di else "N" for x in input_dna]
                return("".join([reverse_compliment_di[x] for x in input_dna])[::-1])
        else:
            return("".join([reverse_compliment_di[x] for x in input_dna])[::-1])
    
    
    def translate_nucleotides(input_dna, reading_frame=0):
        '''Translates nucleotide sequence into corresponding polypeptide'''
        #input_dna in string format, reading_frame enter 0, 1, or 2
        #translates nucleotides sequence into protein (stop codons = "x"); able to infer ambiguous resides
        reading_frame = reading_frame%3
        frame_input = input_dna[reading_frame:].upper()
        frame_input = frame_input.replace('U', 'T')
        frame_input = frame_input.replace('-', '')
        translate_str = ''
        if len(frame_input)%3 != 0:
            frame_input = frame_input[:-(len(frame_input)%3)]
        for i in range(0, len(frame_input), 3):
            if frame_input[i:i+3] in dna_convert_aa_di:
                translate_str += dna_convert_aa_di[frame_input[i:i+3]]
            else:
                ambiguous_combination_li = ['']
                for char in frame_input[i:i+3]:
                    if char in ambiguous_dna_di:
                        ambiguous_combination_li = ambiguous_combination_li*len(ambiguous_dna_di[char])
                        z = 0
                        for y in range(len(ambiguous_dna_di[char])):
                            for x in range(len(ambiguous_combination_li)//len(ambiguous_dna_di[char])):
                                ambiguous_combination_li[z] += ambiguous_dna_di[char][y]
                                z += 1
                    else:
                        for c in range(len(ambiguous_combination_li)):
                            ambiguous_combination_li[c] += char
                translate_check = dna_convert_aa_di[ambiguous_combination_li[0]]
                for codon in ambiguous_combination_li[1:]:
                    if dna_convert_aa_di[codon] != translate_check:
                        translate_check = "fail"
                if translate_check == "fail":
                    translate_str += "X"
                else:
                    translate_str += translate_check
        return(translate_str)
    
    
    def unique_directory(base_path):
        '''Build random, non-existent directory name'''
        while True:
            temp_id = ''.join(random.choice(ascii_uppercase + ascii_lowercase + digits) for r in range(10))
            if not os.path.isdir(base_path+"temp"+temp_id):
                break
        return(base_path+"temp"+temp_id)


    def convert_bool(v):
        # https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
        if isinstance(v, bool):
           return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        
################################################################################
            
class validate:

    def __init__():
        '''The validate class is initialized'''
    
    def date(date_str, pattern="%d-%b-%Y"):
        try:
            if date_str != datetime.strptime(date_str, pattern).strftime(pattern):
                raise ValueError
            return True
        except ValueError:
            return False
        
################################################################################
        
class message:
    
    def __init__():
        '''The message class is initialized'''
        
    def proceed(input_value="Additional setup required. Continue? (y/n) ", \
                response_y="You have selected to continue. Please wait.", \
                response_n="Process terminated."):        
        old_stdout = sys.stdout 
        try:
            sys.stdout = sys.stderr
            while True:
                output = input(input_value)
                if output.lower() == "y" or output.lower() == "yes":
                    if len(response_y) > 0: 
                        print(response_y)
                    return(True)
                elif output.lower() == "n" or output.lower() == "no":
                    if len(response_n) > 0: 
                        print(response_n)
                    return()
                else:
                    print("Invalid response.")  
        finally:
            sys.stdout = old_stdout
        
################################################################################
        
class nested:
    
    def __init__():
        '''The nested class is initialized'''
        
    def temporary(file_path, file_str, function_di):
        '''Creates temporary files from each string within file_li; runs designated dictionary of functions/parameters; removes temporary files'''
        open_file = open(file_path, "w")
        open_file.write(file_str)
        open_file.close()
        output_li = []
        for func in function_di:
            output_li.append(func(*function_di[func]))
        os.remove(file_path)
        return output_li
    
    def delay(function_di, delay, start_time=None):
        '''Ensures that functions take at least a set amount of time to run'''
        output_li = []
        if start_time == None:
            start_time = time.time()
        for func in function_di:
            output_li.append(func(*function_di[func]))
            calc_delay = delay + start_time - time.time()
            if calc_delay > 0:
                time.sleep(calc_delay)
        return output_li

################################################################################

class alignment:
    
    def __init__():
        '''The align class is initialized'''
    
    class mafft:
        
        def __init__():
            '''The mafft subclass is initialized'''
    
        def align(input_fasta, output_fasta):
            '''Performs mafft alignment on provided input fasta and saves to output fasta path (bash required)'''
            #uses bash to call mafft to perform a residue (nucleotide or amino acid) alignment
            mafft_output = re.sub("\.[^\.]*?$", "-mafft.fasta", input_fasta) #allows output_fasta to be the same as input_fasta
            os_call = "mafft "+input_fasta+" > "+mafft_output+" 2>/dev/null"
            subprocess.call(os_call, shell=True)
            read_msa_li2 = fasta.read(mafft_output)
            read_msa_li2[1] = [x.upper() for x in read_msa_li2[1]] #force mafft sequence output to be uppercase (expected mafft format)
            #MAFFT output is not case sensitive
            #AA translation function produces "x" as stop codon to avoid breaking scripts
            read_fasta_li2 = fasta.read(input_fasta)
            for msa_index, head in enumerate(read_msa_li2[0]):
                fasta_index = read_fasta_li2[0].index(head)
                for index1, character1 in enumerate(read_fasta_li2[1][fasta_index]):
                    if character1.islower():
                        for index2, charater2 in enumerate(read_msa_li2[1][msa_index][index1:]):
                            if len(read_msa_li2[1][msa_index][:index1+index2].replace("-", "")) == index1:
                                read_msa_li2[1][msa_index] = read_msa_li2[1][msa_index][:index1+index2] + \
                                read_msa_li2[1][msa_index][index1+index2].lower() + \
                                read_msa_li2[1][msa_index][index1+index2+1:]
            #MAFFT alignment can produce entirely gapped columns          
            read_msa_li2[1] = fasta.remove_column(read_msa_li2[1], 1)[0]
            fasta.write(output_fasta, read_msa_li2)
            os.remove(mafft_output)
            
        
    class hmm:
        
        def __init__():
            '''The hmm subclass is initialized'''
    
        def build(msa_input, hmm_profile):
            '''Constructs hmm profile using from reference input multiple sequence alignment using hmmbuild (bash required)'''
            os_call = ["hmmbuild", hmm_profile, msa_input]
            subprocess.call(os_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
        
        def align(fasta_input, hmm_profile, n=0):
            '''Performs hmmalign on input fasta using provided hmm profile (bash required)'''
            # force align "n" number of unmatched terminal residues
            # reference: https://github.com/EddyRivasLab/hmmer/issues/254
            
            os_call = ["hmmalign", hmm_profile, fasta_input]
            read_str = subprocess.check_output(os_call, stderr=subprocess.DEVNULL).decode("utf-8")
            read_sequence_str = "".join([x[8:].strip() for x in read_str.split("\n") if "Sequence" in x and "#" not in x])
            #treatment for terminal ends
            #force fit condition: non-matching (by hmmalign) ends are re-defined as matching if enough nucleotides/gaps are available
            end_li = [re.findall("^[-a-z]*", read_sequence_str)[0], re.findall("[-a-z]*$", read_sequence_str)[0]]
            trim_sequence_str = read_sequence_str[len(end_li[0]):len(read_sequence_str)-len(end_li[1])]
            if n > 0:
                for i, end in enumerate(end_li):
                    lower_li = [x for x in end if x.islower()]
                    lower_count, gap_count = len(lower_li), end.count("-")
                    # must be more gaps than unmatched
                    if lower_count <= gap_count and lower_count <= n:
                        upper_str = "".join(lower_li).upper()
                        trim_sequence_str  = upper_str*(-i+1) + trim_sequence_str + upper_str*(-i+1)
                        end_li[i] = (gap_count-lower_count)*"-"
            return(trim_sequence_str, end_li)
                

    class blast:
        
        def __init__():
            '''The blast subclass is initialized'''
        
        def build(fasta_input):
            '''Constructs BLAST database under same path/name as input multiple sequence alignment fasta'''
            os_call = ["makeblastdb", "-dbtype=nucl", "-in="+fasta_input]
            subprocess.call(os_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            