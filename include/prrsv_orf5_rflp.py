
import os
import pandas as pd
import re


if __name__ == "__main__":
    from misc import alignment, nested
else:
    from .misc import alignment, nested
    
################################################################################
    
spar_path = re.sub("/include$", "", os.path.dirname(os.path.realpath(__file__)))
    
################################################################################
################################################################################
            
class prrsv_orf5_rflp:
    def __init__(self, temp_dir=""):
        '''The prrsv_orf5_rflp class is initialized'''
        
        self.temporary_path = spar_path + "/temporary"
        if os.path.isdir(temp_dir):
            self.temporary_path = temp_dir
        df = pd.read_excel(spar_path + "/resources/prrsv_orf5_cut_pattern_determination.xlsx")
        self.crop_head_li = [i[i.index("= ")+2:].replace("^", "").replace("(", "[").replace(")", "]") for i in list(df) if "Unnamed" not in i]
        self.cut_site_li = [1, 3, 4]
        self.rflp_di = {i:[[]] for i in self.crop_head_li}
        df = df.rename(columns=df.iloc[0]).drop(df.index[0])
        cut_sites_li = df["Cut Sites"].values.tolist()
        for x in range(len(cut_sites_li)):
            for y in range(len(cut_sites_li[x])):
                if str(cut_sites_li[x][y]).lower() != "nan" and str(cut_sites_li[x][y]).lower() != "-":
                    self.rflp_di[self.crop_head_li[y]].append([int(i) for i in str(cut_sites_li[x][y]).split(",")])
    
            
    def main(self, nucleotide_sequence, annotate_li2, require_full=True):
        '''Assigns restriction fragment length polymorphism (RFLP) pattern values to type 2 PRRSV ORF5 genes'''
        
        def match_rflp(dna_sequence, rflp_value):
            for r in range(len(self.crop_head_li)):
                regex = self.crop_head_li[r]
                #find restriction sites
                finditer_li = re.finditer(r"(?=("+regex+"))", dna_sequence.replace("U", "T"))
                result = [i.start(1) for i in finditer_li]
                fix_result = [i+self.cut_site_li[r]+1 for i in result]
                if fix_result in self.rflp_di[regex]:
                    rflp_value[r] = str(self.rflp_di[regex].index(fix_result)+1)
            return(rflp_value)
            
        rflp_value = ["null"] #does not contain full length ORF5
        orf5_li = [x for x in annotate_li2 if "ORF5" in x[0]]
        if orf5_li != [] and (require_full != True or "full" in orf5_li[0][0]):
            start, stop = [int(x) for x in orf5_li[0][1].split(",")]
            dna_sequence = nucleotide_sequence[start-1:stop]
            rflp_value = match_rflp(dna_sequence, ["X" for i in range(3)])
            if "X" in "".join(rflp_value) or len(nucleotide_sequence) != 603:
                hmmalign_input = self.temporary_path+"/hmmalign_input.fasta"
                hmm_build = spar_path+"/required/PRRSV2/hmm_profiles/hmm_build/ORF5.hmm"
                read_fasta_str = ">Sequence\n"+dna_sequence
                function_di = {
                    alignment.hmm.align: [hmmalign_input, hmm_build]     
                }
                hmmalign_simple = nested.temporary(hmmalign_input, read_fasta_str, function_di)[0]
                hmmalign_simple = re.sub("[a-z]+", "", hmmalign_simple)
                rflp_value = match_rflp(hmmalign_simple, ["X" for i in range(3)])
        return("-".join(rflp_value))
