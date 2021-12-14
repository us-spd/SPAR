
import os
import re
import shutil
import subprocess
import sys
import traceback


if __name__ == "__main__":
    from settings import bash_path
    from misc import alignment, fasta, string #, message
else:
    from .settings import bash_path
    from .misc import alignment, fasta, string #, message

################################################################################

spar_path = re.sub("/include$", "", os.path.dirname(os.path.realpath(__file__)))
valid_organism_li = [d for d in os.listdir(spar_path+"/required/") if os.path.isdir(spar_path+"/required/"+d) and d != "blastdb"]
if bash_path != None:
    os.environ["PATH"] = bash_path
    
################################################################################
################################################################################
    
class species_identify:
    
    def __init__(self, input_organism_li=valid_organism_li, temp_dir=""):
        '''The species_identify class is initialized'''
        
        # check required bash utilities
        required_program_li = ["blastn", "makeblastdb", "blastdbcmd", "awk"]
        missing_programs = string.readable_list([x for x in required_program_li if shutil.which(x) == None])
        if len(missing_programs) > 0:
            raise Exception("The following dependencies could not be found: "+missing_programs+". Please install and re-run.")
        
        self.temporary_path = spar_path + "/temporary"
        if os.path.isdir(temp_dir):
            self.temporary_path = temp_dir
        
        # check if species identification is comprehensive or specific
        self.input_organism_li = [x for x in input_organism_li if x in valid_organism_li] # validate
        if len(self.input_organism_li) > 1:
            if len(self.input_organism_li) == len(valid_organism_li):
                blastdb_dir = spar_path+"/required/blastdb" # check/build for comprehensive blast database
                os_call = "blastdbcmd -info -db "+blastdb_dir+"/blastdb.fasta > /dev/null 2>/dev/null"
                self.check_blastdb_li = [subprocess.call(os_call.replace(".fasta", ""), shell=True), subprocess.call(os_call, shell=True)]
            else:
                blastdb_dir = self.temporary_path # build specific, temporary blastdb
                os_call = "blastdbcmd -info -db "+blastdb_dir+"/blastdb.fasta > /dev/null 2>/dev/null"
                self.check_blastdb_li = []
            # assign blastdb
            self.blastdb = blastdb_dir+"/blastdb.fasta"
            if 0 not in self.check_blastdb_li: #and message.proceed(input_value="BLAST database must be constructed. Continue? (y/n) ", response_n="Not an optional step. Process terminated."):            
                blast_rf_li3 = []
                for organism in self.input_organism_li:
                    fasta_li = [spar_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f for f \
                                in os.listdir(spar_path+"/required/"+organism+"/hmm_profiles/msa_save") \
                                if os.path.isfile(spar_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f) and ".fa" in f]
                    len_li = []
                    for fasta_path in fasta_li:
                        blast_rf_li2 = fasta.read(fasta_path)
                        blast_rf_li2[0] = [x+"--"+organism for x in blast_rf_li2[0]]
                        blast_rf_li2[1] = [x.replace("-", "") for x in blast_rf_li2[1]]
                        len_li.append(len(blast_rf_li2[0]))
                        blast_rf_li3.append(blast_rf_li2)
                blast_rf_li2 = [[head for blast_rf_li2 in blast_rf_li3 for head in blast_rf_li2[0]], \
                                [seq for blast_rf_li2 in blast_rf_li3 for seq in blast_rf_li2[1]]]
                fasta.write(self.blastdb, blast_rf_li2)
                alignment.blast.build(self.blastdb)
                print("BLAST database successfully constructed.", file=sys.stderr)
            # validate blastdb
            self.check_blastdb_li = [subprocess.call(os_call.replace(".fasta", ""), shell=True), subprocess.call(os_call, shell=True)]
            if self.check_blastdb_li[0] == 0:
                self.blastdb = self.blastdb.replace(".fasta", "")
        
        
    def main(self, read_fasta_li2):

        fasta_di = {} #track original FASTA order
        if 0 in self.check_blastdb_li:
            try:
                # identify organism with BLAST
                if len(self.input_organism_li) > 1:
                    #print("Identifying species for "+str(len(read_fasta_li2[1]))+" sequences", file=sys.stderr)
                    for index, sequence in enumerate(read_fasta_li2[1]):
                        fasta.write(self.temporary_path+"/blast_in.fasta", [[str(index)], [sequence.replace("-", "")]])
                        # blastn call/interpretation is largely reliant upon default behaviors; may require customized options
                        os_call = "blastn -db "+self.blastdb +" -query "+self.temporary_path+"/blast_in.fasta -outfmt 6 2>/dev/null | awk -F'\t' '{print $2}' | awk -F'--' '{print $NF}' > "+self.temporary_path+"/blast_out.txt"
                        subprocess.call(os_call, shell=True)
                        with open(self.temporary_path+"/blast_out.txt") as blast_out: # use file to refer to the file object
                            read_out = blast_out.read().split("\n")[:-1] # ignore last newline character
                            if len(read_out) > 0:
                                organism_count_di = {x:read_out.count(x) for x in set(read_out)}
                                max_key_count = max(organism_count_di.items(), key=lambda x : x[1])
                                total_key_count = sum(organism_count_di.values())
                                if max_key_count[0] in self.input_organism_li and max_key_count[1]/total_key_count > 0.5:
                                    if max_key_count[0] not in fasta_di:
                                        fasta_di[max_key_count[0]] = []
                                    fasta_di[max_key_count[0]].append([index, sequence])
                elif len(self.input_organism_li) == 1:
                    fasta_di[self.input_organism_li[0]] = [[index, sequence] for index, sequence in enumerate(read_fasta_li2[1])] 
            except:
                traceback.print_exc(file=sys.stderr)
                print(read_fasta_li2[0], file=sys.stderr) #?
        return(fasta_di)
