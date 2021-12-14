#!/usr/bin/env python3

import argparse
import os
from progress.bar import Bar
import subprocess
import sys
import traceback


from version import __version__
from include.misc import fasta, gff3, string
from include.prrsv_orf5_rflp import prrsv_orf5_rflp
from include.sequence_annotate import sequence_annotate
from include.settings import bash_path, minimal_stderr
from include.species_identify import species_identify  

################################################################################

spar_path = os.path.dirname(os.path.realpath(__file__))
valid_organism_li = [d for d in os.listdir(spar_path+"/required/") if os.path.isdir(spar_path+"/required/"+d) and d != "blastdb"]
if bash_path != None:
    os.environ["PATH"] = bash_path

################################################################################
################################################################################
    
def main():
    '''Main method : parse input arguments and run appropriate operation'''
    
    def str2bool(v):
        out = string.convert_bool(v)
        if out == None:
            raise argparse.ArgumentTypeError('Boolean value expected.')
        return(out)
        
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="process", required=True, help="Type of input (rflp or annotate). See specific help for more options based on input type chosen e.g. python run.py rflp --help")
    
    # version
    parser.add_argument("-v", "--version", action="version", version="SPAR "+__version__)
    
    # annotate
    sp_annotate = subparsers.add_parser("annotate", help="Annotate genomic features. Currently supports "+string.readable_list(sorted(valid_organism_li)))
    sp_annotate.add_argument("input", type=str, help="Input FASTA file name")
    sp_annotate.add_argument("--limit", type=str, help="Annotates only genomic features provided in comma delimited list")
    sp_annotate.add_argument("--organism", type=str, help="Possible organism identity of sequences in input FASTA file provided in comma delimited list")
    sp_annotate.add_argument("--output", type=str, help="Output gff3 file name")
    
    # rflp
    sp_rflp = subparsers.add_parser("rflp", help="Assign RFLP value to PRRSV-2.") #!
    sp_rflp.add_argument("input", type=str, help="Input FASTA file name")
    sp_rflp.add_argument("--full", type=str2bool, default=True, help="Only assign RFLP pattern to input sequences that are complete")
    sp_rflp.add_argument("--output", type=str, help="Output FASTA file name")
    
    args = parser.parse_args()

#   logging.info("Input Arguments : %s", args) #?
    
    read_fasta_li2 = fasta.read(args.input)
    # build random, non-existent directory name
    temp_dir = string.unique_directory(spar_path+"/temporary/")
    subprocess.call(["mkdir", "-p", temp_dir])
    
    # variable assignment based on invoked subparser option
    if args.process == "annotate":
        input_organism_li, choose_match_li = [], []
        if args.organism != None:
            input_organism_li = args.organism.split(",")
        if args.limit != None:
            choose_match_li = args.limit.split(",")
    if args.process == "rflp":
        input_organism_li, choose_match_li = ["PRRSV-1", "PRRSV-2"], ["ORF5"]
        # remove gaps from MSA prior to determining RFLP
        copy_sequence_li = read_fasta_li2[1][:]
        read_fasta_li2[1] = [x.replace("-", "") for x in read_fasta_li2[1]]
    
    # Run

    # processes input fasta based on provided list of possible organism/s and feature/s
    annotate_li4 = [[] for x in read_fasta_li2[0]]
    input_organism_li = [x for x in input_organism_li if x in valid_organism_li]
    if len(input_organism_li) == 0:
        input_organism_li = valid_organism_li[:]

    if minimal_stderr: # ignore affirmative sequence processing and error traceback
        original_stderr = sys.stderr
        sys.stderr = open(os.devnull, "w")
    # check selected organism dependencies
    prepare = sequence_annotate.prepare()
    for organism in input_organism_li:
        if not prepare.main(organism):
            print("Removing "+organism+" from searchable organisms.", file=sys.stderr)
            input_organism_li.remove(organism)
    fasta_di, sum_len = {}, len(read_fasta_li2[1])
    if len(input_organism_li) != 1:
        print("Identifying species for "+str(len(read_fasta_li2[1]))+" sequences", file=sys.stderr)
        identify = species_identify(input_organism_li, temp_dir)
        fasta_di = identify.main(read_fasta_li2)
        sum_len = sum([len(fasta_di[x]) for x in fasta_di])
        blast_complete_str = "Species identification complete"
        if len(read_fasta_li2[1]) > sum_len:
            blast_complete_str += " ("+str(len(read_fasta_li2[1])-sum_len)+" sequences could not be identified)"
        print(blast_complete_str, file=sys.stderr) 
    else:
        fasta_di[input_organism_li[0]] = [[index, sequence] for index, sequence in enumerate(read_fasta_li2[1])]
    print("Annotating "+str(sum_len)+" sequences...", file=sys.stderr)
    bar = Bar("Progress", fill='#', suffix='%(percent)d%%', max=sum_len) # track progress with bar visual
    try:
        # annotate by organism
        for organism in fasta_di:
            process = sequence_annotate.process(organism, temp_dir)
            choose_match_li = [x for x in choose_match_li if x in process.total_cds]
            if choose_match_li == []:
                choose_match_li = process.total_cds        
            for di_li in fasta_di[organism]: # use original FASTA order
                total_annotate_li2 = process.main(di_li[1], choose_match_li)
                if len(total_annotate_li2) > 0:
                    annotate_li4[di_li[0]] = [organism, total_annotate_li2]
                bar.next()
    except:
        traceback.print_exc(file=sys.stderr)
        print(read_fasta_li2[0], file=sys.stderr) #?
    bar.finish()
    print("Annotation process complete", file=sys.stderr)
    if minimal_stderr:
        sys.stderr.close()
        sys.stderr = original_stderr

    # annotate argument
    if args.process == "annotate":
        output_gff3 = gff3.write(read_fasta_li2, annotate_li4, spar_path+"/required/")
        if args.output != None:
            write_file = open(args.output, "w")
            write_file.write(output_gff3)
            write_file.close()
        else:
            print(output_gff3)
            
    # rflp argument
    if args.process == "rflp":
        rflp = prrsv_orf5_rflp(temp_dir)
        print("Assigning RFLP patterns to "+str(len(annotate_li4))+" sequences...", file=sys.stderr)
        bar = Bar("Progress", fill='#', suffix='%(percent)d%%', max=len(annotate_li4)) # track progress with bar visual
        for fasta_index, annotate_li3 in enumerate(annotate_li4):
            if len(annotate_li3) > 0 and annotate_li3[0] == "PRRSV-2":
                read_fasta_li2[0][fasta_index] += "/"+rflp.main(read_fasta_li2[1][fasta_index], annotate_li3[1], args.full)
            else:
                read_fasta_li2[0][fasta_index] += "/na" # invalid organism
            bar.next()
        bar.finish()
        print("RFLP pattern assignment complete", file=sys.stderr)
        # return original sequences
        read_fasta_li2[1] = copy_sequence_li
        if args.output != None and os.path.isfile(args.output):
            fasta.write(args.output, read_fasta_li2)
        else:
            for fasta_index, head in enumerate(read_fasta_li2[0]):
                print(">"+head)
                print(read_fasta_li2[1][fasta_index])
                
    # remove temporary directory
    subprocess.call(["rm", "-rf", temp_dir])


if __name__ == "__main__":
    main()
    
    

