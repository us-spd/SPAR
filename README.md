# Swine Pathogen Analysis Resource (SPAR)
`SPAR` is a python module useful for the analysis of swine pathogen genomes. 

`SPAR` offers high accuracy sequence annotating based on custom built multiple sequence alignment (MSA) templates. Templates can be representative of virtually any kind of continuous genomic feature (ie coding region, gene, peptide) that can be translated into an amino acid sequence. Annotating uses regular expressions derived from template MSAs to rapidly identify the location of genomic features. Conditional further treatment involves alignment to a hidden Markov model (HMM) profile. The HMM profile alignment is used to accurately determine terminal endpoints, identify translational frameshifts, and locate indels.

Comparable annotation pipelines consist of using BLAST to identify highly similar parent sequences. Following pairwise alignment, annotations from the parent sequence/s are transfered to the query. By comparison, these methods are generally faster and sufficiently accurate when a large reference database is available. In the absence of  similar comparison sequences, HMM profile alignment is a good alternative to pairwise alignment that pools available information to maximize alignment accuracy.

`SPAR` also provides a classification method unique to Porcine reproductive and respiratory syndrome virus-2 (PRRSV-2) that assigns restriction fragment length polymorphism (RFLP) patterns. RFLP assignment is flexible; HMM profile alignment is used to compensate for indels in the query sequence that would impede pattern determination. No other comparable public resource is known.


## Setup

OS X & Linux:

Install [Python3.7][python3]<br />
Unzip `required.zip`<br />
Install necessary python packages:<br />
```sh
pip install -r requirements.txt
```
Install [MAFFT][mafft] (version >= 7.310) dependency<br />
Install [BLAST][blast] (version >= 2.10.0+) dependency<br />
Install [HMMER3][hmmer3] (release >= 3.1b2) dependency<br />

Rename include/settings_template.py to include/settings.py and make the following changes:
1) Open bash command line terminal
2) type `echo $PATH` into terminal and copy output
3) replace "None" in line `bash_path = None` with copied string (quoted)


## Documentation

Usage information is available via command line:

```
$ python3 run.py -h
usage: run.py [-h] [-v] {annotate,rflp} ...

positional arguments:
  {annotate,rflp}  Type of input (rflp or annotate). See specific help for more options based on input type chosen e.g. python run.py rflp --help
    annotate       Annotate genomic features. Currently supports CSFV, FMDV, PDCoV, PEDV, PRRSV-1, PRRSV-2, and SVA
    rflp           Assign RFLP value to PRRSV-2.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show program's version number and exit
```

Detailed information for each command can be requested:

<details>
  <summary>Genome annotation</summary>
  <p>
<pre>
$ python3 run.py annotate -h
usage: run.py annotate [-h] [--limit LIMIT] [--organism ORGANISM] [--output OUTPUT] input

positional arguments:
  input                Input FASTA file name

optional arguments:
  -h, --help           show this help message and exit
  --limit LIMIT        Annotates only genomic features provided in comma delimited list
  --organism ORGANISM  Possible organism identity of sequences in input FASTA file provided in comma delimited list
  --output OUTPUT      Output gff3 file name
</pre>

The <i>--organism</i> argument is used to specify the possible species identity for any sequence within the input FASTA file. Identification is based solely on BLAST similarity score. Additional pathogens may be added by creating a pathway under the requirements directory with the following structure: <i>organism_abbreviation/hmm_profiles/msa_save</i>. Add nucleotide MSA templates in FASTA file format to the msa_save directory. Templates can be representative of virtually any kind of continuous genomic feature that can be translated into an amino acid sequence. If a genomic feature contains a translational frameshift, adjust <code>translational_frameshift_di</code> variable as described in include/settings_template.py. Once the necessary required templates have been provided, running the script for the first time will initiate a build process.
  </p>
</details>


<details>
  <summary>PRRSV-2 RFLP classification</summary>
  <p>
<pre>
$ python3 run.py rflp -h
usage: run.py rflp [-h] [--full FULL] [--output OUTPUT] input

positional arguments:
  input            Input FASTA file name

optional arguments:
  -h, --help       show this help message and exit
  --full FULL      Only assign RFLP pattern to input sequences that are complete
  --output OUTPUT  Output FASTA file name
</pre>

Output file format will consist of the input FASTA file contents with RFLP pattern values appended at the end of each header, separated by a forward slash ("/") delimiter. If a sequence is not identified as being PRRSV-2, "na" will be appended. If a sequence does not contain a complete ORF5 gene, "null" will be appended. RFLP determination is only performed on complete genes be default. Use caution when enabling RFLP assignment to partial ORF5 sequences. Restriction sites cannot be inferred from missing residues.
  </p>
</details>


## Contact information

Please contact the developer with any comments, concerns, or questions: blake.inderski@usda.gov<br>
Alternatively, post an issue in this GitHub repository for assistance.<br>

[python3]: https://docs.python-guide.org/starting/install3/linux/
[mafft]: https://mafft.cbrc.jp/alignment/software/
[blast]: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
[hmmer3]: http://hmmer.org/download.html