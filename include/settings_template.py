'''Template for settings file. Original should not be commited to GitHub.'''

# when running run.py in python console, set environment path to match bash
# run the following code in unix shell: echo $PATH

bash_path = None # copy/paste terminal output here

minimal_stderr = False # ignore affirmative sequence processing and error traceback


# describe values necessary for translational frameshift treatment (maximum one translational frameshift site per organism)
translational_frameshift_di = {
  "PEDV": {
    "fs_regex": "AGCACTGATATGGCTTATTTAAACGAGTACGGGGCTCTA",  # regular expression exactly matching surrounding frameshift location (based on consensus of msa_save)
    "fs_index": 9,                                          # index of frameshift location relative to fs_regex match
    "fs_type": -1,                                          # direction of ribosomal slippage
    "fs_product": ["ORF1ab", "nsp12"]                       # affected translation product
  },
  "PDCoV": {
    "fs_regex": "AATTCGGCTTATTTAAACG.GTAACGGGTTCTAGTGA",
    "fs_index": 16,
    "fs_type": -1,
    "fs_product": ["ORF1ab"]
  },
  "PRRSV1": {
    "fs_regex": "TTTAAACTG.TAGCCGCCAGCGGCTTGACCCGCTGTGG",
    "fs_index": 10,
    "fs_type": -1,
    "fs_product": ["ORF1ab", "nsp9"]
  },
  "PRRSV2": {
    "fs_regex": "TTTAAACTG.TAGCCGCCAGCGGCTTGACCCGCTGTGG",
    "fs_index": 10,
    "fs_type": -1,
    "fs_product": ["ORF1ab", "nsp9"]
  }
}