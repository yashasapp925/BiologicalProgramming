"""
Assignment 5: Automatic Primer Identification

In this assignment, you will practice using python's subprocess library, the pyfasta library, and the primer3
command-line-interface to automatically identify primers. This template provides a potential route to a completed
script, but the implementation details are entirely up to you. The checker script, as you will see, does little more
than confirm that your script prints the correct primers to the console for a handful of argument combinations.
Author: Tucker J Lancaster
"""

import argparse
import subprocess
from pyfasta import Fasta

"""
If you have not already, install pyfasta and primer3 into the current environment using the following commands:
conda install -c bioconda primer3
conda install -c bioconda pyfasta
"""

"""
Create your parser and arguments. Your script should expect three positional arguments: First, the name of a
fasta file, second, the chromosome you want to amplify, and third, the position you want amplified. Enforce that
the first two arguments are of type str, and the third is of type int.
"""

parser = argparse.ArgumentParser()

parser.add_argument("fasta_file", type=str)
parser.add_argument("chromosome", type=str)
parser.add_argument("position", type=int)

args = parser.parse_args()

fasta_file = args.fasta_file
chromosome = args.chromosome
position = args.position

"""
Open the fasta file and read in the entirety of the sequence. You can use whatever data structure you prefer.
"""

fasta = Fasta(fasta_file)
entire_sequence = fasta[chromosome]

"""
Identify the sequence 500 bp upstream and downstream of the requested position and store it as a string.
"""

upstream = max(0, position - 500)
downstream = min(len(entire_sequence), position + 500)
sequence_around_position = entire_sequence[upstream:downstream]

"""
Add braces to the DNA sequence 100 bp upstream and downstream of the requested position. Alternatively, you can use 
an argument in the input file to specify this
"""

# upstream = max(0, position - 100)
# downstream = min(len(entire_sequence), position + 100)


"""
Create an input file containing the DNA sequence and specify a product length of 600 - 800 base pairs
"""

input_file = f"""SEQUENCE_ID=sequence
SEQUENCE_TEMPLATE={sequence_around_position}
PRIMER_PRODUCT_SIZE_RANGE=600-800
=
"""
with open('input.txt', 'w') as file:
    file.write(input_file)
"""   
Execute the primer3_core command on the input file you created.
"""

subprocess.run(['primer3_core', '--output', 'output.txt', 'input.txt'])

""" 
Read in and parse the output to identify the sequence of best two primers. Print these out to the user
"""

with open('output.txt', 'r') as f:
    for line in f:
        if 'PRIMER_LEFT_0_SEQUENCE' in line:
            print(line[23:len(line) - 1])
        elif 'PRIMER_RIGHT_0_SEQUENCE' in line:
            print(line[24:len(line) - 1])