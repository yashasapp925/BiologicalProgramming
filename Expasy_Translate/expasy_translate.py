# This asssignment was created by Tucker J Lancaster

def valid_DNA_sequence(DNA):
    # This function takes in a string containing possible DNA sequence and returns True if all nucleotides are valid (A,a,C,c,G,g,T,t) and False if any nucleotide is invalid
    valid_nucleotides = ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']
    for letter in DNA:
        if letter not in valid_nucleotides:
            return False
    return True

def print_DNA_sequence(DNA, mode):
    # This function takes in a string containing valid DNA sequence and a mode that specifies the output format. Prints info to a screen, maximum characters per line
    DNA = DNA.upper()
    for i in range(0,3):
        #This loop is used to create the three possible frames on the forward strand
        num_codons = (len(DNA) - i) // 3
        DNA_sequence = DNA[i:(i + num_codons * 3)]  # Create DNA sequence for current frame. Ensure it is divisible by 3
        translated_sequence = translate(DNA_sequence, mode)
        # Print out direction (5' to 3') and frame to screen
        framenum = str(i)
        print('5\' to 3\' Frame: ' + framenum)
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
        if mode == 'DNA':
            for i in range(0, len(DNA_sequence), 60):
                print(DNA_sequence[i:(i + 60)])
                print(translated_sequence[i:(i + 60)])
        else:
            print(translated_sequence)
    rev_DNA_sequence = reverse_complement(DNA)
    for i in range(0,3):
        #This loop is used to create the three possible frames on the forward strand
        num_codons = (len(rev_DNA_sequence) - i) // 3
        DNA_sequence = rev_DNA_sequence[i:(i + num_codons * 3)]# Create DNA sequence for current frame. Ensure it is divisible by 3
        translated_sequence = translate(DNA_sequence, mode)
        # Print out direction (5' to 3') and frame to screen
        framenum = str(i)
        print('3\' to 5\' Frame: ' + framenum)
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
        if mode == 'DNA':
            for i in range(0, len(DNA_sequence), 60):
                print(DNA_sequence[i:(i + 60)])
                print(translated_sequence[i:(i + 60)])
        else: 
            print(translated_sequence)
def translate(DNA_sequence, mode):
    # Create dictionaries that translate codons into amino acids of appropriate format
    codontable = {
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
    'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    # If VERBOSE, modify the start and stop codons and modify all codons to add space after
    mode = mode.upper()
    if mode == 'VERBOSE':
        codontable['ATG'] = 'Met'
        codontable['TAA'] = 'Stop'
        codontable['TAG'] = 'Stop'
        codontable['TGA'] = 'Stop'
        for key in codontable.keys():
            val = codontable[key]
            codontable[key] = val + ' '
    # If DNA, modify all codons to add space before and after
    if mode == 'DNA':
        for key in codontable.keys():
            val = codontable[key]
            codontable[key] = ' ' + val + ' '
    # Loop through DNA sequence codons
    out_seq = ''
    for i in range(0, len(DNA_sequence), 3):
        out_seq += codontable[(DNA_sequence[i:(i + 3)])]#Determine new amino acid with appropriate format

    return out_seq

def reverse_complement(DNA_sequence):
    # Returns string containing reverse complement of DNA sequence
    reverse_comp = ''
    for n in reversed(DNA_sequence):
        if n == 'A':
            reverse_comp += 'T'
        elif n == 'T':
            reverse_comp += 'A'
        elif n == 'G':
            reverse_comp += 'C'
        elif n == 'C':
            reverse_comp += 'G'
    return reverse_comp

    
# Part I: Determine if the user has entered the appropriate number of argments when they called the script (one). 
# Determine if the user entered one of three valid options for the mode. If there is an error in either of these, 
# print out informative error messages indicating which error was made, what the three valid options are, and then quit the program.
import sys

if len(sys.argv) != 2:
    print("Mode can be one of the following options:\n\tCOMPACT\n\tVERBOSE\n\tDNA")
    sys.exit()

mode = sys.argv[1]
if mode.upper() not in ["COMPACT", "VERBOSE", "DNA"]:
    print("Mode can be one of the following options:\n\tCOMPACT\n\tVERBOSE\n\tDNA")
    sys.exit()
    
# Part II: Create loop to query user for DNA sequence
dna = input("Enter DNA sequence (or Exit to quit the program): ")
if dna.upper() == 'EXIT':
    exit()
while not valid_DNA_sequence(dna):        
    print('Invalid DNA sequence. Characters must be one of A, a, C, c, G, g, T, or t')          # Part IV: Determine if user input is valid DNA sequence. If DNA sequence is not valid, print error message and allow user to enter new DNA sequence.
    dna = input("Enter DNA sequence (or Exit to quit the program): ")
    # Part III: Determine if user wants to exit program
    if dna.upper() == 'EXIT':
        exit()

# Part V: Print out 6 translated frames to the screen in appropriate format
print_DNA_sequence(dna, mode)
    
