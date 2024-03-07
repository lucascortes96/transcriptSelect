import argparse

# Create the parser
parser = argparse.ArgumentParser(description='Process some integers.')

# Add the arguments
parser.add_argument('--number', type=int, help='The number in the filename')

# Parse the arguments
args = parser.parse_args()

# Open the file with your data
filename = f'start_diff_{args.number}bp.txt'
with open(filename, 'r') as f:
    data = f.readlines()

# Extract the genes from the first column and count the unique genes
unique_genes = set(line.split('\t')[0] for line in data)
count_genes = len(unique_genes)
# Extract the transcripts from the second column and count the unique transcripts 
unique_transcripts = set(line.split('\t')[1] for line in data)
count_transcripts = len(unique_transcripts)

# Write the count to a new file
with open('generalStats.txt', 'w') as f:
    f.write("NUMBER OF UNIQUE GENES IN FILE\n")
    f.write(str(count_genes)+"\n")
    f.write("NUMBER OF UNIQUE TRANSCRIPTS IN FILE\n")
    f.write(str(count_transcripts))

