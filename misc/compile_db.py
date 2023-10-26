from collections import Counter
import Bio.SeqIO as bio

records = list(bio.parse("/home/gnampfelix/Documents/term5/thesis/data/sequences.fasta", "fasta"))
with open("/home/gnampfelix/Documents/term5/thesis/data/filtered_genomes.txt") as f:
    assembly_names = f.read().splitlines();

for assembly in assembly_names:
    record = [r for r in records if assembly in r.description]
    if len(record) != 1:
        print(f"skipping {assembly}")
    
