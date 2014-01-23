
peptides = []
with open("pepsin_fragments.txt", 'r') as f:
    for line in f:
        peptides.append(line.rstrip())


sequence = "MQGQKVFTNTWAVRIPGGPAVANSVARKHGFLNLGQIFGDYYHFWHRGVTKRSLSPHRPRHSRL" \
    "QREPQVQWLEQQVAKRRTKR"

his_positions = [[28],[42],[45],[56],[58],[42,45],[56,58]]

for his_pos in his_positions:
    print his_pos
    for peptide in peptides:
        pos = sequence.index(peptide)
        if pos <= his_pos[0] and pos+len(peptide) >= his_pos[0]:
            print peptide
