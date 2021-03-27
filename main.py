import os
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO

from Bio.Align.Applications import ClustalOmegaCommandline


# get df of each species sheet
def read_data(sheet_name, exel_file='protein_tables.xlsx'):
    df = pd.read_excel(exel_file, sheet_name)
    return df


# STEP 1 : get common bacteria set
def get_common_bacteria_set(protein_set, species):
    print("\n                       COMMON BACTERIA SET\n")
    if os.path.exists('out/common_bacteria_set.txt'):
        with open('out/common_bacteria_set.txt', 'r') as f:
            common_bacteria_set = f.read().strip().split("\n")
        for i in common_bacteria_set:
            print(i)

    else:
        common_bacteria_set = []
        for _ in species:
            specie_df = read_data(_)
            proteins_set_in_specie = set(specie_df['Protein name'])
            set_intersection = set(protein_set).intersection(proteins_set_in_specie)
            if set_intersection == set(protein_set):
                common_bacteria_set.append(_)
        with open("out/common_bacteria_set.txt", "w") as f:
            for i in common_bacteria_set:
                f.write(i + "\n")
                print(i)
    print("\n")
    return common_bacteria_set


# STEP 3 : extract gene sequence for each protein in protein_set and for each species in common_bacteria_set
def get_homologous_gene_sequences(protein_set, common_bacteria_set):
    print("                     HOMOLOGOUS GENE SEQUENCE")
    if not os.path.exists('out/genomes.txt'):
        genomes = {}
        for p in protein_set:
            new_p = p.replace(" ", "_")
            genomes_of_p = []
            species = []
            print("             ", p, "                        \n")
            for c_b in list(common_bacteria_set):
                species.append(c_b)

                df = read_data(c_b)
                df_p = df.loc[df['Protein name'] == p]
                first_index = df_p.first_valid_index()
                start, stop = df_p.loc[first_index]['Start'], df_p.loc[first_index]['Stop']

                for seq_record in SeqIO.parse(open('fasta/' + c_b + '.fasta'), "fasta"):
                    c_b_seq = seq_record.seq[start:stop + 1]

                genomes_of_p.append(SeqRecord(c_b_seq, c_b, new_p, ""))
            genomes['p'] = genomes_of_p
            SeqIO.write(genomes_of_p, "out/homologous_gene_sequences/" + new_p + ".fasta", "fasta")

            with open('out/genomes.txt', 'a') as genomes_file:
                genomes_file.write(p + " - " + str(genomes_of_p))


# STEP 4 : build phylogenetic trees
def build_phylogeny_trees():
    path = "out/homologous_gene_sequences/"
    out_path = "out/aligned_homologous_gene_sequences/"
    cmd_file_path = 'out/clustal_commonds.txt'

    for homologous_gene_sequence in os.listdir(path):
        in_file = path + homologous_gene_sequence
        out_file = out_path + homologous_gene_sequence

        if not os.path.exists(out_file):
            clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
            os.system(str(clustalomega_cline))

        msa = AlignIO.read(out_file, 'fasta')

        # Calculate the distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(msa)

        # Print the distance Matrix
        print('\nDistance Matrix\n===================')
        print(dm)

        # Construct the phylogenetic tree using UPGMA algorithm
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        # Draw the phylogenetic tree
        Phylo.draw(tree)

        # Print the phylogenetic tree in the terminal
        print('\nPhylogenetic Tree\n', homologous_gene_sequence)
        Phylo.draw_ascii(tree)


if __name__ == '__main__':
    protein_set = ["site-specific DNA-methyltransferase",
                   "LysR family transcriptional regulator",
                   "helix-turn-helix domain-containing protein",
                   "efflux transporter outer membrane subunit"]

    excel_file = pd.ExcelFile('protein_tables.xlsx')
    species = excel_file.sheet_names

    # STEP 1
    common_bacteria_set = get_common_bacteria_set(protein_set, species)

    # STEP 3
    get_homologous_gene_sequences(protein_set, common_bacteria_set)

    # STEP 4
    build_phylogeny_trees()
