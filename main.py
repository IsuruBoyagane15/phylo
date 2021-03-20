import os
import pandas as pd
import Bio.SeqIO as SeqIO
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from skbio import DNA
import skbio.alignment as a
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO


# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', 1)


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
            set_intersection = set(protein_set) & proteins_set_in_specie
            if set_intersection == protein_set:
                common_bacteria_set.append(_)
        print("common_bacteria_set\n")
        with open("out/common_bacteria_set.txt", "w") as f:
            for i in common_bacteria_set:
                f.write(i + "\n")
                print(i)
    print("\n")
    return common_bacteria_set


# STEP 3 : extract gene sequence for each protein in protein_set and for each species in common_bacteria_set
def get_homologous_gene_sequences(protein_set, common_bacteria_set):
    print("                     HOMOLOGOUS GENE SEQUENCE")
    if not os.path.exists('out/genomes.json'):
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
            SeqIO.write(genomes_of_p, "out/homologous_gene_sequences/" + new_p +".fasta", "fasta")

            # writing .phy files

            # for i in range(len(species)):
            #     species[i] = species[i].replace(".", "_")

            # genomes[p] = {'species': species,
            #               'genomes': genomes_of_p}

            # hls_df = pd.DataFrame.from_dict(genomes[p])
            # print(hls_df)

            # with open("out/homologous_gene_sequences/" + p + ".phy", 'w') as f:
            #     f.write("{:10} {}\n".format(hls_df.shape[0], hls_df.genomes.str.len()[0]))
            #     for row in hls_df.iterrows():
            #         f.write("{:10} {}\n".format(*row[1].tolist()))
            # print("\n")

            # with open("out/homologous_gene_sequences_fsa/" + p + ".fsa", 'w') as f:
            #     for i in range(len(species)):
            #         f.write("> "+ species[i] + "\n")
            #         f.write(genomes_of_p[i] + "\n")

            # short_seq_iterator = (record for record in genomes_of_p if len(record) < 300)
            #
            # SeqIO.write(short_seq_iterator, "short_seqs.fasta", "fasta")
            # break
        with open('out/genomes.json', 'w', encoding='utf8') as json_file:
            json.dump(genomes, json_file, ensure_ascii=False)


def build_phylogeny_trees(file: str):

        aln = AlignIO.read(file, 'fasta')

        # Calculate the distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)

        # Print the distance Matrix
        print('\nDistance Matrix\n===================')
        print(dm)

        # Construct the phylogenetic tree using UPGMA algorithm
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        # Draw the phylogenetic tree
        Phylo.draw(tree)

        # Print the phylogenetic tree in the terminal
        print('\nPhylogenetic Tree\n===================')
        Phylo.draw_ascii(tree)


if __name__ == '__main__':
    protein_set = ["site-specific DNA-methyltransferase",
                   "LysR family transcriptional regulator",
                   "helix-turn-helix domain-containing protein",
                   "efflux transporter outer membrane subunit"]

    # protein_set = {"FUSC family protein", "MFS transporter", "YraN family protein"}

    excel_file = pd.ExcelFile('protein_tables.xlsx')
    species = excel_file.sheet_names

    # STEP 1
    common_bacteria_set = get_common_bacteria_set(protein_set, species)

    # STEP 3
    get_homologous_gene_sequences(protein_set, common_bacteria_set)

    # STEP 4
    for file in os.listdir('out/homologous_gene_sequences'):
        build_phylogeny_trees('out/homologous_gene_sequences/' +file.replace(" ", "_"))
