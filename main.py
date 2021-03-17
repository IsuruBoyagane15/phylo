import os
import pandas as pd
import Bio.SeqIO as SeqIO

# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', 1)


# get df of each species sheet
def read_data(sheet_name, exel_file='protein_tables.xlsx'):
    df = pd.read_excel(exel_file, sheet_name)
    return df


# STEP 1 : get common bacteria set
def get_common_bacteria_set(protein_set, species):
    print("\n=============================        COMMON BACTERIA SET         =============================\n")
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
def extract_gene_sequences(protein_set, common_bacteria_set):
    genomes = {}
    for p in list(protein_set):
        genomes_of_p = {}
        print("=============================        ", p, "         =============================\n")
        for c_b in list(common_bacteria_set):

            for seq_record in SeqIO.parse(open('out/fasta/' + c_b + '.fasta'), "fasta"):
                c_b_seq = seq_record.seq

            df = read_data(c_b)
            df_p = df.loc[df['Protein name'] == p]
            first_index = df_p.first_valid_index()
            start, stop = df_p.loc[first_index]['Start'], df_p.loc[first_index]['Stop']
            genomes_of_p[c_b] = c_b_seq[start:stop + 1]
            print(genomes_of_p)
        print("\n")
        break
        # print(genomes_of_p)


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
    # extract_gene_sequences(protein_set, common_bacteria_set)
