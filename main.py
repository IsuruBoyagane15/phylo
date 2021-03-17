import pandas as pd


# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', 1)

# get df of each species sheet
def read_data(sheet_name, exel_file='protein_tables.xlsx'):
    df = pd.read_excel(exel_file, sheet_name)
    return df


# STEP 1 - get common bacteria set
def get_common_bacteria_set(protein_set, species):
    common_bacteria_set = set()
    for _ in species:
        specie_df = read_data(_)
        proteins_set_in_specie = set(specie_df['Protein name'])
        # if "YraN family protein" in proteins_set_in_specie:
        #     print("bio")
        set_intersection = protein_set & proteins_set_in_specie
        if set_intersection == protein_set:
            common_bacteria_set.add(_)
    return common_bacteria_set


if __name__ == '__main__':
    protein_set = {"site-specific DNA-methyltransferase",
                   "LysR family transcriptional regulator",
                   "helix-turn-helix domain-containing protein",
                   "efflux transporter outer membrane subunit"}

    # protein_set = {"FUSC family protein", "MFS transporter", "YraN family protein"}

    excel_file = pd.ExcelFile('protein_tables.xlsx')
    species = excel_file.sheet_names

    common_bacteria_set = get_common_bacteria_set(protein_set, species)
    print("common_bacteria_set\n")
    print(common_bacteria_set)
