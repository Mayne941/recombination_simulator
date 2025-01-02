import pandas as pd
import random as r
from app.utils.cli_utils import read_fa, save_fa

def get_positions(acc_ids, data_dir, virus, recombine_on, non_recombined_genes):
    data = {}
    for acc_id in acc_ids:
        data[acc_id] = {}
        tmp = pd.read_csv(f"{data_dir}/{virus}/{acc_id}.tsv",sep="\t",header=0,index_col=0).T
        try:
            data[acc_id]["start"] = tmp[f"{recombine_on}"]["Start"]
        except: 
            tmp = tmp.T
            data[acc_id]["start"] = tmp[f"{recombine_on}"]["Start"]
        try:
            data[acc_id]["len"] = tmp[f"{recombine_on}"]["Length"]
        except KeyError:
            data[acc_id]["len"] = tmp[f"{recombine_on}"]["End"] - tmp[f"{recombine_on}"]["Start"] 
    non_recombined_genes.remove(recombine_on)
    return data, non_recombined_genes


def main(fname, recombine_on, expdir, virus, data_dir, common_genes, n_genes=2):
    raw_seqs = read_fa(f"{expdir}/{fname.replace('.fasta','')}/{fname}")
    recombine_on = recombine_on.replace(" ","")
    fasta = {}
    for seq in raw_seqs:
        fasta[seq[0]] = seq[1]

    acc_ids = [i for i in fname.replace(".fasta","").split("_") if not i.replace("-","_") == virus]
    data, non_recombined_genes = get_positions(acc_ids, data_dir, virus, recombine_on, common_genes)

    # if n_genes == 0:
    #     master_genome = donor_genome = acc_ids[0]
    # else:
    master_genome, donor_genome = acc_ids[0], acc_ids[1]
    virus_with_kebab = virus.replace("_","-")

    if n_genes == 0:
        recombinant_sequence = fasta[f">{virus_with_kebab}_{donor_genome}"]
        recombinant_name = f">NOT_RECOMBINED_CONTROL_base_{master_genome}_donor_{recombine_on}_{donor_genome}"
        save_fname = f"{expdir}/{fname.replace('.fasta','')}/{fname.replace('.fasta','_recombined.fasta')}"
        save_fa(save_fname, 
                f"{recombinant_name}\n{recombinant_sequence}")

    else:
        for i in range(0, n_genes):
            print(f"Recombining genomes {donor_genome}, {master_genome}, gene {i+1}")
            if i > 0: 
                recombine_on = r.choice(non_recombined_genes)
                data, non_recombined_genes = get_positions(acc_ids, data_dir, virus, recombine_on, non_recombined_genes)
            donor_sequence = fasta[f">{virus_with_kebab}_{donor_genome}"][data[donor_genome]["start"] + 1 : data[donor_genome]["start"] + 1 + data[donor_genome]["len"]]
            recombinant_sequence = fasta[f">{virus_with_kebab}_{master_genome}"][0:data[master_genome]["start"]] + donor_sequence + fasta[f">{virus_with_kebab}_{master_genome}"][data[master_genome]["start"] + 1 + data[master_genome]["len"]:]
            recombinant_name = f">recombinant_base_{master_genome}_donor_{recombine_on}_{donor_genome}"
            save_fname = f"{expdir}/{fname.replace('.fasta','')}/{fname.replace('.fasta','_recombined.fasta')}"
            save_fa(save_fname, 
                    f"{recombinant_name}\n{recombinant_sequence}")
    print(f"Saved recombinant base {master_genome} vs donor {recombine_on} {donor_genome} to {save_fname}")
    return save_fname


