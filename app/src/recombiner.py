import pandas as pd
import random as r
from app.utils.cli_utils import read_fa, save_fa


def main(fname, recombine_on, expdir, virus):
    raw_seqs = read_fa(f"{expdir}/{fname.replace('.fasta','')}/{fname}")
    fasta = {}
    for seq in raw_seqs:
        fasta[seq[0]] = seq[1]

    acc_ids = [i for i in fname.replace(".fasta","").split("_") if not i == virus]
    data = {}
    for acc_id in acc_ids:
        data[acc_id] = {}
        tmp = pd.read_csv(f"{virus}/{acc_id}.tsv",sep="\t",header=0,index_col=0).T
        data[acc_id]["start"] = tmp[f"{recombine_on} "]["Start "]
        data[acc_id]["len"] = tmp[f"{recombine_on} "]["Length"]

    master_genome  = acc_ids[0] # RM < TODO Randomise? If so, seq_names will need to be re-sorted.
    donor_genome   = acc_ids[1]

    donor_sequence = fasta[f">{virus}_{donor_genome}"][data[donor_genome]["start"] + 1 : data[donor_genome]["start"] + 1 + data[donor_genome]["len"]]
    recombinant_sequence = fasta[f">{virus}_{master_genome}"][0:data[master_genome]["start"]] + donor_sequence + fasta[f">{virus}_{master_genome}"][data[master_genome]["start"] + 1 + data[master_genome]["len"]:]
    recombinant_name = f">recombinant_base_{master_genome}_donor_{recombine_on}_{donor_genome}"
    save_fname = f"{expdir}/{fname.replace('.fasta','')}/{fname.replace('.fasta','_recombined.fasta')}"
    save_fa(save_fname, 
            f"{recombinant_name}\n{recombinant_sequence}")
    print(f"Saved recombinant base {master_genome} vs donor {recombine_on} {donor_genome} to {save_fname}")
    return save_fname

if __name__ == "__main__":
    fname = "OL133743_NC038878.fasta"
    recombine_on = "2A"
    expdir = "./"
    main(fname,recombine_on,expdir)
