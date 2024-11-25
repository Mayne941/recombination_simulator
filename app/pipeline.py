from app.utils.cli_utils import read_fa, save_fa, shell
from app.src.recombiner import main as recombine
from app.utils.api_objects import castanet_req_body

import os
import random as r
import numpy as np
import pandas as pd
import datetime as dt
import requests as req

class RecombinationPipeline:
    def __init__(self, virus, n, endpoint, dwgsim_path):
        self.virus = virus
        self.n = n
        self.seqs = read_fa(f"{self.virus}/{self.virus}.fasta")
        self.expname = f'experiment_{self.virus}_{dt.datetime.now().strftime("%Y_%M_%d_%H_%M_%S")}'
        if not os.path.exists(self.expname):
            os.mkdir(self.expname)
        self.virus_genes = {
            "rhinovirus-a": ["5'UTR ", "1A", "1B", "1C", "1D", "2A", "2B", "2C", "3D"]
        }
        self.dwg_settings = {
            "n_reads": 10000,
            "error_rate": 0.025
        }
        self.castanet_endpoint = endpoint
        self.dwgsim_path = dwgsim_path
        self.results = {}
        self.fails = 0

    def generate_seqs(self):
        '''Get 2 random seqs from database'''
        avail_seqs = np.arange(0,10)
        r.shuffle(avail_seqs)
        replicate_seqs = [self.seqs[avail_seqs[0]], self.seqs[avail_seqs[1]]]
        seq_names = [i[0].replace(">","") for i in replicate_seqs]
        exp_dir = f"{self.expname}/{'_'.join(seq_names)}"
        os.mkdir(exp_dir)
        save_fa(f"{self.expname}/{'_'.join(seq_names)}/{'_'.join(seq_names)}.fasta",
                ''.join([f"{i[0]}\n{i[1]}\n" for i in replicate_seqs]))
        return seq_names, exp_dir

    def generate_syn_reads(self, syn_read_dir, recombined_fasta_fname):
        '''Call DWGsim to make synthetic reads, kill unneeded files'''
        os.mkdir(f"{syn_read_dir}/{syn_read_dir.split('/')[1]}_syn_reads")
        shell(f"{self.dwgsim_path} -N {self.dwg_settings['n_reads']} -1 150 -2 150 -y {self.dwg_settings['error_rate']} {recombined_fasta_fname} {syn_read_dir}/{syn_read_dir.split('/')[1]}_syn_reads/{syn_read_dir.split('/')[1]}")
        for file in os.listdir(f"{syn_read_dir}/{syn_read_dir.split('/')[1]}_syn_reads/"):
            if not ".bwa.read" in file:
                os.remove(f"{syn_read_dir}/{syn_read_dir.split('/')[1]}_syn_reads/{file}")
        print(f"Completed generating synthetic reads for {syn_read_dir}")

    def hit_castanet(self, body):
        print(f"Calling CASTANET")
        resp = req.post(self.castanet_endpoint, json=body)
        if not resp.ok:
            raise Exception(f"Castanet returned a non-200, details:\n{resp.text}")
        else:
            print(f"CASTANET call complete")

    def check_castanet_output(self, body, seq_names, recombined_fasta_fname):
        '''Interrogate Castanet output and compile statistics'''
        self.results["_".join(seq_names)] = {}

        '''Check consensus'''
        probename = self.virus.split("-")[0]
        os.listdir(f'{body["SaveDir"]}/{body["ExpName"]}')
    
        if not probename in os.listdir(f'{body["SaveDir"]}/{body["ExpName"]}/consensus_data/'):
            '''No consensus - error state'''
            self.results["_".join(seq_names)]["consensus_identity"] = 0
        else:
            '''Get consensus stats'''
            cons = f'{body["SaveDir"]}/{body["ExpName"]}/consensus_data/{probename}/{probename}_remapped_consensus_sequence.fasta'
            blast_fname = f'{body["SaveDir"]}/blast_results.tab' 
            outfmat = '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'
            shell(f'makeblastdb -in {recombined_fasta_fname} -dbtype nucl')# RM < TODO test blast active in env
            shell(f'blastn -query {cons} -db {recombined_fasta_fname} -out {blast_fname} -evalue 1E-6 -outfmt {outfmat} -num_alignments 100000 -num_threads 1')
            blast_df = pd.read_csv(blast_fname, sep="\t", names=["qseqid", "sseqid", "pident", "qcovs", "qlen", "slen", "evalue", "bitscore"])
            self.results["_".join(seq_names)]["cons_coverage"] = blast_df["qcovs"][0]
            self.results["_".join(seq_names)]["cons_identity"] = blast_df["pident"][0]

        '''Check read n'''
        depth_df = pd.read_csv(f'{body["SaveDir"]}/{body["ExpName"]}/{body["ExpName"]}_depth.csv')
        depth_df = depth_df[depth_df["probetype"] == probename]
        self.results["_".join(seq_names)]["reads_on_target_ratio"] = depth_df["n_reads_all"].item() / self.dwg_settings["n_reads"]

    def save_results(self, out_fname):
        df = pd.DataFrame(self.results).T
        df.to_csv(out_fname)
        df.describe().round(2).to_csv(out_fname.replace("ALL_STATS","SUMMARY_STATS"))
        print(f"Finished. Data saved to {out_fname}")

    def main(self):
        for replicate in range(self.n):
            try:
                seq_names, exp_dir = self.generate_seqs()
                '''Pick a random gene to recombine seqs on'''
                gene_to_recombine = r.choice(self.virus_genes[self.virus])
                print(f"Round {replicate+1}/{self.n}. Processing genomes: {seq_names}, on gene {gene_to_recombine}")

                '''Call recombiner script to make recombinant genome'''
                recombined_fasta_fname = recombine(f"{'_'.join(seq_names)}.fasta", gene_to_recombine, self.expname, self.virus)

                '''Generate synthetic reads on newly made recombinant'''
                self.generate_syn_reads(exp_dir, recombined_fasta_fname)

                '''Populate Castanet command script and hit it'''
                body = castanet_req_body(exp_dir, self.virus)            
                self.hit_castanet(body)
                
                '''Check output and compile statistics'''
                self.check_castanet_output(body, seq_names, recombined_fasta_fname)
                print(f"Finished processing genomes {seq_names}")
            except: 
                print("Couldn't process this sample")
                self.fails += 1
                continue

        self.save_results(f'{self.expname}/{self.expname}_ALL_STATS.csv')
        print(f"...with {self.fails} errors")

if __name__ == "__main__":
    '''
    Assumes directory containing 
      {virus-name}.fasta, with all ref seqs
      for each virus in fasta, an {accession}.tsv with annotations
    AND Castanet running at endpoint listed in 
    '''
    virus = "rhinovirus-a" # ["rhinovirus-a", "enterovirus-a", "hiv-1", "sars-cov2", "cmv"] 
    n = 100
    endpoint = "http://127.0.0.1:8001/end_to_end/"
    dwgsim_path = "../playground/DWGSIM/dwgsim"

    rp = RecombinationPipeline(virus, n, endpoint, dwgsim_path)
    rp.main()