from utils import read_fa, save_fa, shell
from recombiner import main as recombine
from api_objects import castanet_req_body

import os
import random as r
import numpy as np
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
            "n_reads": 100,
            "error_rate": 0.025
        }
        self.castanet_endpoint = endpoint
        self.dwgsim_path = dwgsim_path
        self.results = {}

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
        shell(f".{self.dwgsim_path} -N {self.dwg_settings['n_reads']} -1 150 -2 150 -y {self.dwg_settings['error_rate']} {recombined_fasta_fname} {syn_read_dir}/{syn_read_dir.split('/')[1]}_syn_reads/{syn_read_dir.split('/')[1]}")
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

    def check_castanet_output(self, body, seq_names):
        '''Interrogate Castanet output and compile statistics'''
        dir_root = body["ExpDir"]
        self.results["_".join(seq_names)] = {}
        '''Check consensus created'''
        breakpoint()

    def main(self):
        for replicate in range(self.n):
            seq_names, exp_dir = self.generate_seqs()
            '''Pick a random gene to recombine seqs on'''
            gene_to_recombine = r.choice(self.virus_genes[self.virus])
            print(f"Round {replicate}/{self.n}. Processing genomes: {seq_names}, on gene {gene_to_recombine}")

            '''Call recombiner script to make recombinant genome'''
            recombined_fasta_fname = recombine(f"{'_'.join(seq_names)}.fasta", gene_to_recombine, self.expname, self.virus)

            '''Generate synthetic reads on newly made recombinant'''
            self.generate_syn_reads(exp_dir, recombined_fasta_fname)

            '''Populate Castanet command script and hit it'''
            body = castanet_req_body(exp_dir, self.virus)            
            self.hit_castanet(body)
            
            '''Check output and compile statistics'''
            self.check_castanet_output(body, seq_names)
            print(f"Finished processing genomes {seq_names}")

if __name__ == "__main__":
    '''
    Assumes directory containing 
      {virus-name}.fasta, with all ref seqs
      for each virus in fasta, an {accession}.tsv with annotations
    AND Castanet running at endpoint listed in 
    '''
    virus = "rhinovirus-a" # ["rhinovirus-a", "enterovirus-a", "hiv-1", "sars-cov2", "cmv"] 
    n = 10
    endpoint = "http://127.0.0.1:8001/end_to_end/"
    dwgsim_path = "../DWGSIM/dwgsim"

    rp = RecombinationPipeline(virus, n, endpoint, dwgsim_path)
    rp.main()