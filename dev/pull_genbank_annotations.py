from Bio import SeqIO, AlignIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sp
import pandas as pd
import os

def shell(args, calling_fn="Misc shell function", ret_output=False):
    '''Call Bash shell with input string as argument'''
    whitelist = ["mafft"]
    _ = sp.Popen(args, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = _.communicate()
    if not any(x.lower() in whitelist for x in whitelist):
        sp_error_handler(
            out, err, f"ERROR: Something went wrong with CLI call: {calling_fn}, dumping STDOUT/STDERR to shell.")
    if ret_output:
        return out + err

def sp_error_handler(out, err, name):
    '''Kill program if error found; used in combo with POpen commands.'''
    if err != b"":
        raise SystemExit(f"***\n{name}:\nOut: {out}\nErr: {err}\n***")


def DownloadGenBankFile(GenomeSeqFile, SeqIDLists, email):
    '''Hit GenBank to get data. Needs user to provide email authentication'''
    if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
        os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))

    Entrez.email = email
    try:
        handle = Entrez.efetch(db="nucleotide", id=", ".join([i for i in SeqIDLists]), rettype="gb", retmode="text")
    except Exception as e:
        raise "Failed to pull Genbank Data from NCBI Entrez with exception"
    
    print("Writing contents of genbank object to file")
    with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
        GenomeSeqFile_handle.write(handle.read())

    handle.close()

def get_genbank(GenomeSeqFile, genomes):
    '''3/10: Check if GenBank files exist; dl if needed. Index & transform.'''
    if not os.path.isfile(GenomeSeqFile):
        '''If no gb file provided, attempt to pull seqs'''
        DownloadGenBankFile(GenomeSeqFile, # RM < TODO We should give users option to provide their own gb file if not on genbank
                            genomes, "someone@provider.com")

    '''Parse gb file'''
    GenBankDict = SeqIO.index(GenomeSeqFile, "genbank")
    CleanGenbankDict = {}
    for i in GenBankDict.items():
        if "u" in str(i[1].seq).lower():
            i[1].seq = i[1].seq.back_transcribe()
        CleanGenbankDict[i[0].split(".")[0]] = i[1]

    return CleanGenbankDict


def get_gene_annotations(genbankies):
    return [(i.qualifiers["gene"][0], int(i.location.start), int(i.location.end), (int(i.location.end)-int(i.location.start))) for i in genbankies.features if i.type == "gene"]

if __name__ == "__main__":
    DATA = {
        # "hiv_1": [
        #     "PP313274",
        #     "OQ551963",
        #     "OR736078",
        #     "OR661087",
        #     "OR483970",
        #     "LC735411",
        #     "PQ046777",
        #     "OR876493",
        #     "PP074169",
        #     "NC_001802"   
        # ],
        # "hcmv": [
        #     "NC_006273",
        #     "PP803491",
        #     "OU342900",
        #     "OQ466311",
        #     "OK000909",
        #     "MT649468",
        #     "MW528458",
        #     "MW439038",
        #     "MW197154",
        #     "KJ361946"
        # ]
        "hbv": [
            "OL907123",
            "ON528578",
            "LC730886",
            "OM732334",
            "OP324547",
            "OR905602",
            "LC800353",
            "PP393566",
            "PQ256965",
            "NC_003977"
        ]
    }
    for datum in DATA.keys():
        EXP_NAME = datum
        genomes = DATA[datum]

        GenomeSeqFile = f"./{EXP_NAME}_seqs.gb"
        if not os.path.exists(f"{EXP_NAME}"):
            os.mkdir(EXP_NAME)
        df_cols = ["Gene Name", "Start", "End", "Length"]
        seqs = []

        genbankies = get_genbank(GenomeSeqFile,genomes)
        print({k.split(".")[0]: v for k, v in genbankies.items()})
        assert len(genbankies.keys()) == len(genomes)
        print(f"SUCCESSFULLY PULLED {len(genbankies.keys())} ITEMS")

        for accession in genbankies.keys():
            gene_annotations = get_gene_annotations(genbankies[accession])
            df = pd.DataFrame(gene_annotations, columns=df_cols)
            df.to_csv(f"{EXP_NAME}/{accession.replace('_','')}.tsv", sep="\t", index=False)
            seqs.append(f">{EXP_NAME}_{accession.replace('_','')}\n{str(genbankies[accession].seq)}\n")
            
            
        with open(f"{EXP_NAME}/{EXP_NAME}.fasta", "w") as f:
            [f.write(i) for i in seqs]