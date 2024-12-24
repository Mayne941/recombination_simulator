
import os 

def castanet_req_body(exp_dir, data_dir, virus):
    return {
        "ExpDir": f"{os.getcwd()}/{exp_dir}/{exp_dir.split('/')[-2]}_syn_reads", 
        "ExpName": f"{exp_dir.split('/')[0]}_{exp_dir.split('/')[-2]}_castanet",
        "SaveDir": f'{os.getcwd()}/{"/".join(exp_dir.split("/")[0:3])}/', 
        "RefStem": f"{os.getcwd()}/{data_dir}/{virus}/{virus}.fasta", 
        "Mapper": "bowtie2",
        "SingleEndedReads": False,
        "MatchLength": 50,
        "DoTrimming": True,
        "TrimMinLen": 36,
        "DoKrakenPrefilter": False,
        "LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
        "ExcludeIds": "9606",
        "RetainIds": "",
        "RetainNames": "",
        "ExcludeNames": "Homo",
        "ConsensusMinD": 5,
        "ConsensusCoverage": 1,
        "ConsensusMapQ": 1,
        "ConsensusCleanFiles": True,
        "GtFile": "",
        "GtOrg": "",
        "KrakenDbDir": f"{os.getcwd()}", 
        "KeepDups": True,
        "Clin": "",
        "DepthInf": "",
        "SamplesFile": "",
        "PostFilt": False,
        "AdaptP": "data/all_adapters.fa",
        "NThreads": "auto"
        }