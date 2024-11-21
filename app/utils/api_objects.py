
import os 

def castanet_req_body(exp_dir, virus):
    return {
        "ExpDir": f"{os.getcwd()}/{exp_dir}/{exp_dir.split('/')[-1]}_syn_reads", 
        "ExpName": exp_dir.split("/")[1],
        "SaveDir": f'{os.getcwd()}/{exp_dir.split("/")[0]}/{exp_dir.split("/")[0]}_castanet', 
        "RefStem": f"{os.getcwd()}/{virus}/{virus}.fasta", 
        "SingleEndedReads": False,
        "MatchLength": 40,
        "DoTrimming": True,
        "TrimMinLen": 36,
        "DoKrakenPrefilter": False,
        "LineageFile": "data/ncbi_lineages_2023-06-15.csv.gz",
        "ExcludeIds": "9606",
        "RetainIds": "",
        "RetainNames": "",
        "ExcludeNames": "Homo",
        "ConsensusMinD": 10,
        "ConsensusCoverage": 5,
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