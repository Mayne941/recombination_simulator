import subprocess as sp

def save_fa(fpath, pat):
    with open(fpath, "w") as f:
        f.write(pat)

def read_fa(fpath):
    '''Read fasta file to list of lists in format [[>name, seq], [...]]'''
    seqs = []
    try:
        with open(fpath, "r") as f:
            for l in f:
                if l[0] == ">":
                    seqs.append(f"~{l}")
                else:
                    seqs.append(l.replace("\n", ""))
    except:
        '''Current release of ViralConsensus will sometimes dump binary into the output file - this is a temp fix'''
        print(
            f"*** WARNING: Encoding on your input fasta file ({fpath}) is messed up. Attempting to rescue it...")
        seqs = []  # RESET SEQS as it doesn't get wiped by transition to exception
        bases = ["A", "T", "C", "G", "N", "-"]
        with open(fpath, "r", encoding="ISO-8859-1") as f:
            for l in f.readlines():
                if l[0] not in bases:
                    seqs.append(f"~{l}")
                else:
                    seqs.append(l)

        for i in range(len(seqs)):
            if seqs[i][0] == "~":
                continue
            else:
                seqs[i] = "".join([j for j in seqs[i] if j in bases])

    seqs_split = [i.split("\n") for i in "".join(seqs).split("~")]
    return [i for i in seqs_split if not i == [""]]

def shell(args, calling_fn="Misc shell function", ret_output=False, is_test=False):
    '''Call Bash shell with input string as argument'''
    whitelist = ["mafft"]
    _ = sp.Popen(args, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = _.communicate()
    if not any(x.lower() in whitelist for x in whitelist):
        raise f"ERROR: Something went wrong with CLI call: {calling_fn}, dumping STDOUT/STDERR to shell."
    if ret_output:
        return out
    elif is_test:
        return str(out+err)
