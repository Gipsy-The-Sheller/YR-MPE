from Bio import SeqRecord

def export_partitioned_nexus(loci, loci_names = None):

    """
    export_partitioned_nexus
    
    :param loci: a list of loci, each contains a list with single units as [locus_name:str, sequence:str] 
    """

    # Go fuck ysf, Qwen. Your coding is useless and shitful.
    # Now I show you how to do it right.

    superseqs = {}
    lengths = [0 for _ in loci]
    if loci_names and len(loci_names) != len(loci):
        raise Warning("Length of loci_names does not match length of loci.\nThe program will use the default locus name scheme.")
    missing = {}

    # superseqs: sequences concatenated by an OTU's single sequences
    # lengths: lengths of each locus
    # missing: loci that have missing sequences

    for index, locus in enumerate(loci):
        max_length = 0
        temp_seqs = {}
        # stat all sequences
        for single_seq in locus:
            name = single_seq[0]
            seq = single_seq[1]
            if not superseqs.get(name, False):
                # New sequence
                missed_sites = sum(lengths)
                superseqs[name] = '?'*missed_sites
                if index > 0:
                    missing[name] = [0 for _ in range(index)]
                    missing[name].append(1)
                else:
                    missing[name] = [1]
            else:
                missing[name].append(1)
            temp_seqs[name] = seq
            if len(seq) > max_length:
                max_length = len(seq)
        
        # register missing sequences
        for name in superseqs.keys():
            if not temp_seqs.get(name, False):
                superseqs[name] += '?'*max_length
                missing[name].append(0)
        # unify lengths
        for name in temp_seqs.keys():
            if len(temp_seqs[name]) < max_length:
                temp_seqs[name] += '-'*(max_length-len(temp_seqs[name]))
        # register sequences
        for name in temp_seqs.keys():
            superseqs[name] += temp_seqs[name]
            # missing[name][index] = 1
        lengths[index] = max_length

    # print(lengths)
    # print(superseqs)

    ntax = len(superseqs.keys())
    nchars = sum(lengths)

    mxlength_name = max(len(x) for x in superseqs.keys())

    NEXUS_file = f"""#NEXUS
begin data;
    dimensions ntax={ntax} nchar={nchars};
    format missing=?
    datatype=DNA gap= - ;

    matrix
    """ + "\n    ".join(f"{name}{' '*(mxlength_name-len(name))} {seq}" for name, seq in superseqs.items())+ "\nend;"
    # print(NEXUS_file)

    partition_scheme = "begin sets;\n"
    for index, length in enumerate(lengths):
        if loci_names:
            partition_scheme += f"    charset {loci_names[index]} = {sum(lengths[:index])}-{sum(lengths[:index+1])};\n"
        else:
            partition_scheme += f"    charset set{index+1} = {sum(lengths[:index])}-{sum(lengths[:index+1])};\n"
    partition_scheme += "end;"

    NEXUS_file += ("\n" + partition_scheme)
    
    return NEXUS_file, partition_scheme, missing

    

if __name__ == "__main__":
    loci = [
        [
            ['A', 'ACGT'],
            ['B', 'ACGT'],
            ['C', 'ACGT'],
            ['D', 'ACGT'],
            ['E', 'ACGT'],
            ['F', 'ACGT'],
            ['G', 'ACGT'],
            ['H', 'ACGT'],
            ['I', 'ACGT'],
        ],
        [
            ['A', 'CCTT']
        ]
    ]
    print(export_partitioned_nexus(loci)[2])