#!/usr/bin/env python

import os

import pandas as pd

from sys import argv

from Bio import AlignIO

multiple_alignments = AlignIO.read(argv[1], "fasta")
fname = os.path.splitext(argv[1])[0].split("/")[1]


seqs = {}

for alignment in multiple_alignments:
    seqs[alignment.id] = str(alignment.seq)

ref = list(filter(lambda x: not x.startswith("E_"), seqs.keys()))[0]
print(ref)

for k, v in seqs.items():
    print("{:50}: {} bp".format(k, len(v)))

result = []
for i in range(len(seqs[ref])):
    for k, v in seqs.items():
        if k != ref:
            sample_base = seqs[k][i]
            ref_base = seqs[ref][i]
            if sample_base != ref_base:
                result.append([k, "{}/{}".format(ref_base, sample_base), i])

df = pd.DataFrame(
    result, columns=["Sample", "Mutation (Ref / Sample)", "Position in Reference"]
)

print(df)

df.to_csv("/tmp/{}.csv".format(fname), index=False)

# def main(args):
#     pass


# if __name__ == "__main__":
#     main(sys.argv[1:])
