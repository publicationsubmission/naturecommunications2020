#!/usr/bin/env python

import sys

import pandas as pd

import utils


def main(args, **kwargs):
    mapping = kwargs["mapping"]
    alignments = utils.read_alignment(args[0], mapping)
    df = pd.DataFrame(
        columns=[
            "Position",
            "Reference",
            "EP-1",
            "EP-2",
            "EP-3",
            "EM-1",
            "EM-2",
            "EM-3",
        ]
    )
    df = utils.tabulate_mutations(alignments, mapping, df)


if __name__ == "__main__":
    mapping = {
        "GCF_000960175.1_ASM96017v1_genomic": "s__Clostridium_sporogenes",
        "GCA_002745415.1_ASM274541v1_genomic": "s__Enterobacteria_phage_lambda",
        "GCF_902161805.1_25426_7_320_genomic": "s__Enterococcus_faecalis",
        "GCF_003316915.1_ASM331691v1_genomic": "s__Lactobacillus_johnsonii",
        "E_plus_M1_S1": "EP-1",
        "E_plus_M2_S2": "EP-2",
        "E_plus_M3_S3": "EP-3",
        "E_minus_M1_S4": "EM-1",
        "E_minus_M2_S5": "EM-2",
        "E_minus_M3_S6": "EM-3",
    }
    main(sys.argv[1:], mapping=mapping)
