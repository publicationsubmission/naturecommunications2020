#!/usr/bin/env python3

import sys

import pandas as pd


def summary(fname, group):
    print(fname)
    df = pd.read_csv(fname, sep="\t", header=None)
    df = df[df[1] != "REF"]
    df = df[df[0] == group]

    return len(set(df[2]))


print(summary(sys.argv[1], "treatment"))
print(summary(sys.argv[2], "treatment"))
