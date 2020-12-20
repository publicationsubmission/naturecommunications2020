import pathlib
import sys

import pandas as pd


def classify(fname, synonymous_out, nonsynonymous_out, empty_out):
    cols = {
        0: "QueryID",
        1: "SubjectID",
        2: "PercIdentity",
        3: "AlignmentLength",
        4: "Mismatches",
        5: "GapOpens",
        6: "QueryStart",
        7: "QueryEnd",
        8: "SubjectStart",
        9: "SubjectEnd",
        10: "Evalue",
        11: "BitScore",
    }
    try:
        df = pd.read_csv(fname.as_posix(), comment="#", sep="\t", header=None)
        df.rename(columns=cols, inplace=True)
        df.sort_values(by=["PercIdentity"], inplace=True, ascending=False)
        label = df.iloc[0, 0]
        sample, pos = label.split(",")
        sample = sample.split(":")[-1]
        pos = pos.split(":")[-1]

        effect = fname.name.split("_")[0]
        if df.iloc[0, 2] != 100:
            with open(nonsynonymous_out, "a+") as handle:
                print(
                    "{}\t{}\t{}".format(effect, sample, pos, "non-synonymous"),
                    file=handle,
                )
        else:
            with open(synonymous_out, "a+") as handle:
                print(
                    "{}\t{}\t{}".format(effect, sample, pos, "synonymous"), file=handle
                )
    except pd.errors.EmptyDataError:
        with open(empty_out, "a+") as handle:
            print(fname, "0 hits", file=handle)


def main():
    files = list(pathlib.Path(sys.argv[1]).rglob("*blast.txt"))
    synonymous_out = str(pathlib.Path(sys.argv[1]) / pathlib.Path("synonymous.tab"))
    nonsynonymous_out = str(
        pathlib.Path(sys.argv[1]) / pathlib.Path("nonsynonymous.tab")
    )
    empty_out = str(pathlib.Path(sys.argv[1]) / pathlib.Path("empty.tab"))
    for i, fname in enumerate(files):
        print("{}/{}".format(i + 1, len(files)))
        classify(fname, synonymous_out, nonsynonymous_out, empty_out)


if __name__ == "__main__":
    main()
