#!/usr/bin/env python
import functools
import itertools
import json
import sys
from collections import Counter, defaultdict, namedtuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import venn


def get_sets(sets):
    """Finds the elements in each set

    :param graphs: list of sets
    :returns: elements in each set
    :rtype: dictionary

    """
    N = len(sets)
    union = functools.reduce(lambda x, y: x | y, sets)
    result = {}
    for n in range(1, 2 ** N):
        key = bin(n).split("0b")[-1].zfill(N)
        value = union
        sets_for_intersection = [sets[i] for i in range(N) if key[i] == "1"]
        sets_for_difference = [sets[i] for i in range(N) if key[i] == "0"]

        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s

        result[key] = len(value)

    return result


def is_treatment(sets):
    return sets["110"] > 0 or sets["101"] > 0


def is_other(sets):
    return (sets["100"] > 0 and sets["010"] > 0 and sets["001"] > 0) or sets["011"] > 0


# def venn_mutations(df, counts, countpos, title):
#     """
#     counts: file name
#     countpos: file name
#     """
#     reference = df.filter(regex="Ref*")
#     positive = df.filter(regex="EP.*")
#     negative = df.filter(regex="EM.*")

#     # initialization
#     # first row bases
#     ref_init = set(list(reference.iloc[0].values))
#     pos_init = set(list(positive.iloc[0].values))
#     neg_init = set(list(negative.iloc[0].values))

#     result = Counter(get_sets([ref_init, pos_init, neg_init]))
#     resultpos = {k: [] for k in result.keys()}
#     for k in result.keys():
#         if result[k] > 0:
#             resultpos[k].append(1)

#     # next row bases
#     for i in range(1, df.shape[0]):
#         if i == 52643:
#             import pdb; pdb.set_trace()
#         ref = set(list(reference.iloc[i].values))
#         pos = set(list(positive.iloc[i].values))
#         if len(pos) > 1 and len(pos & ref) > 0:
#             pos = pos - ref
#         neg = set(list(negative.iloc[i].values))
#         if len(neg) > 1 and len(neg & ref) > 0:
#             neg = neg - ref

#         sets = Counter(get_sets([ref, pos, neg]))
#         for k in sets.keys():
#             if sets[k] > 0:
#                 resultpos[k].append(i + 1)

#         result.update(sets)

#     # we also exclude 100 set because it contains the positions of every row.
#     # Specifically, we update pos = pos-ref and neg=neg-ref, so positive or negative sets
#     # don't contain ref
#     number_of_mutations = sum(
#         [v for k, v in result.items() if k != "111" and k != "100"]
#     )
#     print(f"{number_of_mutations}/{df.shape[0]} are mutations")
#     result = {k: f"{k}:{result[k]}" for k in result.keys()}
#     venn.venn3(
#         result, names=["Reference", "Positive Group", "Negative Group"], setfontsize=8
#     )
#     plt.title(f"{title} ({number_of_mutations} mutations)")
#     plt.show()

#     with open(counts, "w+") as countfile:
#         json.dump(result, countfile, indent=4)

#     with open(countpos, "w+") as countposfile:
#         json.dump(resultpos, countposfile, indent=4)


def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def venn_mutations(df, counts):
    """
    counts: file name
    countpos: file name
    """
    reference = df.filter(regex="Ref*")
    positive = df.filter(regex="EP.*")
    negative = df.filter(regex="EM.*")

    # Ref-Pos-Neg
    result = {"treatment": set(), "other": set()}
    for i in range(df.shape[0]):
        ref = set(list(reference.iloc[i].values))
        pos = set(list(positive.iloc[i].values))
        neg = set(list(negative.iloc[i].values))

        for bases in itertools.product(ref, pos, neg):
            R, P, N = bases
            sets = get_sets([set(R), set(P), set(N)])
            if is_treatment(sets):
                result["treatment"].add(i + 1)
            elif is_other(sets):
                result["other"].add(i + 1)

    with open(counts, "w+") as countfile:
        json.dump(result, countfile, indent=4, default=set_default)

    for k, v in result.items():
        print("{}:{}".format(k, len(v)))


def merge_intervals(intervals):
    """ Merges the intervals

    :param intervals: intervals
    :returns: merged intervals
    :rtype: list

    """

    intervals_ = sorted(intervals, key=lambda x: x[0])
    merged = []
    for interval in intervals_:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])

    merged = [(x[0], x[1]) for x in merged]

    return merged


def coding_region(df, blastoutput):
    columns = [
        "QueryId",
        "SubjectId",
        "PercIdentity",
        "AlignmentLen",
        "Mismatches",
        "GapOpens",
        "QueryStart",
        "QueryEnd",
        "SubjectStart",
        "SubjectEnd",
        "Evalue",
        "BitScore",
    ]
    blast_output = pd.read_csv(blastoutput, names=columns, sep="\t", comment="#")
    blast_ref_hits = blast_output[
        blast_output.QueryId.str.contains("^(?!E_).*", regex=True)
    ]

    blast_ref_hits = blast_ref_hits.sort_values("QueryStart")

    query_start_end = sorted(
        set([(x[0], x[1]) for x in blast_ref_hits.iloc[:, 6:8].values]),
        key=lambda k: k[0],
    )
    query_start_end = [[x[0], x[1]] for x in query_start_end]
    query_start_end_merged = merge_intervals(query_start_end)

    gene_coding = [
        True
        if len(list(filter(lambda x: x[0] <= i <= x[1], query_start_end_merged))) > 0
        else False
        for i in df.Position.values
    ]

    df["GeneCoding"] = gene_coding

    # return df[df["GeneCoding"] == True], query_start_end_merged
    return df

    # occurrence = defaultdict(list)
    # for row in df[df["GeneCoding"] == True].iterrows():
    #     for interval in query_start_end_merged:
    #         if interval[0] <= row[1]["Position"] <= interval[1]:
    #             occurrence["{}-{}".format(interval[0], interval[1])].append(
    #                 row[1]["Position"]
    #             )
    #             break

    # with open(fname, "w+") as outfile:
    #     json.dump(occurrence, outfile, indent=4)

    # intervals = [f"{interval[0]}-{interval[1]}" for interval in occurrence.keys()]
    # x_pos = np.arange(len(intervals))
    # mutations = list(occurrence.values())

    # plt.bar(
    #     np.arange(len(intervals)), list(occurrence.values()), align="center", alpha=0.5
    # )
    # plt.xticks(x_pos, intervals, rotation="vertical")
    # plt.ylabel("Number of Mutations")

    # for a, b in zip(x_pos, mutations):
    #     plt.text(a, b, str(b))

    # plt.title(sys.argv[3])

    # plt.show()


def is_genecoding(df, pos, label):
    """
    df is a pd.DataFrame
    pos is a list
    """
    gene_coding = df[df["Position"].isin(pos)]
    gene_coding = gene_coding[gene_coding["GeneCoding"] == True]
    print(f"{label} => {gene_coding.shape[0]}/{len(pos)} are in coding region")

    return gene_coding


Strain = namedtuple("Strain", ["mutations", "blast", "title"])

s1 = Strain(
    mutations="output/s__Clostridium_sporogenes.csv",
    blast="blast/blast_Clostridium_sporogenes.tab",
    title="Clostridium sporogenes",
)
s2 = Strain(
    mutations="output/s__Enterobacteria_phage_lambda.csv",
    blast="blast/blast_Enterobacteria_phage_lambda.tab",
    title="Enterobacteria phage lambda",
)

s3 = Strain(
    mutations="output/s__Enterococcus_faecalis.csv",
    blast="blast/blast_Enterococcus_faecalis.tab",
    title="Enterococcus faecalis",
)

s4 = Strain(
    mutations="output/s__Lactobacillus_johnsonii.csv",
    blast="blast/blast_Lactobacillus_johnsonii.tab",
    title="Lactobacillus johnsonii",
)

for s in [s1, s2, s3, s4]:
    print(s.title)
    print("=" * 80)
    df = pd.read_csv(s.mutations)
    df.dropna(axis="columns", inplace=True)

    fname = s.title.lower().replace(" ", "_")
    venn_mutations(df, f"/tmp/output/{fname}.json")

    df_coding = coding_region(df, s.blast)

    with open(f"/tmp/output/{fname}.json", "r") as f:
        sets = json.load(f)
        # sets = json.load(open(f"output/resultpos.json", "r"))
        # treatment_effects = ["110", "101"]
        # other_effects = ["011", "001", "010"]

        # genecoding_pos = {}
        # for effect in treatment_effects:
        #     pos = sets[effect]
        #     set_df = is_genecoding(df_coding, pos, effect)
        #     genecoding_pos[effect] = set_df["Position"].values.tolist()

        # for effect in other_effects:
        #     pos = sets[effect]
        #     is_genecoding(df_coding, pos, effect)
        #     genecoding_pos[effect] = set_df["Position"].values.tolist()

        # with open(f"/tmp/output/{fname}_genecoding_pos.json", "w") as handle:
        #     json.dump(genecoding_pos, handle, indent=4)

        genecoding_pos = {}
        for effect in sets.keys():
            pos = sets[effect]
            set_df = is_genecoding(df_coding, pos, effect)
            genecoding_pos[effect] = set_df["Position"].values.tolist()

        with open(f"/tmp/output/{fname}_genecoding_pos.json", "w") as handle:
            json.dump(genecoding_pos, handle, indent=4)

# outputmutations = "output/s__Clostridium_sporogenes.csv"
# blastoutput = "blast/blast_Clostridium_sporogenes.tab"
# title = "Clostridium sporogenes"

# outputmutations = "output/s__Enterobacteria_phage_lambda.csv"
# blastoutput = "blast/blast_Enterobacteria_phage_lambda.tab"
# title = "Enterobacteria phage lambda"

# outputmutations = "output/s__Enterococcus_faecalis.csv"
# blastoutput = "blast/blast_Enterococcus_faecalis.tab"
# title = "Enterococcus faecalis"

# outputmutations = "output/s__Lactobacillus_johnsonii.csv"
# blastoutput = "blast/blast_Lactobacillus_johnsonii.tab"
# title = "Lactobacillus johnsonii"


# df = pd.read_csv(outputmutations)
# df.dropna(axis="columns", inplace=True)

# fname = title.lower().replace(" ", "_")
# venn_mutations(df, f"output/{fname}.json", "output/{fname}pos.json", title)

# df_coding = coding_region(df, blastoutput)

# sets = json.load(open("/tmp/resultpos.json", "r"))
# treatment_effects = ["110", "101"]
# other_effects = ["011", "001", "010"]

# for effect in treatment_effects:
#     pos = sets[effect]
#     is_genecoding(df_coding, pos, effect)

# for effect in other_effects:
#     pos = sets[effect]
#     is_genecoding(df_coding, pos, effect)
