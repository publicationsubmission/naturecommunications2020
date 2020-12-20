import glob
import itertools
import json
import pathlib
from collections import defaultdict, namedtuple

from Bio import SeqIO

# TREATMENT_MUTATIONS = ["110", "101"]
# OTHER_MUTATIONS = ["011", "001", "010"]

MAPPING = {
    "E_plus_M1_S1": "P1",
    "E_plus_M2_S2": "P2",
    "E_plus_M3_S3": "P3",
    "E_minus_M1_S4": "N1",
    "E_minus_M2_S5": "N2",
    "E_minus_M3_S6": "N3",
}


class Strain:
    def __init__(self, fasta, pos, label):
        """Constructor of a Strain class

        :param fasta: fasta file of a strain where reads are samples
        :param pos: json file of a mutation positions for both treatment and other effects

        """
        self.fasta = fasta
        self.pos = pos
        self.label = label
        self.treatment_mutations, self.other_mutations = self.merge_group_mutations()

    def merge_group_mutations(self):
        """ merges treatment and other effects mutations respectively

        :returns: treatment and other mutations
        :rtype: tuple

        """
        treatment = []
        other = []
        with open(self.pos, "r") as handle:
            positions = json.load(handle)
            for k, v in positions.items():
                # if k in TREATMENT_MUTATIONS:
                if k == "treatment":
                    treatment.append(v)
                # elif k in OTHER_MUTATIONS:
                elif k == "other":
                    other.append(v)

        self.treatment_mutations = sorted(itertools.chain(*treatment))
        self.other_mutations = sorted(itertools.chain(*other))

        return self.treatment_mutations, self.other_mutations

    def adjust_boundaries(self, pos, size, max_size):
        """calculates and adjusts (pos-size, pos+size) if they exceed limits (0, max_size)

        :param pos: int
        :param size: size is a window length
        :param max_size: int
        :returns: adjusted start and end positions
        :rtype: tuple

        """
        s, e = pos - size, pos + size
        start = 0 if s < 0 else s
        end = max_size if e > max_size else e

        return start, end

    def extract_seqs(self, read, size, pos):
        """extract DNA sequences for given position

        :param read: fasta read
        :param pos: indexing starts at 1
        :returns: seqs
        :rtype: dictionary where key is pos and value is substring within window

        """
        seqs = dict()
        max_size = len(read)

        for p in pos:
            start, end = self.adjust_boundaries(p - 1, size, max_size)
            seqs[p] = read[start : end + 1]

        return seqs

    def get_mutation_seqs(self, size):
        treatment_seqs = dict()
        other_seqs = dict()
        with open(self.fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                name = MAPPING.get(record.id, "REF")
                treatment_seqs[name] = self.extract_seqs(
                    str(record.seq), size, self.treatment_mutations
                )
                other_seqs[name] = self.extract_seqs(
                    str(record.seq), size, self.other_mutations
                )

        return treatment_seqs, other_seqs


def export_seqs(seqs, fname):
    with open(fname, "w+") as handle:
        json.dump(seqs, handle, sort_keys=True, indent=4)


def main():
    window = 100
    fasta_files = sorted(glob.glob("fasta/*.fa"))
    pos_files = sorted(glob.glob("/tmp/output/*_genecoding_pos.json"))
    for fasta, pos in zip(fasta_files, pos_files):
        print(fasta, pos)
        s = Strain(fasta, pos, pathlib.Path(pos).stem)
        treatment_seq, other_seq = s.get_mutation_seqs(window)

        export_seqs(treatment_seq, f"/tmp/output/{s.label}_treatment.json")
        export_seqs(other_seq, f"/tmp/output/{s.label}_other.json")


if __name__ == "__main__":
    main()
