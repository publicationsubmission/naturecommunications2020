import json
import pathlib
import sys
from io import StringIO

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

NUM_THREADS = 36

subject = sys.argv[1]
with open(sys.argv[2], "r") as handle:
    content = json.load(handle)
    group = sys.argv[3]

    for key in content.keys():
        path = pathlib.Path(key)
        path.mkdir(exist_ok=True)
        for pos in content[key].keys():
            print(pos)
            seq = content[key][pos]
            query = SeqRecord(Seq(seq), id="Sample:{},Pos:{}".format(key, pos))
            SeqIO.write(query, "{}/query.fa".format(key), "fasta")
            output = NcbiblastxCommandline(
                query="{}/query.fa".format(key),
                subject=subject,
                out="{}/{}_{}_blast.txt".format(key, group, pos),
                evalue=10,
                outfmt=7,
                num_threads=NUM_THREADS,
            )()
