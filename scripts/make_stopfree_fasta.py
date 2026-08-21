#!/usr/bin/env python3
# ABOUTME: Builds the stop-codon-free search database by dropping the stop-gain variant entries
# ABOUTME: from proteome_combined.fasta, which the canonical filter can never let through anyway.

# A stop-gain entry is a truncated N-terminal prefix of its wild-type protein, so every tryptic
# peptide it yields is a substring of the canonical proteome and filter_canonical_peptides.sh
# discards it in every run. Those entries can therefore never produce a reported call, while still
# competing in the DIA-NN search, forming GENE;GENE_var protein groups, and sitting in the recall
# denominator as guaranteed misses. Removing them at the FASTA is what makes the whole pipeline
# stop-codon free: both DIA-NN passes, the accession->gene map and the fragment-geometry
# wild-type lookup all read this one file.

import argparse
import os
import re
import sys

# Stop-gain accession, as it appears in '>sp|P60484_R130_*:38|PTEN_HUMAN': base accession, the
# reference residue and position, then '*' for the new stop and the recurrence count. Anchored on
# the whole accession so a '*' anywhere else in the header cannot trigger a removal. Every one of
# the 242 stop-gain entries matches this and no indel or delins entry introduces a stop, so the
# rule needs no second pattern.
STOP_GAIN = re.compile(r'_\*:\d+$')

# No default input: the combined library with the stop-gain entries still in it is an input to
# this repo, not an output of it, so there is no path here that would be right for a reader.


def accession(header):
    """'>sp|P60484_R130_*:38|PTEN_HUMAN Phosphatidyl...' -> 'P60484_R130_*:38'."""
    fields = header[1:].split('|')
    return fields[1] if len(fields) > 1 else header[1:].split()[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True,
                    help="Combined proteome to filter, with the stop-gain variant entries in it")
    ap.add_argument("-o", "--output",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         os.pardir, "data", "fasta", "proteome.fasta"),
                    help="Stop-codon-free FASTA to write (default: %(default)s)")
    ap.add_argument("--expect-dropped", type=int, default=242,
                    help="Stop-gain entries the input is known to hold; the script fails if the "
                         "count differs, so a changed library cannot pass unnoticed "
                         "(default: %(default)s)")
    args = ap.parse_args()

    kept = dropped = 0
    tmp = args.output + ".tmp"
    with open(args.input) as fin, open(tmp, "w") as fout:
        keeping = False
        for line in fin:
            if line.startswith(">"):
                keeping = not STOP_GAIN.search(accession(line.rstrip()))
                kept += keeping
                dropped += not keeping
            if keeping:
                fout.write(line)

    if dropped != args.expect_dropped:
        os.remove(tmp)
        sys.exit(f"ERROR: dropped {dropped} stop-gain entries, expected {args.expect_dropped}. "
                 f"The input library has changed -- rerun with --expect-dropped {dropped} once "
                 f"you have confirmed that is right.")

    os.replace(tmp, args.output)

    # Read back rather than trust the counters: this file is the input to ~9 h of search time.
    with open(args.output) as fh:
        headers = [l.rstrip() for l in fh if l.startswith(">")]
    remaining = [h for h in headers if STOP_GAIN.search(accession(h))]
    if remaining:
        sys.exit(f"ERROR: {len(remaining)} stop-gain entries survived, e.g. {remaining[0]}")

    print(f"{args.input}\n  -> {args.output}")
    print(f"  kept    {len(headers):,} entries")
    print(f"  dropped {dropped:,} stop-gain entries")


if __name__ == "__main__":
    main()
