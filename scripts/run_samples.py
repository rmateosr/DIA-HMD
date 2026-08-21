#!/usr/bin/env python3
# ABOUTME: Cohort description read from optional TSVs -- which runs belong to which sample, which
# ABOUTME: truth-table names are spelled differently, and which samples are pools of others.

# DIA-NN names one run in three forms: the pr_matrix column is a full path ending .raw.dia, the
# parquet Run column ends .raw, and reports print the bare stem. Everything here keys on the bare
# stem so the gate and the classifier cannot disagree about which run a cell belongs to.
#
# All three files are optional, and absent means the generic case rather than an error:
#   no sample map -> every run is its own sample, so nothing is grouped that the user did not ask
#                    to group. An injection structure cannot be guessed from file names: a cohort
#                    naming its replicates S1_1/S1_2 and one naming distinct samples P_1/P_2 are
#                    indistinguishable, and guessing wrong either merges two samples or splits one.
#   no aliases    -> truth-table names are used as they are written
#   no pools      -> no sample is treated as a mixture
#
# This is what replaced the two hardcoded cohort forks of the classifier, which between them knew
# the string "24f201_DIA_", the regex _[12]$, one cell-line spelling and three pool names.

import os

# Longest first: a pr_matrix column ends .raw.dia and must not be read as ending .dia.
RUN_SUFFIXES = (".raw.dia", ".raw", ".dia", ".mzML", ".mzml", ".d")

# A pool whose members are given as this accepts a mutation known in ANY sample of the truth
# table. Used when the mixture's composition is not recorded: it is the conservative reading,
# since it can only turn a would-be false positive into an unscored pool hit.
ANY_MEMBER = "*"


def run_key(name):
    """Any of DIA-NN's three spellings of one run -> the bare stem.

    '/data/24f201_DIA_NCI_01.raw.dia' -> '24f201_DIA_NCI_01'
    'J-PDX0009_1.raw'                 -> 'J-PDX0009_1'
    """
    base = os.path.basename(name)
    for suffix in RUN_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def _read_tsv(path, columns):
    """Yield one tuple per data row, or nothing at all if the file is not there.

    Written without pandas so the loaders stay usable from a script that has not imported it,
    and so a missing optional column is reported by name rather than as a KeyError.
    """
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                missing = [c for c in columns if c not in header]
                if missing:
                    raise SystemExit(
                        f"ERROR: {path} is missing the column(s) {', '.join(missing)}. "
                        f"Found: {', '.join(header)}"
                    )
                continue
            row = dict(zip(header, fields))
            yield tuple(row.get(c, "") for c in columns)


class SampleMap:
    """Run -> sample grouping, plus the label a run is reported under.

    An absent map is a valid map: every run is its own sample and is reported by its stem. That
    is the right default for a user whose replicate structure we do not know, and it makes the
    replicate requirement (--min-replicates) a no-op rather than a silent rejection.
    """

    def __init__(self, rows=None):
        self._sample = {}
        self._label = {}
        for run, sample, label in rows or []:
            key = run_key(run)
            self._sample[key] = sample
            self._label[key] = label or key

    def __bool__(self):
        return bool(self._sample)

    def sample(self, run):
        key = run_key(run)
        return self._sample.get(key, key)

    def label(self, run):
        key = run_key(run)
        return self._label.get(key, key)


def load_sample_map(path):
    """TSV: run, sample, and optionally run_label for how the run is printed in reports."""
    rows = []
    for run, sample, label in _read_tsv(path, ["run", "sample", "run_label"]):
        if not run:
            continue
        rows.append((run, sample or run_key(run), label))
    return SampleMap(rows)


def load_aliases(path):
    """TSV: truth_name, run_name. Maps a truth-table sample name onto the name the runs use."""
    return {t: r for t, r in _read_tsv(path, ["truth_name", "run_name"]) if t}


def load_pools(path):
    """TSV: pool, members. '*' = any truth mutation counts, empty = none does, else 'A;B'."""
    pools = {}
    for pool, members in _read_tsv(path, ["pool", "members"]):
        if not pool:
            continue
        members = members.strip()
        if members == ANY_MEMBER:
            pools[pool] = None
        else:
            pools[pool] = [m for m in members.split(";") if m]
    return pools
