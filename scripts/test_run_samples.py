#!/usr/bin/env python3
# ABOUTME: Self-contained checks for the cohort description loaders -- no search output needed.
# ABOUTME: Run with: python3 scripts/test_run_samples.py

# These loaders decide which runs are injections of one sample, so getting them wrong either
# merges two samples or splits one, and the replicate requirement then rejects or admits the
# wrong calls. Every file they read is optional, and "absent" has to keep meaning the generic
# case rather than becoming an error -- which is what these checks pin.

import os
import sys
import tempfile

import run_samples as rs


def _write(tmpdir, name, text):
    path = os.path.join(tmpdir, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def test_run_key_accepts_every_spelling():
    """DIA-NN spells one run three ways; all must reduce to the same stem."""
    assert rs.run_key("/data/24f201_DIA_NCI_01.raw.dia") == "24f201_DIA_NCI_01"
    assert rs.run_key("24f201_DIA_NCI_01.raw") == "24f201_DIA_NCI_01"
    assert rs.run_key("24f201_DIA_NCI_01") == "24f201_DIA_NCI_01"
    # A sample name containing a dot must survive -- the example files are named this way.
    assert rs.run_key("/d/24f201_DIA_COLO205.1.raw.dia") == "24f201_DIA_COLO205.1"
    return 4


def test_absent_map_is_a_valid_map():
    """No sample map: every run is its own sample, reported under its own stem."""
    m = rs.load_sample_map("")
    assert not m
    assert m.sample("/d/J-PDX0009_1.raw.dia") == "J-PDX0009_1"
    assert m.label("J-PDX0009_1.raw") == "J-PDX0009_1"
    assert rs.load_aliases("") == {}
    assert rs.load_pools("") == {}
    # A configured-but-missing path is the caller's problem to validate, not a crash here.
    assert not rs.load_sample_map("/nonexistent/map.tsv")
    return 6


def test_two_column_map_is_accepted():
    """run_label is optional, so the (run, sample) map the documentation asks for must work."""
    with tempfile.TemporaryDirectory() as d:
        p = _write(d, "map.tsv", "run\tsample\n"
                                 "24f201_DIA_COLO205.1\tCOLO205\n"
                                 "24f201_DIA_COLO205.2\tCOLO205\n")
        m = rs.load_sample_map(p)
        assert m
        assert m.sample("24f201_DIA_COLO205.1.raw.dia") == "COLO205"
        assert m.sample("24f201_DIA_COLO205.2.raw") == "COLO205"
        # With no run_label the run is reported under its own stem.
        assert m.label("24f201_DIA_COLO205.1.raw.dia") == "24f201_DIA_COLO205.1"
    return 4


def test_run_label_overrides_the_reported_name():
    with tempfile.TemporaryDirectory() as d:
        p = _write(d, "map.tsv", "run\tsample\trun_label\n"
                                 "24f201_DIA_NCI_05\tCOLO205\tNCI_05\n")
        m = rs.load_sample_map(p)
        assert m.sample("/x/24f201_DIA_NCI_05.raw.dia") == "COLO205"
        assert m.label("/x/24f201_DIA_NCI_05.raw.dia") == "NCI_05"
        # A run the map does not mention falls back to being its own sample.
        assert m.sample("something_else.raw") == "something_else"
    return 3


def test_missing_required_column_is_an_error():
    """A required column absent is a typo in the file, and must not pass silently."""
    with tempfile.TemporaryDirectory() as d:
        p = _write(d, "map.tsv", "run\tcell\nx\ty\n")
        try:
            rs.load_sample_map(p)
        except SystemExit as exc:
            assert "sample" in str(exc), exc
        else:
            raise AssertionError("a map with no sample column was accepted")
    return 1


def test_pool_members_forms():
    """'*' = any truth mutation counts, empty = none does, else an explicit member list."""
    with tempfile.TemporaryDirectory() as d:
        p = _write(d, "pools.tsv", "pool\tmembers\n"
                                   "NCI7ref\t*\n"
                                   "HeLa\t\n"
                                   "mix\tA549;COLO205\n")
        pools = rs.load_pools(p)
        assert pools["NCI7ref"] is None          # any
        assert pools["HeLa"] == []               # none
        assert pools["mix"] == ["A549", "COLO205"]
        assert "A549" not in pools               # a member is not itself a pool
    return 4


def test_aliases_map_truth_names_onto_run_names():
    with tempfile.TemporaryDirectory() as d:
        p = _write(d, "aliases.tsv", "truth_name\trun_name\nNCIH23\tH23\n")
        assert rs.load_aliases(p) == {"NCIH23": "H23"}
    return 1


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_"):
            continue
        try:
            n = fn()
            print(f"PASS  {name}  ({n} checks)")
        except (AssertionError, SystemExit) as exc:
            # SystemExit too: the loaders exit on a malformed file, and an unexpected one must be
            # reported as this test failing rather than silently ending the whole run.
            failed += 1
            print(f"FAIL  {name}\n      {exc}")
    print("\nall checks passed" if not failed else f"\n{failed} test(s) failed")
    sys.exit(1 if failed else 0)
