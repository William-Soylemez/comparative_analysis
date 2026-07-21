#!/usr/bin/env python3
"""
Minimal GO-slim mapper — parses go-basic.obo + a slim subset and rolls specific
GO terms up to their high-level slim categories. No external dependencies.

A term maps to every slim category among its ancestors (via is_a + part_of),
which is the standard "all covered" GO-slim binning. Terms whose lineage never
reaches a slim category (e.g. many metazoan-specific terms under a fungal slim)
map to nothing and are naturally dropped.
"""

import re


def _parse_obo(path):
    """Return dicts: name, namespace, parents(is_a+part_of), obsolete set, alt->primary."""
    name, namespace, parents, obsolete, alt2primary = {}, {}, {}, set(), {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {}
                continue
            if line.startswith("["):  # a non-Term stanza (e.g. [Typedef])
                cur = None
                continue
            if cur is None or not line:
                continue
            if line.startswith("id: GO:"):
                cur["id"] = line[4:].strip()
                parents.setdefault(cur["id"], set())
            elif line.startswith("name: "):
                name[cur["id"]] = line[6:].strip()
            elif line.startswith("namespace: "):
                namespace[cur["id"]] = line[11:].strip()
            elif line.startswith("is_a: GO:"):
                parents[cur["id"]].add(line.split()[1])
            elif line.startswith("relationship: part_of GO:"):
                parents[cur["id"]].add(line.split()[2])
            elif line.startswith("alt_id: GO:"):
                alt2primary[line[8:].strip()] = cur["id"]
            elif line == "is_obsolete: true":
                obsolete.add(cur["id"])
    return name, namespace, parents, obsolete, alt2primary


def _slim_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith("id: GO:"):
                ids.add(line[4:].strip())
    return ids


class GoSlim:
    def __init__(self, obo_path, slim_path):
        (self.name, self.namespace, self._parents,
         self.obsolete, self._alt) = _parse_obo(obo_path)
        self.slim = _slim_ids(slim_path) & set(self._parents)  # slim ids present in ontology
        self._anc_cache = {}

    def _ancestors(self, term):
        if term in self._anc_cache:
            return self._anc_cache[term]
        seen, stack = set(), [term]
        while stack:
            t = stack.pop()
            for p in self._parents.get(t, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        self._anc_cache[term] = seen
        return seen

    def map_term(self, term):
        """Return the set of slim categories a specific GO term rolls up to."""
        term = self._alt.get(term, term)          # resolve alt ids
        if term not in self._parents:
            return set()                          # unknown / not in this ontology
        covered = ({term} | self._ancestors(term)) & self.slim
        return covered

    def slim_name(self, sid):
        return self.name.get(sid, sid)


if __name__ == "__main__":
    # quick validation on a few known terms
    gs = GoSlim("go_data/go-basic.obo", "go_data/goslim_yeast.obo")
    print(f"slim categories loaded: {len(gs.slim)}")
    tests = {
        "GO:0006260": "DNA replication",
        "GO:0042254": "ribosome biogenesis",
        "GO:0007507": "heart development (metazoan artifact)",
        "GO:0015706": "nitrate transmembrane transport",
        "GO:0006508": "proteolysis",
    }
    for gid, desc in tests.items():
        slims = gs.map_term(gid)
        print(f"\n{gid}  ({desc})")
        for s in sorted(slims):
            print(f"    -> {s}  {gs.slim_name(s)} [{gs.namespace.get(s,'?')}]")
        if not slims:
            print("    -> (no fungal-slim category — dropped)")
