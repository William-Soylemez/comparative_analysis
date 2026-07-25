#!/usr/bin/env python3
"""
Restrict each species' network to its largest connected component (LCC) and
write it out under a "_lcc" alias, so run_full_alignment.py can be pointed
at it with no other changes (it just reads net/<name>.tsv by name).

Requested as a check on whether running IsoRank on a disconnected graph
causes problems -- most of these networks turn out to already be a single
connected component (see printed stats), so this mainly matters for
ncrassa, which has one tiny 2-node island next to the giant component.

Usage:
    python3 make_lcc_networks.py <species> [<species> ...]
Expects net/<species>.tsv (2-column edge list, no header).
Writes net/<species>_lcc.tsv.
"""

import sys
import pandas as pd
import networkx as nx

for sp in sys.argv[1:]:
    df = pd.read_csv(f"net/{sp}.tsv", sep="\t", header=None)
    g = nx.Graph()
    g.add_edges_from(df.values)

    components = list(nx.connected_components(g))
    lcc = max(components, key=len)

    kept = df[df[0].isin(lcc) & df[1].isin(lcc)]
    kept.to_csv(f"net/{sp}_lcc.tsv", sep="\t", header=False, index=False)

    print(f"{sp}: {g.number_of_nodes()} nodes / {g.number_of_edges()} edges, "
          f"{len(components)} component(s) -> LCC {len(lcc)} nodes / {len(kept)} edges "
          f"({g.number_of_nodes() - len(lcc)} nodes dropped)")
