import networkx as nx
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Interaction:
    source: str
    target: str
    enzyme: str
    effect: str
    severity: str
    source_reference: str


class ClinicalRippleEngine:
    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.pair_registry: Dict[Tuple[str, str], Interaction] = {}

    def add_interaction(self, interaction: Interaction) -> None:
        source = interaction.source.strip()
        enzyme = interaction.enzyme.strip()
        target = interaction.target.strip()

        # Register direct pair evidence
        pair_key = (source.lower(), target.lower())
        self.pair_registry[pair_key] = interaction

        # Build metabolic pathway topology
        self.graph.add_node(source, kind="entity")
        self.graph.add_node(enzyme, kind="enzyme")
        self.graph.add_node(target, kind="drug_or_herb")

        # Edge from herb to enzyme: modulation of metabolic machinery
        self.graph.add_edge(
            source,
            enzyme,
            effect="modulates",
            enzyme=enzyme,
        )

        # Edge from enzyme to drug: clearance/metabolism
        self.graph.add_edge(
            enzyme,
            target,
            effect="metabolizes",
            enzyme=enzyme,
        )

    def add_many(self, interactions: Iterable[Interaction]) -> None:
        for interaction in interactions:
            self.add_interaction(interaction)

    def screen_pair(self, herb: str, drug: str) -> dict:
        """
        Screen a specific (herb, drug) pair with strict evidence attribution.
        Distinguishes direct empirical pair evidence from indirect enzyme hypothesis.
        """
        h_clean = herb.strip()
        d_clean = drug.strip()
        pair_key = (h_clean.lower(), d_clean.lower())

        # 1. Direct Clinical Pair Evidence
        if pair_key in self.pair_registry:
            direct_entry = self.pair_registry[pair_key]
            return {
                "status": "direct_evidence",
                "herb": direct_entry.source,
                "drug": direct_entry.target,
                "enzyme": direct_entry.enzyme,
                "effect": direct_entry.effect,
                "severity": direct_entry.severity,
                "source_reference": direct_entry.source_reference,
                "is_direct_evidence": True,
            }

        # 2. Check for indirect mechanistic pathway via shared enzyme
        node_map = {n.lower(): n for n in self.graph.nodes}
        h_node = node_map.get(h_clean.lower(), h_clean)
        d_node = node_map.get(d_clean.lower(), d_clean)

        if h_node in self.graph and d_node in self.graph:
            try:
                path = nx.shortest_path(self.graph, h_node, d_node)
                # Path format: Herb -> Enzyme -> Drug
                intermediate_enzymes = [
                    node for node in path[1:-1]
                    if self.graph.nodes[node].get("kind") == "enzyme"
                ]
                enzyme_name = intermediate_enzymes[0] if intermediate_enzymes else "shared metabolic pathway"
                return {
                    "status": "indirect_mechanism_only",
                    "herb": h_node,
                    "drug": d_node,
                    "enzyme": enzyme_name,
                    "effect": f"Shared substrate/metabolic pathway via {enzyme_name}",
                    "severity": "unknown",
                    "source_reference": (
                        f"No direct clinical study found for {h_clean} + {d_clean}. "
                        f"Mechanism hypothesized based on shared {enzyme_name} pathway."
                    ),
                    "is_direct_evidence": False,
                }
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                pass

        # 3. No direct or indirect connection found
        return {
            "status": "no_interaction_found",
            "herb": h_clean,
            "drug": d_clean,
            "enzyme": None,
            "effect": "No direct or indirect metabolic interaction found in seed database",
            "severity": "none",
            "source_reference": "No empirical or mechanistic entry",
            "is_direct_evidence": False,
        }

    def ripple_frontier(self, sources: Iterable[str]) -> dict[str, list[dict]]:
        """
        Traverse graph to find reachable downstream targets from sources,
        attaching direct pair evidence when present, or marking as mechanism hypothesis.
        """
        frontier: dict[str, list[dict]] = {}
        for source in sources:
            s_clean = source.strip()
            if s_clean not in self.graph:
                frontier[s_clean] = []
                continue

            frontier[s_clean] = []
            for target_node in nx.descendants(self.graph, s_clean):
                if self.graph.nodes[target_node].get("kind") == "enzyme":
                    continue  # Only report clinical entity targets

                screen_res = self.screen_pair(s_clean, target_node)
                record = {
                    "from": s_clean,
                    "node": target_node,
                    "enzyme": screen_res["enzyme"],
                    "effect": screen_res["effect"],
                    "severity": screen_res["severity"],
                    "source_reference": screen_res["source_reference"],
                    "is_direct_evidence": screen_res["is_direct_evidence"],
                }
                if record not in frontier[s_clean]:
                    frontier[s_clean].append(record)
        return frontier
