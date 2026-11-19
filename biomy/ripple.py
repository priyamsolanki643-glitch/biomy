import networkx as nx
from dataclasses import dataclass
from typing import Iterable


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

    def add_interaction(self, interaction: Interaction) -> None:
        source = interaction.source.strip()
        enzyme = interaction.enzyme.strip()
        target = interaction.target.strip()
        self.graph.add_node(source, kind="entity")
        self.graph.add_node(enzyme, kind="enzyme")
        self.graph.add_node(target, kind="drug_or_herb")
        self.graph.add_edge(
            source,
            enzyme,
            effect="modulates",
            severity=interaction.severity,
            source_reference=interaction.source_reference,
            enzyme=enzyme,
        )
        self.graph.add_edge(
            enzyme,
            target,
            effect=interaction.effect,
            severity=interaction.severity,
            source_reference=interaction.source_reference,
            enzyme=enzyme,
        )

    def add_many(self, interactions: Iterable[Interaction]) -> None:
        for interaction in interactions:
            self.add_interaction(interaction)

    def ripple_frontier(self, sources: Iterable[str]) -> dict[str, list[dict]]:
        frontier: dict[str, list[dict]] = {}
        for source in sources:
            if source not in self.graph:
                frontier[source] = []
                continue
            frontier[source] = []
            for node in nx.descendants(self.graph, source):
                path = nx.shortest_path(self.graph, source, node)
                for parent, child in zip(path, path[1:]):
                    edge_data = self.graph.get_edge_data(parent, child)
                    record = {"from": parent, "node": child, **edge_data}
                    if record not in frontier[source]:
                        frontier[source].append(record)
        return frontier
