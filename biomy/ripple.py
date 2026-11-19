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
        self.graph.add_edge(
            interaction.source,
            interaction.target,
            enzyme=interaction.enzyme,
            effect=interaction.effect,
            severity=interaction.severity,
            source_reference=interaction.source_reference,
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
