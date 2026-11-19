import networkx as nx
from typing import Iterable


class MetabolicPageRank:
    def rank(self, graph: nx.DiGraph, weight: str = "weight") -> dict[str, float]:
        if graph.number_of_nodes() == 0:
            return {}
        return nx.pagerank(graph, alpha=0.85, weight=weight)

    def bottlenecks(self, graph: nx.DiGraph, top_k: int = 3) -> list[tuple[str, float]]:
        scores = self.rank(graph)
        return sorted(scores.items(), key=lambda item: item[1], reverse=True)[:max(top_k, 0)]
