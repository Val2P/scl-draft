import networkx as nx
from pyvis.network import Network

# assumption
# supplied file has a format ff:
"""
(n0): p0 p1 p2 ...
(n1): q0 q1 q2 ...
"""


def visualize_file(x: str) -> None:
    clusters = []
    with open(x, "r") as f:
        for line in f.readlines():
            cluster = line.strip().split()[1:]
            clusters.append(cluster)

    G = nx.Graph()

    for cluster in clusters:
        for i in range(len(cluster)):
            for j in range(i + 1, len(cluster)):
                G.add_edge(cluster[i], cluster[j])

    n = Network()
    n.from_nx(G)
    n.show("index.html", notebook=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize clusters produced by p5comp algo"
    )
    parser.add_argument("file", help="path for the cluster file")
    args = parser.parse_args()

    visualize_file(args.file)
