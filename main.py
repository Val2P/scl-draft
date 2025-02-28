from graph import Graph

TARGET = "KroganCore"
trials = 5


def ppath(dataset, trialno = None):
    p = "./datasets/"
    if trialno is None:
        return p + dataset + ".txt"
    else:
        return p + dataset + f"v{trialno}.txt"
PATH = "./datasets"
G = Graph(ppath(TARGET))


for i in range(trials):
    run = i+1
    print("Starting run", run)
    x = lambda u,v: G.Depth_FS(u, v, run)
    G.reweight2(x, ppath(TARGET, run),True)
    print("Finished run", run)



"""
scratch
from graph import Graph


G = Graph("./dataset.txt")

print(G.N("YAL033W"))

G.visualize()
"""
