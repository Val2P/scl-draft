from graph import Graph

TARGET = "Collins"
trials = 4


def ppath(dataset, trialno = None):
    p = "./datasets/"
    if trialno is None:
        return p + dataset + ".txt"
    else:
        return p + dataset + f"v{trialno}.txt"

G = Graph(ppath(TARGET), database_path="./database-files/BIOGRID.txt")


for i in range(trials):
    run = i+1
    print("Starting run", run)
    G.set_depth(run)
    G.reweight_with_reliability(ppath(TARGET, run),True)
    print("Finished run", run)



"""
scratch
from graph import Graph


G = Graph("./dataset.txt")

print(G.N("YAL033W"))

G.visualize()
"""
