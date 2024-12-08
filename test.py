from graph import Graph


G = Graph("./dataset.txt", "./database/biogrid.txt")
G.reweight(True, "./output.txt")
