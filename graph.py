from pyvis.network import Network
import networkx as nx
import pandas as pd

from database import Database

from functools import lru_cache


"""
TODO, FSWeight w/ Depth + Reliability
1. compute the Neighbor Sets up to depth k
2. compute the r1(a,b) for each edge, incorporating Neighbor Sets
3. set up LUT for R1 (normalized r1)
4. do Chua FS Weight
"""

class Graph:
    def __init__(self, dataset_path: str, initial_depth:int = 0, database_path: str|None = None, sep:str = '\t'):

        edges = list()
        nodes = set()

        with open(dataset_path, 'r') as f:
            for line in f.readlines():
                l = line.strip()
                x = l.split('\t')
                if len(x) != 3:
                    x = l.split(" ")
                n1, n2, edge = x

                nodes.add(n1)
                nodes.add(n2)
                edges.append((n1, n2, {"weight": float(edge), "label": edge}))


        self.edges = edges
        self.nodes = nodes
        self.network = nx.Graph()


        # getting neighbors
        # G.adj[n].items()

        self.network.add_nodes_from(nodes)
        self.network.add_edges_from(edges)

        assert initial_depth >= 0
        self.depth = initial_depth


        # reliability computation
        self._n_avg = len(edges) / len(nodes)


        self.database: Database|None = None

        if database_path is not None:
            self.database = Database(dataset_path, sep)


    def set_depth(self, n: int):
        assert isinstance(n, int)
        assert n >= 0
        self.depth = n

        self.N_depth.cache_clear()
        self.Functions.cache_clear()
        self.r1.cache_clear()


    def visualize(self):
        _nt = Network()
        _nt.from_nx(self.network)


        _nt.show('index.html', notebook=False)


    def edgelist_to_file(self, l: list, save_path: str|None):

        def mymap(x):
            return f"{x[0]}\t{x[1]}\t{x[2]}"

        if save_path is not None:
            with open(save_path, "w") as f:
                r = map(mymap, l)
                f.write("\n".join(r))

    #####################
    # Functions that involve the neighbor counting
    ####################

    @lru_cache(maxsize=None)
    def N(self, u: str) -> set[str]:
        """
        Returns the names of the proteins adjacent to protein u
        """

        ret = {u}
        for nbr, _ in self.network.adj[u].items():
            ret.add(nbr)


        return ret

    @lru_cache(maxsize=None)
    def N_depth(self, u: str) -> set[str]:
        """
        Returns the names of the proteins w/ at most `self.depth` levels away from u
        """
        k = self.depth
        assert k > -1


        ret = {u}

        just_added = {u}
        
        # wont run on k = 0, for range(1,1) doesnt run
        for _ in range(1, k+1):
            new_neighbors = set()
            for n in list(just_added):
                new_neighbors  = new_neighbors | self.N(n)

            just_added = new_neighbors - ret
            ret =  ret | just_added


        return ret


    ########################
    # Functional Reliability
    ########################
    @lru_cache(maxsize=None)
    def Functions(self, u: str) -> set[str]:

        if self.database is None:
            return set()

        neighbors = self.N_depth(u) - {u}
        ret: set[str] = set()

        for n in neighbors:
            # set union
            ret = ret | self.database.filter_functions(u, n)

        return ret

    @lru_cache(maxsize=None)
    def r1(self, a: str, b: str):
        if self.database is None: return 1

        Fa = self.Functions(a)
        Fb = self.Functions(b)

        a_and_b = Fa & Fb

        # denominator is zero
        if len(a_and_b) == 0:
            return 0

        a_or_b = Fa | Fb

        return len(a_or_b - a_and_b) / len(a_and_b)



    #####################
    # Functions that involve FS Weight
    ####################

    def lmbda(self, u: set[str], v: set[str]) -> float:
        """
        lambda function, Bioinformatics equation 2
        TODO modify further to include experimental sources

        lambda = max(0, n_avg * r_int - ...)
        were r_int is the fraction of all interaction pairs that share some function
        """

        diff = u - v
        inter = u & v

        return max(0, self._n_avg - (len(diff) + len(inter))  )




    def S_FS(self, u: str, v: str) -> float:
        """
        Simple FS Weight function
        """

        N = self.N

        n_u = N(u)
        n_v = N(v)


        u_and_v = n_u & n_v
        u_min_v = n_u - n_v
        v_min_u = n_v - n_u




        numerator = 2 * len(u_and_v)

        deno_1 = len(u_min_v) + numerator + self.lmbda(n_u, n_v)
        deno_2 = len(v_min_u) + numerator + self.lmbda(n_v, n_u)


        return float(numerator/deno_1) * float(numerator/deno_2)



    def Depth_FS(self, u: str, v: str) -> float:
        """
        FS Weight function with Depth
        `N` is the neighbor counting function
        neighbor counting function is a func that accepts a string and returns a set of neighbor
        """
        d = self.depth

        N = lambda x: self.N_depth(x, d)


        n_u = N(u)
        n_v = N(v)


        u_and_v = n_u & n_v
        u_min_v = n_u - n_v
        v_min_u = n_v - n_u




        numerator = 2 * len(u_and_v)

        deno_1 = len(u_min_v) + numerator + self.lmbda(n_u, n_v)
        deno_2 = len(v_min_u) + numerator + self.lmbda(n_v, n_u)


        return float(numerator/deno_1) * float(numerator/deno_2)





    def reweight(self, func, save_path = None, verbose = False) -> None:
        '''
        reweights graph based on function
        '''

        # func will return an edge weight based on its supplied vertices
        # store the edge-weights
        # put it in a list for saving jic
        new_graph = []

        l = len(self.network.edges)
        i = 1
        if verbose:
            print(f"Running on {l} edges")

        for edge in self.network.edges:
            u = edge[0]
            v = edge[1]

            if verbose:
                print(f"{i}: running func on: `{u}` and `{v}`; ", end='')


            weight = func(u,v)

            if verbose:
                print(f"{weight = }")

            new_graph.append((u,v,weight))

            i += 1


        self.edgelist_to_file(new_graph, save_path)


    def reweight_with_reliability(self, save_path:str|None=None):

        assert self.database is not None

        new_graph = []

        # if greater than this value, considered as infinity
        EPSILON = 100000

        min_r1 = EPSILON
        max_r1 = 0

        # calc each r1 for normalization
        # the costly subroutines *should* be cached on computing r1
        for src, dst, _ in self.edges:
            x = self.r1(src, dst)
            if x > EPSILON: continue

            min_r1 = min(min_r1, x)
            max_r1 = max(max_r1, x)


        def R1(a: str, b: str) -> float|int:
            x = self.r1(a, b)
            if x > EPSILON: return 1

            numerator = x - min_r1
            denominator = max_r1 - min_r1
            return numerator / denominator





        self.edgelist_to_file(new_graph, save_path)
