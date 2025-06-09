from pyvis.network import Network
import networkx as nx
import pandas as pd
from tqdm import tqdm

from database import Database

from functools import lru_cache
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

EXP = "Experimental System"


class Graph:
    def __init__(
        self,
        dataset_path: str,
        initial_depth: int = 0,
        database_path: str | None = None,
        sep: str = "\t",
    ):
        edges = list()
        nodes = set()

        logger.info(f"Reading Dataset: {dataset_path}")

        with open(dataset_path, "r") as f:
            for line in f.readlines():
                l = line.strip()
                x = l.split("\t")
                if len(x) != 3:
                    x = l.split(" ")
                n1, n2, edge = x

                nodes.add(n1)
                nodes.add(n2)
                edges.append((n1, n2, {"weight": float(edge), "label": edge}))

        logger.info(f"Finished loading Dataset")

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

        self.database: Database | None = None

        self._name = dataset_path

        if database_path is not None:
            logger.info(f"Loading Database: {database_path}")
            self.database = Database(database_path, sep)
            logger.info(f"Finished Loading Database")

            self.reliability_db = dict()
            n = self.database.df.shape[0]
            experiments = self.database.df[EXP]

            for e in tqdm(
                experiments.unique(),
                "Computing reliability of experiments based on database",
            ):
                self.reliability_db[e] = (experiments == e).sum() / n

            self.database.trim_rows(nodes)

            # r_int
            _r = 1
            for u, v, _ in tqdm(
                self.edges, "Computing Fraction of pairs that share functions"
            ):
                Fx = self.database.filter_functions(u, v)
                Fy = self.database.filter_functions(v, u)
                if len(Fx & Fy) > 0:
                    _r += 1

            self.r_int = float(_r) / float(len(self.edges))
            """
            self.r_int = 0.45228124311218865
            """

            logger.info(f"{self.r_int = }")

    def set_depth(self, n: int):
        assert isinstance(n, int)
        assert n >= 0
        self.depth = n

        # self.N_depth.cache_clear()
        self.N_depth_recursive.cache_clear()
        self.N_depth_nx.cache_clear()
        self.Functions.cache_clear()
        self.r1.cache_clear()
        # Cache.clear()

    def visualize(self):
        _nt = Network()
        _nt.from_nx(self.network)

        _nt.show("index.html", notebook=False)

    def edgelist_to_file(self, l: list, save_path: str | None):
        def mymap(x):
            return f"{x[0]}\t{x[1]}\t{x[2]:.16f}"

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

    # @lru_cache(maxsize=None)
    def N_depth(self, u: str) -> set[str]:
        """
        Returns the names of the proteins w/ at most `self.depth` levels away from u
        """
        k = self.depth
        assert k > -1

        ret = {u}

        just_added = {u}

        # wont run on k = 0, for range(1,1) doesnt run
        for _ in range(1, k + 1):
            new_neighbors = set()
            for n in list(just_added):
                new_neighbors = new_neighbors | self.N(n)

            just_added = new_neighbors - ret
            ret = ret | just_added

        return ret

    # recursive impl with caching for N_depth
    @lru_cache(maxsize=None)
    def N_depth_diff(self, u: str, d: int) -> set[str]:
        if d < 0:
            return set()
        if d == 0:
            return {u}
        else:
            diff = self.N_depth_diff(u, d - 1) - self.N_depth_diff(u, d - 2)
            ret = set()
            for n in diff:
                ret = ret | self.N(n)

            return ret

    @lru_cache(maxsize=None)
    def N_depth_recursive(self, u: str) -> set[str]:
        ret = set()
        for i in range(0, self.depth + 1):
            ret = ret | self.N_depth_diff(u, i)
        return ret

    # Using nx native functions
    @lru_cache(maxsize=None)
    def N_depth_nx(self, u: str):
        return set(nx.ego_graph(self.network, u, radius=self.depth).nodes)

    ########################
    # Functional Reliability
    ########################
    @lru_cache(maxsize=None)
    def Functions(self, u: str) -> set[str]:
        if self.database is None:
            return set()

        # neighbors = self.N_depth(u) - {u}
        neighbors = self.N_depth_recursive(u) - {u}
        # neighbors = self.N_depth_nx(u) - {u}
        ret: set[str] = set()

        for n in neighbors:
            # set union
            ret = ret | self.database.filter_functions(u, n)

        return ret

    @lru_cache(maxsize=None)
    def r1(self, a: str, b: str):
        if self.database is None:
            return 1

        Fa = self.Functions(a)
        Fb = self.Functions(b)

        a_and_b = Fa & Fb

        # denominator is zero
        if len(a_and_b) == 0:
            return 0

        a_or_b = Fa | Fb

        ret = len(a_or_b - a_and_b) / len(a_and_b)

        logger.debug(f"Computed r1 of {a} -> {b}: {ret}")
        return ret

    #####################
    # Functions that involve FS Weight
    ####################
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

        return float(numerator / deno_1) * float(numerator / deno_2)

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

        return float(numerator / deno_1) * float(numerator / deno_2)

    def reweight_chua(
        self, save_path: str, verbose: bool = False, save_cache: bool = False
    ) -> None:
        """
        chua FS weighting algo
        """
        assert self.database is not None

        logger.info(f"Starting reweighting using Chua's algo with depth = {self.depth}")

        N = self.N_depth_nx
        r_int = self.r_int

        @lru_cache(maxsize=None)
        def r(u: str, v: str) -> float:
            exps = self.database.filter_interactions(u, v)[EXP]
            n = 1.0
            for e in exps.unique():
                n *= (1 - self.reliability_db[e]) ** (exps == e).sum()

            return 1.0 - n

        @lru_cache(maxsize=None)
        def SR_term(u: str, v: str) -> float:
            Nu = N(u)
            Nv = N(v)
            Nu_and_Nv = Nu & Nv

            _term_uv = self._n_avg * r_int - (len(Nu - Nv) + len(Nu & Nv))

            lmbda_uv = max(0, _term_uv)

            numerator = 0.0
            for w in Nu_and_Nv:
                numerator += r(u, w) * r(v, w)
            numerator *= 2.0

            denominator_a = 0.0
            for w in Nu:
                denominator_a += r(u, w)

            for w in Nu_and_Nv:
                denominator_a += r(u, w) * (1 - r(v, w))

            denominator_b = 0.0
            for w in Nu_and_Nv:
                denominator_b += r(u, w) * r(v, w)

            denominator = denominator_a + (2 * denominator_b) + lmbda_uv

            try:
                return float(numerator) / float(denominator)
            except:
                return 0.0

        @lru_cache(maxsize=None)
        def SR_edge(u, v) -> tuple:
            weight = SR_term(u, v) * SR_term(v, u)
            ret = (u, v, weight)

            logger.debug(f"Ran on {u = } | {v = } ; {weight = }")
            return ret

        # caching
        cached_graph = []
        cache_filename = save_path + ".cache"
        if save_cache:
            logger.info("Loading cache")
            with open(cache_filename, "r") as f:
                for line in f.readlines():
                    l = line.strip()
                    x = l.split("\t")
                    if len(x) != 3:
                        x = l.split(" ")
                    n1, n2, edge = x
                    cached_graph.append((n1, n2, edge))
            logger.info("Finished loading cache")

        skip = len(cached_graph)
        new_graph = cached_graph

        for data in tqdm(self.edges, "Running reweighting on graph"):
            if skip > 0:
                skip -= 1
                continue
            else:
                computed_edge = SR_edge(data[0], data[1])
                new_graph.append(computed_edge)

                if save_cache:
                    with open(cache_filename, "a") as f:
                        f.write(
                            f"{computed_edge[0]}\t{computed_edge[1]}\t{computed_edge[2]}\n"
                        )

        r.cache_clear()
        SR_term.cache_clear()
        SR_edge.cache_clear()

        self.edgelist_to_file(new_graph, save_path)
