"""
unused stuff/scratch stuff
"""


class K2Table:
    def __init__(self):
        self._t = dict()

    def _key(self, a1, a2):
        k1 = a1
        k2 = a2

        # comment this out if keys arent symmetric
        # ie: table[a,b] != table[b,a]
        if k1 > k2:
            k1, k2 = k2, k1

        k = (k1, k2)

        return k



    def __getitem__(self, args):

        k = self._key(*args)
        
        if k in self._t.keys():
            return self._t[k]
        else:
            return None

    def __setitem__(self, keys, value):
        k = self._key(*keys)

        self._t[k] = value







    # def _reliability(self, u: str, v: str) -> float:
    #     """
    #     NOTE: assuming interactions are symmetric, that is, A interacts with B is the same as B interacts with A
    #     NOTE: assuming same "weight" of High Throughput and Low Throughput
    #     """
    #
    #     if self.df is not None:
    #         ret = self._mem[u, v]
    #
    #         if ret is None:
    #             col = "Systematic Name Interactor A"
    #             u_total = self.df[self.df[col] == u]
    #             v_total = self.df[self.df[col] == v]
    #             total = u_total.shape[0] + v_total.shape[0]
    #
    #             col = "Systematic Name Interactor B"
    #             uv_total = u_total[self.df[col] == v].shape[0] + v_total[self.df[col] == u].shape[0]
    #
    #             ret = 1.0 - (1.0 - uv_total/total)**uv_total
    #             self._mem[u,v] = ret
    #
    #
    #         return ret
    #     else:
    #         return 1
    #
    #
    # def S_R(self, u: str, v: str) -> float:
    #     n_u = self.N(u)
    #     n_v = self.N(v)
    #
    #     u_and_v = n_u & n_v
    #
    #     # fraction one
    #     n = 0.0
    #     for w in u_and_v:
    #         n += self._reliability(u, w) * self._reliability(v, w)
    #     n *= 2
    #
    #
    #     deno = 0.0
    #
    #     for w in n_u:
    #         deno += self._reliability(u, w)
    #
    #     for w in u_and_v:
    #         deno += self._reliability(u, w) * (1 - self._reliability(v, w))
    #
    #     for w in u_and_v:
    #         deno += 2 * self._reliability(u, w) * self._reliability(v, w)
    #
    #     deno += self.lmbda(n_u, n_v)
    #
    #     ret = float(n/deno)
    #
    #
    #     # fraction two
    #     deno = 0.0
    #
    #     for w in n_v:
    #         deno += self._reliability(v, w)
    #
    #     for w in u_and_v:
    #         deno += self._reliability(v, w) * (1 - self._reliability(u, w))
    #
    #     for w in u_and_v:
    #         deno += 2 * self._reliability(u, w) * self._reliability(v, w)
    #
    #     deno += self.lmbda(n_v, n_u)
    #
    #
    #     return ret * float(n/deno)
    #
    #
    # def S_TR(self, u: str, v: str) -> float:
    #
    #     ret = self.S_R(u, v)
    #
    #     n_u = self.N(u)
    #     for w in n_u:
    #         ret = max(ret, self.S_R(u, w) * self.S_R(v, w))
    #
    #
    #     return ret
    #
    #
    # def Z(self) -> float:
    #     # Z = 1 + Sum_(v in N_u) [S_TR(u,v) + Sum_(w in N_v) [S_TR(u,w)]]
    #     return 1
    #
    #
    # def f_x(self, u: str) -> float:
    #     # f_x(u) = [1/Z] * [lambda * r_int * pi_x + 
    #     #                   Sum_(v in N_u) [S_TR(u,v) * delta(v,x)] +
    #     #                   Sum_(w in N_v) [S_TR(u,w) * delta(w,x)]]
    #     # lambda [0,1]: weight representing the contribution of background frequency to the score
    #     # r_int: fraction of all interaction pairs that share some function defined in S_R
    #     # delta(p,x): 1 if p has function x, 0 otherwise
    #     # pi_x: frequency of function x in annotated proteins
    #     return 1
    #
    #
    # def reweight(self, verbose=False, output_path=None) -> None:
    #     edges = self.edges
    #     ret = []
    #
    #     for i in range(len(edges)):
    #         edge = edges[i]
    #
    #         if verbose: 
    #             print(f"Edge {i+1}/{len(edges)}")
    #             print("Computing Weight")
    #
    #         r = self.S_R(edge[0], edge[1])
    #         ret.append((edge[0], edge[1], r))
    #
    #         if verbose: 
    #             print("Done computing")
    #             print("-"*10)
    #
    #     if output_path is not None:
    #
    #         if verbose:
    #             print(f"Writing output to \'{output_path}\'")
    #
    #         with open(output_path, "a") as f:
    #             for v in ret:
    #                 _v = map(str, v)
    #                 line = '\t'.join(_v)
    #                 f.write(f"{line}\n")
    #
    #     else:
    #         if verbose:
    #             print("no path specified")
    #
    #     if verbose:
    #         print("Done")
