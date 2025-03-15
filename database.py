import pandas as pd
from functools import lru_cache

NAME_A = "Systematic Name Interactor A"
NAME_B = "Systematic Name Interactor B"
EXPERIMENT = "Experimental System"
EXPERIMENT_TYPE = "Experimental System Type"
SCORE = "Score"

class Database:
    def __init__(self, path:str, sep:str='\t'):
        self.df = pd.read_csv(path, sep=sep)

        to_retain = [
            NAME_A,
            NAME_B,
            EXPERIMENT,
            EXPERIMENT_TYPE,
            SCORE
        ]

        self.df = self.df[to_retain]

    @lru_cache(maxsize=None)
    def filter_interactions(self, a: str, b: str) -> pd.DataFrame:
        # if a < b:
        #     a,b = b,a

        ret = self.df[self.df[NAME_A] == a]
        ret = ret[ret[NAME_B] == b]
        return ret

    @lru_cache(maxsize=None)
    def filter_functions(self, a: str, b: str) -> set[str]:
        rows = self.filter_interactions(a, b)
        ret = set(rows[EXPERIMENT].unique())
        return ret

