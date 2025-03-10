import pandas as pd

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

    def filter_interactions(self, a: str, b: str) -> pd.DataFrame:
        # if a < b:
        #     a,b = b,a

        return self.df[(self.df[NAME_A] == a) & (self.df[NAME_B == b])]
