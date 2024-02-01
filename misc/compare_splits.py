import pandas as pd
import argparse
import numpy
import itertools

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-s', '--splits',
                   help="The list of files to compare the splits", nargs="+")

    return (p.parse_args())


class Split:
    def __init__(self, a: set, b: set, weight: float) -> None:
        self.a = a
        self.b = b
        self.weight = weight

    def get_canonical(self) -> frozenset: 
        if (len(self.a) < len(self.b)):
            return frozenset(self.a)
        elif(len(self.a) > len(self.b)):
            return frozenset(self.b)
        else:
            for s0, s1 in zip(self.a, self.b):
                if s0 < s1:
                    return frozenset(self.a)
                elif s0 > s1:
                    return frozenset(self.b)            
        return frozenset(self.a)
    
    def get_weight(self) -> float:
        return self.weight
    
    def is_in(self, other: set) -> bool:
        for o in other:
            if o.get_canonical() == self.get_canonical():
                return True
        return False
    
    """
    Print the split such that each split is a valid search term in SplitsTree6
    """
    def __str__(self) -> str:
        return f"{self.weight}: {'|'.join(self.get_canonical())}"
    
    def __repr__(self) -> str:
        return str(self)

    def __gt__(self, other) -> bool:
        return self.weight > other.weight
"""
    Calculate difference A - B while keeping the class structure
"""
def difference(A: set, B: set) -> set:
    result = set()
    B_canonical = set([b.get_canonical() for b in B])
    for a in A:
        if a.get_canonical() not in B_canonical:
            result.add(a)
    return result

"""
    Calculate the intersection A ∩ B while keeping the class structure
"""
def intersect(A: set, B: set) -> set:
    result = set()
    B_canonical = set([b.get_canonical() for b in B])
    for a in A:
        if a.get_canonical() in B_canonical:
            result.add(a)
    return result

def main():
    args = create_parser()
    splits = []
    for split in args.splits:
        df = pd.read_csv(split, sep="\t")
        taxa = list(df)
        taxa.remove("Splits")
        taxa.remove("Weights")
        df["size"] = df[taxa].sum(axis=1)
        # remove trivial splits
        df = df[(df["size"] > 1) & (df["size"] < len(taxa) - 1)]
        df_reduced = df[taxa]

        # create a canonical listing
        def extract_split(row):
            a = sorted(list(df.columns[(row == 0) & (df.columns.isin(taxa))]))
            b = sorted(list(df.columns[(row == 1) & (df.columns.isin(taxa))]))   
            return Split(a, b, row["Weights"])
        
        splits.append((split, taxa, set(df.apply(extract_split, axis=1).to_list())))

    for a, b in itertools.combinations(splits, 2):
        a_diff_b = sorted(list(difference(a[2], b[2])))
        b_diff_a = sorted(list(difference(b[2], a[2])))
        new_line = "\n"
        print(f"a: {a[0]} vs b: {b[0]}")
        print(f"{len(a[1])} taxa in a, {len(b[1])} taxa in b")
        print(f"total weight in a: {sum([s.weight for s in a[2]])}")
        print(f"total weight in b: {sum([s.weight for s in b[2]])}")
        print(f"splits a \\ b: {len(a_diff_b)}")
        print(f"total weight of difference: {sum([s.weight for s in a_diff_b])}")
        print(f"splits a \\ b:\n{new_line.join([str(e) for e in a_diff_b])}")
        print()
        print(f"splits b \\ a: {len(b_diff_a)}")
        print(f"total weight of difference: {sum([s.weight for s in b_diff_a])}")
        print(f"splits b \\ a:\n{new_line.join([str(e) for e in b_diff_a])}")
        print()
        print()

if __name__ == "__main__":
    main()