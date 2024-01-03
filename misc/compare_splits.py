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



def main():
    args = create_parser()
    splits = []
    for split in args.splits:
        df = pd.read_csv(split, sep="\t")
        taxa = list(df)
        taxa.remove("Splits")
        taxa.remove("Weights")
        df["size"] = df[taxa].sum(axis=1)
        df = df[(df["size"] > 1) & (df["size"] < len(taxa) - 1)]
        df = df[taxa]

        # create a canonical listing
        def extract_split(row):
            split0 = sorted(list(df.columns[row == 0]))
            split1 = sorted(list(df.columns[row == 1]))
            if (len(split0) < len(split1)):
                return frozenset(split0)
            elif(len(split0) > len(split1)):
                return frozenset(split1)
            else:
                for s0, s1 in zip(split0, split1):
                    if s0 < s1:
                        return frozenset(split0)
                    elif s0 > s1:
                        return frozenset(split1)            
            return frozenset(split0)
        splits.append((split, taxa, set(df.apply(extract_split, axis=1).to_list())))

    for a, b in itertools.combinations(splits, 2):
        a_diff_b = a[2] - b[2]
        b_diff_a = b[2] - a[2]
        a_intersect_b = a[2].intersection(b[2])
        print(f"a: {a[0]} vs b: {b[0]}")
        print(f"{len(a[1])} taxa in a, {len(b[1])} taxa in b")
        print(f"splits a \\ b: {len(a_diff_b)}")
        print(f"splits a \\ b: {a_diff_b}")
        print(f"splits b \\ a: {len(b_diff_a)}")
        print(f"splits b \\ a: {b_diff_a}")
        print(f"splits a ∩ b: {len(a_intersect_b)}")
        print(f"splits a ∩ b: {a_intersect_b}")

if __name__ == "__main__":
    main()