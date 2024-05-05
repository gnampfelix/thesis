import pandas as pd
import argparse
import json

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-i', '--input',
                   help="Path to the json exports of the mash sketches", nargs="+")

    return (p.parse_args())

def main():
    sketch_size=10000
    args = create_parser()
    sketches = []
    for path in args.input:
        with open(path) as f:
            sketch = json.load(f)
            sequence_name = sketch["sketches"][0]["name"].split("/")[-1]
            hashes = set(sketch["sketches"][0]["hashes"])
            sketches.append((sequence_name, hashes))

    for i in range(len(sketches)):
        for j in range(i, len(sketches)):
            name_a, a = sketches[i]
            name_b, b = sketches[j]
            union = a.union(b)
            s_union = set(sorted(list(union))[:sketch_size])
            interesction = s_union.intersection(a).intersection(b)
            if len(interesction) == 0:
                print(f"{name_a} vs {name_b} is empty")


if __name__ == "__main__":
    main()
