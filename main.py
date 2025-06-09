# from graph import Graph

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="SCL 2024-2025 Thesis",
        description="PPIN Reweighter using Topology and Experimenal Source",
    )

    parser.add_argument(
        "--input_graph", "-I", required=True, help="path of input graph"
    )
    parser.add_argument(
        "--database_path", "-D", required=True, help="path to BioGRID database"
    )
    parser.add_argument(
        "--depth_from",
        "-df",
        default=0,
        type=int,
        help="include starting from to depth n (must be geq than 0)",
    )
    parser.add_argument(
        "--depth_to",
        "-dt",
        default=0,
        type=int,
        help="include up to depth n (must be geq than 0)",
    )
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--path", "-p", help="path to save reweighted graph")
    parser.add_argument(
        "--cache", "-c", action="store_true", help="cache computations?"
    )

    args = parser.parse_args()

    from graph import Graph

    G = Graph(dataset_path=args.input_graph, database_path=args.database_path)

    for depth in range(args.depth_from, args.depth_to + 1):
        path = args.path
        if path:
            inp = args.input_graph
            sep = ""
            if "\\\\" in inp:
                sep = "\\\\"
            elif "\\" in inp:
                sep = "\\"
            else:
                sep = "/"

            file = inp.split(sep)[-1].split(".")
            name = file[0]
            ext = file[1]

            path += f"{sep}{name}v{depth}.{ext}"

        print(f"Running Algoritm on {args.input_graph} with {depth = }")
        G.set_depth(depth)
        G.reweight_chua(path, verbose=args.verbose, save_cache=args.cache)
        print("done")
