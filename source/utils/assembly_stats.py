import pandas as pd, sys
from Bio.SeqIO import parse
from argparse import ArgumentParser


def store_lengths(files, names):
    r = {}
    for i, f in enumerate(files, start=0):
        r[names[i]] = {}
        for record in parse(f, "fasta"):
            r[names[i]][record.id] = len(record.seq)
    return r


def size_distribute(r,
                    lengths=[0, 100, 250, 500, 1000, 2500, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000,
                             50000, 75000, 100000, 125000, 150000, 200000, 250000, 500000]):
    size_dist = {}
    size_dist_df = pd.DataFrame(columns=["%", "min_length", "num_contigs", "total_length", "assembly"])
    for name, contig_lengths in r.items():
        df = pd.DataFrame(contig_lengths, index=["length"]).T
        size_dist[name] = {}
        for i, l in enumerate(lengths):
            if len(df.loc[df.length >= l]) == 0:
                break
            n, s, p = len(df.loc[df.length >= l]), int(df.loc[df.length >= l].sum()), int(
                df.loc[df.length >= l].sum()) / float(df.sum()) * 100
            size_dist[name][i] = {"min_length": l, "num_contigs": n, "total_length": s, "%": p}
        _ = pd.DataFrame(size_dist[name]).T
        _ = _.assign(assembly=pd.Series([name]*len(_), index=_.index))
        _ = _[size_dist_df.columns]
        size_dist_df = pd.concat([size_dist_df, _], sort=True)
    return size_dist_df


def calculate_n_stats(df):
    df.sort_values("length", inplace=True, ascending=True)
    size = int(df.sum())
    N50_length = N90_length = 0
    cumulative = 0
    for contig in df.index:
        l = df.loc[contig, "length"]
        cumulative += l
        if float(cumulative) >= 0.5 * size and not N50_length:
            N50_length = l
        elif float(cumulative) >= 0.1 * size and not N90_length:
            N90_length = l
    return N50_length, N90_length


def calculate_length_stats(contig_lengths):
    df = pd.DataFrame(contig_lengths, index=["length"]).T
    N50_length, N90_length = calculate_n_stats(df)
    stats = {"contigs": len(df), "total_size_bp": int(df.sum()), "min_length": int(df["length"].min()),
             "max_length": int(df["length"].max()), "avg_length": float(df["length"].mean()),
             "median_length": float(df["length"].median()), "N50_length": N50_length, "N90_length": N90_length}
    return stats


def generate_stat_df(r):
    index = ["contigs", "total_size_bp", "min_length", "max_length", "avg_length", "median_length",
             "N50_length", "N90_length"]
    stats = {}
    for name, contig_lengths in r.items():
        stats[name] = calculate_length_stats(contig_lengths)
    stat_df = pd.DataFrame(stats).T[index]
    return stat_df


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, nargs="+", help="Fasta file of contigs", required=True)
    parser.add_argument("-n", "--names", type=str, nargs="+", help="Assembly names", required=True)
    parser.add_argument("--size-dist-file", type=str,
                        help="Write table of size distributions for different contig lengths to file")
    parser.add_argument("--stat-file", type=str,
                        help="Write table of general statistics (size, lengths etc) to file")
    args = parser.parse_args()
    assert (len(args.names) == len(args.infile))

    r = store_lengths(args.infile, args.names)
    stat_df = generate_stat_df(r)

    if not args.stat_file:
        fh_out = sys.stdout
    else:
        fh_out = open(args.stat_file, 'w')
    stat_df.to_csv(fh_out, sep="\t", index=True)

    if args.size_dist_file:
        size_dist_df = size_distribute(r)
        size_dist_df.to_csv(args.size_dist_file, sep="\t", index=False)


if __name__ == '__main__':
    main()
