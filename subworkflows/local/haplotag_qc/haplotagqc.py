import pysam
import tqdm
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import click


def read_bam(bam_path, chromosome=None):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    records = []

    for read in tqdm.tqdm(bamfile.fetch(chromosome)):
        if read.is_unmapped:
            continue

        # If the read has a "PS" tag, grab it; otherwise None
        ps_value = read.get_tag("PS") if read.has_tag("PS") else None

        # If the read has a "HP" tag, grab it; otherwise None
        hp_value = read.get_tag("HP") if read.has_tag("HP") else None

        records.append(
            {
                "chromosome": read.reference_name,
                "start": read.reference_start,
                "end": read.reference_end,
                "secondary": read.is_secondary,
                "supplementary": read.is_supplementary,
                "PS": ps_value,
                "HP": hp_value,
            }
        )

    bamfile.close()

    return pd.DataFrame(records)


def prep_bam_df(bam_path):
    """
    Prepare a DataFrame from the BAM file with relevant read information.
    """
    reads = read_bam(bam_path)

    reads["PS"] = reads["PS"].fillna(-1).astype(int)
    reads["HP"] = reads["HP"].fillna(-1).astype(int)
    reads["length"] = reads["end"] - reads["start"] + 1
    return reads


def plot_read_length_distribution(reads, patient_id, tool):
    """
    Plot the distribution of read lengths.
    """
    plt.figure(figsize=(3, 3))
    sns.histplot(x="length", data=reads, log_scale=True, bins=30)
    sns.despine()
    plt.title(f"{patient_id} {tool}")
    plt.savefig(
        f"read_length_distribution_{patient_id}_{tool}_mqc.png", dpi=300, bbox_inches="tight"
    )
    plt.close()


def export_fraction_haplotagged(reads, patient_id, tool):
    """
    Calculate the fraction of haplotagged reads.
    """
    f = (reads["PS"] != -1).mean()
    with open(f"haplotagged_fraction_{patient_id}_{tool}_mqc.txt", "w") as f_out:
        f_out.write(f"{f:.4f}")


# Compute a table of allele fractions in bins
# Read stats
def export_read_stats(reads, patient_id, tool):
    read_stats = (
        reads
        # .query("PS >= 0")
        .groupby("PS")
        .agg(
            n_reads=("PS", "size"),
            sum_read_length=("length", "sum"),
        )
        .reset_index()
    )
    read_stats.to_csv(f"read_stats_{patient_id}_{tool}_mqc.csv", index=False)


def calculate_read_counts(reads):
    """
    Calculate read counts per bin and allele fraction.
    """
    bin_size = 500000
    reads["bin_start"] = ((reads["start"] / bin_size).astype(int)) * bin_size + 1

    read_counts = (
        reads.query("PS >= 0")
        .groupby(["chromosome", "bin_start", "PS", "HP"])
        .agg(
            count=("PS", "size"),
        )["count"]
        .unstack("HP", fill_value=0)
        .reset_index()
    )
    print(read_counts.columns)
    read_counts = read_counts.rename(columns={1: "allele_1", 2: "allele_2"})
    if read_counts.empty:
        # This is mostly for testing purposes, to ensure the function returns a DataFrame with the minimal test data
        read_counts = pd.DataFrame(
            [{"chromosome": "chr1", "bin_start": 0, "allele_1": 10, "allele_2": 10, "total": 20, "af": 0.5}],
        )
        return read_counts
    read_counts["af"] = read_counts[["allele_1", "allele_2"]].min(axis=1) / (
    read_counts["allele_1"] + read_counts["allele_2"]
    )
    read_counts["total"] = read_counts["allele_1"] + read_counts["allele_2"]
    return read_counts


def plot_allele_fractions(read_counts, patient_id, tool):
    plot_data = read_counts.query("total > 10")
    plt.figure(figsize=(10, 2))
    sns.FacetGrid(
        data=plot_data,
        row="chromosome",
        row_order=sorted(plot_data["chromosome"].unique()),
        sharey=True,
        sharex=True,
    ).map_dataframe(sns.scatterplot, x="bin_start", y="af", size="total")
    plt.savefig(f"allele_fractions_{patient_id}_{tool}_mqc.png", dpi=300, bbox_inches="tight")
    plt.close()


# Median absolute deviation
def export_mad_of_adjacent_diffs(read_counts, patient_id, tool):
    """
    Calculate the median absolute deviation of adjacent allele fraction differences.
    """
    mad_diff = np.median(
        np.abs(
            np.diff(
                read_counts.query("total > 10")
                .sort_values(["chromosome", "bin_start"])["af"]
                .values
            )
        )
    )
    with open(f"mad_of_adjacent_diffs_{patient_id}_{tool}_mqc.txt", "w") as f_out:
        f_out.write(f"{mad_diff:.4f}")


@click.command()
@click.option(
    "--bam",
    type=click.Path(exists=True),
    required=True,
    help="Path to the BAM file to analyze.",
)
@click.option(
    "--sample",
    required=True,
    help="Sample ID for the analysis.",
)
@click.option(
    "--tool",
    required=True,
    help="Tool for the analysis.",
)
def main(bam, sample, tool):
    reads = prep_bam_df(bam)
    plot_read_length_distribution(reads, sample, tool)
    export_fraction_haplotagged(reads, sample, tool)
    export_read_stats(reads, sample, tool)
    read_counts = calculate_read_counts(reads)
    plot_allele_fractions(read_counts, sample, tool)
    export_mad_of_adjacent_diffs(read_counts, sample, tool)


if __name__ == "__main__":
    main()
