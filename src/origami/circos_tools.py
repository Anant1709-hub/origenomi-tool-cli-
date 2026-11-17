#!/usr/bin/env python3

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


def fast_gc_skew(seq, win=1000):
    """Calculate GC skew in non-overlapping windows."""
    seq = np.frombuffer(seq.upper().encode(), dtype="S1")
    g = (seq == b"G").astype(int)
    c = (seq == b"C").astype(int)

    bins = (len(seq) + win - 1) // win
    skew = np.zeros(bins, dtype=float)

    for i in range(bins):
        start = i * win
        end = min(start + win, len(seq))
        g_win = int(np.sum(g[start:end]))
        c_win = int(np.sum(c[start:end]))
        denom = g_win + c_win
        skew[i] = (g_win - c_win) / denom if denom else 0.0

    return skew


def _find_record(fna_path, accession):
    """Return the FASTA record whose ID exactly matches accession."""
    for record in SeqIO.parse(fna_path, "fasta"):
        if record.id == accession:
            return record
    raise ValueError(f"Accession {accession} not found in FASTA: {fna_path}")


def _normalise_origin(oriC_list):
    """Return the first oriC start coordinate, if an oriC hit exists."""
    if not oriC_list:
        return None
    first = oriC_list[0]
    if isinstance(first, (tuple, list)):
        return int(first[0])
    return int(first)


def _validate_coordinate(name, coord, genome_len, accession):
    if coord is not None and not (0 <= int(coord) <= genome_len):
        raise ValueError(
            f"{name} coordinate {coord} for {accession} is outside "
            f"the sequence length ({genome_len} bp)"
        )


def generate_circos(
    record,
    genome_len,
    skew,
    replicon_type,
    origin,
    origin_name,
    output_dir,
    dnaA=None,
    rotated=False,
    shift_bp=0,
):
    """Create one original or rotated circular GC-skew plot."""
    if not len(skew):
        raise ValueError(f"Cannot calculate GC skew for {record.id}: empty sequence")

    if rotated:
        bin_size = genome_len / len(skew)
        shift_bins = int(round(shift_bp / bin_size))
        skew = np.roll(skew, -shift_bins)

    bins = len(skew)
    angles = np.linspace(0, 2 * np.pi, bins, endpoint=False)

    fig, ax = plt.subplots(
        figsize=(10, 10),
        subplot_kw={"projection": "polar"},
    )

    try:
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi / 2)
        ax.set_axis_off()

        outer_r = 1.0
        ax.plot(
            np.linspace(0, 2 * np.pi, 2000),
            [outer_r] * 2000,
            color="navy",
            lw=4,
        )

        bar_r = 0.80
        bar_scale = 0.25
        positive = skew > 0
        negative = skew < 0

        ax.vlines(
            angles[positive],
            bar_r,
            bar_r + skew[positive] * bar_scale,
            color="lightblue",
            linewidth=0.45,
            alpha=0.9,
        )
        ax.vlines(
            angles[negative],
            bar_r - np.abs(skew[negative]) * bar_scale,
            bar_r,
            color="orange",
            linewidth=0.45,
            alpha=0.9,
        )

        num_ticks = 28
        for i in range(num_ticks):
            pos = genome_len * i / num_ticks
            theta = 2 * np.pi * i / num_ticks

            ax.plot(
                [theta, theta],
                [outer_r, outer_r + 0.03],
                color="navy",
                lw=1.5,
            )
            ax.text(
                theta,
                outer_r + 0.10,
                f"{pos / 1e6:.2f} Mb",
                ha="center",
                va="center",
                fontsize=9,
            )

        def get_angle(coord):
            if rotated:
                coord = (coord - shift_bp) % genome_len
            return 2 * np.pi * (coord / genome_len)

        # dnaA is a chromosome marker only.
        if replicon_type == "chromosome" and dnaA is not None:
            ax.plot(
                [get_angle(dnaA)],
                [outer_r],
                marker="s",
                color="blue",
                markersize=10,
            )

        # oriC for chromosomes; oriV for plasmids.
        ax.plot(
            [get_angle(origin)],
            [outer_r],
            marker="o",
            color="red",
            markersize=10,
        )

        if replicon_type == "chromosome" and dnaA is not None:
            fig.text(
                0.05,
                0.93,
                "\u25A0 dnaA",
                color="blue",
                fontsize=13,
                fontweight="bold",
            )

        fig.text(
            0.05,
            0.89,
            "\u25CF " + origin_name,
            color="red",
            fontsize=13,
            fontweight="bold",
        )

        desc = record.description.split(",")[0]
        parts = desc.split(" ", 1)
        line1 = parts[0]
        line2 = parts[1] if len(parts) > 1 else ""

        ax.text(
            0,
            0,
            f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

        suffix = (
            "_rotated_circular_gcskew"
            if rotated
            else "_final_circos_v9"
        )

        os.makedirs(output_dir, exist_ok=True)
        png = os.path.join(output_dir, f"{record.id}{suffix}.png")
        svg = os.path.join(output_dir, f"{record.id}{suffix}.svg")

        fig.savefig(png, dpi=400, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        return png

    finally:
        plt.close(fig)


def _prepare_plot_inputs(
    fna_path,
    accession,
    dnaA_start,
    oriC_list,
    output_dir,
    oriV_start=None,
    replicon_type="chromosome",
    topology="circular",
):
    """Load one replicon and select the correct origin marker."""
    if topology != "circular":
        print(f"Skipping {accession}: topology is {topology!r}; Circos plots require circular replicons.")
        return None

    record = _find_record(fna_path, accession)
    seq = str(record.seq)
    genome_len = len(seq)

    if replicon_type == "plasmid":
        origin = oriV_start
        origin_name = "oriV"
        dnaA = None
        gc_window = 100
        if origin is None:
            print(f"Skipping {accession}: circular plasmid has no oriV hit.")
            return None
    else:
        origin = _normalise_origin(oriC_list)
        origin_name = "oriC"
        dnaA = dnaA_start
        gc_window = 1000
        if dnaA is None or origin is None:
            print(f"Skipping {accession}: chromosome requires both dnaA and oriC hits.")
            return None

    _validate_coordinate(origin_name, origin, genome_len, accession)
    if dnaA is not None:
        _validate_coordinate("dnaA", dnaA, genome_len, accession)

    skew = fast_gc_skew(seq, win=gc_window)
    os.makedirs(output_dir, exist_ok=True)

    return record, genome_len, skew, replicon_type, dnaA, origin, origin_name


def plot_circos(
    fna_path,
    accession,
    dnaA_start,
    oriC_list,
    output_dir,
    oriV_start=None,
    replicon_type="chromosome",
    topology="circular",
):
    """
    Create the non-rotated Circos plot.

    The first five parameters match main.py. The final three allow
    main.py to pass plasmid oriV and replicon metadata.
    """
    prepared = _prepare_plot_inputs(
        fna_path,
        accession,
        dnaA_start,
        oriC_list,
        output_dir,
        oriV_start=oriV_start,
        replicon_type=replicon_type,
        topology=topology,
    )
    if prepared is None:
        return None

    record, genome_len, skew, rep_type, dnaA, origin, origin_name = prepared
    return generate_circos(
        record,
        genome_len,
        skew,
        rep_type,
        origin,
        origin_name,
        output_dir,
        dnaA=dnaA,
        rotated=False,
    )


def plot_rotated_circos(
    fna_path,
    accession,
    dnaA_start,
    oriC_list,
    shift_bp,
    output_dir,
    oriV_start=None,
    replicon_type="chromosome",
    topology="circular",
):
    """
    Create the rotated Circos plot.

    The first six parameters match main.py. For plasmids, main.py should
    pass shift_bp=oriV_start so oriV is placed at the rotation origin.
    """
    prepared = _prepare_plot_inputs(
        fna_path,
        accession,
        dnaA_start,
        oriC_list,
        output_dir,
        oriV_start=oriV_start,
        replicon_type=replicon_type,
        topology=topology,
    )
    if prepared is None:
        return None

    record, genome_len, skew, rep_type, dnaA, origin, origin_name = prepared

    if rep_type == "plasmid":
        if shift_bp is None:
            shift_bp = origin
    elif shift_bp is None:
        # Keep the main.py-supplied shift when present; otherwise derive it.
        oriC = origin
        a, b = sorted([dnaA, oriC])
        shift_bp = b if (b - a) > genome_len / 2 else a

    _validate_coordinate("rotation", shift_bp, genome_len, record.id)

    return generate_circos(
        record,
        genome_len,
        skew,
        rep_type,
        origin,
        origin_name,
        output_dir,
        dnaA=dnaA,
        rotated=True,
        shift_bp=shift_bp,
    )


def plot_all_replicons(fna_path, report_path, output_dir=None):
    """
    Compatibility helper for standalone use.

    The main workflow should call plot_circos/plot_rotated_circos with
    parsed report metadata, as shown in main.py.
    """
    raise NotImplementedError(
        "Use main.py's generate_all_plots() workflow so oriV, topology, "
        "and replicon type are passed explicitly."
    )
