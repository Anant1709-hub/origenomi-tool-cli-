#!/usr/bin/env python3

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


# ============================================================
# 1. READ REPLICON INFORMATION FROM ORIGAMI REPORT
# ============================================================

def read_replicon_report(report_path, accession):
    """
    Extract topology, replicon type, dnaA, oriC and oriV
    for a particular replicon from the Origami report.
    """

    with open(report_path) as f:
        content = f.read()

    # Split report into individual replicon sections
    sections = re.split(
        r'(?=^>>)',
        content,
        flags=re.MULTILINE
    )

    for section in sections:

        if not section.strip().startswith(">>"):
            continue

        lines = section.strip().splitlines()
        header = lines[0]

        # Match accession ID
        report_accession = header[2:].split()[0]

        if report_accession != accession:
            continue

        # Identify replicon type
        if re.search(r'\bplasmid\b', header, re.IGNORECASE):
            replicon_type = "plasmid"
        else:
            replicon_type = "chromosome"

        # Identify topology
        topology_match = re.search(
            r'\b(circular|linear)\b',
            section,
            re.IGNORECASE
        )

        topology = (
            topology_match.group(1).lower()
            if topology_match else None
        )

        # Extract coordinates
        dnaA = None
        oriC = None
        oriV = None

        match = re.search(
            r'Chosen dnaA:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            dnaA = int(match.group(1))

        match = re.search(
            r'Chosen oriC:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            oriC = int(match.group(1))

        match = re.search(
            r'Chosen oriV:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            oriV = int(match.group(1))

        return {
            "accession": accession,
            "replicon_type": replicon_type,
            "topology": topology,
            "dnaA": dnaA,
            "oriC": oriC,
            "oriV": oriV
        }

    raise ValueError(
        f"Accession {accession} not found in Origami report"
    )


# ============================================================
# 2. GC SKEW CALCULATION
# ============================================================

def fast_gc_skew(seq, win=1000):
    """
    Calculate GC skew using NumPy cumulative operations.
    """

    seq = np.frombuffer(
        seq.upper().encode(),
        dtype='S1'
    )

    g = (seq == b'G').astype(int)
    c = (seq == b'C').astype(int)

    cum_g = np.cumsum(g)
    cum_c = np.cumsum(c)

    bins = len(seq) // win
    skew = np.zeros(bins)

    for i in range(bins):

        start = i * win
        end = start + win

        g_win = cum_g[end - 1] - (
            cum_g[start - 1] if start > 0 else 0
        )

        c_win = cum_c[end - 1] - (
            cum_c[start - 1] if start > 0 else 0
        )

        denom = g_win + c_win

        skew[i] = (
            (g_win - c_win) / denom
            if denom > 0 else 0
        )

    return skew


# ============================================================
# 3. LOAD AND VALIDATE INPUTS
# ============================================================

def load_replicon(fna_path, report_path):

    record = next(SeqIO.parse(fna_path, "fasta"))

    seq = str(record.seq)
    genome_len = len(seq)

    accession = record.id

    info = read_replicon_report(
        report_path,
        accession
    )

    print(f"\nAccession: {accession}")
    print(f"Replicon type: {info['replicon_type']}")
    print(f"Topology: {info['topology']}")
    print(f"Genome length: {genome_len} bp")

    # Skip linear replicons
    if info["topology"] == "linear":

        print(
            f"Skipping {accession}: "
            "Linear replicon - no Circos plot generated."
        )

        return None

    if info["topology"] != "circular":

        print(
            f"Skipping {accession}: "
            "Topology missing or unrecognized."
        )

        return None

    # Select the appropriate replication origin
    if info["replicon_type"] == "chromosome":

        origin = info["oriC"]
        origin_name = "oriC"

    else:

        origin = info["oriV"]
        origin_name = "oriV"

    if origin is None:

        print(
            f"Skipping {accession}: "
            f"No {origin_name} found."
        )

        return None

    # Validate coordinates
    for name in ["dnaA", "oriC", "oriV"]:

        coord = info[name]

        if coord is not None and not (0 <= coord <= genome_len):

            raise ValueError(
                f"{name} coordinate {coord} exceeds "
                f"genome length {genome_len}"
            )

    return record, seq, genome_len, info, origin_name


# ============================================================
# 4. SHARED CIRCOS PLOT FUNCTION
# ============================================================

def generate_circos(
    record,
    genome_len,
    skew,
    info,
    rotated=False,
    shift_bp=0
):

    # Marker information
    dnaA = info["dnaA"]

    if info["replicon_type"] == "chromosome":

        origin = info["oriC"]
        origin_name = "oriC"

    else:

        origin = info["oriV"]
        origin_name = "oriV"

    # Rotate GC skew when required
    if rotated:

        bins = len(skew)
        bin_size = genome_len / bins

        shift_bins = int(round(shift_bp / bin_size))

        skew = np.roll(skew, -shift_bins)

    bins = len(skew)

    angles = np.linspace(
        0,
        2 * np.pi,
        bins,
        endpoint=False
    )

    # --------------------------------------------------------
    # POLAR PLOT SETUP
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10, 10),
        subplot_kw={"projection": "polar"}
    )

    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()

    # --------------------------------------------------------
    # OUTER GENOME RING
    # --------------------------------------------------------

    outer_r = 1.0

    ax.plot(
        np.linspace(0, 2 * np.pi, 2000),
        [outer_r] * 2000,
        color="navy",
        lw=4
    )

    # --------------------------------------------------------
    # GC SKEW SPIKES
    # --------------------------------------------------------

    bar_r = 0.80
    bar_scale = 0.25

    positive = skew > 0
    negative = skew < 0

    # Positive GC skew - Light blue outward spikes
    ax.vlines(
        angles[positive],
        bar_r,
        bar_r + skew[positive] * bar_scale,
        color="lightblue",
        linewidth=0.45,
        alpha=0.9
    )

    # Negative GC skew - Orange inward spikes
    ax.vlines(
        angles[negative],
        bar_r - np.abs(skew[negative]) * bar_scale,
        bar_r,
        color="orange",
        linewidth=0.45,
        alpha=0.9
    )

    # --------------------------------------------------------
    # GENOME COORDINATE TICKS
    # --------------------------------------------------------

    num_ticks = 28

    for i in range(num_ticks):

        pos = genome_len * i / num_ticks

        theta = 2 * np.pi * i / num_ticks

        ax.plot(
            [theta, theta],
            [outer_r, outer_r + 0.03],
            color="navy",
            lw=1.5
        )

        ax.text(
            theta,
            outer_r + 0.10,
            f"{pos / 1e6:.2f} Mb",
            ha="center",
            va="center",
            fontsize=9
        )

    # --------------------------------------------------------
    # COORDINATE TRANSFORMATION FOR MARKERS
    # --------------------------------------------------------

    def get_angle(coord):

        if rotated:
            coord = (coord - shift_bp) % genome_len

        return 2 * np.pi * (coord / genome_len)

    # dnaA marker: chromosome only
    if dnaA is not None:

        theta_dnaA = get_angle(dnaA)

        ax.plot(
            [theta_dnaA],
            [outer_r],
            marker="s",
            color="blue",
            markersize=10
        )

    # oriC or oriV marker
    if origin is not None:

        theta_origin = get_angle(origin)

        ax.plot(
            [theta_origin],
            [outer_r],
            marker="o",
            color="red",
            markersize=10
        )

    # --------------------------------------------------------
    # LEGEND
    # --------------------------------------------------------

    if dnaA is not None:

        fig.text(
            0.05,
            0.93,
            u'\u25A0 dnaA',
            color="blue",
            fontsize=13,
            fontweight="bold"
        )

    fig.text(
        0.05,
        0.89,
        u'\u25CF ' + origin_name,
        color="red",
        fontsize=13,
        fontweight="bold"
    )

    # --------------------------------------------------------
    # CENTER LABEL
    # --------------------------------------------------------

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
        fontweight="bold"
    )

    # --------------------------------------------------------
    # OUTPUT FILENAME
    # --------------------------------------------------------

    base = record.id

    if rotated:

        suffix = "_rotated_circular_gcskew"

    else:

        suffix = "_final_circos_v9"

    png = f"{base}{suffix}.png"
    svg = f"{base}{suffix}.svg"

    plt.savefig(
        png,
        dpi=400,
        bbox_inches="tight"
    )

    plt.savefig(
        svg,
        bbox_inches="tight"
    )

    plt.close()

    return png


# ============================================================
# 5. NON-ROTATED CIRCOS PLOT
# ============================================================

def plot_circos(fna_path, report_path):

    loaded = load_replicon(
        fna_path,
        report_path
    )

    if loaded is None:
        return None

    record, seq, genome_len, info, origin_name = loaded

    # Use higher GC-skew resolution for shorter plasmids
    gc_window = 100 if info["replicon_type"] == "plasmid" else 1000

    skew = fast_gc_skew(seq, win=gc_window)
    png = generate_circos(
        record,
        genome_len,
        skew,
        info,
        rotated=False
    )

    return png


# ============================================================
# 6. ROTATED CIRCOS PLOT
# ============================================================

def plot_rotated_circos(fna_path, report_path):

    loaded = load_replicon(
        fna_path,
        report_path
    )

    if loaded is None:
        return None

    record, seq, genome_len, info, origin_name = loaded

    # Select rotation point
    if info["replicon_type"] == "chromosome":

        dnaA = info["dnaA"]
        origin = info["oriC"]

        if dnaA is not None and origin is not None:

            a, b = sorted([dnaA, origin])

            if b - a > genome_len / 2:

                shift_bp = b

            else:

                shift_bp = a

        elif origin is not None:

            shift_bp = origin

        elif dnaA is not None:

            shift_bp = dnaA

        else:

            print("No rotation reference found.")
            return None

    else:

        # Plasmid: rotate using oriV
        shift_bp = info["oriV"]

    # Use the same GC-skew resolution for original and rotated plots
    gc_window = 100 if info["replicon_type"] == "plasmid" else 1000

    skew = fast_gc_skew(seq, win=gc_window)

    png = generate_circos(
        record,
        genome_len,
        skew,
        info,
        rotated=True,
        shift_bp=shift_bp
    )

    return png, origin_name, shift_bp


# ============================================================
# 7. PROCESS ALL REPLICONS IN A MULTI-FASTA FILE
# ============================================================

def plot_all_replicons(fna_path, report_path):

    results = []

    for record in SeqIO.parse(fna_path, "fasta"):

        # Create a temporary single-record FASTA
        temp_fna = f"{record.id}_temp.fna"

        SeqIO.write(record, temp_fna, "fasta")

        try:

            original = plot_circos(
                temp_fna,
                report_path
            )

            rotated = plot_rotated_circos(
                temp_fna,
                report_path
            )

            results.append({
                "accession": record.id,
                "original_plot": original,
                "rotated_plot": rotated
            })

        finally:

            if os.path.exists(temp_fna):
                os.remove(temp_fna)

    return results


