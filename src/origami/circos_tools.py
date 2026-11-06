import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


# ============================================================
# SHARED FUNCTION
# ============================================================

def fast_gc_skew(seq, win=1000):

    seq = np.frombuffer(seq.upper().encode(), dtype="S1")
    g = (seq == b"G").astype(int)
    c = (seq == b"C").astype(int)

    cum_g = np.cumsum(g)
    cum_c = np.cumsum(c)

    skew = np.empty(len(seq) // win)

    for i in range(len(skew)):
        start = i * win
        end = start + win
        g_win = cum_g[end - 1] - (cum_g[start - 1] if start > 0 else 0)
        c_win = cum_c[end - 1] - (cum_c[start - 1] if start > 0 else 0)
        denom = g_win + c_win
        skew[i] = (g_win - c_win) / denom if denom > 0 else 0

    return skew



# ============================================================
# NON-ROTATED CIRCOS (PER CHROMOSOME)
# ============================================================

def plot_circos(fna_path, seq_id, dnaA_start, oriC_list, out_dir):

    record = None
    for rec in SeqIO.parse(fna_path, "fasta"):
        if rec.id == seq_id:
            record = rec
            break

    if record is None:
        raise ValueError(f"{seq_id} not found in FASTA")

    seq = str(record.seq)
    genome_len = len(seq)

    dnaA = dnaA_start
    oriC = oriC_list[0][0] if oriC_list else None

    skew = fast_gc_skew(seq)
    angles = np.linspace(0, 2 * np.pi, len(skew), endpoint=False)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()

    outer_r = 1.0
    ax.plot(np.linspace(0, 2 * np.pi, 2000), [outer_r] * 2000, color="navy", lw=4)

    bar_r = 0.8
    bar_scale = 0.25

    ax.bar(angles[skew > 0], skew[skew > 0] * bar_scale,
           width=2 * np.pi / len(skew), bottom=bar_r,
           color="lightblue", linewidth=0)

    ax.bar(angles[skew < 0], -skew[skew < 0] * bar_scale,
           width=2 * np.pi / len(skew),
           bottom=bar_r - bar_scale * np.abs(skew[skew < 0]),
           color="orange", linewidth=0)

    # ticks
    num_ticks = 28
    for i in range(num_ticks):
        pos = genome_len / num_ticks * i
        theta = 2 * np.pi * (pos / genome_len)
        ax.plot([theta, theta], [outer_r, outer_r + 0.03], color="navy", lw=1.5)
        ax.text(theta, outer_r + 0.10, f"{pos/1e6:.2f} Mb",
                ha="center", va="center", fontsize=9)

    if dnaA:
        ax.plot([2 * np.pi * (dnaA / genome_len)], [outer_r],
                marker="s", color="blue", markersize=10)

    if oriC:
        ax.plot([2 * np.pi * (oriC / genome_len)], [outer_r],
                marker="o", color="red", markersize=10)

    desc = record.description.split(",")[0]
    parts = desc.split(" ", 1)

    line1 = parts[0]            
    line2 = parts[1] if len(parts) > 1 else "" 

    ax.text(0, 0,
            f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha="center", va="center",
            fontsize=10, fontweight="bold")

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"{seq_id}_original_circos.png")
    plt.savefig(png, dpi=400, bbox_inches="tight")
    plt.close()

    return png


# ============================================================
# ROTATED CIRCOS (PER CHROMOSOME)
# ============================================================

def plot_rotated_circos(fna_path, seq_id, dnaA_start, oriC_list, shift_bp, out_dir):

    record = None
    for rec in SeqIO.parse(fna_path, "fasta"):
        if rec.id == seq_id:
            record = rec
            break

    if record is None:
        raise ValueError(f"{seq_id} not found in FASTA")

    seq = str(record.seq)
    genome_len = len(seq)

    skew = fast_gc_skew(seq)
    bins = len(skew)
    angles = np.linspace(0, 2 * np.pi, bins, endpoint=False)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_axis_off()

    # Outer ring
    outer_r = 1.0
    ax.plot(np.linspace(0, 2 * np.pi, 2000),
            [outer_r] * 2000, color="navy", lw=4)

    # GC skew bars
    bar_scale = 0.25
    r_inner = 0.8

    ax.bar(angles[skew > 0], skew[skew > 0] * bar_scale,
           width=2 * np.pi / bins, bottom=r_inner,
           color="lightblue", linewidth=0)

    ax.bar(angles[skew < 0], -skew[skew < 0] * bar_scale,
           width=2 * np.pi / bins,
           bottom=r_inner - bar_scale * np.abs(skew[skew < 0]),
           color="orange", linewidth=0)

    num_ticks = 28
    tick_spacing = genome_len / num_ticks
    label_offset = 0.08

    for i in range(num_ticks):
        pos = tick_spacing * i
        theta = 2 * np.pi * i / num_ticks

        ax.plot([theta, theta],
                [outer_r, outer_r + 0.03],
                color="navy", lw=1.5)

        ax.text(theta, outer_r + label_offset,
                f"{pos/1e6:.2f} Mb",
                ha="center", va="center", fontsize=9)

    # Rotated markers
    dnaA_rot = (dnaA_start - shift_bp) % genome_len if dnaA_start else None
    oriC_rot = (oriC_list[0][0] - shift_bp) % genome_len if oriC_list else None

    if dnaA_rot is not None:
        ax.plot([2 * np.pi * (dnaA_rot / genome_len)],
                [outer_r],
                marker="s", color="blue", markersize=10)

    if oriC_rot is not None:
        ax.plot([2 * np.pi * (oriC_rot / genome_len)],
                [outer_r],
                marker="o", color="red", markersize=10)

    # Center label
    desc = record.description.split(",")[0]
    parts = desc.split(" ", 1)

    line1 = parts[0]            
    line2 = parts[1] if len(parts) > 1 else ""  

    ax.text(0, 0,
            f"{line1}\n{line2}\n({genome_len:,} bp)",
            ha="center", va="center",
            fontsize=10, fontweight="bold")

    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"{seq_id}_rotated_circos.png")
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.close()

    return png
