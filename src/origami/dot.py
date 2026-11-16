#!/usr/bin/env python3

import os
import re
import tempfile
import subprocess
import matplotlib.pyplot as plt

from Bio import SeqIO


# ============================================================
# 1. EXTRACT MARKERS AND TOPOLOGY FROM REPORT
# ============================================================

def extract_coords(report_path, accession):

    dnaA = None
    oriC = None
    oriV = None

    replicon_type = None
    topology = None

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

        # Extract accession from report header
        report_accession = header[2:].split()[0]

        # Match the requested accession
        if report_accession != accession:
            continue

        # Identify replicon type
        if re.search(r'\bplasmid\b', header, re.IGNORECASE):
            replicon_type = "plasmid"
        else:
            replicon_type = "chromosome"

        # Identify topology from this accession's header only
        topology_match = re.search(
            r'\b(circular|linear)\b',
            header,
            re.IGNORECASE
        )

        topology = (
            topology_match.group(1).lower()
            if topology_match else None
        )

        if topology is None:
            raise ValueError(
                f"Topology missing for {accession} in report."
            )

        # Extract dnaA coordinates
        match = re.search(
            r'Chosen dnaA:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            dnaA = int(match.group(1))

        # Extract oriC coordinates
        match = re.search(
            r'Chosen oriC:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            oriC = int(match.group(1))

        # Extract oriV coordinates
        match = re.search(
            r'Chosen oriV:.*?Coords:\s*(\d+)\s*-\s*(\d+)',
            section,
            re.IGNORECASE | re.DOTALL
        )

        if match:
            oriV = int(match.group(1))

        # Return appropriate markers
        if replicon_type == "plasmid":

            if oriV is None:
                raise ValueError(
                    f"oriV coordinates missing for {accession}."
                )

            return None, oriV, "oriV", replicon_type, topology

        else:

            if dnaA is None or oriC is None:
                raise ValueError(
                    f"dnaA or oriC coordinates missing for {accession}."
                )

            return dnaA, oriC, "oriC", replicon_type, topology

    raise ValueError(
        f"Accession {accession} not found in report."
    )


# ============================================================
# 2. ROTATION CALCULATION
# ============================================================

def compute_rotated_positions(
    dnaA,
    origin,
    genome_len,
    replicon_type,
    topology
):

    # --------------------------------------------------------
    # ANY LINEAR REPLICON: NO ROTATION
    # --------------------------------------------------------

    if topology == "linear":

        origin_rot = origin

        if replicon_type == "chromosome":
            return 0, dnaA, origin_rot

        return 0, None, origin_rot

    # --------------------------------------------------------
    # CIRCULAR PLASMID: ROTATE USING oriV
    # --------------------------------------------------------

    if replicon_type == "plasmid":

        shift_bp = origin

        origin_rot = (origin - shift_bp) % genome_len

        return shift_bp, None, origin_rot

    # --------------------------------------------------------
    # CIRCULAR CHROMOSOME: EXISTING SHORTEST-ARC RULE
    # --------------------------------------------------------

    a, b = sorted([dnaA, origin])

    d = b - a

    shift_bp = b if d > genome_len / 2 else a

    dnaA_rot = (dnaA - shift_bp) % genome_len

    origin_rot = (origin - shift_bp) % genome_len

    return shift_bp, dnaA_rot, origin_rot


# ============================================================
# 3. PARSE NUCMER COORDINATES
# ============================================================

def parse_blocks(coords_path):

    blocks = []

    with open(coords_path) as f:

        for line in f:

            parts = line.split()

            if len(parts) >= 4 and parts[0].isdigit():

                s1, e1, s2, e2 = map(int, parts[:4])

                blocks.append((s1, e1, s2, e2))

    return blocks


# ============================================================
# 4. RUN NUCMER + SHOW-COORDS
# ============================================================

def run_nucmer(
    original_fna,
    comparison_fna,
    work_dir,
    accession
):

    # Unique temporary prefix for every replicon
    prefix = os.path.join(
        work_dir,
        f".{accession}_dotplot_output"
    )

    subprocess.run(
        [
            "nucmer",
            "--mum",
            "-g", "1000",
            "-p", prefix,
            original_fna,
            comparison_fna
        ],
        check=True
    )

    coords_path = os.path.join(
        work_dir,
        f".{accession}_coords.txt"
    )

    with open(coords_path, "w") as out:

        subprocess.run(
            [
                "show-coords",
                "-THrd",
                f"{prefix}.delta"
            ],
            stdout=out,
            check=True
        )

    return coords_path, prefix


# ============================================================
# 5. GENERATE DOT PLOT
# ============================================================

def generate_dotplot(
    blocks,
    dnaA,
    origin,
    dnaA_rot,
    origin_rot,
    origin_name,
    replicon_type,
    topology,
    output_dir,
    accession
):

    plt.figure(figsize=(8, 8))

    # --------------------------------------------------------
    # ALIGNMENT BLOCKS (UNCHANGED)
    # --------------------------------------------------------

    for s1, e1, s2, e2 in blocks:

        same_direction = (e2 - s2) * (e1 - s1) > 0

        color = "blue" if same_direction else "orange"

        plt.plot(
            [s1, e1],
            [s2, e2],
            color=color,
            linewidth=1
        )

    # --------------------------------------------------------
    # CHROMOSOME: dnaA MARKERS (RED)
    # --------------------------------------------------------

    if replicon_type == "chromosome":

        plt.axvline(
            dnaA,
            color="red",
            linewidth=1.8
        )

        plt.axhline(
            dnaA_rot,
            color="red",
            linewidth=1.8
        )

    # --------------------------------------------------------
    # ORIGIN MARKERS: oriC OR oriV (GREEN)
    # --------------------------------------------------------

    plt.axvline(
        origin,
        color="green",
        linestyle=":",
        linewidth=1.8
    )

    plt.axhline(
        origin_rot,
        color="green",
        linestyle=":",
        linewidth=1.8
    )

    # --------------------------------------------------------
    # DYNAMIC TITLE
    # --------------------------------------------------------

    if replicon_type == "chromosome":

        title = "Dot Plot (dnaA = red, oriC = green)"

    else:

        title = "Dot Plot (oriV = green)"

    plt.title(
        title,
        fontsize=11,
        fontweight="bold"
    )

    plt.xlabel("Original Genome Position (bp)")

    if topology == "linear":

        plt.ylabel("Genome Position (bp)")

    else:

        plt.ylabel("Rotated Genome Position (bp)")

    plt.tight_layout()

    # --------------------------------------------------------
    # SAVE ALL PLOTS IN THE SAME OUTPUT DIRECTORY
    # --------------------------------------------------------

    png = os.path.join(
        output_dir,
        f"{accession}_dotplot_highlighted_circos_logic.png"
    )

    svg = os.path.join(
        output_dir,
        f"{accession}_dotplot_highlighted_circos_logic.svg"
    )

    plt.savefig(png, dpi=300)

    plt.savefig(svg)

    plt.close()

    return png


# ============================================================
# 6. MAIN FUNCTION FOR ONE REPLICON
# ============================================================

def plot_dotplot(
    original_fna,
    rotated_fna,
    report_file,
    output_dir,
    plot_output_dir=None
):

    original_path = os.path.join(
        output_dir,
        original_fna
    )

    rotated_path = os.path.join(
        output_dir,
        rotated_fna
    )

    report_path = os.path.join(
        output_dir,
        report_file
    )

    # Keep input FASTA/report files in output_dir, but save
    # generated dot plots in the shared plots directory.
    if plot_output_dir is None:
        plot_output_dir = os.path.join(output_dir, "plots")

    os.makedirs(plot_output_dir, exist_ok=True)

    # --------------------------------------------------------
    # READ ORIGINAL FASTA RECORD
    # --------------------------------------------------------

    original_record = next(
        SeqIO.parse(original_path, "fasta")
    )

    accession = original_record.id

    genome_len = len(original_record.seq)

    # --------------------------------------------------------
    # EXTRACT MARKERS AND TOPOLOGY
    # --------------------------------------------------------

    dnaA, origin, origin_name, replicon_type, topology = extract_coords(
        report_path,
        accession
    )

    print("\nAccession:", accession)
    print("Replicon type:", replicon_type)
    print("Topology:", topology)
    print("Genome length:", genome_len)
    print("Origin:", origin_name, origin)

    # --------------------------------------------------------
    # SELECT FASTA FOR COMPARISON
    # --------------------------------------------------------

    if topology == "linear":

        # Linear chromosome or plasmid:
        # Compare original FASTA against itself
        comparison_path = original_path

        print("Linear replicon detected.")
        print("Rotation disabled.")

    else:

        # Circular chromosome or plasmid:
        # Compare original FASTA against rotated FASTA
        comparison_path = rotated_path

        print("Circular replicon detected.")
        print("Using rotated FASTA.")

    # --------------------------------------------------------
    # COMPUTE MARKER POSITIONS
    # --------------------------------------------------------

    _, dnaA_rot, origin_rot = compute_rotated_positions(
        dnaA,
        origin,
        genome_len,
        replicon_type,
        topology
    )

    print("Original origin position:", origin)
    print("Comparison origin position:", origin_rot)

    if replicon_type == "chromosome":

        print("Original dnaA position:", dnaA)
        print("Comparison dnaA position:", dnaA_rot)

    # --------------------------------------------------------
    # RUN NUCMER AND GENERATE PLOT
    # --------------------------------------------------------

    coords_path = None
    prefix = None

    try:

        coords_path, prefix = run_nucmer(
            original_path,
            comparison_path,
            output_dir,
            accession
        )

        blocks = parse_blocks(coords_path)

        print("Alignment blocks:", len(blocks))

        output_png = generate_dotplot(
            blocks,
            dnaA,
            origin,
            dnaA_rot,
            origin_rot,
            origin_name,
            replicon_type,
            topology,
            plot_output_dir,
            accession
        )

        print("Successfully generated:", output_png)

        return output_png

    finally:

        # Remove temporary NUCmer output files
        if coords_path and os.path.exists(coords_path):
            os.remove(coords_path)

        if prefix:

            for extension in [
                ".delta",
                ".ntref",
                ".ntqry"
            ]:

                temp_file = prefix + extension

                if os.path.exists(temp_file):
                    os.remove(temp_file)


# ============================================================
# 7. PROCESS ALL CHROMOSOMES AND PLASMIDS
# ============================================================

def plot_all_replicons(
    original_fna,
    rotated_fna,
    report_file,
    output_dir
):

    original_path = os.path.join(
        output_dir,
        original_fna
    )

    rotated_path = os.path.join(
        output_dir,
        rotated_fna
    )

    # Shared output folder used by circos_tools.py as well.
    # Both scripts should point to this same parent directory.
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # --------------------------------------------------------
    # LOAD ALL FASTA RECORDS
    # --------------------------------------------------------

    original_records = {
        record.id: record
        for record in SeqIO.parse(original_path, "fasta")
    }

    rotated_records = {
        record.id: record
        for record in SeqIO.parse(rotated_path, "fasta")
    }

    print("\n" + "=" * 65)
    print("ORIGENOMI DOT PLOT GENERATION")
    print("=" * 65)

    print("Original FASTA records:", len(original_records))
    print("Rotated FASTA records:", len(rotated_records))

    # Match original and rotated records by accession
    common_accessions = [
        accession
        for accession in original_records
        if accession in rotated_records
    ]

    missing_accessions = [
        accession
        for accession in original_records
        if accession not in rotated_records
    ]

    print("Matching replicons:", len(common_accessions))

    if missing_accessions:

        print(
            "Replicons missing from rotated FASTA:",
            missing_accessions
        )

    successful = 0
    failed = 0

    # --------------------------------------------------------
    # PROCESS EACH REPLICON
    # --------------------------------------------------------

    for accession in common_accessions:

        print("\n" + "=" * 65)
        print("PROCESSING:", accession)
        print("=" * 65)

        original_record = original_records[accession]
        rotated_record = rotated_records[accession]

        original_temp = None
        rotated_temp = None

        try:

            # ------------------------------------------------
            # CREATE TEMPORARY SINGLE-RECORD FASTA FILES
            # ------------------------------------------------

            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".fna",
                prefix=f".{accession}_original_",
                dir=output_dir,
                delete=False
            ) as temp:

                original_temp = temp.name

            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".fna",
                prefix=f".{accession}_rotated_",
                dir=output_dir,
                delete=False
            ) as temp:

                rotated_temp = temp.name

            # Write original record
            SeqIO.write(
                original_record,
                original_temp,
                "fasta"
            )

            # Write corresponding rotated record
            SeqIO.write(
                rotated_record,
                rotated_temp,
                "fasta"
            )

            # ------------------------------------------------
            # GENERATE DOT PLOT
            # ------------------------------------------------

            output_png = plot_dotplot(
                os.path.basename(original_temp),
                os.path.basename(rotated_temp),
                report_file,
                output_dir,
                plots_dir
            )

            print("PNG:", output_png)

            print(
                "SVG:",
                output_png.replace(".png", ".svg")
            )

            successful += 1

        except Exception as e:

            print(
                f"ERROR processing {accession}: {e}"
            )

            failed += 1

        finally:

            # --------------------------------------------
            # REMOVE TEMPORARY FASTA FILES
            # --------------------------------------------

            if original_temp and os.path.exists(original_temp):
                os.remove(original_temp)

            if rotated_temp and os.path.exists(rotated_temp):
                os.remove(rotated_temp)

    # --------------------------------------------------------
    # FINAL SUMMARY
    # --------------------------------------------------------

    print("\n" + "=" * 65)
    print("DOT PLOT GENERATION COMPLETED")
    print("=" * 65)

    print("Total original replicons:", len(original_records))
    print("Matching replicons:", len(common_accessions))
    print("Successfully generated:", successful)
    print("Failed:", failed)

    print("\nOutput directory:", plots_dir)
