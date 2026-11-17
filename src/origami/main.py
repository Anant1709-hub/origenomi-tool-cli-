import os
import shutil
import re
from Bio import SeqIO
from .circos_tools import plot_circos, plot_rotated_circos
from .dot import plot_dotplot
from .linear import run_linear_synteny_plot
from .genome_stat import compute_genome_stats

def _find_record(fna_path, accession):
    for record in SeqIO.parse(fna_path, "fasta"):
        if record.id == accession:
            return record

    raise ValueError(
        f"Accession {accession} not found in {fna_path}"
    )

def parse_report(report_path):
    """Parse chromosome and plasmid marker information from the report."""
    results = {}

    with open(report_path, encoding="utf-8") as f:
        content = f.read()

    sections = re.split(r"(?=^>>)", content, flags=re.MULTILINE)

    for section in sections:
        section = section.strip()
        if not section.startswith(">>"):
            continue

        header = section.splitlines()[0]
        seq_id = header[2:].split()[0]

        replicon_type = (
            "plasmid"
            if re.search(r"\bplasmid\b", header, re.IGNORECASE)
            else "chromosome"
        )

        topology_match = re.search(
            r"\b(circular|linear)\b",
            section,
            re.IGNORECASE,
        )
        topology = topology_match.group(1).lower() if topology_match else None

        def first_coord(pattern):
            match = re.search(
                pattern,
                section,
                re.IGNORECASE | re.DOTALL,
            )
            return int(match.group(1)) if match else None

        dnaA_start = first_coord(
            r"Chosen dnaA:.*?Coords:\s*(\d+)\s*-\s*(\d+)"
        )
        oriV_start = first_coord(
            r"Chosen oriV:.*?Coords:\s*(\d+)\s*-\s*(\d+)"
        )

        oriC_match = re.search(
            r"Chosen oriC:.*?Coords:\s*((?:\d+\s*-\s*\d+"
            r"(?:\s*[,;]\s*\d+\s*-\s*\d+)*)?)",
            section,
            re.IGNORECASE | re.DOTALL,
        )

        oriC_list = []
        if oriC_match:
            oriC_list = [
                (int(a), int(b))
                for a, b in re.findall(
                    r"(\d+)\s*-\s*(\d+)",
                    oriC_match.group(1),
                )
            ]

        results[seq_id] = {
            "replicon_type": replicon_type,
            "topology": topology,
            "dnaA_found": dnaA_start is not None,
            "oriC_found": bool(oriC_list),
            "oriV_found": oriV_start is not None,
            "dnaA_start": dnaA_start,
            "oriC_list": oriC_list,
            "oriV_start": oriV_start,
        }

    return results


def generate_all_plots(
    fasta_output_path,
    report_path,
    rotated_fasta_path,
    work_dir,
    gcf_fasta_path,
    genome_length,
):
    # One shared directory for all generated plot types.
    plot_dir = os.path.join(work_dir, "plot")
    os.makedirs(plot_dir, exist_ok=True)

    report_info = parse_report(report_path)

    for seq_id, data in report_info.items():
        if data["topology"] != "circular":
            print(f"Skipping Circos for {seq_id}: not circular.")
            continue

        
        if data["replicon_type"] == "plasmid":

            # Circular plasmid: oriV is the rotation anchor.
            if not data["oriV_found"]:
                print(f"Skipping {seq_id}: circular plasmid has no oriV hit.")
                continue

            dnaA_start = None
            oriC_list = []
            oriV_start = data["oriV_start"]
            shift_bp = oriV_start

        else:

            # Circular chromosome: both dnaA and oriC are required.
            if not (data["dnaA_found"] and data["oriC_found"]):
                print(f"Skipping {seq_id}: chromosome lacks dnaA or oriC.")
                continue

            dnaA_start = data["dnaA_start"]
            oriC_list = data["oriC_list"]
            oriV_start = None

            oriC_start = oriC_list[0][0]

            # Calculate circular distances between dnaA and oriC.
            genome_len = len(
                _find_record(gcf_fasta_path, seq_id).seq
            )

            dist_dnaa_to_oric = (
                oriC_start - dnaA_start
            ) % genome_len

            dist_oric_to_dnaa = (
                dnaA_start - oriC_start
            ) % genome_len

            # Select the anchor that keeps both markers together.
            if dist_dnaa_to_oric < dist_oric_to_dnaa:
                shift_bp = dnaA_start
                print(f"{seq_id}: Rotating from dnaA at {shift_bp}")

            else:
                shift_bp = oriC_start
                print(f"{seq_id}: Rotating from oriC at {shift_bp}")


        orig_png = plot_circos(
            gcf_fasta_path,
            seq_id,
            dnaA_start,
            oriC_list,
            plot_dir,
            oriV_start=oriV_start,
            replicon_type=data["replicon_type"],
            topology=data["topology"],
        )

        rot_png = plot_rotated_circos(
            rotated_fasta_path,
            seq_id,
            dnaA_start,
            oriC_list,
            shift_bp,
            plot_dir,
            oriV_start=oriV_start,
            replicon_type=data["replicon_type"],
            topology=data["topology"],
        )

        # Dot plot output goes directly into the same plot directory.
        dot_png = plot_dotplot(
            os.path.basename(fasta_output_path),
            os.path.basename(rotated_fasta_path),
            os.path.basename(report_path),
            output_dir=work_dir,
            accession=seq_id,
            plot_output_dir=plot_dir
        )

        # Do not copy a plot onto itself; preserve the dot.py filename.
        # If your downstream naming convention needs a renamed copy, add it here.
        if dot_png and os.path.abspath(os.path.dirname(dot_png)) != os.path.abspath(plot_dir):
            shutil.copy2(dot_png, os.path.join(plot_dir, f"{seq_id}_dotplot.png"))

        synteny_dir = os.path.join(work_dir, "linear_synteny")
        os.makedirs(synteny_dir, exist_ok=True)

        syn_png = run_linear_synteny_plot(synteny_dir)
        if syn_png and os.path.exists(syn_png):
            destination = os.path.join(plot_dir, f"{seq_id}_linear_synteny.png")
            if os.path.abspath(syn_png) != os.path.abspath(destination):
                shutil.copy2(syn_png, destination)

    return plot_dir
