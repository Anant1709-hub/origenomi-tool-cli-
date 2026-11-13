import os
import shutil
import re

from .circos_tools import plot_circos, plot_rotated_circos
from .dot import plot_dotplot
from .linear import run_linear_synteny_plot
from .genome_stat import compute_genome_stats

def parse_report(report_path):
    results = {}
    current_seq = None

    with open(report_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">>"):
                seq_id = line[2:].split()[0]
                results[seq_id] = {
                    "dnaA_found": False,
                    "oriC_found": False,
                    "dnaA_start": None,
                    "oriC_list": []
                }
                current_seq = seq_id
                continue

            if "Chosen dnaA:" in line:
                next(f)
                coords = next(f).strip()
                m = re.search(r"(\d+)\s*-\s*(\d+)", coords)
                if m:
                    results[current_seq]["dnaA_found"] = True
                    results[current_seq]["dnaA_start"] = int(m.group(1))

            if "Chosen oriC:" in line:
                next(f)
                coords = next(f).strip()
                pairs = re.findall(r"(\d+)\s*-\s*(\d+)", coords)
                if pairs:
                    results[current_seq]["oriC_found"] = True
                    results[current_seq]["oriC_list"] = [(int(a), int(b)) for a, b in pairs]

    return results


def generate_all_plots(fasta_output_path, report_path, rotated_fasta_path, work_dir, gcf_fasta_path,genome_length):

    plot_dir = os.path.join(work_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    report_info = parse_report(report_path)

    for seq_id, data in report_info.items():

        if not (data["dnaA_found"] and data["oriC_found"]):
            continue

        orig_png = plot_circos(
            gcf_fasta_path,
            seq_id,
            data["dnaA_start"],
            data["oriC_list"],
            plot_dir
        )

        shift_bp = data["dnaA_start"]

        rot_png = plot_rotated_circos(
            rotated_fasta_path,
            seq_id,
            data["dnaA_start"],
            data["oriC_list"],
            shift_bp,
            plot_dir
        )

        dot_png = plot_dotplot(
            os.path.basename(fasta_output_path),
            os.path.basename(rotated_fasta_path),
            os.path.basename(report_path),
            output_dir=work_dir,
        )

        if dot_png:
            shutil.copy(dot_png, os.path.join(plot_dir, f"{seq_id}_dotplot.png"))

        synteny_dir = os.path.join(work_dir, "linear_synteny")
        os.makedirs(synteny_dir, exist_ok=True)

        syn_png = run_linear_synteny_plot(synteny_dir)
        if syn_png:
            shutil.copy(syn_png, os.path.join(plot_dir, f"{seq_id}_linear_synteny.png"))
