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
