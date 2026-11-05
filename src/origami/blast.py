import subprocess
from subprocess import DEVNULL

# === Annotation (dnaA / oriC) ===
# Keep your exact outfmt for annotation BLASTs
OUTFMT_ANN = (
    "7 qseqid sseqid qstart qend sstart send qlen slen "
    "pident evalue bitscore mismatch gaps length"
)

def run_blast_annotation(
    query_fasta: str,
    db_prefix: str,
    out_path: str,
    evalue: str = "1e-5",
    quiet: bool = True,
) -> None:
    """
    One BLAST per DB (dnaA, oriC) over the concatenated FASTA.
    """
    args = [
        "blastn",
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", out_path,
        "-outfmt", OUTFMT_ANN,
        "-evalue", str(evalue),
    ]
    subprocess.run(
        args,
        check=True,
        stdout=DEVNULL if quiet else None,
        stderr=DEVNULL if quiet else None,
    )

def parse_blast_annotation(path: str):
    """
    Parse annotation outfmt into dict[qseqid] -> list[hit dicts].
    """
    hits = {}
    num_of_hits = {}
    with open(path) as f:
        lines = iter(f)

        for line in lines:
            if not line.strip() or line.startswith("#"):
                if line.startswith("# Query:"):
                    NZ_id = line.split()[2]
                    num = None

                    for next_line in lines:
                        if "hits" in next_line:
                            num = int(next_line.split()[1])
                            break
                    
                    num_of_hits[NZ_id] = num
                continue
            
            (qid, sid, qstart, qend, sstart, send,
             qlen, slen, pident, evalue, bitscore,
             mismatch, gaps, length) = line.strip().split()

            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)

            entry = {
                "sseqid": sid,
                "qstart": min(qstart, qend),
                "qend":   max(qstart, qend),
                "sstart": min(sstart, send),
                "send":   max(sstart, send),
                "qlen": int(qlen),
                "slen": int(slen),
                "length": int(length),
                "pident": float(pident),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "mismatch": int(mismatch),
                "gaps": int(gaps),
            }
            hits.setdefault(qid, []).append(entry)
    return hits,num_of_hits

def parse_oriv_blast_annotation(path: str,oric_hits: list):
    """
    Parse annotation outfmt into dict[qseqid] -> list[hit dicts].
    """
    replicon_to_be_processed = []
    for key,value in oric_hits.items():
        if value == 0:
           replicon_to_be_processed.append(key) 

    hits = {}
    with open(path) as f:
        lines = iter(f)
        
        for line in lines:
            if not line.strip() or line.startswith("#"):
                if line.startswith("# Query:"):
                    NZ_id = line.split()[2]
                    if NZ_id not in replicon_to_be_processed:
                        for next_line in lines:
                            if "Query" in next_line:
                                break

                continue
            
            (qid, sid, qstart, qend, sstart, send,
             qlen, slen, pident, evalue, bitscore,
             mismatch, gaps, length) = line.strip().split()

            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)

            entry = {
                "sseqid": sid,
                "qstart": min(qstart, qend),
                "qend":   max(qstart, qend),
                "sstart": min(sstart, send),
                "send":   max(sstart, send),
                "qlen": int(qlen),
                "slen": int(slen),
                "length": int(length),
                "pident": float(pident),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "mismatch": int(mismatch),
                "gaps": int(gaps),
            }
            hits.setdefault(qid, []).append(entry)
    return hits


# === Trim (half2 → DB(half1)) ===
# Minimal columns, strict flags; we sort by length and a couple tie-breakers
OUTFMT_TRIM = "6 qseqid sseqid qstart qend sstart send qlen slen pident length evalue"

def run_blast_trim(
    query_fasta: str,
    db_prefix: str,
    out_path: str,
    evalue: str = "1e-10",
    quiet: bool = True,
) -> None:
    """
    BLAST for terminal overlap detection (half2 vs DB(half1)).
    """
    args = [
        "blastn",
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", out_path,
        "-outfmt", OUTFMT_TRIM,
        "-evalue", str(evalue),
        "-dust", "no",
        "-soft_masking", "false",
        "-max_hsps", "1",
    ]
    subprocess.run(
        args,
        check=True,
        stdout=DEVNULL if quiet else None,
        stderr=DEVNULL if quiet else None,
    )

def parse_blast_trim(path: str):
    """
    Parse trim outfmt into a flat list of hits (dicts).
    """
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            (qid, sid, qstart, qend, sstart, send,
             qlen, slen, pident, length, evalue) = line.strip().split()
            rows.append({
                "qseqid": qid, "sseqid": sid,
                "qstart": int(qstart), "qend": int(qend),
                "sstart": int(sstart), "send": int(send),
                "qlen": int(qlen), "slen": int(slen),
                "pident": float(pident), "length": int(length),
                "evalue": float(evalue),
            })
    return rows
