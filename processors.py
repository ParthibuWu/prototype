# processors.py
from __future__ import annotations

import io
from typing import Dict, List, Optional, Set, TextIO, Tuple

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from compression_io import open_text_handle
from format_detect import detect_file_format, sniff_fasta_fastq
from composition import atgc_content, gc_fraction


def _parse_id_list(filter_ids: Optional[str]) -> Optional[Set[str]]:
    if not filter_ids:
        return None
    # comma or whitespace separated
    raw = filter_ids.replace("\n", ",").replace(" ", ",").split(",")
    ids = {x.strip() for x in raw if x.strip()}
    return ids or None


def process_fasta_like_stream(handle: TextIO, file_format: str, mode: str, max_records: int) -> List[Dict]:
    rows: List[Dict] = []
    for i, record in enumerate(SeqIO.parse(handle, file_format), start=1):
        if i > max_records:
            break

        seq = str(record.seq)
        comp = atgc_content(seq, mode=mode)
        gc = gc_fraction(seq, mode=mode)

        rows.append({
            "ID": record.id,
            "Description": getattr(record, "description", ""),
            "Sequence": seq,  # keep for Details tab; remove if memory becomes an issue
            "Length": len(seq),
            "GC_percent": None if gc is None else 100 * gc,
            "A": comp.counts["A"],
            "T": comp.counts["T"],
            "G": comp.counts["G"],
            "C": comp.counts["C"],
            "N": comp.counts["N"],
            "AMB": comp.counts["AMB"],
            "denom_used": comp.denom,
            "First_base": seq[0] if seq else "",
            "Last_base": seq[-1] if seq else "",
        })
    return rows


def process_fastq_stream(handle: TextIO, mode: str, max_records: int, wanted_ids: Optional[Set[str]] = None) -> List[Dict]:
    rows: List[Dict] = []
    kept = 0

    # STREAMING: no handle.read()
    for title, seq, qual in FastqGeneralIterator(handle):
        record_id = title.split(None, 1)[0]

        if wanted_ids is not None and record_id not in wanted_ids:
            continue

        comp = atgc_content(seq, mode=mode)
        gc = gc_fraction(seq, mode=mode)
        avg_q = (sum((ord(c) - 33) for c in qual) / len(qual)) if qual else 0.0

        rows.append({
            "ID": record_id,
            "Title": title,
            "Sequence": seq,      # keep for Details tab
            "Quality": qual,      # optional; remove if too heavy
            "Length": len(seq),
            "GC_percent": None if gc is None else 100 * gc,
            "Avg_quality": avg_q,
            "A": comp.counts["A"],
            "T": comp.counts["T"],
            "G": comp.counts["G"],
            "C": comp.counts["C"],
            "N": comp.counts["N"],
            "AMB": comp.counts["AMB"],
            "denom_used": comp.denom,
        })

        kept += 1
        if kept >= max_records:
            break

    return rows


def get_sequence_stats(content: bytes, filename: str, mode: str = "canonical") -> Dict:
    """
    Stats-only: memory efficient. Computes totals and averages without storing per-record data.
    """
    compression, handle = open_text_handle(content)

    file_format = detect_file_format(filename)
    if file_format is None:
        # sniff FASTA/FASTQ only
        prefix = handle.read(2048)
        file_format = sniff_fasta_fastq(prefix)
        handle.seek(0)

    if file_format is None:
        raise ValueError("Could not detect file format from filename/sniffing.")

    total_bases = 0
    total_gc_percent = 0.0
    count = 0

    if file_format == "fastq":
        for title, seq, qual in FastqGeneralIterator(handle):
            L = len(seq)
            total_bases += L
            gc = gc_fraction(seq, mode=mode)
            total_gc_percent += (0.0 if gc is None else 100 * gc)
            count += 1
    else:
        for record in SeqIO.parse(handle, file_format):
            seq = str(record.seq)
            L = len(seq)
            total_bases += L
            gc = gc_fraction(seq, mode=mode)
            total_gc_percent += (0.0 if gc is None else 100 * gc)
            count += 1

    avg_len = (total_bases / count) if count else 0.0
    avg_gc = (total_gc_percent / count) if count else 0.0

    return {
        "filename": filename,
        "format": file_format,
        "compression": compression,
        "total_sequences": count,
        "total_bases": total_bases,
        "average_length": round(avg_len, 2),
        "average_gc_content": round(avg_gc, 2),
    }


def filter_fastq(content: bytes, filename: str, filter_ids: Optional[str] = None,
                mode: str = "canonical", max_records: int = 200) -> Dict:
    """
    FASTQ-specific processing with optional ID filtering.
    """
    file_format = detect_file_format(filename)
    if file_format != "fastq":
        raise ValueError("This function only accepts FASTQ files")

    wanted_ids = _parse_id_list(filter_ids)
    compression, handle = open_text_handle(content)

    sequences = process_fastq_stream(handle, mode=mode, max_records=max_records, wanted_ids=wanted_ids)

    return {
        "filename": filename,
        "format": file_format,
        "compression": compression,
        "total_sequences": len(sequences),
        "filtered": wanted_ids is not None,
        "filter_count": len(wanted_ids) if wanted_ids else 0,
        "sequences": sequences,
    }


def process_sequences_universal(content: bytes, filename: str, mode: str = "canonical", max_records: int = 200) -> Dict:
    """
    Universal processing: FASTA/FASTQ/GenBank/EMBL (+ gz/bz2) -> preview records list.
    """
    compression, handle = open_text_handle(content)

    file_format = detect_file_format(filename)
    if file_format is None:
        prefix = handle.read(2048)
        file_format = sniff_fasta_fastq(prefix)
        handle.seek(0)

    if file_format is None:
        raise ValueError("Could not detect file format. Use correct extension.")

    if file_format == "fastq":
        sequences = process_fastq_stream(handle, mode=mode, max_records=max_records)
    else:
        sequences = process_fasta_like_stream(handle, file_format=file_format, mode=mode, max_records=max_records)

    return {
        "filename": filename,
        "format": file_format,
        "compression": compression,
        "total_sequences": len(sequences),
        "total_bases": sum(s["Length"] for s in sequences),
        "sequences": sequences,
    }
