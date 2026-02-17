# format_detect.py
from typing import Optional

FASTA_EXT = {".fa", ".fasta"}
FASTQ_EXT = {".fq", ".fastq"}
GENBANK_EXT = {".gb", ".gbk"}
EMBL_EXT = {".embl"}


def detect_file_format(filename: str) -> Optional[str]:
    """Return 'fasta'/'fastq'/'genbank'/'embl' or None."""
    name = filename.lower()
    for ext in FASTA_EXT:
        if name.endswith(ext):
            return "fasta"
    for ext in FASTQ_EXT:
        if name.endswith(ext):
            return "fastq"
    for ext in GENBANK_EXT:
        if name.endswith(ext):
            return "genbank"
    for ext in EMBL_EXT:
        if name.endswith(ext):
            return "embl"
    return None


def sniff_fasta_fastq(prefix: str) -> Optional[str]:
    """Sniff only FASTA/FASTQ from content prefix."""
    s = prefix.lstrip()
    if s.startswith(">"):
        return "fasta"
    if s.startswith("@"):
        return "fastq"
    return None
