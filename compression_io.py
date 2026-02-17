# compression_io.py
import bz2
import gzip
import io
from typing import Tuple, TextIO

GZIP_MAGIC = b"\x1f\x8b"
BZIP2_MAGIC = b"BZh"


def detect_compression_from_bytes(data: bytes) -> str:
    """Return 'gzip', 'bzip2', or 'none' using magic bytes."""
    if data.startswith(GZIP_MAGIC):
        return "gzip"
    if data.startswith(BZIP2_MAGIC):
        return "bzip2"
    return "none"


def open_text_handle(data: bytes, encoding: str = "utf-8") -> Tuple[str, TextIO]:
    """
    Returns (compression, text_handle).
    text_handle is a decoded stream suitable for Bio.SeqIO and FastqGeneralIterator.
    """
    comp = detect_compression_from_bytes(data)
    raw = io.BytesIO(data)

    if comp == "gzip":
        bio = gzip.GzipFile(fileobj=raw, mode="rb")
        return comp, io.TextIOWrapper(bio, encoding=encoding, errors="replace")

    if comp == "bzip2":
        bio = bz2.BZ2File(raw, mode="rb")
        return comp, io.TextIOWrapper(bio, encoding=encoding, errors="replace")

    return comp, io.TextIOWrapper(raw, encoding=encoding, errors="replace")
