from typing import Iterator, List
import gzip

def open_maybe_gzip(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')

def variant_iter(path: str) -> Iterator[List[str]]:
    """Yield split fields for each variant line; skip headers; be tolerant."""
    with open_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                yield None
                continue
            yield parts

def is_snp(ref: str, alt: str) -> bool:
    alts = alt.split(',') if ',' in alt else [alt]
    return len(ref) == 1 and all(len(a) == 1 for a in alts)

def is_indel(ref: str, alt: str) -> bool:
    alts = alt.split(',') if ',' in alt else [alt]
    return any(len(a) != len(ref) for a in alts)

def transitions_and_transversions(ref: str, alt: str):
    purines = {'A','G'}; pyrimidines = {'C','T'}
    if ',' in alt:
        alt = alt.split(',')[0]
    if len(ref) != 1 or len(alt) != 1:
        return 0, 0
    if (ref in purines and alt in purines) or (ref in pyrimidines and alt in pyrimidines):
        return 1, 0
    return 0, 1
