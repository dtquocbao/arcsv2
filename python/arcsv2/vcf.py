"""VCF writer for ARC-SV2 call tables (U1).

Default call output is ``arcsv_out.vcf.gz`` + ``.tbi`` beside ``arcsv_out.tab``.
VCF 4.2 conventions: simple ``DEL|DUP|INV|INS``, complex ``CPX`` (gnomAD-SV style)
plus optional BND mate pairs per novel adjacency.
"""

from __future__ import annotations

import re
from pathlib import Path

_TOKEN = re.compile(r"([A-Za-z]\d*|_)(')?")


def _svtype(row: dict[str, str]) -> str:
    sv = (row.get("svtype") or "").split(",")[0].upper()
    if sv in ("DEL", "DUP", "INV", "INS"):
        return sv
    return "CPX"


def _gt(row: dict[str, str]) -> str:
    if row.get("gt") == "HOM":
        return "1/1"
    return "0/1"


def _sum_supports(text: str | None) -> int:
    total = 0
    for part in (text or "0").split(","):
        part = part.strip()
        if not part or part.upper() == "NA":
            continue
        try:
            total += int(float(part))
        except ValueError:
            continue
    return total


def _is_imprecise(row: dict[str, str]) -> bool:
    """True when every internal junction is midpoint (no clip/split stamp)."""
    sources = [p.strip() for p in (row.get("bp_point_source") or "").split(",") if p.strip()]
    if not sources:
        return False
    return all(src == "mid" for src in sources)


def _ci_pair(uncertainty: str | None, index: int) -> tuple[int, int]:
    """Return (CIPOS/CIEND minus, plus) from a comma-separated uncertainty list."""
    parts = [p.strip() for p in (uncertainty or "").split(",") if p.strip() and p.strip().upper() != "NA"]
    if not parts:
        return 0, 0
    try:
        val = int(float(parts[min(index, len(parts) - 1)]))
    except ValueError:
        return 0, 0
    # Symmetric CI from the recorded half-width when only one number is present.
    return -abs(val), abs(val)


def vcf_header(sample: str, contigs: list[tuple[str, int]], model: str | None = None) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=arcsv2",
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation (DRP/midpoint breakpoints)">',
        '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="ARC-SV complex class (e.g. ddup.del)">',
        '##INFO=<ID=REARR,Number=1,Type=String,Description="Rearrangement block string">',
        '##INFO=<ID=BLOCKS,Number=1,Type=String,Description="Reference block intervals">',
        '##INFO=<ID=AFFECTED_LEN,Number=1,Type=Integer,Description="Affected length (bp)">',
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">',
        '##INFO=<ID=EVENT,Number=1,Type=String,Description="Event identifier for linked BND/CPX records">',
        '##INFO=<ID=MLPROB,Number=1,Type=Float,Description="Calibrated ML probability">',
        '##INFO=<ID=MLCLASS,Number=1,Type=String,Description="ML class">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Split-read support (sum across junctions)">',
        '##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Paired-end support (sum across junctions)">',
        '##FORMAT=<ID=SVSCORE,Number=1,Type=Float,Description="Log LR versus the reference genotype">',
        '##FORMAT=<ID=SVSCORE_NEXT,Number=1,Type=Float,Description="Log LR versus the next genotype">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=MLFAIL,Description="Failed the ML filter">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=CPX,Description="Complex structural variant">',
    ]
    if model:
        lines.append(f"##arcsv2_ml_model={model}")
    for name, length in contigs:
        if length and length > 0:
            lines.append(f"##contig=<ID={name},length={int(length)}>")
        else:
            lines.append(f"##contig=<ID={name}>")
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")
    return "\n".join(lines) + "\n"


def _info(row: dict[str, str], svtype: str) -> str:
    affected = int(float(row.get("len_affected") or 0))
    if svtype == "DEL":
        svlen = -affected
    elif svtype == "INS":
        inslen = row.get("inslen") or str(affected)
        try:
            svlen = int(str(inslen).lstrip(">"))
        except ValueError:
            svlen = affected
    else:
        svlen = affected
    end = int(row.get("maxbp") or row.get("minbp") or 0)
    cipos = _ci_pair(row.get("bp_uncertainty"), 0)
    ciend = _ci_pair(row.get("bp_uncertainty"), -1)
    fields = [
        f"END={end}",
        f"SVTYPE={svtype}",
        f"SVLEN={svlen}",
        f"CIPOS={cipos[0]},{cipos[1]}",
        f"CIEND={ciend[0]},{ciend[1]}",
        f"AFFECTED_LEN={affected}",
        f"REARR={row.get('rearrangement') or '.'}",
        f"BLOCKS={row.get('blocks') or row.get('blocks_point') or '.'}",
    ]
    if svtype == "CPX":
        fields.append(f"CPX_TYPE={row.get('complextype') or 'other'}")
    if _is_imprecise(row):
        fields.append("IMPRECISE")
    if row.get("ml_prob") not in (None, "", "NA"):
        fields.append(f"MLPROB={row['ml_prob']}")
        fields.append(f"MLCLASS={row.get('ml_class', svtype)}")
    return ";".join(fields)


def _sample_format(row: dict[str, str]) -> str:
    sr = _sum_supports(row.get("sr_support"))
    pe = _sum_supports(row.get("pe_support"))
    score = row.get("score_vs_ref") or "0"
    score_next = row.get("score_vs_next") or "0"
    return f"{_gt(row)}:{sr}:{pe}:{score}:{score_next}"


def _rank(name: str) -> int | None:
    """Inverse of block_idx_to_name for a locus rooted at A."""
    if not name or name == "_":
        return None
    letter = name[0]
    rest = name[1:]
    mult = int(rest) if rest else 0
    mod = ord(letter) - ord("A")
    if letter.islower():
        mod -= 6
    if mod < 0 or mod >= 52:
        return None
    return mult * 52 + mod


def _rearrangement_tokens(text: str) -> list[tuple[str, int]]:
    return [(m.group(1), -1 if m.group(2) else 1) for m in _TOKEN.finditer(text or "")]


def _novel_orients(rearrangement: str) -> list[tuple[str, str]]:
    """Breakend orientations for each novel adjacency (orig:sv_classify.py)."""
    tokens = _rearrangement_tokens(rearrangement)
    orients: list[tuple[str, str]] = []
    for (name1, s1), (name2, s2) in zip(tokens, tokens[1:]):
        r1, r2 = _rank(name1), _rank(name2)
        if s1 > 0 and s2 > 0 and r1 is not None and r2 == r1 + 1:
            continue
        left = "+" if s1 < 0 else "-"
        right = "-" if s2 < 0 else "+"
        orients.append((left, right))
    return orients


def bnd_alt_string(orient: str, other_orient: str, chrom: str, other_pos: int) -> str:
    """orig:sv_output.py::bnd_alt_string. ``other_pos`` is the mate's VCF position."""
    alt_after = orient == "-"
    bracket = "]" if other_orient == "-" else "["
    loc = f"{bracket}{chrom}:{other_pos}{bracket}"
    if alt_after:
        return "N" + loc
    return loc + "N"


def _vcf_pos(bp0: int, orient: str) -> int:
    """1-based VCF position. A '-' breakend stays on the 0-based coordinate."""
    pos = bp0 + 1
    if orient == "-":
        return pos - 1
    return pos


def format_record(row: dict[str, str], event_id: str) -> list[str]:
    """Return the CPX/simple record plus one BND pair per novel adjacency."""
    svtype = _svtype(row)
    chrom = row["chrom"]
    pos = int(row["minbp"]) + 1
    filt = "MLFAIL" if row.get("ml_pass") == "0" else "PASS"
    alt = f"<{svtype}>"
    qual = row.get("score_vs_ref", ".")
    sample = _sample_format(row)
    main = "\t".join(
        [
            chrom,
            str(pos),
            event_id,
            "N",
            alt,
            str(qual),
            filt,
            _info(row, svtype),
            "GT:SR:PE:SVSCORE:SVSCORE_NEXT",
            sample,
        ]
    )
    records = [main]
    # Breakend mates for novel adjacencies (orig-compatible; also covers simple DEL/INV).
    junctions = (row.get("sv_bp") or "").split(";")
    orients = _novel_orients(row.get("rearrangement") or "")
    mate_ids = []
    for i, junction in enumerate(junctions):
        if not junction or junction == "NA" or "," not in junction:
            continue
        left_s, right_s = junction.split(",")[:2]
        orient = orients[i] if i < len(orients) else ("-", "+")
        mate_ids.append((int(left_s), int(right_s), orient, f"{event_id}_bnd{i}"))
    for left, right, (o1, o2), bid in mate_ids:
        pos1 = _vcf_pos(left, o1)
        pos2 = _vcf_pos(right, o2)
        alt1 = bnd_alt_string(o1, o2, chrom, pos2)
        alt2 = bnd_alt_string(o2, o1, chrom, pos1)
        info = f"SVTYPE=BND;MATEID={bid}b;EVENT={event_id}"
        info_b = f"SVTYPE=BND;MATEID={bid};EVENT={event_id}"
        records.append(
            "\t".join(
                [chrom, str(pos1), bid, "N", alt1, ".", filt, info, "GT:SR:PE:SVSCORE:SVSCORE_NEXT", sample]
            )
        )
        records.append(
            "\t".join(
                [
                    chrom,
                    str(pos2),
                    f"{bid}b",
                    "N",
                    alt2,
                    ".",
                    filt,
                    info_b,
                    "GT:SR:PE:SVSCORE:SVSCORE_NEXT",
                    sample,
                ]
            )
        )
    return records


def _vcf_sort_key(line: str) -> tuple[str, int]:
    parts = line.split("\t")
    return parts[0], int(parts[1])


def write_vcf(
    rows: list[dict[str, str]],
    path: str | Path,
    sample: str,
    contigs: list[tuple[str, int]],
    model: str | None = None,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    body: list[str] = []
    for i, row in enumerate(rows):
        body.extend(format_record(row, row.get("id") or f"arcsv2_{i}"))
    body.sort(key=_vcf_sort_key)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(vcf_header(sample, contigs, model))
        for line in body:
            handle.write(line + "\n")


def write_vcf_gz(
    rows: list[dict[str, str]],
    path_gz: str | Path,
    sample: str,
    contigs: list[tuple[str, int]],
    model: str | None = None,
) -> Path:
    """Write ``*.vcf.gz`` + ``*.tbi``. Returns the gzip path."""
    import pysam

    path_gz = Path(path_gz)
    if not str(path_gz).endswith(".gz"):
        path_gz = Path(str(path_gz) + ".gz")
    plain_name = path_gz.name[:-3] if path_gz.name.endswith(".gz") else path_gz.name
    if not plain_name.endswith(".vcf"):
        plain_name = plain_name + ".vcf"
    plain = path_gz.with_name(plain_name)
    write_vcf(rows, plain, sample, contigs, model=model)
    pysam.tabix_compress(str(plain), str(path_gz), force=True)
    pysam.tabix_index(str(path_gz), preset="vcf", force=True)
    plain.unlink(missing_ok=True)
    return path_gz


def contigs_from_fai(reference: str | Path | None, fallback: list[tuple[str, int]] | None = None) -> list[tuple[str, int]]:
    """Read ``reference.fai`` contig lengths when present."""
    if reference:
        fai = Path(str(reference) + ".fai")
        if fai.is_file():
            out: list[tuple[str, int]] = []
            for line in fai.read_text(encoding="utf-8").splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    out.append((parts[0], int(parts[1])))
            if out:
                return out
    return list(fallback or [])


def sample_from_bam(bam_path: str | Path | None, default: str = "SAMPLE") -> str:
    """Best-effort ``@RG SM`` from a BAM/CRAM/SAM header."""
    if not bam_path:
        return default
    path = Path(bam_path)
    if not path.is_file():
        return default
    try:
        import pysam
    except ImportError:
        return default
    try:
        with pysam.AlignmentFile(str(path), "r") as handle:
            for rg in handle.header.to_dict().get("RG", []):
                sm = rg.get("SM")
                if sm:
                    return str(sm)
    except (OSError, ValueError):
        return default
    return default


def rows_from_call_rows(rows: list) -> list[dict[str, str]]:
    """Convert ``CallRow`` objects or dicts into string dicts for the writer."""
    out: list[dict[str, str]] = []
    for row in rows:
        if hasattr(row, "as_dict"):
            out.append({k: str(v) for k, v in row.as_dict().items()})
        else:
            out.append({k: str(v) for k, v in dict(row).items()})
    return out
