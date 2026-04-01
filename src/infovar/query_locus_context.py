#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import pysam


# =========================
# Helpers
# =========================

NCBI_GRCH38_CONTIG_MAP = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "M": "NC_012920.1",
    "MT": "NC_012920.1",
}

def json_default(obj: Any):
    if isinstance(obj, set):
        return sorted(obj)
    if isinstance(obj, bytes):
        return obj.decode()
    if hasattr(obj, "tolist"):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def parse_gff3_attributes(attr_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def normalize_contig_for_dbsnp_ncbi(contig: str, available: list[str]) -> str:
    c = str(contig).strip()
    if c.startswith("chr"):
        c = c[3:]

    c = c.upper()
    target = NCBI_GRCH38_CONTIG_MAP.get(c)
    if target is None:
        raise ValueError(f"Contig no soportado para dbSNP NCBI: {contig}")

    if target not in set(available):
        raise ValueError(
            f"Contig mapeado '{target}' no encontrado en VCF dbSNP"
        )
    return target


def normalize_contig_for_target(contig: str, available: list[str]) -> str:
    avail = set(available)
    if contig in avail:
        return contig

    if contig.startswith("chr"):
        alt = contig[3:]
        if alt in avail:
            return alt
    else:
        alt = f"chr{contig}"
        if alt in avail:
            return alt

    raise ValueError(f"Contig '{contig}' no encontrado en el archivo objetivo")


def make_context_string(seq: str, center_index: int) -> str:
    top = " ".join(seq)
    mid = []
    for i, _ in enumerate(seq):
        mid.append("^" if i == center_index else " ")
    return top + "\n" + " ".join(mid)


def feature_label(attrs: dict[str, str], fallback_type: str) -> str:
    for key in ("gene_name", "gene", "Name", "ID", "transcript_name", "transcript_id"):
        if key in attrs:
            return attrs[key]
    return fallback_type


def variant_to_jsonable(rec: pysam.VariantRecord) -> dict[str, Any]:
    info = {}
    for k in rec.info.keys():
        v = rec.info[k]
        if isinstance(v, tuple):
            info[k] = list(v)
        else:
            info[k] = v

    return {
        "contig": rec.contig,
        "pos": rec.pos,
        "id": rec.id,
        "ref": rec.ref,
        "alts": list(rec.alts or []),
        "qual": rec.qual,
        "filter": list(rec.filter.keys()) if rec.filter else [],
        "info": info,
    }


def pick_frequency_fields(info: dict[str, Any]) -> dict[str, Any]:
    """
    Extrae campos candidatos a frecuencias poblacionales de forma heurística.
    No asume un schema específico del VCF.
    """
    out = {}
    for k, v in info.items():
        ku = k.upper()
        if any(token in ku for token in ("AF", "FREQ", "MAF", "COMMON", "TOPMED", "GNOMAD", "1000G")):
            out[k] = v
    return out


@dataclass
class Feature:
    contig: str
    start: int
    end: int
    strand: str
    type: str
    source: str
    phase: str | None
    attributes: dict[str, str]

    @property
    def label(self) -> str:
        return feature_label(self.attributes, self.type)


# =========================
# Queries
# =========================

def fetch_reference_context(
    fasta_path: str,
    chrom: str,
    pos1: int,
    flank: int,
) -> dict[str, Any]:
    fasta = pysam.FastaFile(fasta_path)
    try:
        contig = normalize_contig_for_target(chrom, list(fasta.references))
        start0 = max(0, pos1 - 1 - flank)
        end0 = pos1 + flank
        seq = fasta.fetch(contig, start0, end0).upper()
        if not seq:
            raise ValueError("No se pudo recuperar secuencia del FASTA")

        center_idx = pos1 - 1 - start0
        if center_idx < 0 or center_idx >= len(seq):
            raise ValueError("La posición solicitada quedó fuera del contexto recuperado")

        ref_base = seq[center_idx]
        return {
            "contig": contig,
            "pos": pos1,
            "window_start_1based": start0 + 1,
            "window_end_1based": end0,
            "flank": flank,
            "sequence": seq,
            "center_index_0based_in_window": center_idx,
            "ref_base": ref_base,
            "context_pretty": make_context_string(seq, center_idx),
        }
    finally:
        fasta.close()


def fetch_gff3_features(
    gff3_path: str,
    chrom: str,
    pos1: int,
    flank: int,
) -> list[Feature]:
    tbx = pysam.TabixFile(gff3_path)
    try:
        contig = normalize_contig_for_target(chrom, list(tbx.contigs))
        start0 = max(0, pos1 - 1 - flank)
        end0 = pos1 + flank

        features: list[Feature] = []
        try:
            for line in tbx.fetch(contig, start0, end0):
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue

                seqid, source, ftype, start, end, score, strand, phase, attrs = fields
                feat = Feature(
                    contig=seqid,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    type=ftype,
                    source=source,
                    phase=None if phase == "." else phase,
                    attributes=parse_gff3_attributes(attrs),
                )
                features.append(feat)
        except ValueError:
            pass

        return sorted(features, key=lambda x: (x.start, x.end, x.type))
    finally:
        tbx.close()


def fetch_dbsnp_records(
    vcf_path: str,
    chrom: str,
    pos1: int,
    flank: int,
) -> list[dict[str, Any]]:
    vcf = pysam.VariantFile(vcf_path)
    try:
        contig = normalize_contig_for_dbsnp_ncbi(chrom, list(vcf.header.contigs))
        start0 = max(0, pos1 - 1 - flank)
        end0 = pos1 + flank

        out = []
        for rec in vcf.fetch(contig, start0, end0):
            j = variant_to_jsonable(rec)
            j["population_frequency_fields"] = pick_frequency_fields(j["info"])
            j["is_target_position"] = (rec.pos == pos1)
            out.append(j)
        return out
    finally:
        vcf.close()


# =========================
# SVG rendering
# =========================

def render_svg(
    result: dict[str, Any],
    width: int = 1400,
    track_height: int = 28,
    margin_left: int = 80,
    margin_right: int = 30,
    margin_top: int = 30,
    margin_bottom: int = 40,
) -> str:
    ref = result["reference"]
    feats = result["features"]
    vars_ = result["dbsnp_variants"]

    seq = ref["sequence"]
    flank = ref["flank"]
    target_pos = ref["pos"]
    window_start = ref["window_start_1based"]
    window_end = ref["window_end_1based"]
    span = window_end - window_start + 1

    feature_tracks = []
    for feat in feats:
        placed = False
        for track in feature_tracks:
            if feat.start > track[-1].end:
                track.append(feat)
                placed = True
                break
        if not placed:
            feature_tracks.append([feat])

    n_tracks = max(1, len(feature_tracks))
    height = margin_top + 120 + n_tracks * track_height + 120 + margin_bottom
    plot_x0 = margin_left
    plot_x1 = width - margin_right
    def xmap(pos1: int) -> float:
        rel = (pos1 - window_start) / max(1, span - 1)
        return plot_x0 + rel * (plot_x1 - plot_x0)

    y_axis = margin_top + 30
    y_seq = margin_top + 70
    y_tracks0 = margin_top + 110
    y_dbsnp = y_tracks0 + n_tracks * track_height + 40

    parts = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="monospace">'
    )
    parts.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>')
    parts.append(
        f'<text x="{margin_left}" y="20" font-size="18" font-weight="bold">'
        f'{ref["contig"]}:{target_pos} ref={ref["ref_base"]}</text>'
    )

    # Axis
    parts.append(f'<line x1="{plot_x0}" y1="{y_axis}" x2="{plot_x1}" y2="{y_axis}" stroke="black" stroke-width="1"/>')
    for p in range(window_start, window_end + 1):
        x = xmap(p)
        tick_h = 10 if p == target_pos else 6
        parts.append(f'<line x1="{x}" y1="{y_axis}" x2="{x}" y2="{y_axis - tick_h}" stroke="black" stroke-width="1"/>')
        if p == window_start or p == window_end or p == target_pos or (p - window_start) % max(1, span // 10) == 0:
            parts.append(f'<text x="{x}" y="{y_axis - 14}" font-size="11" text-anchor="middle">{p}</text>')

    # Sequence
    for i, base in enumerate(seq):
        p = window_start + i
        x = xmap(p)
        color = "crimson" if p == target_pos else "black"
        weight = "bold" if p == target_pos else "normal"
        parts.append(
            f'<text x="{x}" y="{y_seq}" font-size="16" text-anchor="middle" fill="{color}" font-weight="{weight}">{base}</text>'
        )

    # Feature tracks
    for tidx, track in enumerate(feature_tracks):
        y = y_tracks0 + tidx * track_height
        parts.append(f'<text x="10" y="{y+12}" font-size="11">track{tidx+1}</text>')
        for feat in track:
            x0 = xmap(max(feat.start, window_start))
            x1 = xmap(min(feat.end, window_end))
            w = max(2.0, x1 - x0)
            fill = {
                "gene": "#8dd3c7",
                "transcript": "#80b1d3",
                "exon": "#bebada",
                "CDS": "#fb8072",
                "five_prime_UTR": "#fdb462",
                "three_prime_UTR": "#b3de69",
            }.get(feat.type, "#d9d9d9")
            parts.append(
                f'<rect x="{x0}" y="{y}" width="{w}" height="14" fill="{fill}" stroke="black" stroke-width="0.5"/>'
            )
            parts.append(
                f'<text x="{(x0+x1)/2}" y="{y+11}" font-size="10" text-anchor="middle">{feat.label}</text>'
            )

    # Target position marker
    tx = xmap(target_pos)
    parts.append(f'<line x1="{tx}" y1="{margin_top}" x2="{tx}" y2="{height - margin_bottom}" stroke="red" stroke-dasharray="4,4"/>')

    # dbSNP variants
    parts.append(f'<text x="{margin_left}" y="{y_dbsnp - 10}" font-size="13" font-weight="bold">dbSNP variants</text>')
    if vars_:
        y = y_dbsnp + 10
        for rec in vars_[:12]:
            label = f'{rec["pos"]} {rec["ref"]}>{",".join(rec["alts"])}'
            if rec.get("id"):
                label = f'{rec["id"]} {label}'
            freq = rec.get("population_frequency_fields", {})
            if freq:
                label += " " + json.dumps(freq, ensure_ascii=False)
            color = "crimson" if rec["pos"] == target_pos else "black"
            parts.append(f'<text x="{margin_left}" y="{y}" font-size="11" fill="{color}">{label}</text>')
            y += 14
    else:
        parts.append(f'<text x="{margin_left}" y="{y_dbsnp + 10}" font-size="11">No dbSNP variants in window</text>')

    parts.append("</svg>")
    return "\n".join(parts)


# =========================
# Main
# =========================

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Consulta referencia, GFF3 y dbSNP para una posición CHROM:POS"
    )
    ap.add_argument("--chrom", required=True, help="Contig, e.g. chr1 o 1")
    ap.add_argument("--pos", required=True, type=int, help="Posición 1-based")
    ap.add_argument("--fasta", required=True, help="FASTA de referencia")
    ap.add_argument("--gff3", required=True, help="GFF3.gz indexado con tabix")
    ap.add_argument("--dbsnp-vcf", required=True, help="VCF.gz indexado con tabix")
    ap.add_argument("--flank", type=int, default=10, help="Bases adyacentes a cada lado")
    ap.add_argument(
        "--output-format",
        choices=["json", "svg"],
        default="json",
        help="json -> STDOUT, svg -> requiere --output o STDOUT",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="-",
        help="Archivo de salida o '-' para STDOUT",
    )

    args = ap.parse_args()

    ref = fetch_reference_context(
        fasta_path=args.fasta,
        chrom=args.chrom,
        pos1=args.pos,
        flank=args.flank,
    )

    feats = fetch_gff3_features(
        gff3_path=args.gff3,
        chrom=ref["contig"],
        pos1=args.pos,
        flank=args.flank,
    )

    dbsnp = fetch_dbsnp_records(
        vcf_path=args.dbsnp_vcf,
        chrom=ref["contig"],
        pos1=args.pos,
        flank=args.flank,
    )

    result = {
        "query": {
            "chrom": args.chrom,
            "pos": args.pos,
            "flank": args.flank,
        },
        "reference": ref,
        "features": feats,
        "dbsnp_variants": dbsnp,
    }

    if args.output_format == "json":
        json_result = {
            **result,
            "features": [asdict(f) for f in feats],
        }
        payload = json.dumps(json_result, ensure_ascii=False, indent=2, default=json_default)
    else:
        payload = render_svg(result)

    if args.output == "-":
        sys.stdout.write(payload)
        if not payload.endswith("\n"):
            sys.stdout.write("\n")
    else:
        Path(args.output).write_text(payload, encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
