#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field


APP_ROOT = Path(__file__).resolve().parent
QUERY_SCRIPT = APP_ROOT / "query_locus_context.py"

DEFAULT_FASTA = "/mnt/cephfs/hot_nvme/hg38/ref/hg38.fa"
DEFAULT_GFF3 = "/mnt/cephfs/hot_nvme/gencode/GRCh38.14/gencode.v49.basic.annotation.gff3.gz"
DEFAULT_DBSNP_VCF = "/mnt/cephfs/hot_nvme/dbsnp/population_frequency/freq.vcf.gz"

BROWSE_ROOT = Path("/mnt/cephfs/hot_nvme").resolve()
TMP_DIR = Path(tempfile.gettempdir()) / "infovar_api"
TMP_DIR.mkdir(parents=True, exist_ok=True)


app = FastAPI(title="InfoVar API", version="0.1.0")


class LocusRequest(BaseModel):
    chrom: str = Field(..., examples=["chr1", "1", "chrX"])
    pos: int = Field(..., ge=1)
    ref: str | None = None
    alt: str | None = None
    flank: int = Field(10, ge=0, le=500)
    fasta: str = DEFAULT_FASTA
    gff3: str = DEFAULT_GFF3
    dbsnp_vcf: str = DEFAULT_DBSNP_VCF


def safe_path_under_root(path_str: str, root: Path = BROWSE_ROOT) -> Path:
    path = Path(path_str).resolve()
    if not str(path).startswith(str(root)):
        raise HTTPException(status_code=403, detail=f"Path fuera de root permitido: {root}")
    return path


def run_query_locus_context(
    *,
    chrom: str,
    pos: int,
    flank: int,
    fasta: str,
    gff3: str,
    dbsnp_vcf: str,
    output_format: str,
    output_path: Path | None = None,
) -> dict[str, Any]:
    cmd = [
        "python",
        str(QUERY_SCRIPT),
        "--chrom",
        chrom,
        "--pos",
        str(pos),
        "--fasta",
        fasta,
        "--gff3",
        gff3,
        "--dbsnp-vcf",
        dbsnp_vcf,
        "--flank",
        str(flank),
        "--output-format",
        output_format,
    ]

    if output_path is not None:
        cmd.extend(["-o", str(output_path)])

    try:
        cp = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise HTTPException(
            status_code=500,
            detail={
                "cmd": cmd,
                "returncode": e.returncode,
                "stdout": e.stdout,
                "stderr": e.stderr,
            },
        ) from e

    return {
        "cmd": cmd,
        "stdout": cp.stdout,
        "stderr": cp.stderr,
    }


def run_query_json(
    chrom: str,
    pos: int,
    flank: int,
    fasta: str,
    gff3: str,
    dbsnp_vcf: str,
) -> dict[str, Any]:
    result = run_query_locus_context(
        chrom=chrom,
        pos=pos,
        flank=flank,
        fasta=fasta,
        gff3=gff3,
        dbsnp_vcf=dbsnp_vcf,
        output_format="json",
        output_path=None,
    )
    try:
        return json.loads(result["stdout"])
    except json.JSONDecodeError as e:
        raise HTTPException(
            status_code=500,
            detail={"message": "JSON inválido desde query_locus_context.py", **result},
        ) from e


def build_svg_path(chrom: str, pos: int) -> Path:
    safe_chrom = "".join(c for c in chrom if c.isalnum() or c in {"_", "-"})
    return TMP_DIR / f"{safe_chrom}_{pos}.svg"


def parse_positions_upload(file_path: Path) -> list[dict[str, Any]]:
    def _open(p: Path):
        if p.suffix == ".gz":
            return gzip.open(p, "rt", encoding="utf-8", newline="")
        return open(p, "rt", encoding="utf-8", newline="")

    rows: list[dict[str, Any]] = []
    with _open(file_path) as fh:
        sample = fh.read(4096)
        fh.seek(0)

        lines = []
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            lines.append(line)
            if len(lines) >= 5:
                break

        if not lines:
            return []

        first = lines[0]
        delimiter = "\t" if first.count("\t") >= first.count(",") else ","

    with _open(file_path) as fh:
        reader = csv.DictReader(
            (line for line in fh if line.strip() and not line.startswith("#")),
            delimiter=delimiter,
        )

        norm_map = {name.lower(): name for name in (reader.fieldnames or [])}

        chrom_col = None
        pos_col = None
        ref_col = None
        alt_col = None

        for cand in ("chrom", "chr", "chromosome", "chr_name"):
            if cand in norm_map:
                chrom_col = norm_map[cand]
                break

        for cand in ("pos", "position", "chr_position", "bp"):
            if cand in norm_map:
                pos_col = norm_map[cand]
                break

        for cand in ("ref", "reference_allele", "other_allele"):
            if cand in norm_map:
                ref_col = norm_map[cand]
                break

        for cand in ("alt", "alternate_allele", "effect_allele"):
            if cand in norm_map:
                alt_col = norm_map[cand]
                break

        if chrom_col is None or pos_col is None:
            raise HTTPException(
                status_code=400,
                detail={
                    "message": "No se pudieron identificar columnas CHROM/POS",
                    "columns": reader.fieldnames,
                },
            )

        for i, row in enumerate(reader, start=1):
            chrom = str(row[chrom_col]).strip()
            pos_raw = str(row[pos_col]).strip()
            if not chrom or not pos_raw:
                continue
            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            rows.append(
                {
                    "rownum": i,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": row.get(ref_col) if ref_col else None,
                    "alt": row.get(alt_col) if alt_col else None,
                }
            )

    return rows


@app.get("/health")
def health():
    return {"ok": True}


@app.get("/chrompos/{chrom}/{pos}")
def chrompos(
    chrom: str,
    pos: int,
    ref: str | None = Query(default=None),
    alt: str | None = Query(default=None),
    flank: int = Query(default=10, ge=0, le=500),
    fasta: str = Query(default=DEFAULT_FASTA),
    gff3: str = Query(default=DEFAULT_GFF3),
    dbsnp_vcf: str = Query(default=DEFAULT_DBSNP_VCF),
):
    data = run_query_json(
        chrom=chrom,
        pos=pos,
        flank=flank,
        fasta=fasta,
        gff3=gff3,
        dbsnp_vcf=dbsnp_vcf,
    )
    data["request"] = {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "flank": flank,
        "fasta": fasta,
        "gff3": gff3,
        "dbsnp_vcf": dbsnp_vcf,
    }
    return JSONResponse(data)


@app.get("/chrompos/{chrom}/{pos}/svg")
def chrompos_svg(
    chrom: str,
    pos: int,
    flank: int = Query(default=10, ge=0, le=500),
    fasta: str = Query(default=DEFAULT_FASTA),
    gff3: str = Query(default=DEFAULT_GFF3),
    dbsnp_vcf: str = Query(default=DEFAULT_DBSNP_VCF),
):
    svg_path = build_svg_path(chrom, pos)

    run_query_locus_context(
        chrom=chrom,
        pos=pos,
        flank=flank,
        fasta=fasta,
        gff3=gff3,
        dbsnp_vcf=dbsnp_vcf,
        output_format="svg",
        output_path=svg_path,
    )

    return FileResponse(svg_path, media_type="image/svg+xml", filename=svg_path.name)


@app.post("/positions/upload")
async def positions_upload(
    file: UploadFile = File(...),
):
    suffix = Path(file.filename or "upload.tsv").suffix
    tmp_path = TMP_DIR / f"upload_{os.getpid()}_{file.filename or 'positions.tsv'}"
    with open(tmp_path, "wb") as out:
        out.write(await file.read())

    rows = parse_positions_upload(tmp_path)
    return {
        "filename": file.filename,
        "n_rows": len(rows),
        "rows": rows,
    }


@app.get("/fs/list")
def fs_list(path: str = Query(default=str(BROWSE_ROOT))):
    p = safe_path_under_root(path)

    if not p.exists():
        raise HTTPException(status_code=404, detail="Path no existe")
    if not p.is_dir():
        raise HTTPException(status_code=400, detail="Path no es directorio")

    items = []
    for child in sorted(p.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
        items.append(
            {
                "name": child.name,
                "path": str(child),
                "is_dir": child.is_dir(),
                "size": child.stat().st_size if child.is_file() else None,
            }
        )

    return {
        "root": str(BROWSE_ROOT),
        "path": str(p),
        "items": items,
    }
