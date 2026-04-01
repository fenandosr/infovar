#!/usr/bin/env python3
from __future__ import annotations

import html
import json
from urllib.parse import urlencode

import requests
from bottle import Bottle, request, response, run


API_BASE = "http://127.0.0.1:8042"
app = Bottle()


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def page(title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="es">
<head>
  <meta charset="utf-8">
  <title>{esc(title)}</title>
  <style>
    body {{
      font-family: sans-serif;
      margin: 1.2rem 2rem;
    }}
    input, select, button, textarea {{
      font-family: inherit;
      font-size: 14px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: 240px 1fr;
      gap: 0.5rem 1rem;
      align-items: center;
      max-width: 1000px;
    }}
    .row {{
      margin-bottom: 0.75rem;
    }}
    .box {{
      border: 1px solid #ccc;
      padding: 0.8rem;
      margin-top: 1rem;
      border-radius: 6px;
    }}
    .mono {{
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      white-space: pre-wrap;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      margin-top: 1rem;
    }}
    th, td {{
      border: 1px solid #ccc;
      padding: 0.4rem 0.5rem;
      text-align: left;
      vertical-align: top;
    }}
    th {{
      background: #f5f5f5;
    }}
    .two {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 1rem;
    }}
    .sidebar {{
      max-height: 500px;
      overflow: auto;
    }}
    .svgbox {{
      border: 1px solid #ccc;
      padding: 0.5rem;
      overflow: auto;
      background: #fff;
    }}
  </style>
</head>
<body>
{body}
</body>
</html>"""


@app.get("/")
def index():
    body = f"""
    <h1>InfoVar UI</h1>

    <div class="box">
      <h2>Consulta locus</h2>
      <form method="get" action="/chrompos">
        <div class="grid">
          <label>CHROM</label><input name="chrom" value="chr1">
          <label>POS</label><input name="pos" value="768448">
          <label>REF</label><input name="ref" value="">
          <label>ALT</label><input name="alt" value="">
          <label>flank</label><input name="flank" value="10">
        </div>
        <div class="row" style="margin-top: 1rem;">
          <button type="submit">Consultar</button>
        </div>
      </form>
    </div>

    <div class="box">
      <h2>Cargar archivo de posiciones</h2>
      <form method="post" action="/upload" enctype="multipart/form-data">
        <input type="file" name="file">
        <button type="submit">Subir</button>
      </form>
      <p>Se esperan columnas tipo CHROM/POS y opcionalmente REF/ALT.</p>
    </div>

    <div class="box">
      <h2>Explorar sistema de archivos</h2>
      <form method="get" action="/fs">
        <input name="path" style="width: 700px;" value="/mnt/cephfs/hot_nvme">
        <button type="submit">Listar</button>
      </form>
    </div>
    """
    return page("InfoVar UI", body)


@app.get("/chrompos")
def chrompos_form_redirect():
    chrom = request.query.get("chrom", "").strip()
    pos = request.query.get("pos", "").strip()
    ref = request.query.get("ref", "").strip()
    alt = request.query.get("alt", "").strip()
    flank = request.query.get("flank", "10").strip()

    if not chrom or not pos:
        return page("Error", "<p>Faltan chrom o pos.</p>")

    qs = {}
    if ref:
        qs["ref"] = ref
    if alt:
        qs["alt"] = alt
    if flank:
        qs["flank"] = flank

    query = urlencode(qs)
    url = f"/chrompos/{chrom}/{pos}"
    if query:
        url += "?" + query

    response.status = 302
    response.set_header("Location", url)
    return ""


@app.get("/chrompos/<chrom>/<pos:int>")
def chrompos_view(chrom: str, pos: int):
    ref = request.query.get("ref")
    alt = request.query.get("alt")
    flank = request.query.get("flank", "10")

    params = {"flank": flank}
    if ref:
        params["ref"] = ref
    if alt:
        params["alt"] = alt

    api_json = requests.get(f"{API_BASE}/chrompos/{chrom}/{pos}", params=params, timeout=300)
    api_json.raise_for_status()
    data = api_json.json()

    api_svg_url = f"{API_BASE}/chrompos/{chrom}/{pos}/svg?{urlencode({'flank': flank})}"
    svg = requests.get(api_svg_url, timeout=300)
    svg.raise_for_status()

    ref_info = data.get("reference", {})
    features = data.get("features", [])
    dbsnp = data.get("dbsnp_variants", [])

    show_table = (
        len(features) <= 30 and len(dbsnp) <= 50
    )

    features_rows = ""
    if show_table:
        for feat in features:
            attrs = feat.get("attributes", {})
            features_rows += (
                "<tr>"
                f"<td>{esc(feat.get('type'))}</td>"
                f"<td>{esc(feat.get('start'))}</td>"
                f"<td>{esc(feat.get('end'))}</td>"
                f"<td>{esc(feat.get('strand'))}</td>"
                f"<td>{esc(attrs.get('gene_name') or attrs.get('Name') or attrs.get('ID'))}</td>"
                "</tr>"
            )

    dbsnp_rows = ""
    if show_table:
        for rec in dbsnp[:50]:
            dbsnp_rows += (
                "<tr>"
                f"<td>{esc(rec.get('id'))}</td>"
                f"<td>{esc(rec.get('pos'))}</td>"
                f"<td>{esc(rec.get('ref'))}</td>"
                f"<td>{esc(','.join(rec.get('alts', [])))}</td>"
                f"<td><div class='mono'>{esc(json.dumps(rec.get('population_frequency_fields', {}), ensure_ascii=False))}</div></td>"
                "</tr>"
            )

    body = f"""
    <h1>Locus {esc(chrom)}:{pos}</h1>
    <p><a href="/">Inicio</a></p>

    <div class="two">
      <div class="box">
        <h2>Parámetros</h2>
        <div class="grid">
          <label>CHROM</label><input value="{esc(chrom)}" readonly>
          <label>POS</label><input value="{esc(pos)}" readonly>
          <label>REF</label><input value="{esc(ref)}" readonly>
          <label>ALT</label><input value="{esc(alt)}" readonly>
          <label>flank</label><input value="{esc(flank)}" readonly>
          <label>ref_base</label><input value="{esc(ref_info.get('ref_base'))}" readonly>
          <label>window</label><input value="{esc(ref_info.get('window_start_1based'))}-{esc(ref_info.get('window_end_1based'))}" readonly>
        </div>
      </div>

      <div class="box">
        <h2>Resumen</h2>
        <div class="grid">
          <label># features</label><input value="{len(features)}" readonly>
          <label># dbSNP records</label><input value="{len(dbsnp)}" readonly>
        </div>
      </div>
    </div>

    <div class="box">
      <h2>SVG</h2>
      <div class="svgbox">
        {svg.text}
      </div>
    </div>

    <div class="two">
      <div class="box">
        <h2>Contexto de referencia</h2>
        <div class="mono">{esc(ref_info.get("context_pretty", ""))}</div>
      </div>

      <div class="box">
        <h2>JSON</h2>
        <div class="mono">{esc(json.dumps(data, ensure_ascii=False, indent=2))}</div>
      </div>
    </div>
    """

    if show_table:
        body += f"""
        <div class="box">
          <h2>Features</h2>
          <table>
            <thead>
              <tr><th>type</th><th>start</th><th>end</th><th>strand</th><th>label</th></tr>
            </thead>
            <tbody>
              {features_rows}
            </tbody>
          </table>
        </div>

        <div class="box">
          <h2>dbSNP</h2>
          <table>
            <thead>
              <tr><th>id</th><th>pos</th><th>ref</th><th>alts</th><th>freq</th></tr>
            </thead>
            <tbody>
              {dbsnp_rows}
            </tbody>
          </table>
        </div>
        """
    else:
        body += """
        <div class="box">
          <h2>Tablas omitidas</h2>
          <p>Se omitieron porque el número de features o variantes dbSNP excede el umbral de visualización.</p>
        </div>
        """

    return page(f"{chrom}:{pos}", body)


@app.post("/upload")
def upload_positions():
    up = request.files.get("file")
    if up is None:
        return page("Error", "<p>No se recibió archivo.</p>")

    files = {"file": (up.filename, up.file.read(), up.content_type or "application/octet-stream")}
    r = requests.post(f"{API_BASE}/positions/upload", files=files, timeout=300)
    r.raise_for_status()
    data = r.json()

    rows_html = []
    for row in data.get("rows", [])[:500]:
        chrom = row["chrom"]
        pos = row["pos"]
        ref = row.get("ref") or ""
        alt = row.get("alt") or ""
        q = {}
        if ref:
            q["ref"] = ref
        if alt:
            q["alt"] = alt
        href = f"/chrompos/{chrom}/{pos}"
        if q:
            href += "?" + urlencode(q)

        rows_html.append(
            "<tr>"
            f"<td>{row['rownum']}</td>"
            f"<td><a href='{esc(href)}'>{esc(chrom)}</a></td>"
            f"<td>{esc(pos)}</td>"
            f"<td>{esc(ref)}</td>"
            f"<td>{esc(alt)}</td>"
            "</tr>"
        )

    body = f"""
    <h1>Archivo cargado</h1>
    <p><a href="/">Inicio</a></p>
    <div class="box">
      <p><b>filename:</b> {esc(data.get("filename"))}</p>
      <p><b>n_rows:</b> {esc(data.get("n_rows"))}</p>
    </div>

    <div class="box">
      <h2>Filas</h2>
      <table>
        <thead>
          <tr><th>rownum</th><th>chrom</th><th>pos</th><th>ref</th><th>alt</th></tr>
        </thead>
        <tbody>
          {"".join(rows_html)}
        </tbody>
      </table>
    </div>
    """
    return page("Upload", body)


@app.get("/fs")
def fs_view():
    path = request.query.get("path", "/mnt/cephfs/hot_nvme")
    r = requests.get(f"{API_BASE}/fs/list", params={"path": path}, timeout=300)
    r.raise_for_status()
    data = r.json()

    rows = []
    parent = None
    if data["path"] != data["root"]:
        parent = str(__import__("pathlib").Path(data["path"]).parent)

    if parent:
        rows.append(
            f"<tr><td>..</td><td></td><td></td><td><a href='/fs?{urlencode({'path': parent})}'>subir</a></td></tr>"
        )

    for item in data["items"]:
        link = f"/fs?{urlencode({'path': item['path']})}" if item["is_dir"] else ""
        rows.append(
            "<tr>"
            f"<td>{esc(item['name'])}</td>"
            f"<td>{'dir' if item['is_dir'] else 'file'}</td>"
            f"<td>{esc(item['size'])}</td>"
            f"<td>{f'<a href=\"{esc(link)}\">abrir</a>' if link else ''}</td>"
            "</tr>"
        )

    body = f"""
    <h1>Filesystem</h1>
    <p><a href="/">Inicio</a></p>
    <div class="box">
      <form method="get" action="/fs">
        <input name="path" style="width: 900px;" value="{esc(data['path'])}">
        <button type="submit">Ir</button>
      </form>
    </div>

    <div class="box sidebar">
      <table>
        <thead>
          <tr><th>name</th><th>type</th><th>size</th><th>action</th></tr>
        </thead>
        <tbody>
          {''.join(rows)}
        </tbody>
      </table>
    </div>
    """
    return page("Filesystem", body)


if __name__ == "__main__":
    run(app, host="127.0.0.1", port=8043, debug=True, reloader=True)
