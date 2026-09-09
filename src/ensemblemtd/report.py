"""Output writers: TSV tables, the standalone HTML summary, and the
ReacNet-style interactive report.

The interactive report reuses the HTML shell ReacNetGenerator emits for a single
run, with its embedded data replaced by the aggregate and two of its row
templates patched so run support shows up next to every species and reaction.
The TSV files are the primary artefacts; the HTML is for browsing.
"""

from __future__ import annotations

import html
import json
import os
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

from .prevalence import ReactionKey

RNGDATA_SCRIPT_ID = "rngdata"


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def link_or_copy(src: Path, dst: Path) -> None:
    """Hard link, symlink, or copy -- whichever the filesystem allows.

    ReacNetGenerator writes its outputs next to its input, so the Li-free file
    has to appear inside the per-run temp directory.  On the cluster's shared
    filesystems hard links are not always available.
    """
    if dst.exists():
        dst.unlink()
    for attempt in (os.link, os.symlink):
        try:
            attempt(src, dst)
            return
        except OSError:
            continue
    shutil.copy2(src, dst)


def write_simple_html(
    path: Path,
    n_runs: int,
    input_files: Sequence[str],
    top_reactions: Sequence[Tuple[str, int, float, int]],
    top_species: Sequence[Tuple[str, int, float]],
) -> None:
    """A dependency-free summary page: run list, top reactions, top species."""
    files_html = "\n".join(
        f"<li>{html.escape(os.path.basename(x))}</li>" for x in input_files
    )
    rxn_rows = "\n".join(
        f"<tr><td>{i + 1}</td><td>{html.escape(k)}</td><td>{rs}</td>"
        f"<td>{pct:.1f}</td><td>{n}</td></tr>"
        for i, (k, rs, pct, n) in enumerate(top_reactions)
    )
    sp_rows = "\n".join(
        f"<tr><td>{i + 1}</td><td>{html.escape(k)}</td><td>{rs}</td>"
        f"<td>{pct:.1f}</td></tr>"
        for i, (k, rs, pct) in enumerate(top_species)
    )
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    path.write_text(
        f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Aggregated ReacNet (Li-removed)</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;margin:18px;}}
.card{{border:1px solid #d0d7de;background:#f6f8fa;padding:10px;border-radius:8px;}}
table{{border-collapse:collapse;width:100%;table-layout:fixed;}}
th,td{{border:1px solid #d0d7de;padding:6px;vertical-align:top;word-wrap:break-word;}}
th{{background:#f6f8fa;}}
.grid{{display:grid;grid-template-columns:1fr;gap:20px;}}
@media(min-width:1200px){{.grid{{grid-template-columns:1fr 1fr;}}}}
</style></head><body>
<h1>Aggregated ReacNet Report (Li-removed)</h1>
<div class="card">
<div><b>Generated:</b> {html.escape(timestamp)}</div>
<div><b>Successful runs:</b> {n_runs}</div>
<div><b>Input runs:</b></div><ul>{files_html}</ul>
</div>
<div class="grid">
<section><h2>Top Reactions</h2>
<table><thead><tr><th>#</th><th>Reaction</th><th>Run support</th><th>% runs</th><th>Total events</th></tr></thead><tbody>{rxn_rows}</tbody></table>
</section>
<section><h2>Top Species</h2>
<table><thead><tr><th>#</th><th>Species</th><th>Run support</th><th>% runs</th></tr></thead><tbody>{sp_rows}</tbody></table>
</section>
</div>
</body></html>
""",
        encoding="utf-8",
    )


def build_reacnet_rngdata(
    species_id_order: Sequence[str],
    species_stats: Dict[str, Dict[str, object]],
    reaction_rows_full: Sequence[Tuple[ReactionKey, int, float, int]],
    reaction_rows_abcd_full: Sequence[Tuple[ReactionKey, int, float, int]],
    speciesshownum: int,
    reactionsshownum: int,
    network_html: str,
    linkreac: Dict[str, List[str]],
) -> dict:
    """Assemble the payload the ReacNet HTML shell expects.

    The short keys are ReacNet's own; ``rs``/``rs_obs``/``m`` are the extra
    ensemble fields the patched templates display.
    """
    species_entries = [
        {
            "s": s,
            "i": i + 1,
            "inchi": str(species_stats[s].get("collapsed_inchi", "")),
            "rep": str(species_stats[s].get("representative_smiles", s)),
            "m": int(species_stats[s].get("member_raw_species_count", 1)),
            "rs": int(species_stats[s]["display_run_support"]),
            "rs_obs": int(species_stats[s]["observed_run_support"]),
            "rs_json": int(species_stats[s]["observed_run_support_json"]),
            "rs_delta": int(species_stats[s]["support_delta_json"]),
            "src": str(species_stats[s]["support_source"]),
            "ev": int(species_stats[s]["involvement_events"]),
            "init": 1 if bool(species_stats[s]["is_initial"]) else 0,
            "init_n": int(species_stats[s]["initial_frame_count"]),
        }
        for i, s in enumerate(species_id_order)
    ]

    def reaction_entries(
        rows: Sequence[Tuple[ReactionKey, int, float, int]]
    ) -> List[dict]:
        return [
            {
                "i": index,
                "l": list(key[0]),
                "r": list(key[1]),
                "n": int(total_events),
                "rs": int(support),
            }
            for index, (key, support, _pct, total_events) in enumerate(rows, start=1)
        ]

    return {
        "speciesshownum": int(speciesshownum),
        "reactionsshownum": int(reactionsshownum),
        "network": [network_html],
        "species": [species_entries],
        "reactions": [reaction_entries(reaction_rows_full)],
        "reactionsabcd": reaction_entries(reaction_rows_abcd_full),
        "linkreac": linkreac,
    }


def inject_rngdata_into_template(template_html: str, rngdata: dict) -> str:
    """Swap the single-run payload in the HTML shell for the aggregate."""
    marker = re.search(
        rf"<script[^>]*id\s*=\s*['\"]?{RNGDATA_SCRIPT_ID}['\"]?[^>]*>",
        template_html,
        flags=re.IGNORECASE,
    )
    if not marker:
        raise RuntimeError("Could not find <script id=rngdata ...> block in template HTML")
    start = marker.end()
    close = template_html.find("</script>", start)
    if close < 0:
        raise RuntimeError("Could not find closing </script> for rngdata block")
    payload = json.dumps(rngdata, ensure_ascii=False)
    return template_html[:start] + payload + template_html[close:]


SPECIES_TEMPLATE_BODY = (
    '<div class="mx-auto col-sm-auto"> '
    '{{include s tmpl="#svgTmpl"/}} '
    "<div>{{:i}}</div>"
    '<div class="small text-muted">runs={{:rs}}</div>'
    '{{if init}}<div class="small text-muted">obs={{:rs_obs}}</div>{{/if}}'
    '<div class="small text-muted">ev={{:ev}}</div>'
    "</div>"
)

REACTION_TEMPLATE_BODY = (
    '<div class="reacid col-sm-auto">{{:i}}</div> '
    '{{include l tmpl="#rsideTmpl"/}} '
    "<div class=col-sm-auto>"
    "<div class=reacnum>{{:n}}</div>"
    "<div><svg height=14.33 width=25>"
    '<path class=narrow d="M24.35 7.613c-3.035 1.11-5.407 2.908-7.113 5.395h-1.299c.585-1.743 1.567-3.39 2.945-4.938H.649V6.278h18.234c-1.378-1.548-2.36-3.2-2.945-4.956h1.299c1.706 2.487 4.078 4.285 7.114 5.395v.896z" font-size=39.506 font-weight=400 /></svg></div>'
    '<div class="small text-muted text-center">runs={{:rs}}</div>'
    "</div> "
    '{{include r tmpl="#rsideTmpl"/}}'
)


def _patch_script_body(src: str, script_id: str, new_body: str) -> Tuple[str, bool]:
    pattern = re.compile(
        rf"(<script[^>]*id\s*=\s*['\"]?{re.escape(script_id)}['\"]?[^>]*>)(.*?)(</script>)",
        flags=re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(src)
    if not match:
        return src, False
    return (
        src[: match.start()] + match.group(1) + new_body + match.group(3) + src[match.end():],
        True,
    )


def patch_reacnet_templates_for_support_labels(
    template_html: str,
) -> Tuple[str, List[str]]:
    """Add run-support labels to the species and reaction row templates.

    A missing template is reported as a warning, not an error: the report still
    renders, just without the extra labels.
    """
    warnings: List[str] = []
    out = template_html
    for script_id, body in (
        ("specTmpl", SPECIES_TEMPLATE_BODY),
        ("rTmpl", REACTION_TEMPLATE_BODY),
    ):
        out, ok = _patch_script_body(out, script_id, body)
        if not ok:
            warnings.append(
                f"template patch skipped: could not find script id={script_id}"
            )
    return out, warnings


def load_top_rows_from_tsv(
    tsv_path: Path,
    required_fields: Sequence[str],
    converters: Dict[str, Callable[[str], object]],
    top_n: int,
) -> List[Tuple[object, ...]]:
    """Read the first ``top_n`` rows of a TSV, skipping any that will not parse.

    Used only by the rebuild path, where the inputs are this program's own
    earlier output and a partially written file is possible.
    """
    out: List[Tuple[object, ...]] = []
    if not tsv_path.exists():
        return out
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        index = {k: i for i, k in enumerate(header)}
        if any(field not in index for field in required_fields):
            return out
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            try:
                row = tuple(
                    converters[field](parts[index[field]]) for field in required_fields
                )
            except (ValueError, KeyError, IndexError):
                continue
            out.append(row)
            if len(out) >= top_n:
                break
    return out
