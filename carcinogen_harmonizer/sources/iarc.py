from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse, urljoin
import hashlib
import re
from io import StringIO, BytesIO

import pandas as pd
import requests

try:
    from ..normalize import clean_text, infer_carcinogen_class, normalize_iarc_group, search_name
except ImportError:
    from normalize import clean_text, infer_carcinogen_class, normalize_iarc_group, search_name
from .local_tables import normalize_key


DEFAULT_IARC_URL = "https://monographs.iarc.who.int/list-of-classifications"
IARC_CLASSIFICATIONS_PAGE = "https://monographs.iarc.who.int/agents-classified-by-the-iarc/26/"
IARC_OFFICIAL_DATA_JS = "https://webapi.iarc.who.int/loc/loc.app.js"
GROUP_LABELS = {
    "Group 1": "Carcinogenic to humans",
    "Group 2A": "Probably carcinogenic to humans",
    "Group 2B": "Possibly carcinogenic to humans",
    "Group 3": "Not classifiable as to its carcinogenicity to humans",
}
CHEMICAL_HINTS = re.compile(
    r"\b(acid|alcohol|aldehyde|amine|aniline|arsenic|benzene|benzidine|benzyl|biphenyl|cadmium|carbon|chloride|chloro|chromium|compound|dibenz|dichloro|dimethyl|dinitro|dioxin|ether|ethylene|formaldehyde|furan|glycid|hydrazine|hydrocarbon|ketone|lead|metal|methyl|nickel|nitro|nitroso|oxide|pesticide|phenol|phosphate|sulfate|toluene|vinyl|benzo|pyrene|fural|furfural|amide|phenyl|pyridine|quinoline|anthracene|naphthyl|nitrosourea)\b",
    re.I,
)
CAS_RE = re.compile(r"\b[0-9]{2,7}-[0-9]{2}-[0-9]\b")


COLUMN_ALIASES = {
    "agent": ["Agent", "agent", "Substance", "substance", "Name", "name", "IARC agent", "iarc_agent_name"],
    "group": ["Group", "group", "IARC group", "iarc_group", "Evaluation", "evaluation"],
    "volume": ["Volume", "volume", "Vol.", "vol", "Monograph", "monograph"],
    "year": ["Year", "year", "Year1", "year1", "Publication year", "publication_year"],
    "cas": ["CAS No.", "CAS", "CASRN", "cas", "casrn", "CAS Number", "cas_number"],
}


def load_iarc_agents(
    source: str | Path | None = None,
    include_groups: Iterable[str] | None = None,
    include_group3: bool = False,
    timeout: int = 30,
) -> pd.DataFrame:
    source = str(source or DEFAULT_IARC_URL)
    include_groups = list(include_groups or ["Group 1", "Group 2A", "Group 2B"])
    if include_group3 and "Group 3" not in include_groups:
        include_groups.append("Group 3")
    raw = _read_iarc_source(source, timeout=timeout)
    out = normalize_iarc_table(raw, source=source)
    out = out[out["iarc_group_normalized"].isin(include_groups)].copy()
    out = out.drop_duplicates(subset=["iarc_agent_key", "iarc_group_normalized"], keep="first").reset_index(drop=True)
    return out


def normalize_iarc_table(df: pd.DataFrame, source: str) -> pd.DataFrame:
    if df.empty:
        raise ValueError("IARC source did not contain any rows")
    cols = _resolve_columns(df)
    missing = [name for name in ["agent", "group"] if name not in cols]
    if missing:
        raise ValueError(f"IARC source missing required column(s): {missing}; found columns: {list(df.columns)}")

    out = pd.DataFrame()
    out["iarc_agent_name"] = df[cols["agent"]].astype(str).map(_clean_text)
    out["iarc_group_raw"] = df[cols["group"]].astype(str).map(_clean_text)
    out["iarc_group_normalized"] = out["iarc_group_raw"].map(normalize_iarc_group)
    out["iarc_group_label"] = out["iarc_group_normalized"].map(GROUP_LABELS).fillna("")
    out["iarc_volume"] = df[cols["volume"]].astype(str).map(_clean_text) if "volume" in cols else ""
    out["iarc_year"] = df[cols["year"]].astype(str).map(_clean_text) if "year" in cols else ""
    out["iarc_casrn"] = df[cols["cas"]].astype(str).map(_extract_first_cas) if "cas" in cols else ""
    out["iarc_agent_key"] = out["iarc_agent_name"].map(normalize_key)
    out["iarc_source_url"] = _source_url(source)
    out["iarc_source_last_checked"] = date.today().isoformat()
    out["agent_entity_type"] = out["iarc_agent_name"].map(infer_agent_entity_type)
    # Treat "Unknown" as chemical-like so PubChem/ClassyFire still get a
    # chance to identify the agent. The alternative is skipping, which
    # leaves these rows permanently unclassifiable.
    out["is_chemical_like"] = out["agent_entity_type"].isin(
        ["Chemical", "MetalOrCompound", "Drug", "Pesticide", "Unknown"]
    )
    out["mixture_or_exposure_flag"] = out["agent_entity_type"].isin(["Mixture", "Exposure", "Occupation", "Process"])
    # For chemical-like agents we leave the class empty so the downstream
    # ``normalize_seed`` step (and later ClassyFire enrichment) can derive
    # it from the agent name / structure with proper source tracking. For
    # non-chemical agents (mixtures, exposures, radiation, biological)
    # the entity type *is* the class.
    out["carcinogen_class_normalized"] = out.apply(
        lambda r: "" if r["is_chemical_like"] else r["agent_entity_type"],
        axis=1,
    )
    out = out[(out["iarc_agent_name"] != "") & (out["iarc_group_normalized"].str.startswith("Group"))].copy()
    return out


def iarc_to_seed_rows(iarc_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in iarc_df.iterrows():
        agent = str(row.get("iarc_agent_name", "")).strip()
        node_id = f"IARC_{_slug(agent)}"
        seed_row = {
            "node_id": node_id,
            "node_label": agent,
            "node_type": "Carcinogen",
            "carcinogen_class": row.get("carcinogen_class_normalized", ""),
            "iarc": row.get("iarc_group_normalized", ""),
            "Status": "1" if row.get("iarc_group_normalized") in {"Group 1", "Group 2A", "Group 2B"} else "0",
            "source_seed": "IARC Monographs",
            "agent_entity_type": row.get("agent_entity_type", ""),
            "is_chemical_like": row.get("is_chemical_like", ""),
            "iarc_agent_name": agent,
            "iarc_group_raw": row.get("iarc_group_raw", ""),
            "iarc_group_label": row.get("iarc_group_label", ""),
            "iarc_volume": row.get("iarc_volume", ""),
            "iarc_year": row.get("iarc_year", ""),
            "iarc_casrn": row.get("iarc_casrn", ""),
            "iarc_source_url": row.get("iarc_source_url", ""),
            "iarc_source_last_checked": row.get("iarc_source_last_checked", ""),
            "expansion_source": "iarc",
        }
        rows.append(seed_row)
    return pd.DataFrame(rows)


def merge_seed_with_iarc_expansion(seed: pd.DataFrame, iarc_agents: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    iarc_seed = iarc_to_seed_rows(iarc_agents)
    seed_norm = _seed_match_keys(seed)
    existing_keys = set(seed_norm["match_key"].dropna().astype(str))
    additions = []
    skipped = []
    for _, row in iarc_seed.iterrows():
        keys = _row_match_keys(row)
        if keys & existing_keys:
            skipped.append({
                "iarc_agent_name": row.get("iarc_agent_name", ""),
                "node_id": row.get("node_id", ""),
                "reason": "already_present_in_seed",
                "matched_keys": ";".join(sorted(keys & existing_keys)),
            })
            continue
        additions.append(row)
        existing_keys.update(keys)
    additions_df = pd.DataFrame(additions)
    skipped_df = pd.DataFrame(skipped)
    combined = pd.concat([seed.copy(), additions_df], ignore_index=True, sort=False).fillna("")
    return combined, skipped_df


def infer_agent_entity_type(agent: str) -> str:
    text = str(agent or "").lower()
    if not text.strip():
        return "Unknown"
    if any(k in text for k in ["virus", "helicobacter", "schistosoma", "opisthorchis", "clonorchis", "infection with", "malaria"]):
        return "BiologicalAgent"
    if any(k in text for k in ["radiation", "radionuclide", "x-radiation", "gamma-radiation", "ultraviolet", "neutron", "radioiodines"]):
        return "Radiation"
    if any(k in text for k in ["occupation", "occupational", "worker", "painter", "hairdresser", "barber", "manufacture", "production", "mining", "founding", "process"]):
        return "Occupation"
    if any(k in text for k in ["smoke", "exhaust", "emissions", "pollution", "beverages", "mate", "diet", "consumption", "dust", "mists", "fumes"]):
        return "Exposure"
    if any(k in text for k in ["mixture", "extract", "oil", "oils", "fuel", "tar", "pitch", "soot"]):
        return "Mixture"
    if any(k in text for k in ["therapy", "treated with", "drug", "cyclosporine", "azathioprine", "tamoxifen", "etoposide", "melphalan", "cisplatin"]):
        return "Drug"
    if any(k in text for k in ["pesticide", "herbicide", "insecticide", "fungicide", "glyphosate", "chlordane", "ddt", "lindane"]):
        return "Pesticide"
    if any(k in text for k in ["arsenic", "beryllium", "cadmium", "chromium", "cobalt", "lead", "nickel", "uranium", "metal", "compounds"]):
        return "MetalOrCompound"
    if CAS_RE.search(text) or CHEMICAL_HINTS.search(text):
        return "Chemical"
    return "Unknown"


def _read_iarc_source(source: str, timeout: int = 30) -> pd.DataFrame:
    parsed = urlparse(source)
    if parsed.scheme in {"http", "https"}:
        host = (parsed.netloc or "").lower()
        # Primary path: pull the official IARC dataset embedded in loc.app.js.
        # Used when the user passes the default URL or any monographs.iarc.who.int URL.
        if source == DEFAULT_IARC_URL or "iarc.who.int" in host:
            try:
                return _read_iarc_official_dataset(timeout=timeout)
            except Exception:
                pass
        try:
            return _read_iarc_url(source, timeout=timeout)
        except RuntimeError:
            if source == DEFAULT_IARC_URL:
                return _read_iarc_url(IARC_CLASSIFICATIONS_PAGE, timeout=timeout)
            raise
    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"IARC source not found: {source}")
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, dtype=str, keep_default_na=False)
    if suffixes.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _read_iarc_official_dataset(timeout: int = 30) -> pd.DataFrame:
    """Fetch the official IARC Monographs classifications dataset.

    The IARC web team publishes the full agent list as a JS module at
    ``loc.app.js``. We extract the ``agents`` array and return it as a
    DataFrame with the column names expected by ``normalize_iarc_table``.
    """
    response = requests.get(
        IARC_OFFICIAL_DATA_JS,
        timeout=timeout,
        headers={"User-Agent": "carcinogen-pipeline/0.2"},
    )
    response.raise_for_status()
    text = response.text
    match = re.search(r"agents\s*:\s*\[", text)
    if not match:
        raise RuntimeError("IARC loc.app.js did not contain an 'agents' array")
    start = match.end() - 1
    array_text = _slice_balanced_array(text, start)
    if array_text is None:
        raise RuntimeError("Could not extract balanced 'agents' array from IARC loc.app.js")
    json_text = _iarc_js_to_json(array_text)
    try:
        import json

        agents = json.loads(json_text)
    except Exception as exc:
        raise RuntimeError(f"Could not parse IARC loc.app.js agents JSON: {exc}") from exc
    rows = []
    for agent in agents:
        rows.append(
            {
                "Agent": str(agent.get("name", "") or ""),
                "Group": f"Group {agent.get('group', '')}".strip() if agent.get("group") else "",
                "Volume": ", ".join(str(v) for v in (agent.get("volume") or [])),
                "Year": str(agent.get("year", "") or ""),
                "CAS": ",".join(str(c) for c in (agent.get("cas") or [])),
            }
        )
    return pd.DataFrame(rows)


def _slice_balanced_array(text: str, start: int) -> str | None:
    """Return the substring starting at ``start`` that contains a balanced JS array."""
    if start >= len(text) or text[start] != "[":
        return None
    depth = 0
    in_string = False
    quote_char = ""
    escaped = False
    for i in range(start, len(text)):
        ch = text[i]
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == quote_char:
                in_string = False
            continue
        if ch in ("'", '"'):
            in_string = True
            quote_char = ch
            continue
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                return text[start : i + 1]
    return None


_IARC_KEY_RE = re.compile(
    r"(?<![A-Za-z0-9_])(name|group|cas|volume|year|yeareval|comment|in_prep)\s*:"
)


def _iarc_js_to_json(js_text: str) -> str:
    """Convert the JS object-literal subset used by IARC loc.app.js to JSON."""
    text = _IARC_KEY_RE.sub(r'"\1":', js_text)
    text = re.sub(r"(?<![A-Za-z0-9_])!0(?![A-Za-z0-9_])", "true", text)
    text = re.sub(r"(?<![A-Za-z0-9_])!1(?![A-Za-z0-9_])", "false", text)
    return _normalize_js_strings(text)


def _normalize_js_strings(text: str) -> str:
    """Convert single-quoted JS strings to JSON-safe double-quoted strings."""
    out: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == '"':
            j = i + 1
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == '"':
                    j += 1
                    break
                j += 1
            out.append(text[i:j])
            i = j
        elif ch == "'":
            j = i + 1
            buf: list[str] = []
            while j < n:
                c = text[j]
                if c == "\\":
                    if j + 1 < n:
                        nxt = text[j + 1]
                        if nxt == "'":
                            buf.append("'")
                        else:
                            buf.append(c)
                            buf.append(nxt)
                        j += 2
                        continue
                    j += 1
                    continue
                if c == "'":
                    j += 1
                    break
                if c == '"':
                    buf.append('\\"')
                else:
                    buf.append(c)
                j += 1
            out.append('"' + "".join(buf) + '"')
            i = j
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _read_iarc_url(url: str, timeout: int = 30) -> pd.DataFrame:
    try:
        tables = pd.read_html(url)
    except Exception:
        tables = []
    for table in tables:
        table = table.astype(str).fillna("")
        cols = _resolve_columns(table)
        if "agent" in cols and "group" in cols:
            return table
    try:
        response = requests.get(url, timeout=timeout, headers={"User-Agent": "carcinogen-pipeline/0.2"})
        response.raise_for_status()
    except Exception as exc:
        raise RuntimeError(f"Could not retrieve IARC URL {url}: {exc}") from exc
    # Fallback: extract downloadable tabular links (csv/xlsx/tsv/txt) from page HTML.
    link_df = _read_iarc_download_links(response.text, base_url=url, timeout=timeout)
    if link_df is not None and not link_df.empty:
        return link_df
    if "Loading" in response.text and "wpDataTables" in response.text or "List of Classifications" in response.text:
        raise RuntimeError(
            "The current IARC classifications table is rendered with browser-side JavaScript. "
            "Could not auto-download a tabular export from the page. "
            "Download/export the IARC classification list as CSV/XLSX and pass it with --iarc-source. "
            f"Page checked: {url}"
        )
    raise RuntimeError(
        "Could not find an Agent/Group table in the IARC page. "
        "Save the IARC list as CSV/XLSX and pass it with --iarc-source."
    )


def _read_iarc_download_links(html: str, base_url: str, timeout: int = 30) -> pd.DataFrame | None:
    hrefs = re.findall(r"""href=["']([^"']+\.(?:csv|xlsx|xls|tsv|txt)(?:\?[^"']*)?)["']""", html, flags=re.I)
    seen = set()
    for href in hrefs:
        full_url = urljoin(base_url, href)
        if full_url in seen:
            continue
        seen.add(full_url)
        try:
            table = _read_iarc_tabular_url(full_url, timeout=timeout)
            if table is not None and not table.empty:
                cols = _resolve_columns(table.astype(str).fillna(""))
                if "agent" in cols and "group" in cols:
                    return table
        except Exception:
            continue
    return None


def _read_iarc_tabular_url(url: str, timeout: int = 30) -> pd.DataFrame | None:
    resp = requests.get(url, timeout=timeout, headers={"User-Agent": "carcinogen-pipeline/0.2"})
    resp.raise_for_status()
    low = url.lower()
    if low.endswith((".xlsx", ".xls")):
        return pd.read_excel(BytesIO(resp.content), dtype=str, keep_default_na=False)
    if low.endswith((".tsv", ".txt")):
        text = resp.content.decode(resp.encoding or "utf-8", errors="replace")
        return pd.read_csv(StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    text = resp.content.decode(resp.encoding or "utf-8", errors="replace")
    return pd.read_csv(StringIO(text), dtype=str, keep_default_na=False)


def _resolve_columns(df: pd.DataFrame) -> dict[str, str]:
    resolved: dict[str, str] = {}
    norm_to_original = {normalize_key(c): c for c in df.columns}
    for target, candidates in COLUMN_ALIASES.items():
        for candidate in candidates:
            key = normalize_key(candidate)
            if key in norm_to_original:
                resolved[target] = norm_to_original[key]
                break
    return resolved


def _clean_text(value: str) -> str:
    return clean_text(value)


def _extract_first_cas(value: str) -> str:
    match = CAS_RE.search(str(value or ""))
    return match.group(0) if match else ""


def _source_url(source: str) -> str:
    parsed = urlparse(str(source))
    if parsed.scheme in {"http", "https"}:
        return source
    return IARC_CLASSIFICATIONS_PAGE


def _slug(value: str, max_len: int = 70) -> str:
    base = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_")
    if not base:
        base = hashlib.sha1(str(value).encode("utf-8")).hexdigest()[:12]
    return base[:max_len]


def _seed_match_keys(seed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in seed.iterrows():
        for key in _row_match_keys(row):
            rows.append({"node_id": row.get("node_id", ""), "match_key": key})
    return pd.DataFrame(rows)


def _row_match_keys(row: pd.Series | dict) -> set[str]:
    row_series = pd.Series(row)
    keys = set()
    for field in ["node_label", "node_id", "iarc_agent_name", "preferred_query_name"]:
        value = str(row_series.get(field, "")).strip()
        if value:
            keys.add(normalize_key(value))
    try:
        query = search_name(row_series)
        if query:
            keys.add(normalize_key(query))
    except Exception:
        pass
    cas = str(row_series.get("iarc_casrn", row_series.get("casrn", ""))).strip()
    if cas:
        keys.add(f"cas:{cas}")
    return {k for k in keys if k}
