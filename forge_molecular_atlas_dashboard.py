"""Standalone, offline-first Forge Molecular Atlas control tower."""
from __future__ import annotations

import json
import re
from typing import Any

import plotly.graph_objects as go
import streamlit as st

from forge_genomics_protein_spine import (
    build_genomics_figure,
    build_genomics_graph,
    build_genomics_protein_spine,
    filter_graph_by_program,
    identity_neighborhood,
    resolve_identity_query,
)
from forge_molecular_atlas import PROGRAMS, assess_consumer_receipt, build_program_projection_manifests, make_consumer_receipt
from forge_molecular_atlas_feed import load_feed_report, load_inbox_candidates
from forge_molecular_atlas_investigation import build_investigation, investigation_is_current, load_report
from forge_mito_next_assay_packet import build_next_assay_packet, packet_summary_rows
from forge_molecular_atlas_adjudication import build_adjudication_desk
from forge_molecular_atlas_operator import load_report as load_operator_report
from forge_mito_evidence_units import load_mito_evidence_delivery
from forge_molecular_assays import build_assay_catalog
from forge_molecular_consumers import CONSUMER_BINDINGS, CONSUMER_PROGRAMS
from forge_molecular_curation import build_molecular_curation_queue
from forge_molecular_sources import SOURCE_NAMES, source_health
from forge_molecular_source_cache import DEFAULT_CACHE_DIR, load_source_statuses
from forge_molecular_source_refresh import refresh_source_review
from forge_molecular_projection_store import DEFAULT_STORE_DIR, read_consumption_status, read_feed_health

st.set_page_config(page_title="Forge Molecular Atlas", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")
st.markdown("""<style>
.stApp,.stSidebar{background:#101820;color:#e7efed}.block-container{padding-top:1.25rem;max-width:1600px}
[data-testid="stMetric"]{background:#17262d;border:1px solid #263841;border-radius:8px;padding:8px}
.stMetric label,[data-testid="stMetricLabel"]{color:#b9cfcd!important}[data-testid="stMetricValue"]{color:#f4faf8!important}
[data-baseweb="tab"]{color:#b9cfcd!important}[data-baseweb="tab"][aria-selected="true"]{color:#6ed8d0!important}
[data-testid="stWidgetLabel"] p,[data-testid="stWidgetLabel"] label,[data-testid="stMarkdownContainer"] label{color:#dce9e6!important;font-weight:600}
[data-baseweb="select"]>div,[data-testid="stTextInput"] input{background:#17262d!important;color:#f4faf8!important;border:1px solid #58716f!important;border-radius:6px!important}
[data-baseweb="select"] span,[data-baseweb="menu"] li,[role="option"]{color:#f4faf8!important}
[data-baseweb="select"]>div:focus-within,[data-testid="stTextInput"] input:focus{border-color:#6ed8d0!important;box-shadow:0 0 0 2px #245b5a!important}
[data-testid="stRadio"] label,[data-testid="stCheckbox"] label{color:#dce9e6!important}
[data-testid="stRadio"] label p,[data-testid="stCheckbox"] label p{color:#dce9e6!important}
.stCaption,[data-testid="stCaptionContainer"]{color:#b9cfcd!important}.atlas-card{background:#17262d;border:1px solid #263841;border-radius:8px;padding:12px;height:100%}
[data-testid="stAlert"]{color:#f4faf8!important}[data-testid="stPlotlyChart"] .textpoint text,[data-testid="stPlotlyChart"] .annotation-text{fill:#F4F7F5!important}
[data-testid="stPlotlyChart"]{min-height:520px}
</style>""", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def load_catalog() -> dict[str, Any]:
    return build_genomics_protein_spine()


def build_feed_statuses() -> list[dict[str, Any]]:
    """Return retained reviewed-feed health without current-catalog validation side effects."""
    health = read_feed_health(store_dir=DEFAULT_STORE_DIR)
    rows = []
    for health_row in health:
        program = health_row["program"]
        configured_names = health_row["configured_consumers"]
        drift = health_row.get("catalog_drift")
        status = "ready" if health_row["reviewed_pin_available"] else "unpinned"
        if drift and status == "ready":
            status = "drifted"
        valid_count = health_row["valid_current_consumption"]
        rows.append({
            "program": program,
            "status": status,
            "boundary": "certified-ready" if health_row["reviewed_pin_available"] else "unbound",
            "connection count": len(configured_names),
            "named consumers": ", ".join(configured_names) if configured_names else "none",
            "configured pin": "retained reviewed pin" if health_row["reviewed_pin_available"] else "unavailable",
            "current catalog": "drift detected; reviewed snapshot retained" if drift else "matches reviewed pin",
            "actual consumption": f"{valid_count} valid current receipt(s)" if valid_count else "no valid current receipt",
            "consumed consumers": ", ".join(health_row.get("valid_consumers", [])) if valid_count else "none",
            "consumption hash": ", ".join(health_row.get("valid_receipt_hashes", [])) if valid_count else "none",
            "latest consumed at": health_row.get("latest_consumed_at") or "—",
            "invalid receipts": health_row["invalid_count"],
            "stale receipts": health_row["stale_count"],
            "unattributed invalid receipts": health_row.get("unattributed_invalid_count", 0),
        })
    return rows


def architecture_figure(entity: dict[str, Any]) -> go.Figure:
    fig = go.Figure()
    for feature in entity.get("architecture") or []:
        fig.add_trace(go.Bar(name=feature["label"], x=[feature["end"] - feature["start"] + 1], y=[feature["type"]], base=[feature["start"]], orientation="h", hovertemplate=f"{feature['label']}: residues {feature['start']}–{feature['end']}<extra></extra>"))
    fig.update_layout(title="Protein architecture · residue-range identity context", barmode="stack", height=270, paper_bgcolor="#101820", plot_bgcolor="#101820", font_color="#e7efed", showlegend=False, margin=dict(l=10, r=10, t=45, b=10), xaxis_title="Residue position")
    return fig


_PDB_ID_RE = re.compile(r"^[0-9][A-Z0-9]{3}$")
_ALPHAFOLD_ACCESSION_RE = re.compile(r"^AF-[A-Z0-9]{6,10}-F1$")
_UNIPROT_RE = re.compile(r"^[A-Z0-9]{6,10}$")
MOLSTAR_VIEWER_BASE_URL = "https://molstar.org/viewer/?pdb="


def approved_pdb_id(identity: dict[str, Any]) -> str | None:
    """Return only a verified four-character PDB ID present on this record."""
    value = (identity.get("identifiers") or {}).get("pdb")
    return value if isinstance(value, str) and _PDB_ID_RE.fullmatch(value) else None


def molstar_viewer_url(identity: dict[str, Any]) -> str | None:
    """Build the fixed Mol* URL from the selected catalog record, never user text."""
    pdb_id = approved_pdb_id(identity)
    return f"{MOLSTAR_VIEWER_BASE_URL}{pdb_id}" if pdb_id else None


def alphafold_accession(identity: dict[str, Any]) -> str | None:
    value = (identity.get("identifiers") or {}).get("alphafold")
    return value if isinstance(value, str) and _ALPHAFOLD_ACCESSION_RE.fullmatch(value) else None


def alphafold_record_url(identity: dict[str, Any]) -> str | None:
    """Return an AlphaFold link only when the record has a validated accession."""
    accession = alphafold_accession(identity)
    uniprot = (identity.get("identifiers") or {}).get("uniprot")
    if not accession or not isinstance(uniprot, str) or not _UNIPROT_RE.fullmatch(uniprot):
        return None
    return f"https://alphafold.ebi.ac.uk/entry/{uniprot}"


def protein_studio_state(identity: dict[str, Any]) -> dict[str, str | None]:
    """Summarize safe, record-bound structure links for the Protein Studio UI."""
    identifiers = identity.get("identifiers") or {}
    uniprot = identifiers.get("uniprot")
    return {
        "uniprot": uniprot if isinstance(uniprot, str) and _UNIPROT_RE.fullmatch(uniprot) else None,
        "pdb": approved_pdb_id(identity),
        "alphafold": alphafold_accession(identity),
        "molstar_url": molstar_viewer_url(identity),
        "alphafold_url": alphafold_record_url(identity),
    }


def source_status_rows() -> list[dict[str, str]]:
    rows = []
    for row in load_source_statuses(DEFAULT_CACHE_DIR):
        rows.append({
            "source system": row["source"],
            "health": row.get("status", "degraded"),
            "receipt state": row.get("receipt_state", "awaiting validated receipt"),
            "release/revision": row.get("release/revision", "—"),
            "last success": row.get("last_success", "—"),
            "last failure": row.get("last_failure", row.get("failure_reason", "—")),
            "age": row.get("age", "—"),
            "data": row.get("data", "not loaded"),
        })
    return rows


def literature_candidate_rows(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Flatten staged literature candidates into bounded review rows."""
    rows = []
    for record in records:
        if record.get("record_type") != "literature_evidence_candidate":
            continue
        rows.append({
            "title": str(record.get("title") or ""),
            "compound": str(record.get("compound") or ""),
            "lane": str(record.get("lane") or ""),
            "evidence IDs": ", ".join(str(item) for item in record.get("evidence_ids") or []),
            "journal": str(record.get("journal") or ""),
            "year": str(record.get("year") or ""),
            "linked identities": ", ".join(str(link.get("entity_id")) for link in record.get("identity_links") or []),
            "status": str(record.get("status") or ""),
        })
    return rows


def edge_candidate_rows(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Flatten staged evidence-edge candidates into bounded review rows."""
    rows = []
    for record in records:
        if record.get("record_type") != "evidence_edge_candidate":
            continue
        rows.append({
            "source": str(record.get("source") or ""),
            "target": str(record.get("target") or ""),
            "relation": str(record.get("relation") or ""),
            "evidence IDs": ", ".join(str(item) for item in record.get("evidence_ids") or []),
            "status": str(record.get("status") or ""),
        })
    return rows


def feed_run_caption(report: dict[str, Any] | None) -> str:
    """Summarize the last feed run without exposing or applying anything."""
    if not report:
        return "No feed run has been recorded yet. Run the Titans Med literature feed to stage review-required studies."
    return (
        f"Last feed run {str(report.get('generated_at') or 'unknown')} ({str(report.get('mode') or 'unknown')}) · "
        f"{report.get('candidates_staged', 0)} literature candidates and {report.get('edges_staged', 0)} edge candidates staged · "
        f"{report.get('duplicates_skipped', 0)} duplicates skipped · {report.get('studies_seen', 0)} studies seen. "
        "Every staged record is review-required; nothing is applied automatically."
    )


def investigation_metric_rows(bundle: dict[str, Any]) -> list[dict[str, int]]:
    metrics = bundle.get("metrics", {})
    return [{
        "candidates": metrics.get("candidates", 0),
        "publication inventory": metrics.get("publication_count", metrics.get("candidates", 0)),
        "evidence units": metrics.get("evidence_unit_count", 0),
        "deduplication": metrics.get("deduplication", 0),
        "provisional units": metrics.get("fallback_provisional_units", 0),
        "clusters": metrics.get("clusters", 0),
        "compounds": metrics.get("compounds", 0),
        "contradiction flags": metrics.get("contradiction_flags", 0),
        "critical": metrics.get("critical", 0),
        "high": metrics.get("high", 0),
        "medium": metrics.get("medium", 0),
        "low": metrics.get("low", 0),
    }]


def ranked_candidate_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    return [{
        "rank": index + 1,
        "title": str(item.get("title") or ""),
        "compound": str(item.get("compound") or ""),
        "lane": str(item.get("lane") or ""),
        "endpoint": str(item.get("primary_endpoint") or ""),
        "year": str(item.get("year") or ""),
        "score": item.get("triage_score", 0),
        "band": str(item.get("priority_band") or ""),
        "evidence unit": str(item.get("evidence_unit_id") or ""),
    } for index, item in enumerate(bundle.get("ranked_candidates", []))]


def ranked_evidence_unit_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    """Render the normalized, ranked evidence-unit view without raw text."""
    return [{
        "rank": index + 1,
        "evidence unit ID": str(item.get("evidence_unit_id") or ""),
        "cohort identity": str(item.get("cohort_identity") or ""),
        "compound": str(item.get("compound") or ""),
        "endpoint": str(item.get("endpoint") or ""),
        "representative title": str(item.get("representative_publication_title") or item.get("title") or ""),
        "publications (inventory)": item.get("publication_count", 0),
        "dependency state": str(item.get("dependency_state") or ""),
        "adjudication state": str(item.get("adjudication_status") or ""),
        "registry bound": bool(item.get("registry_bound", False)),
        "score": item.get("triage_score", 0),
    } for index, item in enumerate(bundle.get("ranked_evidence_units", []))]


def cluster_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    return [{
        "cluster": str(item.get("cluster_id") or ""),
        "compound": str(item.get("compound") or ""),
        "endpoint": str(item.get("endpoint") or ""),
        "studies": item.get("candidate_count", 0),
        "publications (inventory)": item.get("publication_count", item.get("candidate_count", 0)),
        "evidence units": item.get("evidence_unit_count", item.get("candidate_count", 0)),
        "representative": str(item.get("representative_title") or ""),
        "top score": item.get("representative_score", 0),
    } for item in bundle.get("clusters", [])]


def contradiction_rows(bundle: dict[str, Any]) -> list[dict[str, str]]:
    return [{
        "compound": str(item.get("compound") or ""),
        "flag": str(item.get("flag_type") or ""),
        "description": str(item.get("description") or ""),
    } for item in bundle.get("contradictions", [])]


def synthesis_rows(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    return [{
        "compound": str(item.get("compound") or ""),
        "studies": item.get("candidate_count", 0),
        "publications (inventory)": item.get("publication_count", item.get("candidate_count", 0)),
        "evidence units": item.get("evidence_unit_count", item.get("candidate_count", 0)),
        "clusters": item.get("cluster_count", 0),
        "top signal": str(item.get("top_candidate_title") or ""),
        "flags": ", ".join(item.get("flags") or []) or "none",
        "summary": str(item.get("summary") or ""),
    } for item in bundle.get("compound_syntheses", [])]


def adjudication_desk_rows(desk: dict[str, Any], operator_report: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    receipts = {
        str(item.get("queue_id") or ""): item
        for item in (operator_report or {}).get("receipts") or []
        if isinstance(item, dict)
    }
    rows = []
    for item in desk.get("packets", []):
        receipt = receipts.get(str(item.get("queue_id") or ""), {})
        rows.append({
            "priority": item.get("priority_band", ""),
            "score": item.get("priority_score", 0),
            "compound": item.get("compound", ""),
            "endpoint": item.get("endpoint", ""),
            "nct": item.get("nct_id", ""),
            "conflict": item.get("conflict_reason", ""),
            "operator role": receipt.get("recommended_role") or item.get("operator_recommended_role", ""),
            "operator status": receipt.get("operator_status") or item.get("operator_status", ""),
            "claims": ", ".join(item.get("linked_claim_ids") or []) or "none",
            "status": item.get("status", ""),
        })
    return rows


def assay_identity_rows(catalog: dict[str, Any]) -> list[dict[str, str]]:
    return [{
        "assay ID": item.get("id", ""),
        "label": item.get("label", ""),
        "endpoint": item.get("endpoint_id", ""),
        "assay": item.get("assay", ""),
        "matrix": item.get("matrix", ""),
        "population": item.get("population", ""),
    } for item in catalog.get("assays", [])]


def graph_for_program(graph: Any, program: str | None) -> Any:
    """Pure dashboard adapter: return a detached program-filtered graph."""
    return filter_graph_by_program(graph, program)


def graph_for_focus(graph: Any, stable_id: str | None, depth: int = 1) -> Any:
    """Pure dashboard adapter: return a detached identity neighborhood."""
    return identity_neighborhood(graph, stable_id, depth) if stable_id else graph.copy()


def graph_for_search(graph: Any, query: str, bundle: dict[str, Any]) -> Any:
    """Return a deterministic detached graph containing query matches."""
    if not query.strip():
        return graph.copy()
    matches = {hit["id"] for hit in resolve_identity_query(query, bundle)}
    return graph.subgraph(sorted(node for node in graph if node in matches)).copy()


def provenance_inspector_rows(entity: dict[str, Any]) -> list[dict[str, str]]:
    """Expose record-bound provenance without changing catalog state."""
    return [{"source": url, "record": entity["id"]} for url in entity.get("source_urls", [])]


def connected_consumer_rows(entity: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"program": program, "consumers": ", ".join(sorted(consumer for consumer, programs in CONSUMER_PROGRAMS.items() if program in programs)) or "none"}
        for program in entity.get("programs", [])
    ] or [{"program": "none", "consumers": "none"}]


def program_display_name(program: str) -> str:
    """Title-case program names for the UI while retaining canonical values."""
    return "All programs" if program == "All programs" else program.title()


RESEARCH_BOUNDARY_TEXT = (
    "Identity edges locate stable records. Context/Hypothesis edges preserve provenance and boundaries; "
    "they are not mechanistic proof, clinical interpretation, pathogenicity, treatment, or causal conclusions."
)


filter_graph = graph_for_program
build_neighborhood = graph_for_focus


try:
    catalog = load_catalog()
except ValueError as exc:
    st.error(f"Molecular identity spine failed closed: {exc}")
    st.stop()

st.title("Forge Molecular Atlas")
st.caption("Shared molecular identity control tower · research-planning-only · loopback dashboard")
st.warning("Research planning only. Identity, context, and evidence are separate layers. Nothing here is diagnostic, causal, therapeutic, dosing, eligibility, or clinical advice.")

entities = catalog["entities"]
proteins = [item for item in entities if item["kind"] == "protein"]
genes = [item for item in entities if item["kind"] == "gene"]
cards = st.columns(6)
source_rows = source_status_rows()
source_receipt_count = sum(row["receipt state"] == "validated receipt" for row in source_rows)
for card, label, value in zip(cards, ("Genes", "Proteins", "Edges", "Seeded identity variants", "PDB-linked proteins", "Source receipts"), (len(genes), len(proteins), len(catalog["edges"]), len(catalog["seeded_variants"]), sum(bool(p["identifiers"].get("pdb")) for p in proteins), source_receipt_count)):
    card.metric(label, value)

feed_statuses = build_feed_statuses()
ready_feed_count = sum(row["status"] == "ready" for row in feed_statuses)
with st.expander(f"Program feed control tower · {ready_feed_count}/5 feeds ready", expanded=False):
    st.caption("Five program feeds · explicit hash-pinned research consumers · unbound only when no verified binding exists.")
    st.dataframe(feed_statuses, width='stretch', hide_index=True)

identity_tab, protein_tab, genome_tab, evidence_tab, curation_tab, sources_tab = st.tabs(["Molecular graph explorer", "Protein architecture", "Genome & mtDNA identity", "Evidence boundary", "Intelligent curation", "Source health & review"])

with identity_tab:
    st.markdown("**Identity layer** · stable IDs and typed edges only. Context is not evidence.")
    programs = ["All programs", *PROGRAMS]
    control_columns = st.columns([1.05, 1.8, 2.7, 1.35], gap="small")
    with control_columns[0]:
        selected_program = st.selectbox("Program", programs, format_func=program_display_name, key="atlas_program")
    program_value = None if selected_program == "All programs" else selected_program
    full_graph = build_genomics_graph(catalog)
    program_graph = graph_for_program(full_graph, program_value)
    search_graph = program_graph
    choices = sorted(program_graph.nodes, key=lambda node: (program_graph.nodes[node]["display_label"], node))
    with control_columns[1]:
        query = st.text_input("Find identity", placeholder="Search label, stable ID, or official identifier", key="atlas_identity_search")
    if query:
        search_graph = graph_for_search(program_graph, query, catalog)
        choices = sorted(search_graph.nodes, key=lambda node: (program_graph.nodes[node]["display_label"], node))
        if not choices:
            st.info("No matching molecular identities were found for the current search and program filter.")
    focus_options = ["No focus", *choices]
    with control_columns[2]:
        selected = st.selectbox("Focus entity (searchable)", focus_options, format_func=lambda value: "No focus" if value == "No focus" else f"{program_graph.nodes[value]['display_label']} · {value}", key="atlas_identity")
    with control_columns[3]:
        st.markdown("**View**")
        view_columns = st.columns(2, gap="small")
        with view_columns[0]:
            depth = st.radio("Depth", [1, 2], horizontal=True, key="atlas_depth", format_func=lambda value: f"{value}-hop", label_visibility="visible")
        with view_columns[1]:
            show_labels = st.checkbox("Labels", value=False, key="atlas_show_labels")
    focused_graph = graph_for_focus(search_graph, None if selected == "No focus" else selected, depth)
    st.plotly_chart(build_genomics_figure(focused_graph, height=520, show_labels=show_labels), use_container_width=True, config={"displayModeBar": False})
    st.caption("Node symbols encode entity kind; node size is a visual category cue. Edge color/dash encodes epistemic status.")
    st.markdown("**Edge status legend** · solid teal = identity · dashed green = context · long-dash slate = assay context · dotted slate = hypothesis (review-required).")
    st.info(RESEARCH_BOUNDARY_TEXT)
    if selected != "No focus":
        item = next(entity for entity in entities if entity["id"] == selected)
        st.subheader("Identity")
        st.dataframe([{"label": item["label"], "kind": item["kind"], "stable ID": item["id"], "locus/scope": item.get("locus", item.get("scope", "not specified")), "taxonomy": item.get("taxonomy", "not specified")}], width='stretch', hide_index=True)
        st.subheader("Official identifiers")
        st.dataframe([{"system": key, "identifier": value} for key, value in item.get("identifiers", {}).items()], width='stretch', hide_index=True)
        st.subheader("Provenance and source links")
        provenance = provenance_inspector_rows(item)
        st.markdown("<br>".join(f"- [{row['source']}]({row['source']})" for row in provenance) or "No source links recorded", unsafe_allow_html=False)
        st.subheader("Programs")
        st.write(", ".join(item.get("programs") or []) or "No programs recorded")
        st.subheader("Connected-consumer impact")
        st.dataframe(connected_consumer_rows(item), width='stretch', hide_index=True)

with protein_tab:
    st.markdown("**Context layer** · architecture ranges and external records are navigational identity context, not functional claims.")
    selected = st.selectbox("Protein", [p["id"] for p in proteins], format_func=lambda value: next(p["label"] for p in proteins if p["id"] == value), key="atlas_protein")
    item = next(p for p in proteins if p["id"] == selected)
    studio = protein_studio_state(item)
    st.markdown(f"**{item['label']}** · UniProt `{item['identifiers']['uniprot']}`")
    record_links = []
    if studio["pdb"]:
        record_links.append(f"[PDB {studio['pdb']}](https://www.rcsb.org/structure/{studio['pdb']})")
    if studio["alphafold_url"]:
        record_links.append(f"[AlphaFold {studio['alphafold']}](https://alphafold.ebi.ac.uk/entry/{studio['uniprot']})")
    st.markdown(" · ".join(record_links) if record_links else "No approved external structure record links")
    st.caption("Real public record links are provided; this dashboard does not call them.")
    st.subheader("Protein Studio")
    st.caption("Reference visualization only · external structure context is not evidence, diagnosis, or a claim input.")
    st.markdown(
        f"**Selected identity provenance** · UniProt `{studio['uniprot'] or 'not available'}` · "
        f"PDB `{studio['pdb'] or 'none approved'}` · AlphaFold `{studio['alphafold'] or 'not available'}`"
    )
    if studio["molstar_url"]:
        st.markdown("**Experimental PDB reference visualization**")
        st.components.v1.iframe(studio["molstar_url"], height=520, scrolling=False)
        st.caption("Mol* is pinned to the selected identity record's approved PDB ID; it cannot alter catalog evidence or claims.")
    else:
        st.info("No approved experimental PDB structure is available for this selected identity.")
        if studio["alphafold_url"]:
            st.markdown(f"[Open AlphaFold predicted-structure reference]({studio['alphafold_url']})")
        else:
            st.caption("No AlphaFold reference accession is present on this identity record.")
    st.plotly_chart(architecture_figure(item), use_container_width=True, config={"displayModeBar": False})

with genome_tab:
    st.markdown("**Non-coordinate identity view** · labels and assembly/accession context, not a genome browser or variant interpretation.")
    genome_rows = [{"identity": g["label"], "stable ID": g["id"], "locus": g.get("locus", "nuclear"), "assembly/accession": g.get("scope", "GRCh38 / NC_012920.1"), "role": "identity record only"} for g in genes]
    st.dataframe(genome_rows, width='stretch', hide_index=True)
    st.info("mtDNA entries remain catalogued as identity objects. No coordinate rendering, pathogenicity, dosing, risk, or eligibility claims are attached.")

with evidence_tab:
    st.markdown("**Evidence layer boundary**")
    st.warning("Evidence edges are review-required hypotheses. They require explicit evidence IDs and human review; they cannot carry claim status or mutate eligibility.")
    st.dataframe([{"edge state": "review-required", "epistemic status": "hypothesis", "evidence required": True, "linked claim IDs": [], "allowed action": "stage candidate only"}], width='stretch', hide_index=True)
    st.caption("Identity and context provenance resolve what an object is and where it came from. They are not mechanistic proof.")
    st.subheader("Titans Med literature feed · staged review inbox")
    feed_report = load_feed_report()
    staged_records = load_inbox_candidates()
    st.caption(feed_run_caption(feed_report))
    literature_rows = literature_candidate_rows(staged_records)
    edge_rows = edge_candidate_rows(staged_records)
    if literature_rows:
        st.dataframe(literature_rows, width='stretch', hide_index=True)
    if edge_rows:
        st.dataframe(edge_rows, width='stretch', hide_index=True)
    if literature_rows or edge_rows:
        st.caption("Candidates are staged automatically from Titans Med literature pulls. Nothing here is applied to the catalog, promoted to a claim, or mutated without explicit human review.")
    else:
        st.info("No staged literature candidates yet. Run the Titans Med literature feed after a literature pull to stage review-required studies for human review.")
    st.subheader("Investigation & triage")
    st.caption("Deterministic, explainable baseline (not an ML finding). Ranks staged candidates so human review happens once at the end of the pipeline, not per study.")
    stored_investigation = load_report()
    current_report = investigation_is_current(stored_investigation, staged_records)
    investigation = stored_investigation if current_report else build_investigation(staged_records)
    if stored_investigation and not current_report:
        st.caption("Investigation recalculated from the current inbox; the saved report had older or unverified inputs.")
    st.dataframe(investigation_metric_rows(investigation), width='stretch', hide_index=True)
    st.warning(investigation["metrics"].get("raw_count_warning", "Publication counts are inventory-only; evidence units are the downstream counting basis."))
    st.caption("Evidence units are the downstream counting basis; publications are shown as traceability inventory. Research-planning-only.")
    unit_rows = ranked_evidence_unit_rows(investigation)
    if unit_rows:
        st.dataframe(unit_rows, width='stretch', hide_index=True)
    ranked_rows = ranked_candidate_rows(investigation)
    if ranked_rows:
        st.subheader("Publication inventory · ranked candidates")
        st.dataframe(ranked_rows, width='stretch', hide_index=True)
        st.dataframe(cluster_rows(investigation), width='stretch', hide_index=True)
        contradictions = contradiction_rows(investigation)
        if contradictions:
            st.warning("Contradiction flags: " + "; ".join(f"{row['compound']} ({row['flag']})" for row in contradictions))
            st.dataframe(contradictions, width='stretch', hide_index=True)
        st.dataframe(synthesis_rows(investigation), width='stretch', hide_index=True)
        st.caption("Investigation output is review-required and research-planning-only; it does not mutate claims, eligibility, or catalog facts.")
        st.download_button("Download investigation report JSON", json.dumps(investigation, sort_keys=True), "molecular_atlas_investigation.json", "application/json", key="atlas_investigation_download")
    else:
        st.info("No staged candidates to investigate yet. Run the Titans Med literature feed first.")
    st.subheader("Mito next-assay packet")
    st.caption("One screening-only experiment from this feed and investigation. Claim eligibility is not mutated.")
    packet = build_next_assay_packet(investigation, feed_report=feed_report)
    packet_rows = packet_summary_rows(packet)
    if packet_rows:
        st.dataframe(packet_rows, width='stretch', hide_index=True)
        st.caption(packet["claim_boundary"])
    else:
        st.info("No next-assay packet could be ranked yet.")
    st.subheader("Atlas adjudication desk")
    st.caption("Role-conflict packets from the Mito evidence bundle. The local operator resolves study role automatically; labs can later turn human gates back on. Claim status and eligibility stay frozen.")
    try:
        delivery = load_mito_evidence_delivery(limit=200, persist_snapshot=False)
        if delivery.get("source") == "snapshot":
            st.warning(f"Saved evidence snapshot · age {delivery.get('age_seconds', 0):.0f} seconds. Live evidence is unavailable; this view may be out of date.")
        else:
            st.caption("Evidence delivery: live · " + str(delivery.get("generated_utc") or "timestamp unavailable"))
        desk = build_adjudication_desk(delivery["bundle"])
        operator_report = load_operator_report()
        st.write(
            f"Packets: {desk['metrics']['packets']} · "
            + ", ".join(f"{band} {desk['metrics']['bands'].get(band, 0)}" for band in ("critical", "high", "medium", "low"))
        )
        if operator_report:
            metrics = operator_report.get("metrics") or {}
            st.write(
                f"Operator: {metrics.get('auto_resolved', 0)} auto-resolved · "
                f"{metrics.get('needs_full_text', 0)} need full text · "
                f"{metrics.get('needs_human', 0)} held for a lab gate. "
                "No claim mutation."
            )
        rows = adjudication_desk_rows(desk, operator_report)
        if rows:
            st.dataframe(rows, width='stretch', hide_index=True)
        else:
            st.info("No role-conflict packets are queued for Atlas adjudication.")
    except (RuntimeError, ValueError) as exc:
        st.error(f"Atlas adjudication desk failed closed: {exc}")

with curation_tab:
    st.markdown("**Intelligent curation** · bounded review queue from declared catalog gaps")
    st.info("Model boundary: untrained deterministic baseline; not an ML finding. Candidates are hypotheses for human review only.")
    queue = build_molecular_curation_queue(catalog, consumers=CONSUMER_BINDINGS)
    if queue:
        st.dataframe(queue, width='stretch', hide_index=True)
        st.caption("No action is taken: this tab does not persist, apply, or mutate catalog facts, evidence edges, claim status, or eligibility.")
    else:
        st.info("No review-required curation candidates are present. No action is available.")

with sources_tab:
    st.markdown("**Source health** · refresh is explicit, bounded, and research-only")
    st.dataframe(source_rows, width='stretch', hide_index=True)
    st.subheader("Refresh one source record")
    st.caption("Nothing is fetched on page load. A refresh reads and writes only the supplied source cache; it never changes the catalog, consumer pins, or claim eligibility.")
    refresh_source_name = st.selectbox("Source", ["HGNC", "UniProt", "Ensembl"], key="atlas_refresh_source")
    refresh_identifier = st.text_input("Identifier", placeholder="HGNC symbol, UniProt accession, or Ensembl gene stable ID", key="atlas_refresh_identifier")
    if st.button("Refresh source", key="atlas_refresh_submit"):
        with st.spinner("Refreshing source record…"):
            st.session_state["atlas_source_refresh_report"] = refresh_source_review(refresh_source_name, refresh_identifier, cache_dir=DEFAULT_CACHE_DIR)
        st.rerun()
    refresh_report = st.session_state.get("atlas_source_refresh_report")
    if refresh_report:
        classification = refresh_report.get("classification", "unknown")
        if classification == "failed":
            st.error(f"Refresh failed: {refresh_report.get('error', 'unknown error')}")
        else:
            st.success(f"Refresh review: {classification}")
        st.json({key: value for key, value in refresh_report.items() if key not in {"before", "after"}})
        st.download_button("Download source refresh review JSON", json.dumps(refresh_report, sort_keys=True, indent=2), "molecular_atlas_source_refresh.json", "application/json", key="atlas_source_refresh_download")
    st.subheader("Assay identity catalog")
    st.caption("Assay, matrix, population, and endpoint pins. Variance packets cannot leave screening-only unless they match one of these identities.")
    assay_catalog = build_assay_catalog()
    st.caption(f"Catalog hash: {assay_catalog['content_hash']}")
    st.dataframe(assay_identity_rows(assay_catalog), width='stretch', hide_index=True)
    st.subheader("Source-change review")
    st.info("No active source changes are applied. Refresh output is a bounded review report only; additions, accession changes, transcript drift, assembly conflicts, deprecated symbols, provenance changes, and license updates remain review-required until a validated receipt and explicit human decision exist.")
    st.download_button("Download feed status JSON", json.dumps(feed_statuses, sort_keys=True, indent=2), "molecular_atlas_feed_status.json", "application/json", key="atlas_status_download")
