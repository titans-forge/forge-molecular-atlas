"""Research-only data plane for the shared genomics/protein identity spine.

The atlas is a consumer of :mod:`forge_genomics_protein_spine`.  It owns no
identity facts and cannot change claim status or eligibility.  All returned
records are plain standard-library values so callers can serialize receipts
without introducing a persistence or network dependency.
"""
from __future__ import annotations

import copy
import hashlib
import json
from typing import Any, Mapping

from forge_genomics_protein_spine import (
    ALLOWED_PROGRAMS,
    BUNDLE_SCHEMA,
    build_genomics_protein_spine,
    catalog_snapshot,
    program_projection,
    validate_genomics_protein_spine,
)

ATLAS_SCHEMA = "titans-med-molecular-atlas/v1"
PROGRAMS = ("mito", "neuro", "pgx", "cardio", "onco")
RECEIPT_STATUSES = frozenset({"ready", "unpinned", "stale", "incompatible", "blocked"})
REVIEW_REQUIRED = "review-required"
RESEARCH_PLANNING_ONLY = True


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def content_hash(value: Mapping[str, Any]) -> str:
    """Return the SHA-256 of a canonical mapping, excluding its hash field."""
    material = {key: value[key] for key in value if key not in {"content_hash", "manifest_hash", "projection_hash"}}
    return hashlib.sha256(_json(material).encode("utf-8")).hexdigest()


def _snapshot_hash(snapshot: Mapping[str, Any]) -> str:
    """Hash snapshot content without changing the established manifest hash."""
    material = {key: snapshot[key] for key in snapshot if key != "snapshot_hash"}
    return hashlib.sha256(_json(material).encode("utf-8")).hexdigest()


def _checked_catalog(catalog: Mapping[str, Any] | None) -> dict[str, Any]:
    value = build_genomics_protein_spine() if catalog is None else copy.deepcopy(dict(catalog))
    validate_genomics_protein_spine(value)
    return value


def build_catalog_snapshot(catalog: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Publish immutable metadata for one validated identity catalog."""
    value = _checked_catalog(catalog)
    metadata = catalog_snapshot(value)
    result = {
        "schema": ATLAS_SCHEMA,
        "spine_schema": metadata["schema"],
        "accessed": metadata["accessed"],
        "spine_hash": metadata["content_hash"],
        "counts": copy.deepcopy(metadata["counts"]),
        "seeded_variants": metadata["seeded_variants"],
        "seeded_haplotypes": metadata["seeded_haplotypes"],
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
    }
    result["snapshot_hash"] = _snapshot_hash(result)
    return result


def validate_catalog_snapshot(snapshot: Mapping[str, Any]) -> bool:
    required = {"schema", "spine_schema", "spine_hash", "counts", "snapshot_hash"}
    if not isinstance(snapshot, Mapping) or not required <= set(snapshot) or snapshot.get("schema") != ATLAS_SCHEMA:
        raise ValueError("invalid catalog snapshot")
    if snapshot.get("research_planning_only") is not True or snapshot.get("claim_eligibility_mutation") is not False:
        raise ValueError("catalog snapshot must be planning-only")
    if _snapshot_hash(snapshot) != snapshot["snapshot_hash"]:
        raise ValueError("snapshot hash mismatch")
    return True


def build_program_projection_manifest(
    program: str, catalog: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Build a hash-addressed projection manifest for one allowed program."""
    if program not in PROGRAMS:
        raise ValueError(f"unknown program: {program}")
    value = _checked_catalog(catalog)
    projection = program_projection(value, program)
    manifest = {
        "schema": ATLAS_SCHEMA,
        "program": program,
        "spine_schema": BUNDLE_SCHEMA,
        "spine_hash": value["content_hash"],
        "projection": projection,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
    }
    manifest["projection_hash"] = content_hash(projection)
    manifest["manifest_hash"] = content_hash(manifest)
    return manifest


def build_program_projection_manifests(catalog: Mapping[str, Any] | None = None) -> dict[str, dict[str, Any]]:
    value = _checked_catalog(catalog)
    return {program: build_program_projection_manifest(program, value) for program in PROGRAMS}


def validate_program_projection_manifest(manifest: Mapping[str, Any]) -> bool:
    required = {"schema", "program", "spine_schema", "spine_hash", "projection", "projection_hash", "manifest_hash"}
    if not isinstance(manifest, Mapping) or not required <= set(manifest):
        raise ValueError("incomplete projection manifest")
    if manifest["schema"] != ATLAS_SCHEMA or manifest["program"] not in PROGRAMS:
        raise ValueError("incompatible projection manifest")
    projection = manifest["projection"]
    if content_hash(projection) != manifest["projection_hash"]:
        raise ValueError("projection hash mismatch")
    validate_genomics_protein_spine(projection)
    if content_hash(manifest) != manifest["manifest_hash"]:
        raise ValueError("manifest hash mismatch")
    if manifest.get("research_planning_only") is not True or manifest.get("claim_eligibility_mutation") is not False:
        raise ValueError("projection must be planning-only")
    return True


def make_consumer_receipt(
    program: str,
    *,
    projection_hash: str | None = None,
    manifest_hash: str | None = None,
    spine_hash: str | None = None,
    status: str | None = None,
    consumer: str = "molecular-atlas",
) -> dict[str, Any]:
    """Create a receipt; a missing projection pin is intentionally unpinned."""
    if program not in PROGRAMS:
        raise ValueError(f"unknown program: {program}")
    if status is not None and status not in RECEIPT_STATUSES:
        raise ValueError("invalid receipt status")
    return {
        "schema": ATLAS_SCHEMA,
        "consumer": consumer,
        "program": program,
        "projection_hash": projection_hash,
        "manifest_hash": manifest_hash,
        "spine_hash": spine_hash,
        "status": status or ("unpinned" if not projection_hash else "ready"),
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
    }


def assess_consumer_receipt(receipt: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    """Return a fail-closed status without changing either input."""
    if not isinstance(receipt, Mapping) or receipt.get("program") != manifest.get("program"):
        return "incompatible"
    if receipt.get("status") == "blocked":
        return "blocked"
    try:
        validate_program_projection_manifest(manifest)
    except (TypeError, ValueError):
        return "incompatible"
    if not receipt.get("projection_hash"):
        return "unpinned"
    if receipt.get("projection_hash") != manifest["projection_hash"]:
        return "stale"
    if receipt.get("manifest_hash") not in {None, manifest["manifest_hash"]}:
        return "stale"
    if receipt.get("spine_hash") not in {None, manifest["spine_hash"]}:
        return "stale"
    if receipt.get("research_planning_only") is not True or receipt.get("claim_eligibility_mutation") is not False:
        return "blocked"
    return "ready"


def pin_consumer_receipt(receipt: Mapping[str, Any], manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Return a new receipt pinned to a validated manifest."""
    validate_program_projection_manifest(manifest)
    if receipt.get("program") != manifest["program"]:
        raise ValueError("receipt and manifest program mismatch")
    result = dict(receipt)
    result.update({"projection_hash": manifest["projection_hash"], "manifest_hash": manifest["manifest_hash"], "spine_hash": manifest["spine_hash"], "status": "ready"})
    return result


def validate_consumer_receipt(receipt: Mapping[str, Any]) -> bool:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != ATLAS_SCHEMA or receipt.get("program") not in PROGRAMS:
        raise ValueError("invalid consumer receipt")
    if receipt.get("status") not in RECEIPT_STATUSES:
        raise ValueError("invalid consumer receipt status")
    if receipt.get("research_planning_only") is not True or receipt.get("claim_eligibility_mutation") is not False:
        raise ValueError("receipt must be planning-only")
    return True


def stage_candidate_diff(active: Mapping[str, Any], proposed: Mapping[str, Any]) -> dict[str, Any]:
    """Describe a proposed projection change without mutating the active feed."""
    if "projection" in active and "projection" in proposed:
        active = active["projection"]
        proposed = proposed["projection"]
    validate_genomics_protein_spine(active)
    validate_genomics_protein_spine(proposed)
    active_entities = {item["id"]: item for item in active["entities"]}
    proposed_entities = {item["id"]: item for item in proposed["entities"]}
    active_edges = {(item["source"], item["target"], item["relation"]): item for item in active["edges"]}
    proposed_edges = {(item["source"], item["target"], item["relation"]): item for item in proposed["edges"]}
    result = {
        "schema": ATLAS_SCHEMA,
        "status": REVIEW_REQUIRED,
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
        "active_hash": content_hash(active),
        "proposed_hash": content_hash(proposed),
        "added_entity_ids": sorted(set(proposed_entities) - set(active_entities)),
        "removed_entity_ids": sorted(set(active_entities) - set(proposed_entities)),
        "changed_entity_ids": sorted(key for key in set(active_entities) & set(proposed_entities) if active_entities[key] != proposed_entities[key]),
        "added_edges": [copy.deepcopy(proposed_edges[key]) for key in sorted(set(proposed_edges) - set(active_edges))],
        "removed_edges": [copy.deepcopy(active_edges[key]) for key in sorted(set(active_edges) - set(proposed_edges))],
    }
    result["diff_hash"] = content_hash(result)
    return result


def make_evidence_edge_candidate(
    source: str, target: str, relation: str, *, evidence_ids: list[str] | tuple[str, ...] = (), rationale: str = ""
) -> dict[str, Any]:
    """Create a review-required edge proposal; candidates cannot carry claims."""
    if not all(isinstance(value, str) and value.strip() for value in (source, target, relation)):
        raise ValueError("evidence edge identifiers and relation are required")
    if source == target or not isinstance(evidence_ids, (list, tuple)) or any(not isinstance(item, str) or not item for item in evidence_ids):
        raise ValueError("invalid evidence edge")
    return {
        "schema": ATLAS_SCHEMA,
        "record_type": "evidence_edge_candidate",
        "source": source,
        "target": target,
        "relation": relation,
        "evidence_ids": list(evidence_ids),
        "rationale": str(rationale),
        "status": REVIEW_REQUIRED,
        "epistemic_status": "hypothesis",
        "evidence_required": True,
        "linked_claim_ids": [],
        "research_planning_only": True,
        "claim_eligibility_mutation": False,
    }


def validate_evidence_edge_candidate(candidate: Mapping[str, Any], catalog: Mapping[str, Any] | None = None) -> bool:
    if not isinstance(candidate, Mapping) or candidate.get("record_type") != "evidence_edge_candidate":
        raise ValueError("invalid evidence edge candidate")
    forbidden = {"claim_status", "eligibility", "automatic_status_change", "claim_ids"}
    if forbidden & set(candidate) or candidate.get("status") != REVIEW_REQUIRED or candidate.get("epistemic_status") != "hypothesis":
        raise ValueError("evidence edge is not review-required")
    value = _checked_catalog(catalog)
    known = {entity["id"] for entity in value["entities"]}
    if candidate.get("source") not in known or candidate.get("target") not in known or candidate["source"] == candidate["target"]:
        raise ValueError("evidence edge references unknown identity")
    if candidate.get("linked_claim_ids") != [] or candidate.get("claim_eligibility_mutation") is not False:
        raise ValueError("candidate cannot mutate claims")
    return True


# Short aliases for consumers that use the data-plane vocabulary directly.
catalog_manifest = build_catalog_snapshot
projection_manifest = build_program_projection_manifest
projection_manifests = build_program_projection_manifests
consumer_receipt_status = assess_consumer_receipt
candidate_diff = stage_candidate_diff

__all__ = [name for name in globals() if not name.startswith("_")]
