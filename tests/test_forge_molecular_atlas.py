import copy

import pytest

from forge_genomics_protein_spine import finalize_bundle
from forge_molecular_atlas import (
    PROGRAMS,
    assess_consumer_receipt,
    build_catalog_snapshot,
    build_program_projection_manifest,
    build_program_projection_manifests,
    make_consumer_receipt,
    make_evidence_edge_candidate,
    pin_consumer_receipt,
    stage_candidate_diff,
    validate_evidence_edge_candidate,
    validate_program_projection_manifest,
    validate_catalog_snapshot,
)


def test_all_program_feeds_are_deterministic_and_hash_addressed():
    first = build_program_projection_manifests()
    second = build_program_projection_manifests()
    assert first == second
    assert tuple(first) == PROGRAMS
    assert all(validate_program_projection_manifest(item) for item in first.values())
    assert {item["program"] for item in first.values()} == set(PROGRAMS)
    assert all(item["projection_hash"] == item["projection"]["content_hash"] for item in first.values())


def test_catalog_snapshot_publishes_immutable_spine_metadata():
    snapshot = build_catalog_snapshot()
    assert snapshot["spine_hash"]
    assert snapshot["snapshot_hash"]
    assert snapshot["research_planning_only"] is True
    assert snapshot["claim_eligibility_mutation"] is False


def test_catalog_snapshot_roundtrip_and_mutation_detection():
    snapshot = build_catalog_snapshot()
    assert validate_catalog_snapshot(snapshot)
    changed = copy.deepcopy(snapshot)
    changed["accessed"] = changed["accessed"] + "-changed"
    with pytest.raises(ValueError, match="snapshot hash"):
        validate_catalog_snapshot(changed)


def test_projection_hash_drift_is_detected():
    manifest = build_program_projection_manifest("mito")
    changed = copy.deepcopy(manifest)
    changed["projection"]["entities"][0]["label"] += " changed"
    with pytest.raises(ValueError, match="projection hash"):
        validate_program_projection_manifest(changed)


@pytest.mark.parametrize("status", ["unpinned", "stale", "incompatible", "blocked"])
def test_receipt_statuses_fail_closed(status):
    manifest = build_program_projection_manifest("neuro")
    receipt = make_consumer_receipt("neuro")
    if status == "stale":
        receipt = make_consumer_receipt("neuro", projection_hash="0" * 64)
    elif status == "incompatible":
        receipt = make_consumer_receipt("mito")
    elif status == "blocked":
        receipt = make_consumer_receipt("neuro", status="blocked")
    assert assess_consumer_receipt(receipt, manifest) == status


def test_receipt_can_be_pinned_without_mutating_original():
    manifest = build_program_projection_manifest("cardio")
    receipt = make_consumer_receipt("cardio")
    pinned = pin_consumer_receipt(receipt, manifest)
    assert assess_consumer_receipt(pinned, manifest) == "ready"
    assert receipt["projection_hash"] is None


def test_candidate_diff_is_staged_and_does_not_mutate_active_snapshot():
    active = build_program_projection_manifest("onco")["projection"]
    proposed = copy.deepcopy(active)
    proposed["entities"][0]["label"] += " (candidate)"
    finalize_bundle(proposed)
    diff = stage_candidate_diff(active, proposed)
    assert diff["status"] == "review-required"
    assert diff["changed_entity_ids"]
    assert diff["active_hash"] != diff["proposed_hash"]
    assert active == build_program_projection_manifest("onco")["projection"]


def test_evidence_edge_candidate_is_review_only_and_validates_identity_ids():
    candidate = make_evidence_edge_candidate("gene:HGNC:613", "protein:P02649", "candidate_supports", evidence_ids=["EVID-1"])
    assert validate_evidence_edge_candidate(candidate)
    assert candidate["linked_claim_ids"] == []
    invalid = dict(candidate, target="gene:missing")
    with pytest.raises(ValueError, match="unknown identity"):
        validate_evidence_edge_candidate(invalid)
    with pytest.raises(ValueError, match="not review-required"):
        validate_evidence_edge_candidate(dict(candidate, claim_status="eligible"))


def test_candidate_does_not_mutate_claim_eligibility():
    candidate = make_evidence_edge_candidate("gene:HGNC:14581", "protein:Q9BXM7", "candidate_context")
    assert candidate["claim_eligibility_mutation"] is False
    assert "claim_status" not in candidate
