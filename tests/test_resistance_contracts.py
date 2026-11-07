"""Contract checks for ResistanceProfile, AntibioticCandidate results and Certificate.

The contracts are described in docs/contracts.md. These tests check the shape of the
fixtures and the rules that tie their fields together, so a lane that changes a contract
has to change the fixtures and this file in the same PR.
"""

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent / "fixtures"


def load(name):
    return json.loads((ROOT / f"{name}.json").read_text(encoding="utf-8"))


PROFILE = load("resistance_profile")
RESULT = load("antibiotic_candidates")
CERTS = load("certificate")

ARO = r"ARO:[0-9]{7}"
PROFILE_EVIDENCE = r"P[0-9]+"
CANDIDATE_EVIDENCE = r"E[0-9]+"
UUID = r"[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
RELEASE = r"[a-z][a-z0-9_]*-[0-9]+(?:\.[0-9]+)*"
INCHIKEY = r"[A-Z]{14}-[A-Z]{10}-[A-Z]"

INPUT_KINDS = {"species_preset", "gene_list", "amrfinderplus"}
CONFIDENCE = {"low", "medium", "high"}
EVASION = {"unaffected", "partly_affected", "defeated", "unknown"}
SOURCES = {"approved", "generated"}
CERT_STATUS = {"not_run", "certified", "failed", "unavailable"}
RULES = ["mw", "logp", "hbd", "hba", "rot_bonds", "psa"]
FEATURES = {"activity", "target_hit", "docking", "safety"}


def fullmatch(pattern, value):
    assert isinstance(value, str) and re.fullmatch(pattern, value), value


def keys(obj, expected):
    assert set(obj) == set(expected), (sorted(set(obj) ^ set(expected)), obj.get("name", ""))


# ---------- ResistanceProfile ----------


def test_profile_shape():
    keys(PROFILE, [
        "contract_version", "profile_id", "input", "unresolved_inputs", "releases", "genes",
        "mechanisms", "defeated_classes", "affected_targets", "evidence", "confidence", "caveats",
    ])
    assert PROFILE["contract_version"] == "1"
    fullmatch(UUID, PROFILE["profile_id"])
    assert PROFILE["input"]["kind"] in INPUT_KINDS
    assert PROFILE["confidence"]["level"] in CONFIDENCE
    assert PROFILE["confidence"]["basis"]
    assert PROFILE["caveats"], "a profile always carries caveats"
    for release in PROFILE["releases"]:
        fullmatch(RELEASE, release)


def test_profile_ids_and_references():
    evidence_ids = {e["id"] for e in PROFILE["evidence"]}
    assert len(evidence_ids) == len(PROFILE["evidence"])
    for e in PROFILE["evidence"]:
        fullmatch(PROFILE_EVIDENCE, e["id"])
        fullmatch(ARO, e["aro_id"])
        assert e["release"] in PROFILE["releases"]

    gene_ids = {g["aro_id"] for g in PROFILE["genes"]}
    assert len(gene_ids) == len(PROFILE["genes"])
    for g in PROFILE["genes"]:
        fullmatch(ARO, g["aro_id"])
        assert g["mechanisms"], f"{g['name']} resolves to no mechanism"
        assert set(g["evidence_ids"]) <= evidence_ids and g["evidence_ids"]
        for term in g["mechanisms"] + g["defeated_classes"]:
            fullmatch(ARO, term["aro_id"])

    for group in ("mechanisms", "defeated_classes", "affected_targets"):
        for item in PROFILE[group]:
            assert set(item["gene_aro_ids"]) <= gene_ids
            assert set(item.get("evidence_ids", [])) <= evidence_ids


def test_profile_summaries_agree_with_genes():
    """The roll-ups must be derivable from the per-gene rows, with no extra entries."""
    classes = {}
    mechanisms = {}
    for g in PROFILE["genes"]:
        for c in g["defeated_classes"]:
            classes.setdefault(c["aro_id"], set()).add(g["aro_id"])
        for m in g["mechanisms"]:
            mechanisms.setdefault(m["aro_id"], set()).add(g["aro_id"])
    assert {c["aro_id"]: set(c["gene_aro_ids"]) for c in PROFILE["defeated_classes"]} == classes
    assert {m["aro_id"]: set(m["gene_aro_ids"]) for m in PROFILE["mechanisms"]} == mechanisms


def test_unresolved_inputs_are_reported_not_dropped():
    for item in PROFILE["unresolved_inputs"]:
        assert item["text"] in PROFILE["input"]["requested"]
        assert item["reason"]


# ---------- Certificate ----------


def all_certificates():
    yield from CERTS["certificates"]
    for c in RESULT["candidates"]:
        if c["certificate"] is not None:
            yield c["certificate"]


@pytest.mark.parametrize("cert", list(all_certificates()), ids=lambda c: f"{c['inchikey']}-{c['ran']}")
def test_certificate_rules_and_flags(cert):
    keys(cert, [
        "inchikey", "rule_set", "ran", "certified", "rules", "failed_rules", "lean_version",
        "seconds", "detail",
    ])
    fullmatch(INCHIKEY, cert["inchikey"])
    assert [r["rule"] for r in cert["rules"]] == RULES
    for r in cert["rules"]:
        assert r["margin"] == pytest.approx(r["value"] - r["limit"], abs=0.01)
        assert r["passed"] == (r["margin"] <= 0)
    if cert["ran"]:
        assert cert["failed_rules"] == [r["rule"] for r in cert["rules"] if not r["passed"]]
        assert cert["certified"] == (not cert["failed_rules"])
        assert cert["lean_version"]
    else:
        # Never "certified" when Lean did not run, and say why.
        assert cert["certified"] is False
        assert cert["failed_rules"] == []
        assert cert["lean_version"] is None
        assert cert["detail"]


def test_certificate_limits_match_the_python_rule_set():
    """The fixture limits must be the real antibiotic_relaxed numbers, not invented ones."""
    from htb_verify.rules import RULE_SETS

    t = RULE_SETS["antibiotic_relaxed"]
    expected = {
        "mw": t.mw_x10 / 10, "logp": t.logp_x100 / 100, "hbd": t.hbd, "hba": t.hba,
        "rot_bonds": t.rot_bonds, "psa": t.psa_x10 / 10,
    }
    for cert in all_certificates():
        assert cert["rule_set"] == "antibiotic_relaxed"
        assert {r["rule"]: r["limit"] for r in cert["rules"]} == expected


def test_standalone_certificates_cover_each_outcome():
    outcomes = {(c["ran"], c["certified"]) for c in CERTS["certificates"]}
    assert outcomes == {(True, True), (True, False), (False, False)}


# ---------- AntibioticCandidate result ----------


def test_result_shape():
    keys(RESULT, [
        "contract_version", "result_id", "profile_id", "pathogen", "releases", "ranking", "candidates",
    ])
    fullmatch(UUID, RESULT["result_id"])
    assert RESULT["profile_id"] == PROFILE["profile_id"]
    assert set(RESULT["ranking"]["order"]) <= FEATURES | {"evasion"}
    for release in RESULT["releases"]:
        fullmatch(RELEASE, release)


CANDIDATES = RESULT["candidates"]


@pytest.mark.parametrize("c", CANDIDATES, ids=lambda c: c["name"])
def test_candidate_fields(c):
    keys(c, [
        "candidate_id", "name", "source", "smiles", "inchikey", "drug_class", "mechanism_of_action",
        "target", "rank", "evasion", "features", "explanation", "smoke_test", "certificate_status",
        "certificate", "generated", "evidence",
    ])
    assert c["source"] in SOURCES
    assert c["evasion"]["status"] in EVASION
    assert c["certificate_status"] in CERT_STATUS
    assert set(c["features"]) == FEATURES
    fullmatch(INCHIKEY, c["inchikey"])
    for f in c["features"].values():
        assert set(f) == {"available", "value", "unit", "evidence_ids"}
        if not f["available"]:
            assert f["value"] is None and f["evidence_ids"] == []
    assert c["explanation"]["summary"]
    assert "missing" in c["explanation"] and "caveats" in c["explanation"]


@pytest.mark.parametrize("c", CANDIDATES, ids=lambda c: c["name"])
def test_candidate_references_resolve(c):
    own = {e["id"] for e in c["evidence"]}
    assert len(own) == len(c["evidence"])
    for e in c["evidence"]:
        fullmatch(CANDIDATE_EVIDENCE, e["id"])
        assert e["statement"]
    for f in c["features"].values():
        assert set(f["evidence_ids"]) <= own

    profile_evidence = {e["id"] for e in PROFILE["evidence"]}
    assert set(c["evasion"]["profile_evidence_ids"]) <= profile_evidence

    genes = {g["aro_id"]: g for g in PROFILE["genes"]}
    for hit in c["evasion"]["defeated_by"]:
        assert hit["gene_aro_id"] in genes
        assert hit["mechanism_aro_id"] in {m["aro_id"] for m in genes[hit["gene_aro_id"]]["mechanisms"]}


@pytest.mark.parametrize("c", CANDIDATES, ids=lambda c: c["name"])
def test_evasion_status_agrees_with_the_profile(c):
    """A candidate is `defeated` exactly when its class is in the profile's defeated classes."""
    defeated_classes = {d["aro_id"] for d in PROFILE["defeated_classes"]}
    drug_class = c["drug_class"]
    status = c["evasion"]["status"]
    if drug_class is None:
        assert status == "unknown"
    elif drug_class["aro_id"] is None:
        assert status in {"unaffected", "unknown"}
    elif drug_class["aro_id"] in defeated_classes:
        assert status in {"defeated", "partly_affected"}
        assert c["evasion"]["defeated_by"], "a defeated class must name the gene that defeats it"
    else:
        assert status == "unaffected"
    if status in {"unaffected", "unknown"}:
        assert c["evasion"]["defeated_by"] == []


def test_ranks_cover_approved_candidates_only():
    approved = [c for c in CANDIDATES if c["source"] == "approved"]
    generated = [c for c in CANDIDATES if c["source"] == "generated"]
    assert sorted(c["rank"] for c in approved) == list(range(1, len(approved) + 1))
    # Generated molecules are never ranked with approved drugs.
    assert all(c["rank"] is None for c in generated)
    # The ordered rule puts evasion first: no defeated candidate outranks an unaffected one.
    order = {"unaffected": 0, "partly_affected": 1, "defeated": 2}
    by_rank = sorted(approved, key=lambda c: c["rank"])
    levels = [order[c["evasion"]["status"]] for c in by_rank]
    assert levels == sorted(levels)


def test_generated_candidates_carry_provenance():
    for c in CANDIDATES:
        if c["source"] == "generated":
            assert c["generated"] and c["generated"]["model"]
            assert "seed" in c["generated"] and "model_version" in c["generated"]
        else:
            assert c["generated"] is None


@pytest.mark.parametrize("c", CANDIDATES, ids=lambda c: c["name"])
def test_smoke_test_and_certificate_agree(c):
    smoke = c["smoke_test"]
    assert smoke["valid"]
    # A violation list and `passed` must say the same thing.
    assert smoke["passed"] == (not smoke["violations"])
    for v in smoke["violations"]:
        assert v["rule"] in RULES + ["sa", "alert"]
        assert v["margin"] == pytest.approx(v["value"] - v["limit"], abs=0.01)

    status, cert = c["certificate_status"], c["certificate"]
    if status == "not_run":
        assert cert is None
        return
    assert cert["inchikey"] == c["inchikey"]
    assert cert["rule_set"] == smoke["rule_set"]
    if status == "certified":
        assert cert["ran"] and cert["certified"] and smoke["passed"]
    elif status == "failed":
        assert cert["ran"] and not cert["certified"]
        # Python and Lean must name the same failing property rules.
        assert {v["rule"] for v in smoke["violations"] if v["rule"] in RULES} == set(cert["failed_rules"])
    elif status == "unavailable":
        assert not cert["ran"] and not cert["certified"]


def test_embedded_certificates_match_the_standalone_file():
    standalone = {(c["inchikey"], c["ran"]): c for c in CERTS["certificates"]}
    for c in CANDIDATES:
        cert = c["certificate"]
        if cert is not None:
            assert standalone[(cert["inchikey"], cert["ran"])] == cert
