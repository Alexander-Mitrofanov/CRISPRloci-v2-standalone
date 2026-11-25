"""Dependency-light authority for the executed self-discovery settings."""

from dataclasses import dataclass

DISCOVERY_POLICY = "vmatch-source-discovery-v1"
SEED_TEXT_LIMIT = 55
MAX_LOCAL_CHAIN_SPACER_LENGTH = 140
RESULT_LIMIT = 1_000_000
FUZZY_ERROR = "{i<=3,d<=3,s<=3,i+d+s<=6}"
FAST_SEEDS = {"balanced": (23, -10, 100), "legacy": (23, 25, 60),
              "sensitive": (21, -10, 100)}


@dataclass(frozen=True)
class VmatchContract:
    requested_profile: str
    profile_applies: bool
    l_arguments: tuple
    e: int = 1
    evalue: int = 1
    result_limit: int = RESULT_LIMIT


def resolve_vmatch_contract(fast_run, profile="balanced"):
    if type(fast_run) is not bool or profile not in FAST_SEEDS:
        raise ValueError("fast_run must be boolean and profile balanced, legacy, or sensitive")
    return VmatchContract(profile, fast_run, FAST_SEEDS[profile] if fast_run else (21, -10, 100))


def effective_discovery(fast_run, profile, flank, engine):
    contract = resolve_vmatch_contract(fast_run, profile)
    return {
        "policy": DISCOVERY_POLICY,
        "requested_fast_profile": contract.requested_profile,
        "fast_profile_applies": contract.profile_applies,
        "vmatch": {"-l": list(contract.l_arguments), "-e": contract.e,
                   "-evalue": contract.evalue, "result_limit": contract.result_limit},
        "seed_text_truncation": SEED_TEXT_LIMIT,
        "local_chain_max_gap": MAX_LOCAL_CHAIN_SPACER_LENGTH,
        "fuzzy_flank": flank, "fuzzy_engine": engine, "fuzzy_error": FUZZY_ERROR,
    }
