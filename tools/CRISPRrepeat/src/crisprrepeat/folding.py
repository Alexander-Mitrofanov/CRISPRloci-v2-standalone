"""Exact numeric ViennaRNA interfaces; local probabilities are never row-normalized."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import math
from pathlib import Path
import tempfile

from .inputs import MAX_CONTEXT_LENGTH, MAX_REPEAT_LENGTH

MAX_PAIRS = 250000


@dataclass(frozen=True)
class Settings:
    temperature: float = 37.0
    window: int = 150
    span: int = 80
    no_lonely_pairs: bool = False

    def __post_init__(self):
        if type(self.temperature) not in {int, float}:
            raise ValueError("Temperature must be numeric")
        if type(self.window) is not int or type(self.span) is not int:
            raise ValueError("Window and span must be integers")
        if type(self.no_lonely_pairs) is not bool:
            raise ValueError("Lonely-pair policy must be boolean")
        if not math.isfinite(self.temperature) or not 0 <= self.temperature <= 100:
            raise ValueError("Temperature must be between 0 and 100 Celsius")
        if not 10 <= self.window <= 300 or not 5 <= self.span <= self.window:
            raise ValueError("Require 10 <= window <= 300 and 5 <= span <= window")


def probability(value):
    value = float(value)
    if not math.isfinite(value) or not -1e-8 <= value <= 1 + 1e-8:
        raise ValueError("Invalid probability returned by ViennaRNA")
    return min(1.0, max(0.0, value))


def structure_pairs(structure):
    stack, pairs = [], []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            pairs.append([stack.pop(), i])
        elif char != ".":
            raise ValueError("Invalid dot-bracket structure")
    if stack:
        raise ValueError("Unbalanced dot-bracket structure")
    return sorted(pairs)


class Folder:
    def __init__(self, settings=Settings()):
        import RNA

        if RNA.__version__ != "2.7.0":
            raise ValueError("This release requires ViennaRNA 2.7.0")
        self.RNA, self.settings = RNA, settings
        RNA.params_load_RNA_Turner2004()
        model = self.model()
        self._parameters = RNA.param(model)
        self._partition_parameters = RNA.exp_param(model)
        with tempfile.TemporaryDirectory(prefix="crisprrepeat-parameters-") as directory:
            path = Path(directory) / "parameters.par"
            RNA.params_save(str(path))
            parameter_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        self.provenance = {
            "engine": "ViennaRNA",
            "version": RNA.__version__,
            "parameters": "RNA_Turner2004",
            "parameter_sha256": parameter_sha,
            "dangles": 2,
            "salt_molar": 1.021,
            "minimum_loop_size": 3,
            "allow_GU": True,
            "pseudoknots": False,
            "numeric_cutoff": 0.0,
            "pair_coordinates": "zero-based indices; maximum span includes endpoints",
            **asdict(settings),
        }
        self.cache = {}

    def model(self):
        # Constructor keywords rebuild ViennaRNA's internal pair tables, unlike
        # assigning energy_set after construction. Ignore process-global defaults.
        # Isolation is unrestricted; context() applies the requested local bounds.
        return self.RNA.md(
            temperature=self.settings.temperature,
            dangles=2, noLP=int(self.settings.no_lonely_pairs), salt=1.021,
            noGU=0, noGUclosure=0, gquad=0, circ=0,
            min_loop_size=3, compute_bpp=1,
            energy_set=0, logML=0, special_hp=1,
            betaScale=1.0, pf_smooth=1, sfact=1.07, backtrack=1,
            max_bp_span=-1, window_size=-1,
        )

    def _compound(self, sequence, model, options=0):
        """Copy owned Turner2004 energy/Boltzmann tables before native computation."""
        compound = self.RNA.fold_compound(sequence, model, options)
        compound.params_subst(self._parameters)
        compound.exp_params_subst(self._partition_parameters)
        # Adjust copied model details; owned tables remain unchanged between folds.
        for parameters in (compound.params, compound.exp_params):
            parameters.model_details.window_size = model.window_size
            parameters.model_details.max_bp_span = model.max_bp_span
        return compound

    def isolated(self, sequence):
        if not isinstance(sequence, str) or not 1 <= len(sequence) <= MAX_REPEAT_LENGTH:
            raise ValueError("Isolated repeat must contain 1–200 nucleotides")
        if set(sequence) - set("ACGU"):
            return {"status": "unsupported", "reason": "ambiguous_bases"}
        if sequence in self.cache:
            return self.cache[sequence]
        fc = self._compound(sequence, self.model())
        structure, mfe = fc.mfe()
        if type(mfe) not in {int, float} or not math.isfinite(mfe):
            raise ValueError("Nonfinite or invalid folding energy")
        if not isinstance(structure, str) or len(structure) != len(sequence):
            raise ValueError("Incomplete MFE structure")
        motif = structure_pairs(structure)
        if any(j - i < 4 or sequence[i] + sequence[j] not in {"AU", "UA", "GC", "CG", "GU", "UG"}
               for i, j in motif):
            raise ValueError("MFE pair violates folding model")
        fc.exp_params_rescale(mfe)
        _, ensemble = fc.pf()
        if type(ensemble) not in {int, float} or not math.isfinite(ensemble):
            raise ValueError("Nonfinite or invalid folding energy")
        matrix = fc.bpp()
        pairs = []
        for i in range(1, len(sequence) + 1):
            for j in range(i + 1, len(sequence) + 1):
                try:
                    p = probability(matrix[i][j])
                except (IndexError, TypeError) as error:
                    raise ValueError("Incomplete global pair probabilities") from error
                if p > 0:
                    pairs.append([i - 1, j - 1, p])
        unpaired = [1.0] * len(sequence)
        for i, j, p in pairs:
            unpaired[i] -= p
            unpaired[j] -= p
        result = {
            "status": "completed",
            "method": "global_partition_function",
            "structure": structure,
            "mfe_kcal_mol": float(mfe),
            "ensemble_kcal_mol": float(ensemble),
            "pairs": pairs,
            "unpaired": [probability(p) for p in unpaired],
            "motif_pairs": motif,
        }
        self.cache[sequence] = result
        return result

    def context(self, sequence):
        if not isinstance(sequence, str) or not 1 <= len(sequence) <= MAX_CONTEXT_LENGTH:
            raise ValueError("Array context must contain 1–5000 nucleotides")
        n = len(sequence)
        window = min(n, self.settings.window)
        span = min(window, self.settings.span)
        metadata = {
            "method": "local_window_pair_averages",
            "window": window,
            "span": span,
            "requested_window": self.settings.window,
            "requested_span": self.settings.span,
        }
        if set(sequence) - set("ACGU"):
            return {"status": "unsupported", "reason": "ambiguous_bases", **metadata}
        if n < 5:
            return {"status": "unsupported", "reason": "context_too_short", **metadata}
        md = self.model()
        md.window_size, md.max_bp_span = window, span
        fc = self._compound(sequence, md, self.RNA.OPTION_WINDOW)
        pairs, unpaired, rows, accessibility_rows, errors = [], [None] * n, set(), set(), []

        def callback(pr, size, i, maximum, kind, _):
            # SWIG can swallow callback exceptions; collect errors and fail after return.
            if errors:
                return
            try:
                if type(i) is not int or not 1 <= i <= n:
                    raise ValueError("Invalid local probability row")
                if kind & self.RNA.PROBS_WINDOW_BPP:
                    if i in rows or type(size) is not int or not i <= size <= n:
                        raise ValueError("Duplicate or incomplete local pair row")
                    rows.add(i)
                    for j in range(i + 1, size + 1):
                        if pr[j] is None:
                            raise ValueError("Missing local pair probability")
                        p = probability(pr[j])
                        if p > 0:
                            if j - i + 1 > span:
                                raise ValueError("Local pair exceeds requested span")
                            if len(pairs) >= MAX_PAIRS:
                                raise ValueError("Pair output limit exceeded")
                            pairs.append([i - 1, j - 1, p])
                elif kind & self.RNA.PROBS_WINDOW_UP:
                    if i in accessibility_rows or size < 1:
                        raise ValueError("Duplicate or incomplete accessibility row")
                    accessibility_rows.add(i)
                    unpaired[i - 1] = probability(pr[1])
                else:
                    raise ValueError("Unexpected local probability kind")
            except Exception as error:
                errors.append(str(error))

        success = fc.probs_window(1, self.RNA.PROBS_WINDOW_BPP | self.RNA.PROBS_WINDOW_UP, callback, None)
        if not success or errors or len(rows) != n or len(accessibility_rows) != n:
            raise ValueError("Incomplete local probabilities: " + "; ".join(errors[:1]))
        return {"status": "completed", **metadata, "pairs": pairs, "unpaired": unpaired}


def comparison(isolated, context, offset=0, context_pairs=None):
    base = {
        "status": "unsupported",
        "reason": None,
        "motif_origin": "isolated_mfe_hypothesis",
        "motif_pairs": [],
        "isolated_pair_support": None,
        "context_pair_support": None,
        "support_change": None,
    }
    if isolated["status"] != "completed":
        return {**base, "reason": isolated.get("reason", "isolated_unavailable")}
    motif = isolated["motif_pairs"]
    base["motif_pairs"] = motif
    if not motif:
        return {**base, "reason": "no_reference_pairs"}
    p = {(i, j): value for i, j, value in isolated["pairs"]}
    base["isolated_pair_support"] = sum(p.get(tuple(pair), 0.0) for pair in motif) / len(motif)
    if context is None or context["status"] != "completed":
        return {**base, "reason": "context_not_supplied" if context is None else context.get("reason")}
    if any(j - i + 1 > context["span"] for i, j in motif):
        return {**base, "reason": "motif_exceeds_context_span"}
    q = context_pairs if context_pairs is not None else {(i, j): value for i, j, value in context["pairs"]}
    support = sum(q.get((i + offset, j + offset), 0.0) for i, j in motif) / len(motif)
    return {
        **base,
        "status": "completed",
        "context_pair_support": support,
        "support_change": support - base["isolated_pair_support"],
    }
