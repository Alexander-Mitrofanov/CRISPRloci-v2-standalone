"""Cached runtime for the trusted frozen cassette architecture model."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import joblib
import numpy as np

from casandra.cassette_architecture import cassette_features


@dataclass(frozen=True)
class ArchitectureModel:
    """A trusted Joblib estimator loaded once per worker process."""

    directory: Path
    estimator: Any

    @classmethod
    def from_directory(cls, model_dir: Path) -> ArchitectureModel:
        # Joblib is pickle-based. Callers must verify and allowlist the bundle
        # before reaching this method; never load an uploaded web-user model.
        estimator = joblib.load(Path(model_dir) / "model.joblib")
        if hasattr(estimator, "set_params"):
            # The serialized training pipeline carries n_jobs=6 in its final
            # ExtraTrees step. Cassette inference is tiny and must not create
            # hidden per-request pools.
            thread_parameters = {
                key: 1
                for key in estimator.get_params(deep=True)
                if key == "n_jobs" or key.endswith("__n_jobs")
            }
            if thread_parameters:
                estimator.set_params(**thread_parameters)
        return cls(Path(model_dir), estimator)

    def classify(self, cassette, genes_by_id, predictions_by_id) -> tuple[str | None, float]:
        features = cassette_features(cassette, genes_by_id, predictions_by_id)
        cas_type = cassette["classification"]["type"]
        if cas_type is None:
            return None, 0.0
        probabilities = self.estimator.predict_proba([features])[0]
        classes = self.estimator.classes_
        allowed = [
            index for index, label in enumerate(classes) if str(label).startswith(f"{cas_type}-")
        ]
        index = (
            max(allowed, key=lambda value: probabilities[value])
            if allowed
            else int(np.argmax(probabilities))
        )
        return str(classes[index]), float(probabilities[index])


@lru_cache(maxsize=4)
def _cached_architecture_model(resolved_model_dir: str) -> ArchitectureModel:
    return ArchitectureModel.from_directory(Path(resolved_model_dir))


def load_architecture_model(model_dir: Path) -> ArchitectureModel:
    return _cached_architecture_model(str(Path(model_dir).resolve()))


def classify_architecture_local(model, cassette, genes_by_id, predictions_by_id):
    loaded = (
        model
        if isinstance(model, ArchitectureModel) or callable(getattr(model, "classify", None))
        else load_architecture_model(model)
    )
    return loaded.classify(cassette, genes_by_id, predictions_by_id)
