"""Dependency-local runtime for the frozen cassette architecture model."""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np

from casandra.cassette_architecture import cassette_features
from casandra.integrity import verify_artifact_manifest


def classify_architecture_local(model_dir, cassette, genes_by_id, predictions_by_id):
    model_dir = Path(model_dir)
    verify_artifact_manifest(model_dir)
    model = joblib.load(model_dir / "model.joblib")
    features = cassette_features(cassette, genes_by_id, predictions_by_id)
    cas_type = cassette["classification"]["type"]
    if cas_type is None:
        return None, 0.0
    probabilities = model.predict_proba([features])[0]
    classes = model.classes_
    allowed = [
        index for index, label in enumerate(classes)
        if str(label).startswith(f"{cas_type}-")
    ]
    index = (
        max(allowed, key=lambda value: probabilities[value])
        if allowed
        else int(np.argmax(probabilities))
    )
    return str(classes[index]), float(probabilities[index])
