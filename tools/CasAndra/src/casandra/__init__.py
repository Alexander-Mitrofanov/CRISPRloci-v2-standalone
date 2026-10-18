"""Cas gene detection and CRISPR-Cas cassette classification."""

__version__ = "0.3.0.dev2"


def __getattr__(name: str):
    if name in {"Analyzer", "AnalysisOptions", "GenomeInput", "RunResult"}:
        from casandra.api import AnalysisOptions, Analyzer, GenomeInput, RunResult

        return {
            "Analyzer": Analyzer,
            "AnalysisOptions": AnalysisOptions,
            "GenomeInput": GenomeInput,
            "RunResult": RunResult,
        }[name]
    raise AttributeError(name)


__all__ = ("AnalysisOptions", "Analyzer", "GenomeInput", "RunResult", "__version__")
