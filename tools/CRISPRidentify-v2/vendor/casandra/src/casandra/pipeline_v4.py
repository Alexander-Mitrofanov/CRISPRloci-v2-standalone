"""Final genome pipeline using linear-time ORF extraction and fast HMM scan."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from casandra.cassette import construct_cassettes
from casandra.cassette_filter import passes_genome_evidence_gate
from casandra.cassette_hybrid_v2 import classify_hybrid_cassette_v2
from casandra.integrity import validate_model_bundle
from casandra.orf_fast import extract_genome_proteins_fast
from casandra.pipeline import _jsonl, _write_jsonl, _write_outputs
from casandra.production_v4 import predict_proteins_v4
from casandra.supplied import extract_genbank_proteins, extract_gff3_proteins
from casandra.utils import sha256_file


def predict_genome_v4(model_bundle: Path, output_dir: Path, *, genome_fasta=None, genbank=None, gff3=None, gene_mode="auto", translation_table=11, threads=6):
    if output_dir.exists(): raise FileExistsError(f"Refusing to replace run output: {output_dir}")
    if genbank is not None and (genome_fasta is not None or gff3 is not None): raise ValueError("GenBank input is mutually exclusive")
    if genbank is None and genome_fasta is None: raise ValueError("A genome FASTA or GenBank input is required")
    validate_model_bundle(model_bundle)
    config = json.loads((model_bundle / "config.json").read_text(encoding="utf-8"))
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)); started = time.perf_counter()
    try:
        extraction = staging / "extraction"
        if genbank is not None: extract_genbank_proteins(genbank, extraction); source=genbank; input_mode="genbank_supplied_cds"
        elif gff3 is not None: extract_gff3_proteins(genome_fasta, gff3, extraction, translation_table); source=genome_fasta; input_mode="gff3_supplied_cds"
        else: extract_genome_proteins_fast(genome_fasta, extraction, mode=gene_mode, translation_table=translation_table); source=genome_fasta; input_mode="pyrodigal_de_novo"
        predict_proteins_v4(model_bundle / config["protein_model_directory"], extraction / "proteins.faa", staging / "protein_predictions.jsonl", threads=threads)
        genes=list(_jsonl(extraction / "genes.jsonl")); predictions=list(_jsonl(staging / "protein_predictions.jsonl"))
        candidates=construct_cassettes(genes,predictions,**config["cassette_parameters"])
        genes_by_id={x["protein_id"]:x for x in genes}; predictions_by_id={x["sequence_id"]:x for x in predictions}; accepted=[]; rejected=[]
        for cassette in candidates:
            cassette["final_classification"]=classify_hybrid_cassette_v2(cassette,model_bundle/config["cassette_architecture_directory"],genes_by_id,predictions_by_id)
            keep,reason=passes_genome_evidence_gate(cassette,genes_by_id,predictions_by_id); cassette["genome_evidence_gate"]={"accepted":keep,"rule":reason}; (accepted if keep else rejected).append(cassette)
        _write_outputs(staging,genes,predictions,accepted); _write_jsonl(staging/"rejected_cassette_candidates.jsonl",rejected)
        run={"schema_version":4,"program":"CasAndra","bundle_id":config["bundle_id"],"bundle_role":config["bundle_role"],"created_at":datetime.now(timezone.utc).isoformat(),"input":{"path":str(source),"sha256":sha256_file(source),"mode":input_mode},"cpu_threads":threads,"wall_seconds":time.perf_counter()-started,"genes":len(genes),"cas_proteins":sum(x["is_cas"] for x in predictions),"raw_cassette_candidates":len(candidates),"rejected_cassette_candidates":len(rejected),"cas_cassettes":len(accepted),"classified_cassettes":len(accepted),"crispr_array_prediction":False,"offline_inference":True,"genome_evidence_gate":"supported subtype (or Type VII); Class 1 >=2 Cas genes; Class 2 >=2 genes or singleton length>=400 aa and positive profile score>=50"}
        run["hierarchy_conflict_masked_cassettes"] = sum(
            cassette["final_classification"]["evidence"][
                "hierarchy_conflict_masked"] for cassette in accepted
        )
        (staging/"run.json").write_text(json.dumps(run,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        manifest={"schema_version":4,"bundle_manifest_sha256":sha256_file(model_bundle/"manifest.json"),"files":{}}
        for path in sorted(staging.rglob("*")):
            if path.is_file(): manifest["files"][str(path.relative_to(staging))]={"size":path.stat().st_size,"sha256":sha256_file(path)}
        (staging/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8"); staging.replace(output_dir); return run
    except Exception:
        shutil.rmtree(staging,ignore_errors=True); raise

