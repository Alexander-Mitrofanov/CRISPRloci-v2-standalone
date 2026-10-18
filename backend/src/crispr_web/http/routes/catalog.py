"""Read-only public reference catalogue, independent of private job access."""

from __future__ import annotations

from typing import Annotated, Callable, Literal

from fastapi import APIRouter, HTTPException, Path, Query

from ...catalog import (
    ENTITIES,
    Catalog,
    CatalogCursorError,
    CatalogQueryError,
    CatalogUnavailable,
)
from ...catalog_models import CatalogArray, CatalogAssembly, CatalogPage, CatalogSummary
from ...config import Settings


def _read(call: Callable):
    try:
        return call()
    except CatalogUnavailable as error:
        raise HTTPException(
            status_code=503,
            detail="reference catalogue is temporarily unavailable",
            headers={"Retry-After": "30"},
        ) from error
    except CatalogCursorError as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    except CatalogQueryError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


def build_catalog_router(settings: Settings) -> APIRouter:
    router = APIRouter(prefix="/api/v1/catalog", tags=["catalogue"])
    catalog = Catalog(settings.catalog_database_path, settings.catalog_search_index_path)

    @router.get("", response_model=CatalogSummary, operation_id="catalog_summary")
    def summary():
        return _read(catalog.summary)

    def endpoint(name: str):
        def page(
            q: Annotated[str | None, Query(max_length=1000)] = None,
            assembly_accession: Annotated[str | None, Query(
                min_length=1, max_length=80, pattern=r"^GC[AF]_\d+\.\d+$"
            )] = None,
            array_id: Annotated[int | None, Query(ge=1, lt=2**63 - 1)] = None,
            cursor: Annotated[str | None, Query(max_length=2000)] = None,
            limit: Annotated[int, Query(ge=1, le=100)] = 25,
            match: Literal["exact", "similar"] = "exact",
        ):
            return _read(lambda: catalog.page(
                name, q=q, assembly_accession=assembly_accession, array_id=array_id,
                cursor=cursor, limit=limit, match=match,
            ))
        return page

    for name, entity in ENTITIES.items():
        router.add_api_route(
            "/" + name,
            endpoint(name),
            methods=["GET"],
            response_model=CatalogPage[entity.model],
            operation_id="catalog_" + name.replace("-", "_"),
            description=(
                "Bounded keyset pagination. Assembly search is an accession prefix; "
                "sequence, call ID, Cas annotation/profile and subtype searches are exact. "
                "Global exact scientific search uses an offline snapshot-pinned index. "
                "Cas names also match complete gene-name tokens in source annotations/profiles. "
                "Repeat/spacer counts describe occurrences, including explicit spacer deletions."
                " For repeats/spacers, match=similar finds full-length sequences with at most one "
                "literal nucleotide substitution on either strand. IUPAC symbols are not wildcards; "
                "insertions, deletions and local alignments are not searched."
            ),
        )

    @router.get("/assemblies/{accession}", response_model=CatalogAssembly, operation_id="catalog_assembly")
    def assembly(accession: Annotated[str, Path(max_length=80, pattern=r"^GC[AF]_\d+\.\d+$")]):
        result = _read(lambda: catalog.detail("assemblies", accession))
        if result is None:
            raise HTTPException(status_code=404, detail="catalogue assembly not found")
        return result

    @router.get("/arrays/{array_id}", response_model=CatalogArray, operation_id="catalog_array")
    def array(array_id: Annotated[int, Path(ge=1, lt=2**63 - 1)]):
        result = _read(lambda: catalog.detail("arrays", array_id))
        if result is None:
            raise HTTPException(status_code=404, detail="catalogue array not found")
        return result

    return router
