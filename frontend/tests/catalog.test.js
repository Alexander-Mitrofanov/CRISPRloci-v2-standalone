import { describe, expect, it, vi } from "vitest";
import { createApiClient } from "../src/api.js";
import {
  catalogCell,
  catalogCsv,
  catalogFasta,
  catalogInterval,
} from "../src/features/catalog/index.js";
import { catalogArrayFixture, catalogRows } from "./support/catalogFixture.js";

describe("catalog client and scientific formatting", () => {
  it("uses contract routes, encodes searches and sends no job capability", async () => {
    const fetcher = vi.fn().mockResolvedValue({ ok: true, json: async () => ({ items: [] }) });
    const client = createApiClient("https://analysis.example.org", fetcher);
    await client.catalogSummary();
    expect(fetcher.mock.calls[0][0]).toBe("https://analysis.example.org/api/v1/catalog");
    await client.catalogPage("arrays", {
      q: "AC+GT",
      assembly_accession: "GCA_001.1",
      limit: 25,
      cursor: "a/b=",
    });
    const [url, options] = fetcher.mock.calls[1];
    expect(new URL(url).searchParams.get("q")).toBe("AC+GT");
    expect(new URL(url).searchParams.get("cursor")).toBe("a/b=");
    expect(options.headers).toEqual({ Accept: "application/json" });
    expect(options.credentials).toBe("omit");
    await client.catalogArray(101);
    expect(fetcher.mock.calls[2][0]).toBe("https://analysis.example.org/api/v1/catalog/arrays/101");
    await expect(client.catalogPage("private-jobs")).rejects.toMatchObject({
      code: "invalid_catalog_entity",
    });
    expect(fetcher).toHaveBeenCalledTimes(3);
  });

  it("preserves missing values, true zero and deletion boundaries", () => {
    expect(catalogCell({ score: null }, "score")).toBe("Not reported");
    expect(catalogCell({ score: 0 }, "score")).toBe("0");
    expect(catalogCell({ confidence: 0.000028 }, "confidence")).toBe("0.000028");
    expect(catalogCell({ strand: "unknown" }, "strand")).toBe("Unknown");
    expect(catalogInterval({ start: null, end: 10 })).toBe("Not reported");
    expect(catalogInterval({ kind: "deletion", boundary_position: 0 })).toBe(
      "Boundary 0 (0-based)",
    );
  });

  it("exports only actual sequences, identifies consensus and escapes CSV values", () => {
    const fasta = catalogFasta(catalogRows.spacers, "spacers");
    expect(fasta).toContain("GATTACAGATTACAGATTACAGATTA");
    expect(fasta).not.toContain("302");
    expect(catalogFasta([catalogArrayFixture], "arrays")).toContain("repeat_consensus\nACGTACGT");
    expect(
      catalogCsv(
        [{ annotation: '=DANGEROUS("x")' }, { annotation: null }],
        [["annotation", "Annotation"]],
      ),
    ).toBe('"Annotation"\r\n"\'=DANGEROUS(""x"")"\r\n""');
  });
});
