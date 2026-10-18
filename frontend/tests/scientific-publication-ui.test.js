import { render, screen, waitFor } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import { api } from "../src/api.js";
import ResultsView from "../src/components/results/ResultsView.vue";
import { sanitizeAdapterMembership } from "../src/utils/results.js";
import { cloneResultJob, resultCredential } from "./support/resultFixture.js";

afterEach(() => vi.restoreAllMocks());

describe("scientific result publication", () => {
  it("loads exact group membership from a sanitized manifest for older jobs", async () => {
    const job = cloneResultJob();
    const fullGroup = job.summary.adapter.groups[0];
    job.summary.adapter.groups = [
      {
        name: fullGroup.name,
        array_count: fullGroup.array_count,
        repeat_key: fullGroup.repeat_key,
      },
    ];
    job.artifacts = [
      {
        artifact_id: "manifest-1",
        name: "adapter/manifest.json",
        size_bytes: 900,
        media_type: "application/json",
      },
    ];
    const spy = vi
      .spyOn(api, "downloadArtifact")
      .mockResolvedValue(
        new Blob([JSON.stringify({ groups: [fullGroup] })], { type: "application/json" }),
      );
    render(ResultsView, { props: { job, credential: resultCredential } });
    await waitFor(() => expect(document.querySelectorAll(".group-member")).toHaveLength(5));
    expect(spy).toHaveBeenCalledWith(
      resultCredential.jobId,
      "manifest-1",
      resultCredential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
  });

  it("preserves bounded orientation provenance from adapter manifests", () => {
    const groups = sanitizeAdapterMembership({
      groups: [
        {
          name: "group",
          arrays: [
            {
              source_id: "source",
              array_id: "array",
              category: "Bona-fide",
              spacer_count: 4,
              strand: "-",
              input_sequence_orientation: "source",
              ccdb_strand: "+",
            },
          ],
        },
      ],
    });
    expect(groups[0].arrays[0]).toMatchObject({
      input_sequence_orientation: "source",
      ccdb_strand: "+",
    });
    const bounded = sanitizeAdapterMembership({
      groups: [
        {
          name: "group",
          arrays: [{ input_sequence_orientation: "x".repeat(100), ccdb_strand: "y".repeat(100) }],
        },
      ],
    });
    expect(bounded[0].arrays[0].input_sequence_orientation).toHaveLength(32);
    expect(bounded[0].arrays[0].ccdb_strand).toHaveLength(32);
  });

  it("bounds browser buffering and filters empty artifacts", () => {
    const job = cloneResultJob();
    job.artifacts = [
      { artifact_id: "empty", name: "empty.txt", size_bytes: 0 },
      { artifact_id: "report", name: "report.txt", size_bytes: 42 },
    ];
    render(ResultsView, {
      props: { job, credential: resultCredential, maxArchiveBytes: 128 * 1024 * 1024 },
    });
    expect(screen.getByText(/buffered in this browser tab/i)).toHaveTextContent(/128 MiB/i);
    expect(screen.queryByText("empty.txt")).not.toBeInTheDocument();
    expect(screen.getByText("report.txt")).toBeInTheDocument();
  });
});
