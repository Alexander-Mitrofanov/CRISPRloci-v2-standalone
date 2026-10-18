import { describe, expect, it } from "vitest";
import { render, screen, cleanup } from "@testing-library/vue";
import DefenseContext from "../src/components/results/annotations/DefenseContext.vue";

describe("Defense context evidence", () => {
  it("retains older-model and activity limitations when there are no calls", () => {
    render(DefenseContext, {
      props: { available: true, counts: { systems: 0, genes: 42, components: 0, hits: 5 } },
    });
    expect(screen.getByText(/model release 1.2.3/)).toBeInTheDocument();
    expect(screen.getByText(/do not establish activity/)).toBeInTheDocument();
    expect(
      screen.getByText(/No modeled systems detected in the assessed sequence/),
    ).toBeInTheDocument();
    cleanup();
  });
  it("distinguishes truncated previews and unknown counts from absence", () => {
    render(DefenseContext, {
      props: { available: true, systemsTruncated: true, hitsTruncated: true },
    });
    expect(screen.getByText(/No system calls in this preview/)).toBeInTheDocument();
    expect(screen.getByText(/Unknown system calls/)).toBeInTheDocument();
    expect(screen.getByText(/Profile-hit preview is limited/)).toBeInTheDocument();
    cleanup();
  });
});
