import { fireEvent, render, screen } from "@testing-library/vue";
import { flushPromises } from "@vue/test-utils";
import { afterEach, describe, expect, it, vi } from "vitest";

import App from "../src/App.vue";
import { api } from "../src/api.js";

const configured = api.configured;

afterEach(() => {
  api.configured = configured;
  vi.useRealTimers();
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

describe("analysis service recovery", () => {
  it("enables Compute after a transient startup outage without losing the draft", async () => {
    api.configured = true;
    vi.useFakeTimers();
    window.history.replaceState(null, "", "/?method=detection");
    vi.spyOn(api, "health")
      .mockRejectedValueOnce(new Error("Temporary outage"))
      .mockResolvedValue({ version: "2.1.0" });
    vi.spyOn(api, "config").mockResolvedValue({ max_total_bases: 10_000 });
    render(App);
    await flushPromises();
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    expect(screen.getByText("Service unavailable", { exact: true })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();

    await vi.advanceTimersByTimeAsync(10_000);
    await flushPromises();

    expect(screen.getByText("Service ready", { exact: true })).toBeInTheDocument();
    expect(screen.getByRole("textbox")).toHaveValue(">a\nACGT\n");
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeEnabled();
    expect(screen.getByText("1+ record · 10,000 bases max")).toBeInTheDocument();
    await vi.advanceTimersByTimeAsync(60_000);
    expect(api.health).toHaveBeenCalledTimes(2);
  });
});
