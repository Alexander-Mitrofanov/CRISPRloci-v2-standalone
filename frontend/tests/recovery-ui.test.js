import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { defineComponent } from "vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import App from "../src/App.vue";
import { api } from "../src/api.js";
import RecoveryLink from "../src/components/jobs/RecoveryLink.vue";
import { useJobSession } from "../src/composables/useJobSession.js";
import { serializeRecoveryHash } from "../src/jobStore.js";

const credential = {
  jobId: "0123456789abcdef0123456789abcdef",
  accessToken: "a".repeat(43),
  expiresAt: "2099-01-01T00:00:00Z",
};
const nextCredential = {
  jobId: "fedcba9876543210fedcba9876543210",
  accessToken: "z".repeat(43),
  expiresAt: "2099-01-01T00:00:00Z",
};

afterEach(() => {
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

describe("private job recovery links", () => {
  it("copies the fragment link and explains its bearer capability", async () => {
    const writeText = vi.fn().mockResolvedValue(undefined);
    Object.defineProperty(navigator, "clipboard", {
      configurable: true,
      value: { writeText },
    });
    render(RecoveryLink, { props: { credential } });

    const notice = screen.getByRole("complementary", { name: /private job recovery link/i });
    expect(notice).toHaveTextContent(/reopen this job until it expires/i);
    expect(notice).toHaveTextContent(/anyone with the full link can access/i);
    await fireEvent.click(screen.getByRole("button", { name: /copy recovery link/i }));

    expect(writeText).toHaveBeenCalledOnce();
    expect(new URL(writeText.mock.calls[0][0]).hash).toBe(serializeRecoveryHash(credential));
    expect(await screen.findByText(/keep the link private/i)).toBeInTheDocument();
  });

  it("resets copied state when the active job capability changes", async () => {
    const writeText = vi.fn().mockResolvedValue(undefined);
    Object.defineProperty(navigator, "clipboard", {
      configurable: true,
      value: { writeText },
    });
    const view = render(RecoveryLink, { props: { credential } });
    await fireEvent.click(screen.getByRole("button", { name: /copy recovery link/i }));
    expect(await screen.findByRole("button", { name: /link copied/i })).toBeInTheDocument();

    await view.rerender({ credential: nextCredential });
    expect(screen.getByRole("button", { name: /copy recovery link/i })).toBeInTheDocument();
    expect(screen.queryByText(/keep the link private/i)).not.toBeInTheDocument();
  });

  it("ignores clipboard completion from the previous job capability", async () => {
    let finishCopy;
    const writeText = vi.fn().mockImplementation(
      () =>
        new Promise((resolve) => {
          finishCopy = resolve;
        }),
    );
    Object.defineProperty(navigator, "clipboard", {
      configurable: true,
      value: { writeText },
    });
    const view = render(RecoveryLink, { props: { credential } });
    await fireEvent.click(screen.getByRole("button", { name: /copy recovery link/i }));
    await view.rerender({ credential: nextCredential });
    finishCopy();
    await Promise.resolve();
    await Promise.resolve();

    expect(screen.getByRole("button", { name: /copy recovery link/i })).toBeInTheDocument();
    expect(screen.queryByText(/keep the link private/i)).not.toBeInTheDocument();
  });

  it("restores a job directly from a valid hash link", async () => {
    Object.defineProperty(Element.prototype, "scrollIntoView", {
      configurable: true,
      value: vi.fn(),
    });
    window.history.replaceState(null, "", `/${serializeRecoveryHash(credential)}`);
    vi.spyOn(api, "health").mockResolvedValue({ version: "1.0.0" });
    vi.spyOn(api, "config").mockResolvedValue({ api_version: "1.0.0" });
    vi.spyOn(api, "getJob").mockResolvedValue({ status: "queued", mode: "orientation" });

    render(App);
    const heading = await screen.findByRole("heading", { name: "Queued" });
    await waitFor(() =>
      expect(api.getJob).toHaveBeenCalledWith(
        credential.jobId,
        credential.accessToken,
        expect.objectContaining({ signal: expect.any(AbortSignal) }),
      ),
    );
    expect(heading).toBeInTheDocument();
    expect(window.location.hash).toBe(serializeRecoveryHash(credential));
    expect(screen.queryByText(/recovery file/i)).not.toBeInTheDocument();

    const historyLength = window.history.length;
    await fireEvent.click(screen.getByRole("link", { name: "References" }));
    expect(window.location.hash).toBe(serializeRecoveryHash(credential));
    expect(window.history.length).toBe(historyLength);
  });

  it("fails closed and removes a malformed recovery fragment", async () => {
    window.history.replaceState(null, "", "/#job=not-a-capability");
    vi.spyOn(api, "health").mockResolvedValue({ version: "1.0.0" });
    vi.spyOn(api, "config").mockResolvedValue({ api_version: "1.0.0" });
    const getJob = vi.spyOn(api, "getJob");

    render(App);
    expect(await screen.findByText(/recovery link is invalid/i)).toHaveAttribute("role", "alert");
    expect(window.location.hash).toBe("");
    expect(getJob).not.toHaveBeenCalled();
  });

  it("scrubs a malformed recovery fragment introduced after boot", async () => {
    vi.spyOn(api, "health").mockResolvedValue({ version: "1.0.0" });
    vi.spyOn(api, "config").mockResolvedValue({ api_version: "1.0.0" });
    const getJob = vi.spyOn(api, "getJob");
    render(App);

    window.history.pushState(null, "", "/#job=not-a-capability");
    window.dispatchEvent(new HashChangeEvent("hashchange"));

    expect(await screen.findByText(/recovery link is invalid/i)).toHaveAttribute("role", "alert");
    expect(window.location.hash).toBe("");
    expect(getJob).not.toHaveBeenCalled();
  });

  it("ignores stale cancellation results after switching recovery links", async () => {
    Object.defineProperty(Element.prototype, "scrollIntoView", {
      configurable: true,
      value: vi.fn(),
    });
    window.history.replaceState(null, "", `/${serializeRecoveryHash(credential)}`);
    vi.spyOn(api, "getJob").mockImplementation(async (jobId) => ({
      status: "complete",
      marker: jobId === credential.jobId ? "job-a" : "job-b",
    }));
    let resolveFirstCancel;
    let resolveSecondCancel;
    const cancelJob = vi.spyOn(api, "cancelJob").mockImplementation(
      (jobId) =>
        new Promise((resolve) => {
          if (jobId === credential.jobId) resolveFirstCancel = resolve;
          else resolveSecondCancel = resolve;
        }),
    );
    const Harness = defineComponent({
      setup: () => useJobSession(api),
      template: `<p data-testid="job-marker">{{ job && job.marker }}</p><p data-testid="poll-error">{{ pollError }}</p><button type="button" :disabled="cancelling" @click="cancel">Cancel</button>`,
    });
    render(Harness);
    expect(await screen.findByText("job-a")).toBeInTheDocument();

    await fireEvent.click(screen.getByRole("button", { name: "Cancel" }));
    window.history.pushState(null, "", `/${serializeRecoveryHash(nextCredential)}`);
    window.dispatchEvent(new HashChangeEvent("hashchange"));
    expect(await screen.findByText("job-b")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Cancel" })).not.toBeDisabled();

    await fireEvent.click(screen.getByRole("button", { name: "Cancel" }));
    expect(cancelJob).toHaveBeenCalledTimes(2);
    resolveSecondCancel({ job: { status: "cancelled", marker: "job-b-cancelled" } });
    expect(await screen.findByText("job-b-cancelled")).toBeInTheDocument();
    resolveFirstCancel({ job: { status: "cancelled", marker: "stale-job-a" } });
    await Promise.resolve();
    await Promise.resolve();
    expect(screen.getByTestId("job-marker")).toHaveTextContent("job-b-cancelled");
    expect(screen.getByTestId("poll-error")).toBeEmptyDOMElement();
  });
});
