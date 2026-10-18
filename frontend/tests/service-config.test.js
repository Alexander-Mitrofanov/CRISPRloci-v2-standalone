import { defineComponent } from "vue";
import { render } from "@testing-library/vue";
import { flushPromises } from "@vue/test-utils";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import { useServiceConfig } from "../src/composables/useServiceConfig.js";

function mountService(client) {
  let state;
  const view = render(
    defineComponent({
      setup() {
        state = useServiceConfig(client);
        return () => null;
      },
    }),
  );
  return { ...state, unmount: view.unmount };
}

function deferred() {
  let resolve;
  let reject;
  const promise = new Promise((yes, no) => {
    resolve = yes;
    reject = no;
  });
  return { promise, resolve, reject };
}

function client() {
  return {
    configured: true,
    health: vi.fn().mockResolvedValue({ version: "2.1.0" }),
    config: vi.fn().mockResolvedValue({ max_total_bases: 1000 }),
  };
}

beforeEach(() => vi.useFakeTimers());
afterEach(() => vi.useRealTimers());

describe("service readiness retries", () => {
  it("retains explicitly advertised scientific modes", async () => {
    const api = client();
    api.config.mockResolvedValue({ modes: ["detection", "protospacer", null] });
    const state = mountService(api);
    await flushPromises();
    expect(state.service.value.modes).toEqual(["detection", "protospacer"]);
  });

  it("retries continued outages at a bounded rate and stops after unmount", async () => {
    const api = client();
    api.health.mockRejectedValue(new Error("Unavailable"));
    const state = mountService(api);
    await flushPromises();
    await vi.advanceTimersByTimeAsync(9999);
    expect(api.health).toHaveBeenCalledTimes(1);
    await vi.advanceTimersByTimeAsync(20_001);
    expect(api.health).toHaveBeenCalledTimes(4);
    state.unmount();
    await vi.advanceTimersByTimeAsync(60_000);
    window.dispatchEvent(new Event("focus"));
    window.dispatchEvent(new Event("online"));
    expect(api.health).toHaveBeenCalledTimes(4);
  });

  it.each(["focus", "online"])("retries on %s without overlapping requests", async (event) => {
    const api = client();
    const pending = deferred();
    api.health.mockRejectedValueOnce(new Error("Unavailable")).mockReturnValueOnce(pending.promise);
    const state = mountService(api);
    await flushPromises();
    window.dispatchEvent(new Event(event));
    window.dispatchEvent(new Event(event));
    await vi.advanceTimersByTimeAsync(60_000);
    expect(api.health).toHaveBeenCalledTimes(2);
    pending.resolve({ version: "recovered" });
    await flushPromises();
    expect(state.service.value).toMatchObject({ state: "online", version: "recovered" });
    window.dispatchEvent(new Event(event));
    expect(api.health).toHaveBeenCalledTimes(2);
  });

  it.each(["resolve", "reject"])(
    "ignores a superseded request that later calls %s",
    async (settle) => {
      const api = client();
      const pending = deferred();
      api.health.mockReturnValueOnce(pending.promise);
      const state = mountService(api);
      const firstSignal = api.health.mock.calls[0][0].signal;
      await state.refresh();
      expect(firstSignal.aborted).toBe(true);
      pending[settle](
        settle === "resolve" ? { version: "obsolete" } : new Error("Obsolete outage"),
      );
      await flushPromises();
      await vi.advanceTimersByTimeAsync(60_000);
      expect(state.service.value).toMatchObject({ state: "online", version: "2.1.0" });
      expect(api.health).toHaveBeenCalledTimes(2);
    },
  );

  it("aborts the remaining request after config fails and retries both reads", async () => {
    const api = client();
    const pending = deferred();
    api.health.mockReturnValueOnce(pending.promise);
    api.config.mockRejectedValueOnce(new Error("Config unavailable"));
    const state = mountService(api);
    await flushPromises();
    expect(api.health.mock.calls[0][0].signal.aborted).toBe(true);
    await vi.advanceTimersByTimeAsync(10_000);
    expect(state.service.value.state).toBe("online");
    expect(state.limits.value.maxBases).toBe(1000);
    expect(api.config).toHaveBeenCalledTimes(2);
    pending.resolve({ version: "obsolete" });
    await flushPromises();
  });

  it("aborts an in-flight read on unmount and ignores its late result", async () => {
    const api = client();
    const pending = deferred();
    api.health.mockReturnValueOnce(pending.promise);
    const state = mountService(api);
    state.unmount();
    expect(api.health.mock.calls[0][0].signal.aborted).toBe(true);
    pending.resolve({ version: "obsolete" });
    await flushPromises();
    expect(state.service.value.state).toBe("checking");
    expect(state.limits.value.maxBases).toBe(0);
    await vi.advanceTimersByTimeAsync(60_000);
    expect(api.health).toHaveBeenCalledTimes(1);
  });

  it("makes no requests when the API endpoint is not configured", async () => {
    const api = { ...client(), configured: false };
    const state = mountService(api);
    await vi.advanceTimersByTimeAsync(60_000);
    window.dispatchEvent(new Event("focus"));
    window.dispatchEvent(new Event("online"));
    expect(state.service.value.state).toBe("offline");
    expect(api.health).not.toHaveBeenCalled();
    expect(api.config).not.toHaveBeenCalled();
  });
});
