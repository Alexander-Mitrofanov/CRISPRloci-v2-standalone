import { readFileSync } from "node:fs";

import { describe, expect, it } from "vitest";

import { applyFrameBootPolicy } from "../src/frameGuard.js";

function installBootMarkup() {
  document.body.innerHTML = `
    <a id="skip-link" hidden>Skip</a>
    <main id="frame-blocked-message" hidden>
      Blocked
      <a id="frame-blocked-open" href="./" target="_blank" rel="noopener noreferrer">Open</a>
    </main>
    <div id="root" hidden></div>
  `;
}

describe("public-host frame boot policy", () => {
  it("ships every boot target hidden before JavaScript runs", () => {
    const html = readFileSync("index.html", "utf8");
    const parsed = new DOMParser().parseFromString(html, "text/html");

    expect(parsed.getElementById("root")?.hasAttribute("hidden")).toBe(true);
    expect(parsed.getElementById("skip-link")?.hasAttribute("hidden")).toBe(true);
    expect(parsed.getElementById("frame-blocked-message")?.hasAttribute("hidden")).toBe(true);
    expect(parsed.getElementById("frame-blocked-open")?.getAttribute("target")).toBe("_blank");
  });

  it("reveals the application only in a top-level browsing context", () => {
    installBootMarkup();
    const topLevelWindow = {};
    topLevelWindow.self = topLevelWindow;
    topLevelWindow.top = topLevelWindow;

    expect(applyFrameBootPolicy(topLevelWindow, document)).toBe(true);
    expect(document.documentElement).toHaveAttribute("data-frame-boot", "allowed");
    expect(document.getElementById("root")).not.toHaveAttribute("hidden");
    expect(document.getElementById("skip-link")).not.toHaveAttribute("hidden");
    expect(document.getElementById("frame-blocked-message")).toHaveAttribute("hidden");
  });

  it("fails closed in a frame and reveals only the safe notice", () => {
    installBootMarkup();
    const recoveryUrl = `https://example.test/workbench/#job=${"a".repeat(32)}.${"b".repeat(43)}`;
    const framedWindow = { self: {}, top: {}, location: { href: recoveryUrl } };

    expect(applyFrameBootPolicy(framedWindow, document)).toBe(false);
    expect(document.documentElement).toHaveAttribute("data-frame-boot", "denied");
    expect(document.getElementById("root")).toHaveAttribute("hidden");
    expect(document.getElementById("skip-link")).toHaveAttribute("hidden");
    expect(document.getElementById("frame-blocked-message")).not.toHaveAttribute("hidden");
    expect(document.getElementById("frame-blocked-open")).toHaveAttribute("href", recoveryUrl);
  });

  it("fails closed when frame ancestry cannot be inspected", () => {
    installBootMarkup();
    const inaccessibleWindow = { self: {} };
    Object.defineProperty(inaccessibleWindow, "top", {
      get() {
        throw new DOMException("Blocked", "SecurityError");
      },
    });

    expect(applyFrameBootPolicy(inaccessibleWindow, document)).toBe(false);
    expect(document.getElementById("root")).toHaveAttribute("hidden");
    expect(document.getElementById("frame-blocked-message")).not.toHaveAttribute("hidden");
  });

  it("fails closed when required boot markup is missing", () => {
    installBootMarkup();
    document.getElementById("frame-blocked-message").remove();
    const topLevelWindow = {};
    topLevelWindow.self = topLevelWindow;
    topLevelWindow.top = topLevelWindow;

    expect(applyFrameBootPolicy(topLevelWindow, document)).toBe(false);
    expect(document.documentElement).toHaveAttribute("data-frame-boot", "denied");
    expect(document.getElementById("root")).toHaveAttribute("hidden");
  });
});
