import { isRecoveryHash } from "./hash.js";

export function installRecoverySectionNavigation(hasRecovery, browser = window) {
  const documentObject = browser.document;
  if (!documentObject?.addEventListener) return () => {};

  function onClick(event) {
    if (
      !hasRecovery() ||
      event.defaultPrevented ||
      event.button !== 0 ||
      event.metaKey ||
      event.ctrlKey ||
      event.shiftKey ||
      event.altKey
    )
      return;
    const element = event.target?.closest ? event.target : event.target?.parentElement;
    const anchor = element?.closest?.('a[href^="#"]');
    if (!anchor || anchor.target === "_blank") return;
    const hash = anchor.getAttribute("href");
    if (!hash || hash === "#" || isRecoveryHash(hash)) return;

    let id;
    try {
      id = decodeURIComponent(hash.slice(1));
    } catch {
      return;
    }
    const target = documentObject.getElementById(id);
    if (!target) return;
    event.preventDefault();
    const reducedMotion = browser.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    target.scrollIntoView({ behavior: reducedMotion ? "auto" : "smooth", block: "start" });
  }

  documentObject.addEventListener("click", onClick, true);
  return () => documentObject.removeEventListener("click", onClick, true);
}
