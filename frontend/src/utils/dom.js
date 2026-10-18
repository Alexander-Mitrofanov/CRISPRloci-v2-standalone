export function preferredScrollBehavior() {
  return window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth";
}

export function revealSection(id, headingSelector = "h2") {
  const region = document.getElementById(id);
  if (!region) return;
  region.scrollIntoView({
    behavior: preferredScrollBehavior(),
    block: "start",
  });
  region.querySelector(headingSelector)?.focus({ preventScroll: true });
}
