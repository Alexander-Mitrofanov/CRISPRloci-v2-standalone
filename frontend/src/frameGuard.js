export function applyFrameBootPolicy(windowObject, documentObject) {
  const root = documentObject?.getElementById("root");
  const skipLink = documentObject?.getElementById("skip-link");
  const blockedMessage = documentObject?.getElementById("frame-blocked-message");
  const blockedOpen = documentObject?.getElementById("frame-blocked-open");

  if (blockedOpen && typeof windowObject?.location?.href === "string") {
    blockedOpen.href = windowObject.location.href;
  }

  let topLevel = false;
  try {
    topLevel = Boolean(windowObject && windowObject.self && windowObject.self === windowObject.top);
  } catch {
    topLevel = false;
  }

  const allowed = Boolean(topLevel && root && skipLink && blockedMessage);
  if (documentObject?.documentElement) {
    documentObject.documentElement.dataset.frameBoot = allowed ? "allowed" : "denied";
  }

  if (root) root.hidden = !allowed;
  if (skipLink) skipLink.hidden = !allowed;
  if (blockedMessage) blockedMessage.hidden = allowed;
  return allowed;
}
