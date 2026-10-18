import { nextTick, onMounted, onScopeDispose, ref, watch } from "vue";
import { saveBlob } from "../../utils/download.js";
import { createGenomeBrowser } from "./renderer.js";

export function useGenomeBrowser(props, host, onSelect, create = createGenomeBrowser) {
  const ready = ref(false);
  const error = ref("");
  let browser;
  let disposed = false;

  onMounted(async () => {
    try {
      const value = await create(host.value, props.sources, props.tracks, onSelect);
      if (disposed) {
        value.dispose();
        return;
      }
      browser = value;
      ready.value = true;
      if (props.active) browser.visible();
    } catch (_failure) {
      if (!disposed)
        error.value =
          "The genome viewer could not start. Your result tables and downloads remain available.";
    }
  });
  watch(
    () => props.active,
    async (active) => {
      await nextTick();
      if (active && !disposed) browser?.visible();
    },
  );
  onScopeDispose(() => {
    disposed = true;
    browser?.dispose();
    browser = undefined;
  });

  async function navigate(source, start = 1, end = source.length) {
    if (!ready.value) return;
    if (
      !Number.isSafeInteger(start) ||
      !Number.isSafeInteger(end) ||
      start < 1 ||
      end < start ||
      end > source.length
    ) {
      error.value = `Choose an interval between 1 and ${source.length.toLocaleString()} on this record.`;
      return;
    }
    error.value = "";
    try {
      await browser.search(`${source.chr}:${start}-${end}`);
    } catch (_failure) {
      error.value = "This interval could not be displayed.";
    }
  }

  async function zoom(direction) {
    if (!browser) return;
    try {
      await (direction === "in" ? browser.zoomIn() : browser.zoomOut());
    } catch (_failure) {
      error.value = "The view could not be resized.";
    }
  }

  function exportFigure(source) {
    if (!browser) return;
    try {
      const document = new DOMParser().parseFromString(browser.svg(), "image/svg+xml");
      const svg = document.documentElement;
      const ns = "http://www.w3.org/2000/svg";
      const width = Number(svg.getAttribute("width"));
      const height = Number(svg.getAttribute("height"));
      if (!Number.isFinite(width) || !Number.isFinite(height))
        throw new Error("Invalid figure extent");
      const plot = document.createElementNS(ns, "g");
      while (svg.firstChild) plot.appendChild(svg.firstChild);
      plot.setAttribute("transform", "translate(0 56)");
      svg.appendChild(plot);
      svg.setAttribute("height", String(height + 56));
      svg.setAttribute("viewBox", `0 0 ${width} ${height + 56}`);
      for (const [y, text] of [
        [20, `CRISPRloci · ${source.inputId.slice(0, 100)}`],
        [40, `${browser.loci()} · 1-based inclusive · reported evidence`],
      ]) {
        const label = document.createElementNS(ns, "text");
        label.setAttribute("x", "12");
        label.setAttribute("y", String(y));
        label.setAttribute("font-family", "sans-serif");
        label.setAttribute("font-size", "13");
        label.setAttribute("fill", "#292b2f");
        label.textContent = text;
        svg.appendChild(label);
      }
      const title = document.createElementNS(ns, "title");
      title.textContent = `CRISPRloci genomic context: ${source.inputId}`;
      const description = document.createElementNS(ns, "desc");
      description.textContent =
        `${source.chr} = ${source.inputId}; ${browser.loci()}. ` +
        "Coordinates are 1-based inclusive. Reported evidence only; proximity does not establish function. " +
        props.tracks.map((track) => track.label).join("; ");
      svg.prepend(title, description);
      saveBlob(
        new Blob([new XMLSerializer().serializeToString(document)], { type: "image/svg+xml" }),
        "crisprloci-genome.svg",
      );
    } catch (_failure) {
      error.value = "The figure could not be exported.";
    }
  }
  return { ready, error, navigate, zoom, exportFigure };
}
