const MAX_EXAMPLE_BYTES = 5_000_000;

export async function readExampleResponse(response) {
  if (!response.ok || Number(response.headers?.get("content-length")) > MAX_EXAMPLE_BYTES)
    throw new Error("The example demonstration could not be loaded.");
  if (!response.body?.getReader) {
    const text = await response.text();
    if (new TextEncoder().encode(text).byteLength > MAX_EXAMPLE_BYTES)
      throw new Error("The example exceeds the browser read limit.");
    return text;
  }
  const reader = response.body.getReader();
  const chunks = [];
  let bytes = 0;
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      bytes += value.byteLength;
      if (bytes > MAX_EXAMPLE_BYTES) {
        await reader.cancel();
        throw new Error("The example exceeds the browser read limit.");
      }
      chunks.push(value);
    }
    return await new Blob(chunks).text();
  } finally {
    reader.releaseLock();
  }
}
