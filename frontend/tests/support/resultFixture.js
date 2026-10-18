import snapshot from "../../public/example-result.json";

export const resultCredential = Object.freeze({
  jobId: "0123456789abcdef0123456789abcdef",
  accessToken: "a".repeat(43),
  expiresAt: "2099-01-01T00:00:00Z",
});

export const cloneResultJob = () => structuredClone(snapshot.job);
