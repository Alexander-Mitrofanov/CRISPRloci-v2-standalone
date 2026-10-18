// Catalog transport code contains the public field name `accession`. Keep actual
// genomic identities out of the bundle; example/static data remains stricter.
const accessionValue =
  /(?<![A-Za-z0-9])(?:(?:CP|FR|LN|LR|AP)\d{6}|GC[AF]_\d{9}\.\d+|(?:NC|NZ|NW|NT)_\d{6,}\.\d+)/i;
const staticIdentityKey =
  /(?:organism|strain|accession|ncbi_url|region_start_1based|region_end_1based)\s*[=:"]/i;
const bundleIdentityKey =
  /(?:organism|strain|ncbi_url|region_start_1based|region_end_1based)\s*[=:"]/i;

export function hasForbiddenPublicIdentity(text, { compiledBundle = false } = {}) {
  return (
    accessionValue.test(text) || (compiledBundle ? bundleIdentityKey : staticIdentityKey).test(text)
  );
}
