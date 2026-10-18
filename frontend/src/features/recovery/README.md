# Recovery feature

- `credential.js` validates the exact in-memory job credential shape.
- `hash.js` owns the strict `#job=<locator>.<token>` wire format and recovery URL construction.
- `browserHash.js` is the only module that reads or replaces the browser fragment.
- `sectionNavigation.js` scrolls same-page links without replacing an active recovery fragment or
  creating duplicate capability-bearing history entries.
- `index.js` is the feature façade used by components and composables.

The fragment is a bearer secret. It may never be copied into an API path/query, log, public
snapshot, cookie, local storage, or session storage. API calls continue to authenticate only with
the `Authorization` header.
