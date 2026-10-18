# Frontend agent guide

The frontend is a static Vue 3/Vite application. Preserve the existing
GitHub Pages build and the stable public API origin boundary.

## Ownership

- `src/components/submission/`: form presentation and narrow input controls.
- `src/components/jobs/`: recovery-link and progress UI.
- `src/components/results/`: scientific result views.
- `src/composables/`: stateful orchestration; keep pure transformations out.
- `src/features/`: cohesive domain implementations behind stable top-level
  façades such as `src/example.js` and feature-local `index.js` files.
- `src/styles/`: ordered CSS entry point plus narrow domain files.
- `src/utils/`: framework-neutral helpers shared by more than one feature.

## Rules

- Prefer one focused SFC, composable, validator, or formatter per file.
- Keep coupled template/logic/style concepts together by feature; do not create
  generic dumping-ground modules.
- Components import a feature façade, not another feature's internal file.
- Normalize backend aliases once before presentation. New result components
  should accept narrow props rather than a generic raw `Object`.
- Routes and shared enums live in
  `src/contracts/public-api-v1.json`; backend OpenAPI tests enforce drift.
- Keep direct HTTP calls in API modules or orchestration composables. The source
  boundary check rejects transport imports in presentation components.
- Preserve header-only API authentication, strict `#job=` recovery parsing,
  exact example hashing, CSP/frame boot policy, and `null`/missing scientific semantics.
- `src/styles/index.css` import order is part of the cascade. Follow
  `src/styles/README.md`.
- Production builds require an exact HTTPS `VITE_API_BASE_URL`; no secrets may
  be introduced through `VITE_*`.

## Checks

Run these from the directory containing this guide; that is also the root of
the exported public repository.

- Feature/unit test: `npm test -- --run <test-file>`
- Source structure: `npm run check:source`
- Public export scanner: `npm run scan:public`
- Tests and production build: `npm test && npm run build`

Inside the private monorepo, `../scripts/check.sh frontend` is the equivalent
repository wrapper. It is intentionally absent from the exported public repo.

Only `frontend/` is exported publicly. Keep its nested
`.github/workflows/pages.yml` valid for a repository where this directory is
the root.
