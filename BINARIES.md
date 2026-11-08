# Windows binary bundle

The Initializer downloads one archive and unpacks it into `YR_MPE/bin`.

- Archive name: `YR-MPE-bins-win64.zip` (a `.tar.gz` is also accepted).
- Accepted layouts:
  - `bin/mafft/...`, `bin/iqtree3/...`, `bin/mb/...` (wrapped in a top-level folder is fine)
  - `YR_MPE/bin/...`
- Do not include a full `YR_MPE/` source tree in the archive; the Initializer
  copies only the detected `bin/` directory.
- Default mirrors are listed in `binary_sources.json`. A custom URL can be
  pasted into the Initializer, so an institutional or temporary mirror does
  not require a code change.
- Keep `config.template.json` in sync when a bundled executable is added or
  renamed.

The Gitee mirror should expose the same asset name and the same `bin/` layout
as the GitHub release. The `{version}` placeholder in `binary_sources.json`
is replaced with `YR_MPE.__version__` at download time.
