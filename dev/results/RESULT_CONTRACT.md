# Plugin result contract (v1)

> Pilot: One-Step OrthoFinder. The same channel should replace the per-plugin
> `import_dataset_signal` / `import_alignment_signal` handlers over time.

## What is wrong today

1. **One signal per plugin.** `import_dataset_signal`, `import_alignment_signal`,
   `batch_import_alignment_signal`, `exportToDataset` — each with its own payload
   shape and its own host handler (`add_seqdbg_export_to_workspace`,
   `handle_seqmatrix_import`, ...).
2. **Ad-hoc transport classes.** `DatasetObject` is redefined inline in several
   host methods; fields differ; there is no schema and no validation.
3. **Two item models.** `TempDatasetItem` (NEXUS import) and the local
   `DatasetItem` in `dataset_manager.py` are separate from
   `dataset_models.DatasetItem`. Conversion code is duplicated and lossy:
   `dataset_manager` hardcodes `ITEM_TYPE_ALIGNMENT` even for sequence items.
4. **Items are dropped on import.** `SingleGeneWorkspace.add_dataset()` creates
   a `DatasetInfo` but explicitly skips its items (`pass`), so the dataset has
   no children in the selection manager and cannot be used downstream.
5. **Persistence is duplicated.** Selection state lives in
   `dataset_selection_state.json`; display data lives in
   `DatasetInfo.settings['dataset_items']`; sequences are not serialized, so
   file paths are the only bridge.
6. **No provenance.** A dataset does not record which plugin/tool/run produced
   it, nor its input/output directories.
7. **Inconsistent controls.** `add_dataset()` uses `DatasetButton`;
   `refresh_workspace_layout()` rebuilds plain `QToolButton`s; tooltips carry
   only a name; there is no summary/context menu.

## v1 design

### 1. Transport DTOs (stdlib only)

`YR_MPE/runtime/results.py`:

```
ResultSpec
  schema = "yrmpe.result.v1"
  source: str
  datasets: [DatasetSpec]
  artifacts: [path]
  summary: dict
  warnings: [str]

DatasetSpec
  name, description, source
  items: [DatasetItemSpec]
  partition_scheme: str | None
  summary: dict
  metadata: dict
  dataset_id: str | None

DatasetItemSpec
  loci_name
  item_type (sequence/alignment/model/...)
  sequences: [SeqRecord-like]
  file_path, format
  is_aligned, length, sequence_count
  data: dict          # orthogroup, map/protein paths, ...
  metadata: dict
  parent_index: int | None
  selected: bool
```

`DatasetItemSpec.validate()` and `DatasetSpec.validate()` reject incomplete
items early. The DTO does not import PyQt or the host's dataset manager.

### 2. One signal on `BasePlugin`

```python
result_ready = pyqtSignal(object)

def export_result(self, result: ResultSpec) -> bool:
    problems = result.validate()
    if problems:
        ... log and return False
    self.result_ready.emit(result)
    return True
```

Plugins no longer define their own export signal. The host connects once:

```python
window.result_ready.connect(self.receive_plugin_result)
```

### 3. One host receiver

`YR_MPEA_Widget.receive_plugin_result(result)` validates the spec, calls
`SingleGeneWorkspace.add_result_spec(result)`, and reports errors. It does not
know OrthoFinder-specific fields.

`SingleGeneWorkspace.add_result_spec()`:

1. Create/reuse a `DatasetInfo` with a unique name.
2. Convert each `DatasetItemSpec` to `dataset_models.DatasetItem`; set
   `item_type`, sequences, counts, `data`, `metadata`, parent links, ownership
   UUID, and selection state.
3. Register every item with `DatasetSelectionManager.add_item()` so the
   selection tree and downstream tools see them.
4. Persist a canonical `settings['dataset_items']` (with `id`, `item_type`,
   `parent_id`, `data`, `metadata`) plus `settings['source']`,
   `settings['summary']`, `settings['provenance']`.
5. Add a `DatasetButton` control with a summary tooltip and a context menu
   ("Open Dataset Manager", "Show summary").

### 4. Legacy adapter

`add_dataset(dataset)` keeps working for SeqDBG/SeqMatrix. It detects:

- `ResultSpec` / `DatasetSpec` → new importer;
- legacy `DatasetObject` / list of items → convert to `DatasetSpec` and use the
  same importer.

This fixes the existing plugins without touching them first.

### 5. OrthoFinder pilot

The Output tab gets **Export to workspace**, enabled when the plugin was opened
from YR-MPEA and a run completed. It builds:

- one `DatasetItemSpec` per single-copy orthogroup (`OGxxxx.cds.fna`),
- `item_type="sequence"`, `is_aligned=False`,
- `data` with the `.faa`, `.map.tsv`, and OG id,
- a `DatasetSpec` named `OrthoFinder_<output-dir>`,
- provenance: input folder/file, output dir, genetic code, counts.

The button emits one `ResultSpec`; the host creates the dataset and button.

## Deliberately not in v1

- Full replacement of `DatasetManager`'s local item model. It still loads from
  `settings['dataset_items']`; v1 makes that data lossless and type-correct.
- A global result bus / message-bus integration. The Qt signal is the GUI
  transport; `ResultSpec` is serializable, so it can move to the bus later.
- Multi-result routing (model/tree/distance). The DTO already has `item_type`,
  but the host importer currently builds datasets.
- Dataset manifest files. Provenance is stored in `DatasetInfo.settings`; a
  future `dataset_manifest.json` can make output directories self-describing.

## Acceptance

- OrthoFinder → YR-MPEA creates one dataset with N single-copy OG items.
- The dataset button shows name, loci × taxa, source.
- Double-click opens Dataset Manager with N rows; item types are preserved.
- `DatasetSelectionManager` contains the items; selecting the dataset selects
  its children.
- SeqDBG/SeqMatrix legacy exports still work and now register their items.
- Headless tests cover DTO validation, workspace import, and OrthoFinder spec
  construction without running OrthoFinder.

## Implementation (v1)

Files:

- `YR_MPE/runtime/results.py` — DTOs and validation.
- `YR_MPE/templates/base_plugin_ui.py` — `result_ready` signal and
  `export_result()`.
- `YR_MPE/platforms/mpea.py` — `receive_plugin_result()`,
  `SingleGeneWorkspace.add_result_spec()` / `add_dataset_spec()`, legacy
  adapter, `DatasetButton` tooltip/context menu, dataset summary label.
- `YR_MPE/plugins/orthofinder_plugin.py` — `Export to workspace` button and
  `_build_result_spec()`.
- `YR_MPE/platforms/methods/dataset_selection_manager.py` — restores
  sequences from `file_path` after loading a workspace.
- `YR_MPE/platforms/methods/dataset_manager.py` — preserves `item_type`,
  `format`, `data`, `metadata` when loading/saving; no longer forces
  `ITEM_TYPE_ALIGNMENT`.

Tests and renderers:

- `dev/results/test_results.py` — 7 tests (DTO validation, workspace import,
  legacy adapter, unique names, OrthoFinder spec from fake output).
- `dev/results/render_workspace.py` — workspace control + Dataset Manager
  screenshots.
- Screenshots: `dev/results/screenshots/01_workspace_dataset.png`,
  `02_dataset_manager.png`; plugin output: `dev/gbff_to_orthofinder/screenshots/04_output.png`.

Verified with the real 62-OG output: dataset `OrthoFinder_plugin_full`,
62 sequence items × 4 taxa, Dataset Manager shows 62 rows.

## Next refactor (v2, not done)

- Make `DatasetSelectionManager.items` the single source of truth and remove
  the local `DatasetItem` model in `dataset_manager.py` (it still duplicates
  fields and re-implements loading).
- Add `dataset_manifest.json` next to plugin outputs so a dataset can be
  re-imported without the workspace state.
- Route model/tree/distance results through the same `ResultSpec` channel;
  v1 only imports datasets.
- Bulk-add items in `DatasetSelectionManager` to avoid rebuilding the
  selection tree per item; relevant when a dataset has thousands of loci.
- Integrate with the message bus (`kind=result`) once the bus lands.
