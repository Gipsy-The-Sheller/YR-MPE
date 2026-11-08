# YR-MPE plugin development standard (v1)

**Status:** draft for the old-plugin refactor, 2026-11-03.
**Applies to:** every module under `YR_MPE/plugins/`, and every external
entry point that reuses `BasePlugin`.
**Normative words:** **MUST**, **SHOULD**, **MAY** are used in the usual
RFC 2119 sense. A deviation from a MUST needs a written reason in the plugin
module docstring.

This document is the single source of truth for new plugins. The wiki page
`YR-MPE.wiki/development/plugin_development.md` is historical and does not
override this file.

## 0. Two plugin kinds

### 0.1 Analysis / data plugin

An analysis plugin runs a method (external binary or built-in Python) and
produces one or more workspace objects: alignment, model, distance matrix,
tree, tree set, chain, selection table, and so on.

**MUST** follow the full contract:

- `BasePlugin` subclass + `BaseProcessThread` subclass when a process is
  involved;
- four inherited tabs (`Input & Parameters`, `Output Preview`, `Console`,
  `Citation`);
- `ResultSpec` export through `BasePlugin.export_result()` when the result
  belongs in the workspace;
- one `PluginEntry.run(...)` returning the plugin widget.

Reference implementations:

- built-in analysis: `YR_MPE/plugins/codon_ztest_plugin.py`,
  `YR_MPE/plugins/rf_plugin.py`;
- dataset producer: `YR_MPE/plugins/orthofinder_plugin.py`
  (`_build_result_spec()` + `export_to_workspace()`);
- external binary wrapper: `YR_MPE/plugins/mafft_plugin.py`,
  `YR_MPE/plugins/iqtree_plugin.py` (legacy UI, being migrated).

### 0.2 Viewer / utility plugin

A viewer consumes an existing object and does not produce a dataset (IcyTree,
MiniTracer, a preview tool). It **MAY** use a specialised widget instead of
the four-tab layout, but **MUST** still:

- use the `PluginEntry.run(import_from, import_data, workdir)` shape;
- declare `plugin_name`, `tool_name`, `input_types`, `output_types`, and
  citation keys;
- never import `mpea.py` or mutate the workspace directly;
- never create dataset objects; viewers only request/receive data.

If a viewer starts to produce a reusable result, it becomes an analysis
plugin and must follow section 7.

## 1. Module, class, and naming contract

### 1.1 Required artefacts

Every analysis plugin module provides exactly one primary tool wrapper:

```text
YR_MPE/plugins/<tool>_plugin.py
    <Tool>Thread(BaseProcessThread)
    <Tool>Plugin(BasePlugin)
    <Tool>PluginEntry
```

- `<Tool>Thread` is **MUST** when the plugin launches an external process.
- `<Tool>Plugin` **MUST** implement `init_plugin_info()`,
  `setup_input_tab()`, `setup_output_tab()`, and `run_analysis()`.
  `<Tool>Plugin` **MUST** implement `stop_analysis()` when it starts a
  cancellable worker; a fast in-process analysis may hide the inherited Stop
  button instead of pretending to support cancellation.
- `<Tool>PluginEntry.run()` **MUST** accept
  `import_from=None, import_data=None, workdir=None, **kwargs` and return the
  plugin widget.

### 1.2 Naming

| Thing | Rule | Example |
|---|---|---|
| module | lowercase snake_case, `_plugin.py` | `mafft_plugin.py` |
| thread class | `CamelCaseThread` | `MAFFTThread` |
| plugin class | `CamelCasePlugin` | `MAFFTPlugin` |
| entry class | `CamelCasePluginEntry` | `MAFFTPluginEntry` |
| factory method | `get_<id>_plugin()` | `get_mafft_plugin()` |
| result `source` | stable lowercase snake_case id | `"mafft"`, `"one_step_orthofinder"` |
| semantic action | `<plugin_id>.<verb>` | `orthofinder.modify_options` |
| icon | `YR_MPE/icons/software/<tool>.svg` | `software/iqtree.svg` |
| run name | human-readable, no version | `IQ-TREE run 1` |

Do not rename an existing `source` value or action id without a migration;
agents and persisted `workspace.json` records depend on them.

### 1.3 Import boundaries

Plugins **MUST NOT** import:

- `YR_MPE.platforms.mpea` or any host window;
- `DatasetSelectionManager`, `DatasetItem`, `DatasetInfo`;
- `PluginAdapter` / `PluginExecutor` (legacy adapter layer).

A plugin may import:

- `YR_MPE.templates.base_plugin_ui.BasePlugin`;
- `YR_MPE.templates.base_process_thread.BaseProcessThread`;
- `YR_MPE.runtime.results` (`ResultSpec`, `DatasetSpec`,
  `DatasetItemSpec`);
- `YR_MPE.runtime.citations`;
- `YR_MPE.runtime.subprocess_utils` (`run_hidden`, `popen_hidden`);
- its own helper modules under `YR_MPE/plugins/components/`.

The host owns workspace mutation. The plugin only emits a validated
`ResultSpec`; `YR_MPEA_Widget.receive_plugin_result()` converts it.

### 1.4 Third-party dependencies

Use what the runtime already ships: PyQt5, Biopython, matplotlib, numpy, and
the packaged binaries. A new runtime dependency **MUST** be added to the
runtime repository, the initializer, and the relevant validation workflow
before the plugin is merged.

## 2. Plugin metadata, configuration, citations

### 2.1 `init_plugin_info()`

Every plugin **MUST** set:

```python
self.plugin_name = "MAFFT 7"                 # UI title
self.tool_name = "MAFFT"                     # must match config.json exactly
self.citation_keys = ["katoh2013mafft"]      # keys in citations/*.bib|ris
self.input_types = {"Alignment FASTA": ["fas", "fna", "fa", "fasta"]}
self.output_types = {"Alignment FASTA": ["fas"]}
self.plugin_path = os.path.join(os.path.dirname(__file__), "..")
```

`self.citation` is legacy text and **SHOULD** stay empty for new plugins.
The formatted citation comes from `BasePlugin.get_citation()`.

### 2.2 `config.json`

A plugin that wraps an external binary **MUST** have an entry in
`YR_MPE/config.json`:

```json
{
  "name": "MAFFT",
  "path": "/bin/mafft/mafft.bat",
  "citations": ["katoh2013mafft"]
}
```

Rules:

- `tool_name` in the plugin and `name` in `config.json` must match
  case-insensitively; `config()` resolves by this name.
- `path` is relative to `YR_MPE/` and starts with `/`.
- Do not hard-code `E:/...`, `C:/...`, or a user home directory.
- A built-in plugin with no binary **MAY** override `config()` to return
  `True`, but must still expose citation keys. Do not create a fake binary
  entry.
- A plugin that needs a second binary (BLAST+, MPI, BEAGLE, a model tool)
  **MUST** declare and validate it by name before `run_analysis()`; do not
  assume it is on `PATH`.
- If the binary or its supported version changes, update the citation entry
  (if the paper changed), `settings.json` version, runtime `MANIFEST.tsv`,
  and `SHA256SUMS`.

### 2.3 Citations

Citations are data, never inline formatted text:

- store the BibTeX/RIS entry in `YR_MPE/citations/software.bib` (or `.ris`);
- store only the citation key in `config.json` or
  `BasePlugin.citation_keys`;
- use one `@software`/`@article` entry per tool, with DOI when available;
- never paste an APA/MLA string into a plugin constant;
- if a plugin bundles several tools, list all keys in the order they should
  be displayed.

## 3. UI contract

### 3.1 Inherited four tabs

Use `BasePlugin`'s tabs; do not recreate them:

1. **Input & Parameters** — `setup_input_tab()`;
2. **Output Preview** — `setup_output_tab()`;
3. **Console** — inherited;
4. **Citation** — inherited.

The input tab **SHOULD** be a `QScrollArea` containing `QGroupBox`
sections: input files, parameters, runtime/engine. Long parameter forms must
scroll rather than force the dialog taller than the screen.

### 3.2 Controls and text

- UI text is English, sentence case. Use `palette(mid)` for hint text and
  never hard-code a light background.
- Use `QFormLayout` for label/control pairs and a `QHBoxLayout` for
  file-edit + Browse.
- Put the primary action in the inherited bottom control bar; override the
  run button text only when it clarifies the operation
  (`self.run_button.setText("Run Z-test")`).
- Do not add a second top-level `QTabWidget`; output sub-tabs are allowed
  inside the Output Preview tab.
- Use an icon from `YR_MPE/icons/` (prefer `software/<tool>.svg`). Do not
  embed base64 images or add a raster icon when an SVG exists.
- Disabled placeholder actions must say why: `"CodeML (planned)"`, not a
  silent no-op.
- Tooltips and hints must not rely on colour alone.

### 3.3 Run lifecycle

`run_analysis()` **MUST**:

1. validate all inputs and binaries synchronously, before starting a thread;
2. show a `QMessageBox.warning()` with a concrete reason on validation
   failure;
3. set `is_running = True`, disable Run, enable Stop, show the progress bar
   (indeterminate if the total is unknown);
4. switch to the Console tab and log the exact command/parameters;
5. start the thread and connect `progress`, `console_output`, `result_ready`
   (if used), `finished`, and `error`;
6. restore the UI in both success and error handlers.

`finished` means "the worker stopped successfully"; `error` means "the
worker failed". Do not emit `error` and then also present a success state.
Cancellation **MUST** be implemented with `thread.stop()` + `thread.wait()`;
never destroy a running `QThread`.

`closeEvent()` **MUST** ask before closing while `is_running` is true, then
stop and wait. An agent-initiated close must refuse rather than show a modal
dialog (section 9).

### 3.4 Output preview

The output tab **MUST** be useful without opening an external folder:

- summary counts/status at the top;
- a table/tree/text preview of the parsed result;
- explicit empty state before a run;
- an export action for tables and figures when applicable.

Do not parse the output file only inside `paintEvent()` or a button handler;
keep parsing in a pure helper so it can be tested headlessly.

## 4. Process/thread contract

### 4.1 `BaseProcessThread`

All external processes **MUST** run through `BaseProcessThread` (or an
equivalent QThread that uses the same signal names and hidden-process
helpers). The thread:

- emits `progress(str)` for user-facing progress;
- emits `console_output(str, str)` for every useful line;
- emits `result_ready(object)` for the parsed result when the plugin uses
  one;
- emits `finished(output_files, html_files)` exactly once on success;
- emits `error(str)` on failure;
- checks `self.is_running` between loop iterations;
- never touches a `QWidget`, `QMessageBox`, or `QFileDialog`.

`get_tool_name()` **MUST** return the tool name shown in progress/logs.

### 4.2 Hidden subprocesses

`run_hidden()` / `popen_hidden()` from
`YR_MPE/runtime/subprocess_utils.py` are **MANDATORY** on Windows. Direct
`subprocess.run()` / `subprocess.Popen()` calls in plugins are forbidden
because they allocate a console window when YRTools runs under `pythonw`.
The only exception is a plugin whose purpose is an interactive console
(`console_plugin.py`); document that choice and never use it for an analysis
tool.

Use argument lists, not shell strings:

```python
cmd = [str(self.tool_path), "-s", str(alignment), "-T", str(threads)]
result = self.execute_command(cmd, cwd=str(output_dir))
```

Pass `cwd` explicitly. Do not rely on the process working directory, and do
not put user paths inside a shell-quoted string.

### 4.3 stdout, stderr, and exit status

- Classify by the process return code, not by the stream name. A tool may
  write progress to stderr; that is not automatically an error.
- Echo command lines with `console_output.emit(command, "command")`.
- Emit `console_output.emit(line, "output")` for normal stdout/stderr and
  `"warning"` / `"error"` only for actual problems.
- On a non-zero return code, include the last N lines of stdout/stderr in the
  `error` message; do not silently return an empty result.
- If the tool writes a known report file, parse the file as well as stdout;
  do not depend on log formatting when a stable file exists.

### 4.4 Working directories and files

- Use `self.workdir` when the host supplies it; otherwise create a
  `tempfile.mkdtemp(prefix="yr_mpe_<tool>_")`.
- Write only under `workdir` or a user-selected output directory.
- Never write into `YR_MPE/plugins/`, `YR_MPE/bin/`, or the package root.
- Never overwrite an input file. Use `--redo`/`--force` only when the user
  explicitly requested it.
- Clean only temporary files created by the plugin; keep user outputs and
  artifacts for provenance.
- Store output paths in the result (`DatasetItemSpec.file_path`,
  `ResultSpec.artifacts`), not only in the console log.

### 4.5 Cancellation

Long loops must check `self.is_running`. `stop()` should kill the process
tree on Windows (the base class does this), wait for readers to finish, and
emit a console warning. Do not use `QThread.terminate()`.

## 5. Input and import contract

### 5.1 Two entry paths

Every plugin must work in both modes:

1. **Standalone** — no `import_data`; the user selects files/folders in the
   input tab.
2. **From YR-MPEA** — `import_from == "YR_MPEA"`; the host supplies a file
   path or a structured payload. The plugin should preselect the file/edit
   fields and still let the user adjust parameters.

`BasePlugin.handle_import_data()` handles only a list of Biopython
`SeqRecord`s. A plugin that accepts richer input (trees, tree sets, distance
matrices, model reports) **MUST** override it and document the accepted
shape.

### 5.2 Input validation

Before starting a run, check:

- the file/folder exists and is readable;
- the format and suffix are accepted (FASTA, Newick, NEXUS, TSV, ...);
- the minimum number of records/trees/partitions is present;
- no empty sequences, duplicate names, or internal stops when the method
  requires it;
- the output directory is writable.

Warnings that do not invalidate the run go to the Console and, if the plugin
produces a result, to `ResultSpec.warnings`. Errors must name the offending
file and the expected format.

### 5.3 Selected basket payloads

When the host sends several selected objects, the wrapper is responsible for
converting them to a temporary file or a dict payload. The plugin must not
walk the workspace. Recommended payload shapes:

```python
# one alignment
import_data = {"alignment_file": "...", "display_name": "OG0000008"}

# a tree set plus a species tree
import_data = {
    "tree_file": "...",
    "labels": ["gene tree 1", "gene tree 2"],
    "groups": [
        {"group_id": "treeset:...", "group_name": "Gene trees", "kind": "treeset"},
        {"group_id": "tree:...", "group_name": "Species tree", "kind": "tree"},
    ],
}
```

Keep the plugin usable when the same file is opened manually; the payload is
a convenience, not the only input path.

## 6. Plugin result contract

This is the most important migration target. New analysis plugins **MUST NOT**
add a bespoke `export_*_signal`. They build a `ResultSpec` and call
`self.export_result(result)`.

### 6.1 Objects

From `YR_MPE/runtime/results.py`:

```text
ResultSpec
  source: str
  datasets: list[DatasetSpec]
  artifacts: list[str]
  summary: dict
  warnings: list[str]
  schema = "yrmpe.result.v1"

DatasetSpec
  name, description, source, summary, metadata
  items: list[DatasetItemSpec]

DatasetItemSpec
  loci_name, item_type, sequences, file_path, format,
  is_aligned, length, sequence_count, data, metadata,
  parent_index, selected, import_id, run_id
```

`export_result()` validates and emits on `BasePlugin.result_ready`; the host
connects that signal once to `YR_MPEA_Widget.receive_plugin_result()`.

### 6.2 Item types

Use the constants in `runtime.results`; do not invent strings.

| item_type | `file_path` should point to | `data` preview |
|---|---|---|
| `sequence` | FASTA | `orthogroup`, source, counts |
| `alignment` | aligned FASTA | `is_aligned=True`, trim stats |
| `model` | `.iqtree`/JSON/report | model name, parameters, BIC |
| `distance` | TSV/PHYLIP matrix | matrix labels + values |
| `phylogeny` | Newick | `newick`, support scheme |
| `treeset` | multi-tree Newick | `data`: list of tree dicts |
| `chain` | MCMC chain file | chain count, burn-in, tool |
| `variant` | VCF/TSV | variant counts |
| `coalescent` | Newick/summary | method, support |
| `clock` | dating tree/table | clock model, calibration |

Every item **MUST** have either `sequences` or `file_path`. For model,
distance, tree, tree set, and chain items, provide a real file path so the
workspace can persist and reopen the object.

### 6.3 Provenance and grouping

- `DatasetSpec.source` and `ResultSpec.source` use the same stable id.
- `DatasetSpec.metadata` records input paths, parameters, tool version, and
  runtime path; do not put only a human sentence.
- `DatasetSpec.summary` holds counts used by the dataset card.
- `DatasetItemSpec.parent_index` refers to the index of the parent item in
  the same `DatasetSpec.items` list. It **MUST** point backwards (acyclic).
- `DatasetItemSpec.run_id` is filled when the host supplies a run context;
  never invent a run id.
- `ResultSpec.artifacts` lists reports, logs, and files that are not items.
- `ResultSpec.warnings` lists non-fatal issues that a user should see.
- All values in `data` and `metadata` must be JSON-serializable; convert
  `Path`, tuples, and Biopython objects to strings/lists before export.

### 6.4 Minimal result example

```python
result = ResultSpec(
    source="mafft",
    datasets=[DatasetSpec(
        name="MAFFT alignment",
        source="mafft",
        summary={"sequences": len(records), "length": length},
        metadata={"input": str(alignment_path), "mode": "auto"},
        items=[DatasetItemSpec(
            loci_name=Path(alignment_path).stem,
            item_type="alignment",
            sequences=records,
            file_path=str(out_path),
            format="fasta",
            is_aligned=True,
            data={"tool": "MAFFT", "version": version},
        )],
    )],
    artifacts=[str(log_path)],
)
self.export_result(result)
```

For multi-output methods, put all outputs in one dataset or in several
datasets with explicit names; do not emit one dataset per output file by
accident.

### 6.5 Viewers and non-dataset results

A viewer may keep its own result display, but any reusable parsed object
(tree, distance, chain) **SHOULD** be exportable as a `DatasetItemSpec` so the
workspace, inspector, and agents can see it. A selection table or RF matrix
that is only a report may be saved as a TSV and shown in the output tab; do
not emit an empty `ResultSpec` just to carry an artifact path — `ResultSpec`
requires at least one dataset. If the report should be available to agents,
add a preview/export action instead of a fake dataset.

## 7. Run/provenance instrumentation

The workspace model is `item -> run -> output item`; a run is an analysis
edge with `input_ids`, `output_ids`, `method`, `tool`, `parameters`, and
`parameter_schema`. See `dev/storage/WORKSPACE_STORE.md`.

**Implemented today:** `WorkspaceStore.new_run()` / `finish_run()` and
`DatasetSelectionManager.create_run()` / `finish_run()` exist, but the
generic plugin wrapper does not call them for every plugin yet.

**Standard for new plugins:**

- If the host passes a run context (`run_id`, `input_ids`, `dataset_id`),
  fill `DatasetItemSpec.run_id` on all outputs and do not create a second
  run.
- If the plugin starts its own analysis inside YR-MPEA, the host wrapper
  **SHOULD** create the run before opening the plugin and finish it after a
  successful export. Do not fabricate provenance.
- Until the host wrapper is migrated, record method/tool/parameters in
  `DatasetSpec.metadata` and `RunRecord.parameters` when available; a missing
  run edge is a migration item, not a reason to block the plugin.

The target parameter schema for `/inspector.rerun_run` and agents is:

```json
{
  "model": {"type": "string", "required": true, "default": "GTR+F+G4",
            "choices": ["GTR+F+G4", "HKY+F+G4"], "label": "Model"},
  "threads": {"type": "integer", "default": 2, "min": 1, "max": 16}
}
```

The same dict drives the form, validation, command construction, and agent
call. Do not keep a second "agent-only" parameter path.

## 8. Agent / semantic API

A new analysis plugin **MUST** expose the same JSON-friendly surface as the
One-Step OrthoFinder plugin when it has user-configurable options. A
single-action plugin with no options **MUST** still expose `inspect_status`,
`inspect_console`, `inspect_actions`, and `agent_invoke`; its option methods
may return an empty/read-only shape. Viewers may skip the option surface only
if they have no meaningful actions, and must document that decision. The
complete example is `dev/agent/ORTHOFINDER_API.md`; the minimum method set is:

```python
inspect_options(as_json=False)       # current values + constraints
modify_options(options)              # atomic patch; rejects unknown keys
inspect_status()                     # idle/running/completed/failed/stopped
inspect_console(since=0, limit=200, level=None)
preview_results()                    # read-only summary + files
inspect_actions()                    # id/label/enabled/reason
agent_contract()                     # schema, state, actions
agent_invoke(action_id, params)      # {"ok": true/false, ...}
```

Rules:

- All return values are JSON-serializable.
- `modify_options` is atomic: one invalid key rejects the whole patch and
  leaves the UI/state unchanged.
- The agent path never opens a modal dialog. Use `show_dialogs=False` or an
  equivalent guard.
- `inspect_status()` uses the state enum `idle|running|completed|failed|
  stopped`, exposes `can_run`, `can_stop`, `can_export`, `can_close`, and
  reports progress as `{text, visible, indeterminate, value}`.
- `inspect_console()` is a cursor window with `total`, `since`, `next`,
  `lines`; callers pass `next` as the next `since`.
- `agent_invoke` returns `{"ok": false, "error": "..."}` rather than raising
  through the bridge.
- GUI buttons and agent actions call the same handlers. Never synthesize
  mouse/keyboard events for normal automation.
- The bridge only dispatches registered actions. Do not add an `eval`,
  arbitrary-shell, or arbitrary-Python action.

Legacy plugins add this surface in Phase C of the migration policy; a new
plugin should not merge without it. A viewer or a single-button utility may
skip the option methods only when it has no meaningful options; mark that
decision in the module docstring.

## 9. Host registration and menu wiring

### 9.1 Inside YR-MPEA

1. Add `PluginFactory.get_<id>_plugin()` in
   `YR_MPE/factories/plugin_factory.py`.
2. Add a `QAction` to the correct toolbar menu in
   `YR_MPE/platforms/mpea.py`.
3. Add an `open_<tool>_wrapper()` method that:
   - creates a `QDialog` (minimum 800×600, tool icon);
   - calls `entry.run(import_from=..., import_data=..., workdir=self.get_workdir())`;
   - connects `result_ready` to `_connect_plugin_result(wrapper)` for data
     plugins;
   - adds the widget and calls `dialog.exec_()`.
4. For a data plugin, verify the result reaches
   `YR_MPEA_Widget.receive_plugin_result()` and creates a dataset.
5. Add the tool icon under `YR_MPE/icons/software/`.
6. Add a user-manual section and a `docs/tools_inventory.md` row.

### 9.2 Standalone YRTools entry

Add an entry to `settings.json` only if the plugin should also appear
outside YR-MPEA:

```json
{
  "name": "MAFFT",
  "version": "7.520",
  "entry_point": "YR_MPE.plugins.mafft_plugin:MAFFTPluginEntry",
  "runtime": {"type": "polypyplug", "dependencies": []}
}
```

The entry point must be constructible with no arguments (or only optional
`config`/`plugin_path`) and accept `run()` with no arguments. Internal
factory registration does not replace `settings.json`; the two serve
different entry paths.

### 9.3 Do not use the legacy adapter layer

`PluginAdapter` and `PluginExecutor` are a previous integration layer. New
plugins do not register with them and do not add signals they expect. The
supported path is `BasePlugin` + `PluginEntry` + `PluginFactory` (and
`settings.json` for the standalone launcher).

## 10. Error handling, logging, and UX

- Log the exact command and the output directory before execution.
- Use `add_console_message(message, level)` with
  `info|command|output|warning|error`; do not use `error` for expected
  progress.
- A failed validation must not start a thread or write output files.
- A failed run must leave the previous result visible, or show a clear
  "failed" state; never display a partial result as complete.
- Report tool version when available (`tool --version` or the report file)
  in the console and in `DatasetSpec.metadata`.
- Do not use a modal dialog for every warning; the console is the primary
  log, dialogs are for decisions and blocking errors.
- Long jobs need a Stop button and a visible progress state. Do not freeze
  the GUI waiting for `QThread.wait()` without a stop path.
- On Windows, never assume a console exists; see section 4.2.

## 11. Testing standard

Every plugin **MUST** have at least:

1. a pure-helper test (parser, command builder, result builder) that runs
   without a display;
2. an offscreen GUI smoke test that constructs the plugin, checks metadata,
   validates empty input, and exercises the result preview with fixture data;
3. an integration test against the packaged binary when the plugin wraps one,
   skipped when the binary is absent.

Test location and naming:

```text
dev/<area>/test_<plugin>.py
```

Use the standard offscreen preamble:

```python
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
app = QApplication.instance() or QApplication([])
```

Integration tests use:

```python
@unittest.skipUnless(TOOL.is_file(), "runtime not available")
class ToolIntegrationTests(unittest.TestCase):
    ...
```

Required assertions for a data plugin:

- `ResultSpec.validate() == []`;
- every `DatasetItemSpec` has a file path or sequences;
- `data`/`metadata` survive a JSON round trip;
- `parent_index` is acyclic and points backwards;
- the output file named by `file_path` exists at export time;
- `export_result()` returns `True` and a connected receiver gets the object.

Run the whole `dev/` suite serially before merging. Do not commit generated
outputs, caches, or temporary files.

## 12. Security and filesystem rules

- Treat imported files as untrusted input: validate suffixes, sizes, and
  record counts; never execute a file from the input directory.
- Pass user paths as separate command arguments; no `shell=True`.
- Do not write outside the workspace/output directory chosen by the user.
- Do not delete or overwrite a user file without an explicit confirmation.
- Do not embed tokens, absolute machine paths, or credentials in a plugin.
- The local agent bridge exposes only registered actions; do not add an
  action that runs arbitrary commands.
- A plugin that downloads data must make the network access explicit and
  report the URL and checksum in the result metadata.

## 13. Definition of done

A new analysis plugin is ready when all of these are true:

- [ ] module, thread, plugin, and entry classes follow section 1.1;
- [ ] `plugin_name`, `tool_name`, `input_types`, `output_types`, and
      `plugin_path` are set;
- [ ] external binary path is in `config.json`, with citation keys;
- [ ] input validation happens before the thread starts;
- [ ] all child processes use `run_hidden`/`popen_hidden`;
- [ ] the worker emits `progress`, `console_output`, `result_ready` (if
      used), `finished`, and `error` correctly;
- [ ] `stop_analysis()` stops and waits;
- [ ] data results are a validated `ResultSpec`, not a bespoke signal;
- [ ] every result item has sequences or a real file path;
- [ ] `Data`/`metadata` are JSON-serializable;
- [ ] the host imports the result into the workspace;
- [ ] the plugin works both standalone and from Selected;
- [ ] new analysis plugin exposes the semantic status/console/actions/
      `agent_invoke` surface (option methods when it has options);
- [ ] an icon, menu entry, settings entry (if standalone), and manual
      section exist;
- [ ] pure, GUI, and integration tests pass;
- [ ] no console window appears under `pythonw`;
- [ ] no inline formatted citation remains.

## 14. Legacy plugin migration policy

Do not rewrite a plugin in one pass. Migrate in this order, keeping the
plugin runnable after each step:

### Phase A — mechanical safety

- replace direct `subprocess.run`/`Popen` with `run_hidden`/`popen_hidden`;
- move blocking work out of the GUI thread into `BaseProcessThread`;
- replace inline citation strings with citation keys;
- make tool paths come from `config.json`;
- add the missing input validation and Stop handling.

### Phase B — result contract

- replace bespoke `export_*_signal` emission with
  `DatasetItemSpec`/`DatasetSpec`/`ResultSpec` and `export_result()`;
- write real artifact files for model/distance/tree/chain results;
- connect `result_ready` to `_connect_plugin_result()` in the host wrapper;
- keep the old signal temporarily only if another caller still needs it,
  and mark it deprecated.

### Phase C — provenance and agents

- create/finish the `RunRecord` for each analysis edge;
- add `parameter_schema` and `RunRecord.parameters`;
- add `inspect_options` / `modify_options` / `inspect_status` /
  `inspect_console` / `preview_results` / `inspect_actions` /
  `agent_invoke` when the plugin has meaningful options.

### Mapping table

| Legacy pattern | Target |
|---|---|
| `export_model_result_signal` | `DatasetItemSpec(item_type="model")` |
| `export_distance_result_signal` | `DatasetItemSpec(item_type="distance")` |
| `export_phylogeny_result_signal` | `DatasetItemSpec(item_type="phylogeny")` + TreeSet |
| `export_chain_result_signal` | `DatasetItemSpec(item_type="chain")` |
| `DatasetObject` / `TempDatasetItem` | `DatasetSpec` / `DatasetItemSpec` |
| plugin-specific output dir | `workdir` + `RunRecord` outputs |
| formatted citation constant | citation key in `config.json` / `citations/` |
| `QThread` + ad-hoc subprocess | `BaseProcessThread` |
| custom export dialog | output tab + `ResultSpec` |
| `PluginAdapter` registration | `PluginEntry` + `PluginFactory` |

## Appendix A — analysis-plugin contract skeleton

This skeleton shows the required shape, not a complete tool. Copy the
validated lifecycle and result export from `rf_plugin.py` and
`orthofinder_plugin.py`, and keep parsing in pure helpers.

```python
from __future__ import annotations

import os
from pathlib import Path

from Bio import SeqIO
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QFormLayout, QGroupBox, QLabel, QLineEdit, QMessageBox, QVBoxLayout,
)

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread


class ExampleThread(BaseProcessThread):
    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_file, output_file, workdir=None):
        super().__init__(tool_path, [input_file], [], workdir=workdir)
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)
        self.tool_path = Path(tool_path)

    def get_tool_name(self):
        return "Example tool"

    def execute_commands(self):
        cmd = [str(self.tool_path), "-i", str(self.input_file),
               "-o", str(self.output_file)]
        self.progress.emit("Running Example tool...")
        result = self.execute_command(cmd, cwd=str(self.output_file.parent))
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "Example tool failed")
        if not self.output_file.is_file():
            raise RuntimeError("Example tool did not write %s" % self.output_file)
        self.result_ready.emit({"output": str(self.output_file)})
        self.finished.emit([str(self.output_file)], [])


class ExamplePlugin(BasePlugin):
    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._payload = None
        self.thread = None
        super().__init__(import_from, import_data, **kwargs)
        self.run_button.setText("Run Example tool")

    def init_plugin_info(self):
        self.plugin_name = "Example tool"
        self.tool_name = "Example tool"
        self.citation = []
        self.citation_keys = ["example2026"]
        self.input_types = {"Alignment FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"Alignment FASTA": ["fas"]}
        self.plugin_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), ".."
        )

    def setup_input_tab(self):
        layout = QVBoxLayout(self.input_tab)
        group = QGroupBox("Input")
        form = QFormLayout(group)
        self.input_edit = QLineEdit()
        form.addRow("Alignment:", self.input_edit)
        layout.addWidget(group)
        layout.addStretch()

    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        self.summary = QLabel("No result yet.")
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        layout.addStretch()

    def run_analysis(self):
        source = self.input_edit.text().strip()
        if not Path(source).is_file():
            QMessageBox.warning(self, "Input missing", "Select an alignment.")
            return
        if not self.tool_path or not Path(self.tool_path).is_file():
            QMessageBox.warning(self, "Tool missing", "Configure Example tool.")
            return

        base = Path(self.workdir) if self.workdir else Path.cwd()
        output = base / "example_output.fas"
        output.parent.mkdir(parents=True, exist_ok=True)
        self._payload = None
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running Example tool...")
        self.tab_widget.setCurrentIndex(2)
        self.add_console_message("Running Example tool on %s" % source, "info")

        self.thread = ExampleThread(
            self.tool_path, source, output, workdir=self.workdir
        )
        self.thread.progress.connect(self.progress_bar.setFormat)
        self.thread.console_output.connect(self.add_console_message)
        self.thread.result_ready.connect(self._on_result)
        self.thread.finished.connect(self._finished)
        self.thread.error.connect(self._failed)
        self.thread.start()

    def _on_result(self, payload):
        self._payload = payload

    def _finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        if self._payload and self.result_has_consumer():
            try:
                self.export_result(self._build_result_spec())
            except Exception as exc:
                self.add_console_message(
                    "Result export failed: %s" % exc, "error"
                )
        self.summary.setText("Completed: %s" % ", ".join(output_files))
        self.tab_widget.setCurrentIndex(1)

    def _failed(self, message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(message, "error")
        QMessageBox.critical(self, "Example tool failed", message)

    def _build_result_spec(self):
        path = Path(self._payload["output"])
        records = list(SeqIO.parse(str(path), "fasta"))
        return ResultSpec(
            source="example_tool",
            datasets=[DatasetSpec(
                name="Example alignment",
                source="example_tool",
                summary={"sequences": len(records)},
                metadata={"input": self.input_edit.text().strip()},
                items=[DatasetItemSpec(
                    loci_name=path.stem,
                    item_type="alignment",
                    sequences=records,
                    file_path=str(path),
                    format="fasta",
                    is_aligned=True,
                    data={"tool": "Example tool"},
                )],
            )],
            artifacts=[str(path)],
        )

    def stop_analysis(self):
        if self.thread and self.thread.isRunning():
            self.thread.stop()
            self.thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)


class ExamplePluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return ExamplePlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
```

The skeleton still needs `setup_output_tab()` to show parsed records, and a
run/result test. Those are mandatory in section 11.

## Appendix B — reference files

| Concern | Reference |
|---|---|
| four-tab base, console, citations, export | `YR_MPE/templates/base_plugin_ui.py` |
| worker signals, hidden process execution, stop | `YR_MPE/templates/base_process_thread.py` |
| result objects and validation | `YR_MPE/runtime/results.py` |
| citation lookup and BibTeX/RIS parsing | `YR_MPE/runtime/citations.py` |
| hidden subprocess helpers | `YR_MPE/runtime/subprocess_utils.py` |
| ResultSpec producer | `YR_MPE/plugins/orthofinder_plugin.py` |
| built-in analysis plugin | `YR_MPE/plugins/codon_ztest_plugin.py` |
| process + UI + result preview | `YR_MPE/plugins/rf_plugin.py` |
| agent option/status/console/results API | `dev/agent/ORTHOFINDER_API.md` |
| agent bridge transport | `dev/agent/AGENT_BRIDGE.md` |
| workspace model and provenance | `dev/storage/WORKSPACE_STORE.md` |
| result contract tests | `dev/results/test_results.py` |
