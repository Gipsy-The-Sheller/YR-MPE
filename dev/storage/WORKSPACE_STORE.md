# Workspace store (v1 design)

> Goal: replace the current mix of `workspace_state.json`,
> `dataset_selection_state.json`, `history/` and plugin-specific output
> directories with one index plus self-describing per-run folders.

## Layout

```
<workspace>/
  workspace.json
  imports/
    <import_id>/
      import.json
      files/...
  runs/
    <run_id>/
      manifest.json
      console.jsonl
      outputs/...
  datasets/
    <dataset_id>/
      items/
        <item_id>/
          data.<ext>
  recycled/
```

Inspired by PhyloSuite's per-work folders (`GenBank_File/<work>/files`,
`flowchart`, `recycled`) and its per-work settings, but with one JSON index
instead of scattering state across INI files.

## Core model

- **Item = data node.** A sequence set, alignment, model, distance matrix,
  tree, chain, treeset, etc.
- **Run = analysis method edge.** `inputs[] -> outputs[]`. A Bayesian run can
  produce several outputs (tree + chain + model parameters); all are listed.
- **Dataset = grouping of data nodes.** A dataset owns items, not runs.
- **Stream = ordered run chain inside a dataset.** Example:
  `OrthoFinder import -> MAFFT -> IQ-TREE`. Streams are named and persisted.
- **Import = a run with no inputs.** Every import gets its own `import.json`
  with source paths, SHA-256, importer, parameters, and the item IDs it
  created. This is the "precise location" requested: any item can be traced
  to the exact import and source file.

## workspace.json

```json
{
  "schema": "yrmpe.workspace.v1",
  "workspace_id": "...",
  "name": "...",
  "created_at": "...",
  "updated_at": "...",
  "datasets": {
    "<dataset_id>": {
      "name": "...",
      "description": "...",
      "source": "one_step_orthofinder",
      "item_ids": ["..."],
      "stream_ids": ["..."],
      "settings": {}
    }
  },
  "items": {
    "<item_id>": {
      "type": "sequence",
      "name": "OG0000008",
      "dataset_id": "...",
      "format": "fasta",
      "file": "datasets/<dataset_id>/items/<item_id>/data.fna",
      "length": 1059,
      "sequence_count": 4,
      "is_aligned": false,
      "metadata": {},
      "import_id": "...",
      "run_id": "...",
      "parent_id": null
    }
  },
  "runs": {
    "<run_id>": {
      "name": "IQ-TREE run 1",
      "method": "iqtree",
      "tool": "IQ-TREE 3.0.1",
      "parameters": {},
      "input_ids": ["..."],
      "output_ids": ["..."],
      "status": "completed",
      "created_at": "...",
      "completed_at": "...",
      "manifest": "runs/<run_id>/manifest.json",
      "stream_id": "..."
    }
  },
  "streams": {
    "<stream_id>": {
      "name": "OrthoFinder -> MAFFT -> IQ-TREE",
      "dataset_id": "...",
      "run_ids": ["..."],
      "created_at": "..."
    }
  },
  "imports": {
    "<import_id>": {
      "source_paths": ["..."],
      "source_hashes": {"...": "sha256:..."},
      "importer": "one_step_orthofinder",
      "created_at": "...",
      "item_ids": ["..."],
      "manifest": "imports/<import_id>/import.json"
    }
  },
  "settings": {}
}
```

## Run manifest

`runs/<run_id>/manifest.json` is self-describing even if workspace.json is
lost:

```json
{
  "schema": "yrmpe.run.v1",
  "run_id": "...",
  "name": "IQ-TREE run 1",
  "method": "iqtree",
  "tool": "IQ-TREE",
  "tool_version": "3.0.1",
  "command": ["iqtree3", "-s", "..."],
  "parameters": {},
  "inputs": ["<item_id>"],
  "outputs": ["<item_id>", "<item_id>"],
  "status": "completed",
  "created_at": "...",
  "completed_at": "...",
  "cwd": "...",
  "log": "console.jsonl",
  "artifacts": ["..."]
}
```

## Import manifest

```json
{
  "schema": "yrmpe.import.v1",
  "import_id": "...",
  "source_paths": ["..."],
  "source_hashes": {"...": "sha256:..."},
  "importer": "one_step_orthofinder",
  "parameters": {},
  "created_at": "...",
  "dataset_id": "...",
  "items": [
    {"id": "...", "type": "sequence", "name": "OG0000008",
     "file": "datasets/.../data.fna", "source": "..."}
  ]
}
```

## Lifecycle

1. `WorkspaceStore(root)` ensures the layout.
2. New import: `new_import(source_paths, importer, parameters)` copies files
   (or records external paths), computes hashes, writes `import.json`.
3. Imported items: `add_item(..., import_id=import_id, copy=...)`. The item
   file lives under `datasets/<dataset_id>/items/<item_id>/`.
4. Analysis: `new_run(method, tool, inputs, parameters)` creates the run dir
   and manifest; when the tool finishes, `finish_run(output_ids=...)`.
5. Each output item records `run_id`; the inspector can walk
   `item -> run -> inputs/outputs -> stream`.
6. `save()` writes `workspace.json` atomically. Item hydration is lazy from
   `file`; missing files are reported, not silently dropped.

## Migration

`migrate_legacy(dataset_selection_state.json)`:

- preserve dataset/item UUIDs so existing references still work;
- create dataset folders and copy item files when they exist;
- mark migrated items with `metadata.migrated_from = "dataset_selection_state.json"`;
- leave the old file as `dataset_selection_state.json.bak`.

`workspace_state.json` and `history/` are not migrated; they are redundant
with the new store.

## PhyloSuite mapping

| PhyloSuite | YR-MPE store |
|---|---|
| `GenBank_File/<work>/files` | `imports/<import_id>/files` |
| `Other_File/<work>/flowchart` | `runs/<run_id>/outputs` |
| `recycled` | `recycled/` |
| `data_settings.ini` | `workspace.json` + per-run manifest |
| per-work folder | per-import and per-run folders |

## Deliberately out of scope for v1

- SQLite/transactions.
- Content-addressed blobs and deduplication (hash is recorded, but files are
  copied).
- Automatic garbage collection.
- Full plugin instrumentation (plugins must call `new_run` / `finish_run`;
  this will be wired plugin by plugin).
