# OrthoFinder plugin semantic API

The One-Step OrthoFinder plugin exposes a JSON-friendly control surface so an
agent can configure, run, monitor and export without clicking Qt widgets.

All methods return JSON-serializable dicts. Paths are strings.

## Options

```python
plugin.inspect_options()            # dict
plugin.inspect_options(as_json=True) # JSON text
```

Shape:

```json
{
  "schema": "yrmpe.orthofinder.options.v1",
  "tool": "OrthoFinder",
  "options": {
    "source_type": "gbff",
    "input_path": "...",
    "output_dir": "...",
    "genetic_code": null,
    "skip_pseudogenes": true,
    "internal_stop": "skip",
    "keep_terminal_stop": false,
    "run_mode": "full",
    "existing_results": "",
    "orthofinder_path": "...",
    "blast_dir": "...",
    "threads": 2
  },
  "constraints": {
    "source_type": ["gbff", "gbff_records", "proteome"],
    "genetic_code": [null, 1, 2, 4, 5, 11],
    "internal_stop": ["skip", "warn", "keep"],
    "run_mode": ["full", "parse_only", "existing"],
    "threads": {"min": 1, "max": 16}
  }
}
```

`modify_options()` accepts a dict, a JSON string, or the whole payload returned
by `inspect_options()` (the inner `options` key is unwrapped):

```python
plugin.modify_options({
    "source_type": "gbff",
    "input_path": "E:/data/gb",
    "output_dir": "E:/data/orthofinder",
    "genetic_code": 5,
    "skip_pseudogenes": True,
    "internal_stop": "warn",
    "run_mode": "parse_only",
    "threads": 4,
})
```

It applies atomically: one invalid/unknown key rejects the whole patch.

```json
{
  "ok": false,
  "applied": [],
  "rejected": {"threads": "threads must be between 1 and 16"},
  "errors": ["threads must be between 1 and 16"],
  "options": { "...": "unchanged current options" }
}
```

## Status and console

```python
plugin.inspect_status()
plugin.inspect_console(since=0, limit=200, level=None)
```

`inspect_status()`:

```json
{
  "schema": "yrmpe.orthofinder.status.v1",
  "state": "idle|running|completed|failed|stopped",
  "running": false,
  "wizard_step": 0,
  "progress": {"text": "", "visible": false, "indeterminate": true, "value": null},
  "elapsed": null,
  "summary": {},
  "last_error": "",
  "console_lines": 0,
  "can_run": true,
  "can_stop": false,
  "can_export": false,
  "can_close": true
}
```

`inspect_console()` returns a cursor-friendly window:

```json
{
  "total": 42,
  "since": 10,
  "next": 30,
  "lines": [
    {"index": 10, "timestamp": "22:13:01", "type": "info", "message": "..."}
  ]
}
```

Use `next` as the next `since`, so an agent can tail without replaying the log.

## Results

```python
plugin.preview_results()
```

Returns:

```json
{
  "ready": true,
  "output_dir": "...",
  "summary": {"species": 4, "cds": 335, "single_copy_ogs": 62},
  "orthogroups_tsv": ".../Orthogroups.tsv",
  "single_copy_ogs": [
    {"orthogroup": "OG0000008", "cds_fasta": "...", "protein_fasta": "...", "map_tsv": "..."}
  ],
  "files": [
    {"role": "cds_map", "path": "...", "exists": true, "size": 1234}
  ]
}
```

`files` lists the shared outputs (proteomes, CDS, cds_map, summary); each
single-copy OG's three files are in `single_copy_ogs`.

## Actions and next import

```python
plugin.inspect_actions()
```

Actions carry `id`, `label`, `enabled`, `reason` and optional parameter schema.
`next_import` points at the dataset import:

```json
{
  "next_import": {
    "action": "orthofinder.export_to_workspace",
    "target": "YR-MPEA dataset",
    "enabled": true
  }
}
```

`export_to_workspace()` returns `{ok, dataset, loci}` or `{ok: false, error}`.
It only works when the plugin was opened from YR-MPEA and the result channel is
connected.

## Dispatch and closing

```python
plugin.agent_contract()
plugin.agent_invoke("orthofinder.inspect_status")
plugin.agent_invoke("orthofinder.modify_options", {"options": {...}})
plugin.agent_invoke("orthofinder.run", {"show_dialogs": False})
plugin.agent_invoke("orthofinder.stop")
plugin.agent_invoke("orthofinder.preview_results")
plugin.agent_invoke("orthofinder.export_to_workspace")
plugin.agent_invoke("orthofinder.close")
```

`can_close()` returns `{can_close, running, reason}`. Closing while the pipeline
runs is refused for agents; the GUI close event asks the user whether to stop
and close. The host can decide whether to keep the object after `close()`.

## Design rules

- The agent path never opens modal dialogs (`start_pipeline(show_dialogs=False)`).
- `modify_options` is atomic; no half-applied state.
- Status/console/results are read-only and can be polled.
- `next_import` is explicit rather than inferred from button positions.
- The same methods back the GUI buttons; there is no parallel automation path.
