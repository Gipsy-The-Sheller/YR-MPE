# One-Step OrthoFinder skill

## Purpose

Turn GB/GBFF/GFF inputs into a workspace dataset of single-copy orthogroups.
The plugin runs OrthoFinder, extracts CDS for each single-copy orthogroup and
imports the resulting alignments/sequences into YPhylo Workbench.

## When to use it

- You have organelle or bacterial genome annotations and want orthologous
  loci for downstream alignment and phylogenetics.
- You need a dataset that can be sent to MAFFT/MACSE, ModelTest-NG and
  RAxML-NG.

## Semantic actions

This plugin is the most complete agent-facing plugin. Use
`workbench.list_windows` to get its `window_id`, then
`workbench.invoke_plugin` with:

- `orthofinder.inspect_options` - current options and constraints.
- `orthofinder.modify_options` - atomically update options.
- `orthofinder.validate` - validate annotations and the OrthoFinder runtime.
- `orthofinder.run` - start without opening modal dialogs.
- `orthofinder.stop` - stop the running pipeline.
- `orthofinder.inspect_status` - state, progress and capabilities.
- `orthofinder.inspect_console` - cursor window of console lines.
- `orthofinder.inspect_actions` - native action list and `next_import`.
- `orthofinder.preview_results` - generated files and single-copy OGs.
- `orthofinder.export_to_workspace` - idempotent import of the completed run.

## Typical flow

```text
workbench.open_orthofinder
workbench.list_windows -> window_id
orthofinder.inspect_options
orthofinder.modify_options {source_type: "gbff", input_path: "...", output_dir: "...", genetic_code: 11, threads: 2}
orthofinder.validate
orthofinder.run
orthofinder.inspect_status
orthofinder.preview_results
orthofinder.export_to_workspace
workbench.workspace
```

After import, each single-copy orthogroup is a dataset item and can be sent to
the alignment tools.

## Pitfalls

- Validation must pass before `orthofinder.run`; check `inspect_status`.
- The export action is idempotent; call it again after a completed run instead
  of duplicating the dataset manually.
- The plugin cannot close while the pipeline is running.
