# Workflow Runner skill

## Purpose

Execute a compiled `yrmpe.workflow.v1` plan. Each step is an ordinary YR-MPE
plugin opened in headless mode; the runner applies the step parameters, runs
the tool, collects the typed result and feeds the produced files to later
steps through placeholder bindings.

## Input

A JSONL file or an agent transcript. The JSONL form uses one object per line:

```json
{"schema": "yrmpe.workflow.v1", "kind": "workflow", "name": "demo",
 "inputs": {"raw": "raw.fas"}}
{"id": "mafft", "op": "mafft.run",
 "outputs": {"alignment": "aln.fas"},
 "params": {"input": "%input.raw%", "output": "%mafft.alignment.fas%"}}
```

- `inputs` maps a role name to an external file. The key becomes
  `%input.<role><ext>%`.
- A step `outputs` object names each produced artifact. The key becomes
  `%<step>.<role><ext>%` and is available to later steps.
- `params` keys are the plugin control ids from `control_bindings()`.
- Lines without `op`/`step` are treated as headers or pruned discovery
  actions.

## GUI

- **Input & Parameters**: workflow file, output directory and one path field
  per declared input role.
- **Procedures**: the compiled DAG, one square per tool and one circle per
  artifact.
- **Output Preview**: one entry per produced file.
- **Console**: each step's command line and tool output.

## Agent actions

- `workflow.inspect_plan` — compiled steps, inputs and current outputs.
- `workflow.validate` — checks op mapping, input files and output directory.
- `workflow.run` — compile (if needed) and execute.
- `workflow.stop` — stop the current step and the run.
- `workflow.inspect_status` — step count, progress and error.
- `workflow.inspect_console` — console tail.

Generic actions (`plugin.modify_options`, `plugin.run`, `plugin.stop`,
`plugin.inspect_console`) also work on this window.

## Contract

- The runner writes `<output_dir>/workflow_run.json` with the plan, per-step
  status, command line, produced files and error.
- At the end, all item outputs are merged into one `ResultSpec`
  (`source="workflow-runner"`) and sent to the workspace.
- A step without `_build_result_spec()` cannot be executed; the runner
  reports the step as failed instead of inventing an import.

## Pitfalls

- Placeholder paths in `params` must match an input role or a previous step
  output; the compiler rejects unknown placeholders.
- The declared output path is a hint. The actual file comes from the step's
  `ResultSpec`; later steps are rebound to that actual path.
- A step that opens a modal dialog (missing input, missing binary) cannot be
  used headlessly. Validate inputs and tool configuration first.
