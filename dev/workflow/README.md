# Placeholder workflow compiler (prototype)

This is the MEGA-CC-like prototype for turning a flat list of tool calls into
a portable workflow. It is implemented in
`YR_MPE/runtime/workflow_compiler.py` and is not yet wired into the GUI or
agent bridge.

## What it does

Each step can declare output files. The compiler assigns each output a
symbolic placeholder and rewrites every later parameter that refers to the
same path:

```json
{
  "inputs": {"raw": "data/raw.fas"},
  "steps": [
    {
      "id": "mafft",
      "op": "mafft.run",
      "outputs": {"alignment": "runs/mafft/aln.fas"},
      "params": {"input": "data/raw.fas", "output": "runs/mafft/aln.fas"}
    },
    {
      "id": "modeltest",
      "op": "modeltest.run",
      "outputs": {"model": "runs/modeltest/model.out"},
      "params": {"input": "runs/mafft/aln.fas", "output": "runs/modeltest/model.out"}
    }
  ]
}
```

compiles to:

```json
{
  "schema": "yrmpe.workflow.v1",
  "inputs": {"raw": {"placeholder": "%input.raw.fas%", "original": "data/raw.fas"}},
  "steps": [
    {"id": "mafft", "op": "mafft.run", "depends_on": [],
     "params": {"input": "%input.raw.fas%", "output": "%mafft.alignment.fas%"},
     "outputs": {"alignment": "%mafft.alignment.fas%"}},
    {"id": "modeltest", "op": "modeltest.run", "depends_on": ["mafft"],
     "params": {"input": "%mafft.alignment.fas%", "output": "%modeltest.model.out%"},
     "outputs": {"model": "%modeltest.model.out%"}}
  ],
  "artifacts": {
    "%mafft.alignment.fas%": {"step": "mafft", "role": "alignment", "ext": ".fas"},
    "%modeltest.model.out%": {"step": "modeltest", "role": "model", "ext": ".out"}
  }
}
```

The compiler:

- collects all output declarations first, so forward references work;
- performs path alias matching (`normpath`/`normcase`, Windows separators);
- computes `depends_on` and emits a topological order;
- rejects duplicate outputs, unknown placeholders and cycles;
- warns about undeclared external paths.

## API

```python
from YR_MPE.runtime.workflow_compiler import WorkflowPlan, compile_workflow

plan = compile_workflow(spec)
plan.to_dict()

# Bind output placeholders under a run directory.
bindings = plan.bindings(output_base="C:/work/run-001")
for step in plan.resolved_steps(bindings):
    print(step["op"], step["params"])
```

`WorkflowPlan.from_dict()` validates and restores the JSON form.

## Deliberate limits

- Outputs must be declared with `outputs` or inferred from top-level keys such
  as `output`, `output_file`, `prefix`, `-o` or `--output`. It does not guess
  arbitrary tool-specific output names.
- A placeholder names an artifact boundary, not a glob. `%s03.aln.fas%` is
  explicit; `%*.fas%` is intentionally unsupported.
- Execution is not implemented here. The plan is a specification; the
  workspace store / run executor will bind placeholders to
  `runs/<run_id>/outputs/` and record statuses.
- Nondeterministic steps (MCMC, bootstrap) still need seed/version metadata
  before replay can claim reproducibility.

## CLI

```bash
python -m YR_MPE.runtime.workflow_compiler workflow.jsonl
python -m YR_MPE.runtime.workflow_compiler workflow.json --resolve C:/work/run-001
```

`--resolve` prints the step list with all placeholders bound under the given
run directory, which is the hand-off point to a future executor.

## Previewer and step configuration

`YR_MPE/platforms/methods/workflow_preview.py` renders a graph with data
circles, tool squares and orthogonal arrows. Open it from
`WORKSPACE -> Workflow Preview...`, or preview a compiled plan with
`workbench.preview_workflow`.

Tool nodes are clickable. In plan previews the node carries `op`, `inputs`
and `outputs`; clicking it opens the plugin in workflow mode:

- path inputs are locked to `%...%` placeholders and the browse buttons are
  disabled;
- the plugin Run button becomes **Close & Apply**;
- applying emits `workflow_configured` with the step id, placeholders,
  parameters, controls and CLI preview;
- the plan/executor can then store the configuration without re-running the
  tool inside the dialog.

## Importing a workflow

Use `WORKSPACE -> Import Workflow...`, or drop a `.jsonl`/`.json` file onto
the Workbench. `demo_pipeline.jsonl` is a small MAFFT -> ModelTest-NG ->
RAxML-NG example with a header line:

```json
{"schema":"yrmpe.workflow.v1","kind":"workflow","name":"...","inputs":{"raw":"raw.fas"}}
```

The importer compiles JSONL through `compile_jsonl()`, compiles a raw JSON
spec through `compile_workflow()`, and accepts an already compiled
`yrmpe.workflow.v1` plan. The resulting graph opens in the previewer.
