"""Tests for the placeholder workflow compiler prototype."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.runtime.workflow_compiler import (  # noqa: E402
    WorkflowCompileError, WorkflowPlan, compile_jsonl, compile_workflow,
)


def chain_spec():
    return {
        "name": "mafft-modeltest-raxml",
        "inputs": {"raw": "data/raw.fas"},
        "steps": [
            {
                "id": "mafft",
                "op": "mafft.run",
                "outputs": {"alignment": "runs/mafft/aln.fas"},
                "params": {
                    "input": "data/raw.fas",
                    "output": "runs/mafft/aln.fas",
                    "threads": 4,
                },
            },
            {
                "id": "modeltest",
                "op": "modeltest.run",
                "outputs": {"model": "runs/modeltest/model.out"},
                "params": {
                    "input": "runs/mafft/aln.fas",
                    "output": "runs/modeltest/model.out",
                },
            },
            {
                "id": "raxml",
                "op": "raxmlng.run",
                "outputs": {"tree": "runs/raxml/tree.nwk"},
                "params": {
                    "msa": "runs/mafft/aln.fas",
                    "model_file": "runs/modeltest/model.out",
                    "output": "runs/raxml/tree.nwk",
                },
            },
        ],
    }


class WorkflowCompilerTests(unittest.TestCase):
    def test_compiles_out_in_chain(self):
        plan = compile_workflow(chain_spec())
        steps = {step["id"]: step for step in plan.steps}
        self.assertEqual(steps["mafft"]["params"]["input"], "%input.raw.fas%")
        self.assertEqual(
            steps["mafft"]["params"]["output"], "%mafft.alignment.fas%"
        )
        self.assertEqual(
            steps["modeltest"]["params"]["input"], "%mafft.alignment.fas%"
        )
        self.assertEqual(steps["modeltest"]["depends_on"], ["mafft"])
        self.assertEqual(
            steps["raxml"]["params"]["msa"], "%mafft.alignment.fas%"
        )
        self.assertEqual(
            steps["raxml"]["params"]["model_file"], "%modeltest.model.out%"
        )
        self.assertEqual(steps["raxml"]["depends_on"], ["mafft", "modeltest"])
        self.assertEqual(
            [step["id"] for step in plan.steps],
            ["mafft", "modeltest", "raxml"],
        )

    def test_forward_reference_is_topologically_ordered(self):
        spec = {
            "steps": [
                {
                    "id": "raxml",
                    "op": "raxmlng.run",
                    "outputs": {"tree": "tree.nwk"},
                    "params": {"msa": "aln.fas"},
                },
                {
                    "id": "mafft",
                    "op": "mafft.run",
                    "outputs": {"alignment": "aln.fas"},
                    "params": {"input": "raw.fas"},
                },
            ]
        }
        plan = compile_workflow(spec)
        self.assertEqual([step["id"] for step in plan.steps], ["mafft", "raxml"])
        raxml = next(step for step in plan.steps if step["id"] == "raxml")
        self.assertEqual(raxml["params"]["msa"], "%mafft.alignment.fas%")

    def test_infers_output_keys(self):
        plan = compile_workflow({
            "steps": [
                {"op": "tool.run", "params": {"input": "a.fas", "output": "b.fas"}},
                {"op": "tool2.run", "params": {"input": "b.fas"}},
            ]
        })
        self.assertIn("%s01.output.fas%", plan.artifacts)
        self.assertEqual(
            plan.steps[1]["params"]["input"], "%s01.output.fas%"
        )

    def test_duplicate_outputs_are_rejected(self):
        spec = {
            "steps": [
                {"id": "a", "op": "x", "outputs": {"out": "same.fas"}},
                {"id": "b", "op": "y", "outputs": {"out": "same.fas"}},
            ]
        }
        with self.assertRaisesRegex(WorkflowCompileError, "duplicate output"):
            compile_workflow(spec)

    def test_cycles_are_rejected(self):
        spec = {
            "steps": [
                {"id": "a", "op": "x", "outputs": {"out": "a.fas"},
                 "params": {"input": "b.fas"}},
                {"id": "b", "op": "y", "outputs": {"out": "b.fas"},
                 "params": {"input": "a.fas"}},
            ]
        }
        with self.assertRaisesRegex(WorkflowCompileError, "cycle"):
            compile_workflow(spec)

    def test_undeclared_external_path_is_a_warning(self):
        plan = compile_workflow({
            "steps": [
                {"op": "tool.run", "outputs": {"out": "b.fas"},
                 "params": {"input": "external.fas"}},
            ]
        })
        self.assertTrue(any("external.fas" in warning for warning in plan.warnings))

    def test_bindings_resolve_outputs_and_inputs(self):
        plan = compile_workflow(chain_spec())
        bindings = plan.bindings(output_base="C:/work/run1")
        self.assertEqual(bindings["%input.raw.fas%"], "data/raw.fas")
        self.assertTrue(
            bindings["%mafft.alignment.fas%"].replace("\\", "/")
            .endswith("work/run1/mafft/alignment.fas")
        )
        resolved = plan.resolve(
            "prefix=%mafft.alignment.fas%;model=%modeltest.model.out%",
            bindings,
        )
        self.assertIn("alignment.fas", resolved)
        self.assertIn("model.out", resolved)

    def test_compiles_jsonl_commands(self):
        plan = compile_jsonl(
            "# mafft -> raxml\n"
            '{"id":"mafft","op":"mafft.run","outputs":{"alignment":"aln.fas"},'
            '"params":{"input":"raw.fas","output":"aln.fas"}}\n'
            '{"id":"raxml","op":"raxmlng.run","outputs":{"tree":"tree.nwk"},'
            '"params":{"msa":"aln.fas","output":"tree.nwk"}}\n'
        )
        self.assertEqual([step["id"] for step in plan.steps], ["mafft", "raxml"])
        self.assertEqual(plan.steps[1]["params"]["msa"], "%mafft.alignment.fas%")


    def test_prunes_read_only_and_failed_steps(self):
        plan = compile_workflow({
            "steps": [
                {"id": "inspect", "op": "workbench.skills",
                 "params": {}},
                {"id": "menu", "op": "workbench.menu", "params": {}},
                {"id": "bad", "op": "mafft.run",
                 "outputs": {"alignment": "bad.fas"},
                 "params": {"input": "raw.fas"},
                 "result": {"ok": False, "error": "not found"}},
                {"id": "mafft", "op": "mafft.run",
                 "outputs": {"alignment": "aln.fas"},
                 "params": {"input": "raw.fas", "output": "aln.fas"}},
            ]
        })
        self.assertEqual([step["id"] for step in plan.steps], ["mafft"])
        self.assertTrue(any("read-only" in warning for warning in plan.warnings))
        self.assertTrue(any("failed" in warning for warning in plan.warnings))

    def test_pruning_can_be_disabled(self):
        plan = compile_workflow({
            "steps": [
                {"id": "inspect", "op": "workbench.skills", "params": {}},
                {"id": "mafft", "op": "mafft.run",
                 "outputs": {"alignment": "aln.fas"},
                 "params": {"input": "raw.fas", "output": "aln.fas"}},
            ]
        }, prune_non_executable=False)
        self.assertEqual([step["id"] for step in plan.steps], ["inspect", "mafft"])

    def test_json_roundtrip(self):
        plan = compile_workflow(chain_spec())
        restored = WorkflowPlan.from_dict(plan.to_dict())
        self.assertEqual(restored.to_dict()["steps"], plan.to_dict()["steps"])
        self.assertEqual(
            restored.to_dict()["artifacts"], plan.to_dict()["artifacts"]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
