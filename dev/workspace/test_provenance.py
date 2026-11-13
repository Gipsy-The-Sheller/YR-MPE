"""Run provenance: every imported result becomes a persistent RunRecord.

Run from the repository root:

    python dev/workspace/test_provenance.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_models import DatasetItem  # noqa: E402
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402
from YR_MPE.runtime.results import (  # noqa: E402
    DatasetItemSpec,
    DatasetSpec,
    ResultSpec,
    RunSpec,
)

APP = QApplication.instance() or QApplication(sys.argv)

FASTA = ">a\nATGCAA\n>b\nATGCAA\n"


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="yrmpe-prov-")).resolve()
        self.inputs = self.root / "inputs"
        self.outputs = self.root / "outputs"
        self.inputs.mkdir()
        self.outputs.mkdir()
        self.input = self.inputs / "OG1.fna"
        self.input.write_text(FASTA, encoding="utf-8")

        self.manager = DatasetSelectionManager()
        self.manager.set_workspace_path(str(self.root))
        self.workspace = SingleGeneWorkspace(
            dataset_selection_manager=self.manager
        )

        dataset_id = self.manager.create_dataset("input")
        item = DatasetItem(item_type="sequence")
        item.loci_name = "OG1"
        item.dataset_id = dataset_id
        item.file_path = str(self.input)
        item.sequences = [
            SeqRecord(Seq("ATGCAA"), id="a"),
            SeqRecord(Seq("ATGCAA"), id="b"),
        ]
        item.sequence_count = 2
        item.length = 6
        self.assertTrue(self.manager.add_item(item, dataset_id))
        self.input_item_id = item.id
        self.manager._save_state()

    def _result(self):
        output = self.outputs / "OG1.aln.fna"
        output.write_text(
            ">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8"
        )
        return output, ResultSpec(
            source="mafft",
            tool="mafft",
            tool_version="7.520",
            seed="42",
            datasets=[DatasetSpec(
                name="MAFFT alignment",
                source="mafft",
                items=[DatasetItemSpec(
                    loci_name="OG1",
                    item_type="alignment",
                    file_path=str(output),
                    format="fasta",
                    is_aligned=True,
                    sequences=[
                        SeqRecord(Seq("ATGCAA"), id="a"),
                        SeqRecord(Seq("ATGCAA"), id="b"),
                    ],
                )],
                metadata={
                    "input": str(self.input),
                    "parameters": ["--auto", "--thread", "4"],
                },
            )],
            artifacts=[str(output)],
            summary={"command": ["mafft", "--auto", str(self.input)]},
        )

    def test_plugin_result_creates_run(self):
        output, result = self._result()
        dataset_ids = self.workspace.add_result_spec(result)
        self.assertEqual(len(dataset_ids), 1)

        runs = list(self.manager.store.index.runs.values())
        self.assertEqual(len(runs), 1)
        run = runs[0]
        self.assertEqual(run.method, "mafft")
        self.assertEqual(run.tool, "mafft")
        self.assertEqual(run.tool_version, "7.520")
        self.assertEqual(run.seed, "42")
        self.assertEqual(run.parameters["arguments"], ["--auto", "--thread", "4"])
        self.assertEqual(run.command, ["mafft", "--auto", str(self.input)])
        self.assertIn(self.input_item_id, run.input_ids)
        self.assertEqual(len(run.output_ids), 1)
        self.assertEqual(run.status, "completed")
        self.assertIn(str(output), run.artifacts)
        self.assertTrue(run.artifact_hashes[str(output)].startswith("sha256:"))

        item = self.manager.get_items_by_dataset(dataset_ids[0])[0]
        self.assertEqual(item.run_id, run.id)

        # Reopen: the run and the item back-reference both survive.
        reopened = DatasetSelectionManager()
        reopened.set_workspace_path(str(self.root))
        reopened_run = reopened.store.index.runs[run.id]
        self.assertEqual(reopened_run.method, "mafft")
        self.assertEqual(reopened_run.tool_version, "7.520")
        self.assertIn(item.id, reopened_run.output_ids)
        self.assertEqual(reopened.items[item.id].run_id, run.id)

    def test_workflow_runs_are_recorded_per_step(self):
        first = self.outputs / "step1.aln.fna"
        first.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        second = self.outputs / "step2.trim.fna"
        second.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        result = ResultSpec(
            source="workflow-runner",
            tool="workflow-runner",
            datasets=[DatasetSpec(
                name="workflow",
                source="workflow-runner",
                items=[
                    DatasetItemSpec(
                        loci_name="step1:OG1", item_type="alignment",
                        file_path=str(first), format="fasta", is_aligned=True,
                        sequences=[SeqRecord(Seq("ATGCAA"), id="a")],
                        metadata={"workflow_step": "step1"},
                    ),
                    DatasetItemSpec(
                        loci_name="step2:OG1", item_type="alignment",
                        file_path=str(second), format="fasta", is_aligned=True,
                        sequences=[SeqRecord(Seq("ATGCAA"), id="a")],
                        metadata={"workflow_step": "step2"},
                    ),
                ],
            )],
            artifacts=[str(first), str(second)],
            runs=[
                RunSpec(
                    method="mafft", tool="mafft",
                    command=["mafft", "--auto"],
                    parameters={"args": ["--auto"]},
                    input_paths=[str(self.input)],
                    output_paths=[str(first)],
                    metadata={"step_id": "step1", "workflow_id": "wf1"},
                ),
                RunSpec(
                    method="trimal", tool="trimal",
                    command=["trimal", "-in", str(first)],
                    parameters={"args": ["-gappyout"]},
                    input_paths=[str(first)],
                    output_paths=[str(second)],
                    metadata={"step_id": "step2", "workflow_id": "wf1"},
                ),
            ],
        )
        dataset_ids = self.workspace.add_result_spec(result)
        self.assertEqual(len(dataset_ids), 1)
        runs = sorted(
            self.manager.store.index.runs.values(), key=lambda run: run.method
        )
        self.assertEqual([run.method for run in runs], ["mafft", "trimal"])
        by_method = {run.method: run for run in runs}
        self.assertEqual(len(by_method["mafft"].output_ids), 1)
        self.assertEqual(len(by_method["trimal"].output_ids), 1)
        self.assertEqual(by_method["mafft"].metadata["step_id"], "step1")
        self.assertEqual(by_method["trimal"].metadata["workflow_id"], "wf1")
        # Step 2 consumes step 1's output, which is now a workspace item.
        self.assertEqual(
            by_method["trimal"].input_ids, by_method["mafft"].output_ids
        )
        for run in runs:
            self.assertTrue(run.artifact_hashes)
            self.assertEqual(run.status, "completed")

        items = self.manager.get_items_by_dataset(dataset_ids[0])
        self.assertEqual(
            {item.run_id for item in items},
            {by_method["mafft"].id, by_method["trimal"].id},
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
