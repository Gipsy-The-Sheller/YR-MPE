"""Tests for workflow graph construction, layout and preview rendering."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.workflow_preview import (  # noqa: E402
    WorkflowPreviewDialog, WorkflowPreviewWidget,
)
from YR_MPE.runtime.workflow_compiler import compile_workflow  # noqa: E402
from YR_MPE.runtime.workflow_graph import (  # noqa: E402
    WorkflowGraph, graph_from_plan, layered_layout,
)

APP = QApplication.instance() or QApplication(sys.argv)


def sample_plan():
    return compile_workflow({
        "name": "preview-demo",
        "inputs": {"raw": "data/raw.fas"},
        "steps": [
            {"id": "mafft", "op": "mafft.run",
             "outputs": {"alignment": "aligned.fas"},
             "params": {"input": "data/raw.fas", "output": "aligned.fas"}},
            {"id": "modeltest", "op": "modeltest.run",
             "outputs": {"model": "model.out"},
             "params": {"input": "aligned.fas", "output": "model.out"}},
            {"id": "iqtree", "op": "iqtree.run",
             "outputs": {"tree": "iq.tree"},
             "params": {"msa": "aligned.fas", "output": "iq.tree"}},
            {"id": "raxml", "op": "raxmlng.run",
             "outputs": {"tree": "final.tree"},
             "params": {"msa": "iq.tree", "model_file": "model.out",
                        "output": "final.tree"}},
        ],
    })


class WorkflowGraphTests(unittest.TestCase):
    def test_plan_graph_has_branch_and_merge_without_self_edges(self):
        graph = graph_from_plan(sample_plan())
        node_ids = {node.id for node in graph.nodes}
        self.assertEqual(len([n for n in graph.nodes if n.kind == "tool"]), 4)
        self.assertEqual(len([n for n in graph.nodes if n.kind == "data"]), 5)
        edges = {(edge.source, edge.target) for edge in graph.edges}
        self.assertIn(
            ("artifact:%mafft.alignment.fas%", "step:modeltest"), edges
        )
        self.assertIn(
            ("artifact:%iqtree.tree.tree%", "step:raxml"), edges
        )
        self.assertIn(
            ("artifact:%modeltest.model.out%", "step:raxml"), edges
        )
        self.assertFalse(any(source == target for source, target in edges))
        self.assertTrue(edges <= {
            (source, target) for source in node_ids for target in node_ids
        })

    def test_layered_layout_places_tools_right_of_inputs(self):
        graph = graph_from_plan(sample_plan())
        positions = layered_layout(graph)
        self.assertEqual(len(positions), len(graph.nodes))
        by_id = graph.node_map()
        mafft = positions["step:mafft"]
        raw = positions["input:%input.raw.fas%"]
        alignment = positions["artifact:%mafft.alignment.fas%"]
        self.assertLess(raw[0], mafft[0])
        self.assertLess(mafft[0], alignment[0])

    def test_preview_widget_renders_nodes_and_edges(self):
        graph = graph_from_plan(sample_plan())
        widget = WorkflowPreviewWidget(graph)
        widget.resize(1200, 600)
        widget.show()
        for _ in range(4):
            APP.processEvents()
        self.assertEqual(
            len(widget._node_items),
            len([node for node in graph.nodes]),
        )
        pixmap = widget.grab()
        self.assertFalse(pixmap.isNull())
        self.assertGreater(pixmap.width(), 100)
        self.assertGreater(pixmap.height(), 100)
        widget.close()

    def test_dialog_exposes_preview(self):
        graph = graph_from_plan(sample_plan())
        dialog = WorkflowPreviewDialog(graph)
        try:
            self.assertIsInstance(dialog.preview, WorkflowPreviewWidget)
            self.assertIn("Workflow Preview", dialog.windowTitle())
        finally:
            dialog.close()

    def test_clicking_tool_node_calls_step_opener(self):
        graph = graph_from_plan(sample_plan())
        seen = []

        def opener(context):
            seen.append(context)
            return {"id": context.get("step_id"), "op": context.get("op")}

        dialog = WorkflowPreviewDialog(graph, step_opener=opener)
        try:
            dialog.preview.node_clicked.emit("step:mafft")
            self.assertEqual(len(seen), 1)
            self.assertEqual(seen[0]["op"], "mafft.run")
            self.assertIn("input", seen[0]["inputs"])
            self.assertEqual(
                graph.node_map()["step:mafft"].metadata["config"]["op"],
                "mafft.run",
            )
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
