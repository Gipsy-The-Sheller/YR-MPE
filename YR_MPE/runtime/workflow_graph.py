"""Graph model and simple layered layout for workflow previews.

The graph is deliberately Qt-free. The preview widget renders nodes and edges
from this model; agents can also fetch the same structure as JSON through the
host action layer if needed.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from .workflow_compiler import PLACEHOLDER_RE, WorkflowPlan

GRAPH_SCHEMA = "yrmpe.workflow.graph.v1"

DATA_ICONS = {
    "sequence": "file/sequence.svg",
    "alignment": "file/sequence.svg",
    "model": "file/model.svg",
    "distance": "file/distance.svg",
    "phylogeny": "file/phylogeny.svg",
    "treeset": "file/treeset.svg",
    "chain": "file/chain.svg",
    "dataset": "file/dataset.svg",
    "artifact": "file/sequence.svg",
}

TOOL_ICON_HINTS = [
    ("orthofinder", "file/dataset.svg"),
    ("modelfinder", "find_model.svg"),
    ("modeltest", "find_model.svg"),
    ("iq-tree", "software/iqtree.svg"),
    ("iqtree", "software/iqtree.svg"),
    ("raxml", "software/raxml.svg"),
    ("mafft", "software/mafft.svg"),
    ("macse", "software/macse.svg"),
    ("trimal", "software/trimal.svg"),
    ("gblocks", "software/gblocks.svg"),
    ("muscle", "software/muscle.svg"),
    ("clustal", "software/clustalo.svg"),
    ("mrbayes", "software/mrbayes.svg"),
    ("phylobayes", "software/phylobayes.svg"),
    ("mpboot", "mp.svg"),
    ("decenttree", "bionj.svg"),
    ("astral", "software/astral.svg"),
    ("caster", "software/caster.svg"),
    ("treerecs", "coalescent.svg"),
    ("lsd2", "software/lsd2.svg"),
    ("minitracer", "software/minitracer.svg"),
    ("rf", "software/robinsonfoluds.svg"),
    ("concaten", "software/concatenator.svg"),
]


@dataclass
class WorkflowNode:
    id: str
    kind: str  # "data" or "tool"
    label: str
    icon: str = ""
    item_type: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "kind": self.kind,
            "label": self.label,
            "icon": self.icon,
            "item_type": self.item_type,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "WorkflowNode":
        return cls(
            id=str(data.get("id", "")),
            kind=str(data.get("kind", "data")),
            label=str(data.get("label", "")),
            icon=str(data.get("icon", "")),
            item_type=str(data.get("item_type", "")),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass
class WorkflowEdge:
    source: str
    target: str
    label: str = ""

    def to_dict(self) -> Dict[str, str]:
        return {"source": self.source, "target": self.target,
                "label": self.label}


@dataclass
class WorkflowGraph:
    nodes: List[WorkflowNode] = field(default_factory=list)
    edges: List[WorkflowEdge] = field(default_factory=list)
    title: str = "Workflow"
    schema: str = GRAPH_SCHEMA

    def add_node(self, node: WorkflowNode) -> WorkflowNode:
        if not any(existing.id == node.id for existing in self.nodes):
            self.nodes.append(node)
        return node

    def add_edge(self, source: str, target: str, label: str = "") -> None:
        if source == target:
            return
        if any(edge.source == source and edge.target == target
               for edge in self.edges):
            return
        self.edges.append(WorkflowEdge(source, target, label))

    def node_map(self) -> Dict[str, WorkflowNode]:
        return {node.id: node for node in self.nodes}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "title": self.title,
            "nodes": [node.to_dict() for node in self.nodes],
            "edges": [edge.to_dict() for edge in self.edges],
        }


def icon_for_item_type(item_type: str) -> str:
    return DATA_ICONS.get(str(item_type or "").lower(), "file/sequence.svg")


def icon_for_tool(name: str) -> str:
    lowered = str(name or "").lower()
    for needle, icon in TOOL_ICON_HINTS:
        if needle in lowered:
            return icon
    return "find_model.svg"


def _label_from_path(path: str) -> str:
    base = os.path.basename(str(path).rstrip("/\\"))
    return base or str(path)


def graph_from_workspace(manager) -> WorkflowGraph:
    """Build a graph from workspace runs/items, falling back to plain items."""
    graph = WorkflowGraph(title="Workspace workflow")
    store = getattr(manager, "store", None) if manager is not None else None
    if store is None:
        if manager is not None:
            for item in manager.items.values():
                graph.add_node(WorkflowNode(
                    id="item:%s" % item.id,
                    kind="data",
                    label=getattr(item, "loci_name", "")
                          or getattr(item, "name", "") or item.id,
                    icon=icon_for_item_type(getattr(item, "item_type", "")),
                    item_type=getattr(item, "item_type", ""),
                ))
        return graph

    index = store.index
    records = index.items
    for run in index.runs.values():
        tool_id = "run:%s" % run.id
        input_placeholders = {
            ("input" if len(run.input_ids) == 1 else "input_%d" % index):
                "%%item.%s%%" % item_id
            for index, item_id in enumerate(run.input_ids)
        }
        output_placeholders = {
            ("output" if len(run.output_ids) == 1 else "output_%d" % index):
                "%%item.%s%%" % item_id
            for index, item_id in enumerate(run.output_ids)
        }
        graph.add_node(WorkflowNode(
            id=tool_id,
            kind="tool",
            label=run.tool or run.name or run.method,
            icon=icon_for_tool(run.tool or run.method or run.name),
            metadata={
                "run_id": run.id,
                "status": run.status,
                "method": run.method,
                "op": run.method,
                "inputs": input_placeholders,
                "outputs": output_placeholders,
                "command": list(run.command),
            },
        ))
        for item_id in run.input_ids:
            record = records.get(item_id)
            data_id = "item:%s" % item_id
            graph.add_node(WorkflowNode(
                id=data_id,
                kind="data",
                label=(record.name if record is not None else item_id),
                icon=icon_for_item_type(
                    record.type if record is not None else "artifact"
                ),
                item_type=record.type if record is not None else "artifact",
                metadata={"item_id": item_id},
            ))
            graph.add_edge(data_id, tool_id)
        for item_id in run.output_ids:
            record = records.get(item_id)
            data_id = "item:%s" % item_id
            graph.add_node(WorkflowNode(
                id=data_id,
                kind="data",
                label=(record.name if record is not None else item_id),
                icon=icon_for_item_type(
                    record.type if record is not None else "artifact"
                ),
                item_type=record.type if record is not None else "artifact",
                metadata={"item_id": item_id},
            ))
            graph.add_edge(tool_id, data_id)
        if not run.output_ids:
            for index_artifact, artifact in enumerate(run.artifacts):
                path = str(artifact)
                data_id = "artifact:%s:%d" % (run.id, index_artifact)
                graph.add_node(WorkflowNode(
                    id=data_id,
                    kind="data",
                    label=_label_from_path(path),
                    icon="file/sequence.svg",
                    item_type="artifact",
                    metadata={"path": path, "run_id": run.id},
                ))
                graph.add_edge(tool_id, data_id)

    if not graph.nodes and manager is not None:
        for item in manager.items.values():
            graph.add_node(WorkflowNode(
                id="item:%s" % item.id,
                kind="data",
                label=getattr(item, "loci_name", "")
                      or getattr(item, "name", "") or item.id,
                icon=icon_for_item_type(getattr(item, "item_type", "")),
                item_type=getattr(item, "item_type", ""),
            ))
    return graph


def graph_from_plan(plan: WorkflowPlan | Mapping[str, Any]) -> WorkflowGraph:
    """Build a preview graph from a compiled placeholder workflow."""
    if not isinstance(plan, WorkflowPlan):
        plan = WorkflowPlan.from_dict(plan)
    graph = WorkflowGraph(title=plan.name or "Workflow")
    input_nodes: Dict[str, str] = {}
    for name, info in plan.inputs.items():
        placeholder = info.get("placeholder", "")
        node_id = "input:%s" % placeholder
        input_nodes[placeholder] = node_id
        graph.add_node(WorkflowNode(
            id=node_id,
            kind="data",
            label=name,
            icon="file/sequence.svg",
            item_type="input",
            metadata={"placeholder": placeholder,
                      "original": info.get("original", "")},
        ))
    artifact_nodes: Dict[str, str] = {}
    for placeholder, info in plan.artifacts.items():
        node_id = "artifact:%s" % placeholder
        artifact_nodes[placeholder] = node_id
        graph.add_node(WorkflowNode(
            id=node_id,
            kind="data",
            label="%s.%s" % (info.get("role", "out"),
                              info.get("ext", "").lstrip(".")),
            icon=icon_for_item_type(
                "model" if "model" in placeholder.lower()
                else "phylogeny" if "tree" in placeholder.lower()
                else "sequence"
            ),
            item_type="artifact",
            metadata={"placeholder": placeholder,
                      "step": info.get("step", ""),
                      "original": info.get("original", "")},
        ))
    step_nodes: Dict[str, str] = {}
    for step in plan.steps:
        tool_id = "step:%s" % step["id"]
        step_nodes[step["id"]] = tool_id
        params = step.get("params", {}) or {}
        own_outputs = set(_iter_placeholders(step.get("outputs", {})))
        step_inputs = {
            key: value for key, value in params.items()
            if any(True for _ in _iter_placeholders(value))
            and not all(placeholder in own_outputs
                        for placeholder in _iter_placeholders(value))
        }
        graph.add_node(WorkflowNode(
            id=tool_id,
            kind="tool",
            label=step.get("op", step["id"]),
            icon=icon_for_tool(step.get("op", "")),
            metadata={"step_id": step["id"],
                      "op": step.get("op", ""),
                      "inputs": step_inputs,
                      "outputs": dict(step.get("outputs", {})),
                      "params": dict(params),
                      "depends_on": list(step.get("depends_on", []))},
        ))
    for step in plan.steps:
        tool_id = step_nodes[step["id"]]
        own_outputs = set(_iter_placeholders(step.get("outputs", {})))
        for placeholder in _iter_placeholders(step.get("params", {})):
            if placeholder in own_outputs:
                continue
            if placeholder in artifact_nodes:
                graph.add_edge(artifact_nodes[placeholder], tool_id)
            elif placeholder in input_nodes:
                graph.add_edge(input_nodes[placeholder], tool_id)
        for placeholder in own_outputs:
            if placeholder in artifact_nodes:
                graph.add_edge(tool_id, artifact_nodes[placeholder])
    return graph


def _iter_placeholders(value: Any) -> Iterable[str]:
    if isinstance(value, Mapping):
        for item in value.values():
            yield from _iter_placeholders(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _iter_placeholders(item)
    elif isinstance(value, str):
        for match in PLACEHOLDER_RE.finditer(value):
            yield match.group(0)


def layered_layout(graph: WorkflowGraph, x_spacing: int = 240,
                   y_spacing: int = 150) -> Dict[str, Tuple[float, float]]:
    """Assign nodes to alternating data/tool layers with a simple top-down order."""
    nodes = graph.node_map()
    incoming: Dict[str, List[str]] = {node_id: [] for node_id in nodes}
    outgoing: Dict[str, List[str]] = {node_id: [] for node_id in nodes}
    for edge in graph.edges:
        if edge.source in nodes and edge.target in nodes:
            outgoing[edge.source].append(edge.target)
            incoming[edge.target].append(edge.source)

    indegree = {node_id: len(incoming[node_id]) for node_id in nodes}
    ready = [node_id for node_id, degree in indegree.items() if degree == 0]
    order: List[str] = []
    while ready:
        current = ready.pop(0)
        order.append(current)
        for child in outgoing[current]:
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
    if len(order) != len(nodes):
        order = list(nodes)

    layer: Dict[str, int] = {}
    for node_id in order:
        predecessors = [layer[parent] for parent in incoming[node_id]
                        if parent in layer]
        layer[node_id] = (max(predecessors) + 1) if predecessors else 0

    grouped: Dict[int, List[str]] = {}
    for node_id in order:
        grouped.setdefault(layer[node_id], []).append(node_id)

    positions: Dict[str, Tuple[float, float]] = {}
    max_rows = max((len(rows) for rows in grouped.values()), default=1)
    canvas_height = (max_rows - 1) * y_spacing
    for layer_index, row in grouped.items():
        offset = (max_rows - len(row)) * y_spacing / 2.0
        for index, node_id in enumerate(row):
            positions[node_id] = (
                90 + layer_index * x_spacing,
                canvas_height / 2.0 + offset + index * y_spacing,
            )
    return positions
