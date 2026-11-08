"""Compile a command sequence into a portable placeholder workflow.

This is the MEGA-CC-like prototype: every declared output file gets a
symbolic placeholder, and any later parameter that names the same path is
rewritten to that placeholder. The resulting JSON describes both the commands
and their artifact graph, so it can be stored, reviewed and replayed.

Example::

    plan = compile_workflow({
        "inputs": {"raw": "data/raw.fas"},
        "steps": [
            {"op": "mafft.run",
             "params": {"input": "data/raw.fas", "output": "runs/mafft/aln.fas"}},
            {"op": "raxmlng.run",
             "params": {"msa": "runs/mafft/aln.fas", "output": "runs/raxml/tree.nwk"}},
        ],
    })
    plan.to_dict()["steps"][1]["params"]["msa"]
    # "%s01.output.fas%"

The compiler is deliberately stdlib-only and has no Qt dependency. A step is
a plain dict; ``outputs`` may be explicit::

    {"op": "mafft.run", "outputs": {"alignment": "aln.fas"}, "params": {...}}

When ``outputs`` is absent, keys named ``output``/``out``/``output_file``/
``output_path``/``output_dir``/``prefix`` (and the ``-o``/``--output``/
``--prefix`` forms) are used as output declarations. Explicit declarations are
authoritative when both forms are present.

Before compilation, read-only discovery steps (menu, skill, workspace and
plugin inspection calls) and steps whose recorded response is an error are
pruned from the transcript; ``prune_non_executable=False`` preserves them for
auditing.
"""

from __future__ import annotations

import copy
import json
import os
import re
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

WORKFLOW_SCHEMA = "yrmpe.workflow.v1"
PLACEHOLDER_RE = re.compile(r"%([A-Za-z0-9_.-]+)%")
OUTPUT_KEYS = {
    "output", "out", "output_file", "output_path", "output_dir",
    "prefix", "-o", "--output", "--prefix", "-prefix",
}

# Discovery and inspection actions are useful in an agent transcript but do
# not change workspace state, so they should not become DAG nodes.
READ_ONLY_OPS = {
    "workbench.info",
    "workbench.menu",
    "workbench.workspace",
    "workbench.sidebar",
    "workbench.skills",
    "workbench.help_doc",
    "workbench.workflow_graph",
    "workbench.preview_workflow",
    "workbench.show_workflow_preview",
    "workbench.list_windows",
    "workbench.jobs",
    "workbench.job_status",
    "workbench.inspect_job",
    "plugin.help_doc",
    "plugin.inspect_actions",
    "plugin.inspect_status",
    "plugin.inspect_controls",
    "plugin.inspect_console",
    "plugin.preview_results",
}
READ_ONLY_PREFIXES = (
    "workbench.inspect.",
    "workbench.show_",
    "workbench.preview_",
)
FAILED_STATUSES = {"error", "failed", "failure", "cancelled", "canceled"}


class WorkflowCompileError(ValueError):
    """Raised when a workflow cannot be compiled into a valid DAG."""

    def __init__(self, problems: Iterable[str]):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems) or "workflow compile failed")


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _sanitize_token(value: str, fallback: str = "x") -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or ""))
    token = token.strip("._-").replace("..", "_")
    return token or fallback


def _extension(path: str) -> str:
    suffix = os.path.splitext(os.path.basename(str(path).rstrip("/\\")))[1]
    if suffix.startswith(".") and re.fullmatch(r"\.[A-Za-z0-9]+", suffix):
        return suffix
    return ""


def _normalize_path(value: Any) -> str:
    if not isinstance(value, str):
        return ""
    text = value.strip().strip("\"'")
    if not text:
        return ""
    text = os.path.expandvars(os.path.expanduser(text))
    return os.path.normcase(os.path.normpath(text))


def _looks_like_path(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    text = value.strip()
    if not text or text.startswith("-"):
        return False
    if any(separator in text for separator in ("/", "\\")):
        return True
    if text.startswith("."):
        return True
    return bool(_extension(text))


def _infer_outputs(params: Mapping[str, Any]) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    for key, value in params.items():
        normalized_key = str(key).strip().lower()
        if normalized_key not in OUTPUT_KEYS:
            continue
        if isinstance(value, str) and value.strip():
            outputs[_sanitize_token(str(key), "out")] = value
        elif isinstance(value, (list, tuple)):
            for index, item in enumerate(value):
                if isinstance(item, str) and item.strip():
                    outputs["%s_%d" % (_sanitize_token(str(key), "out"), index)] = item
    return outputs


def _rewrite_value(
    value: Any,
    output_index: Mapping[str, Dict[str, Any]],
    input_index: Mapping[str, str],
    dependencies: set,
) -> Any:
    if isinstance(value, dict):
        return {
            key: _rewrite_value(item, output_index, input_index, dependencies)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [
            _rewrite_value(item, output_index, input_index, dependencies)
            for item in value
        ]
    if isinstance(value, tuple):
        return [
            _rewrite_value(item, output_index, input_index, dependencies)
            for item in value
        ]
    if not isinstance(value, str):
        return value
    normalized = _normalize_path(value)
    if normalized and normalized in output_index:
        artifact = output_index[normalized]
        dependencies.add(artifact["step"])
        return artifact["placeholder"]
    if normalized and normalized in input_index:
        return input_index[normalized]
    return value


def _is_read_only_step(step: Mapping[str, Any]) -> bool:
    op = str(step.get("op", "")).strip().lower()
    if op in READ_ONLY_OPS:
        return True
    return any(op.startswith(prefix) for prefix in READ_ONLY_PREFIXES)


def _payload_failed(payload: Any) -> bool:
    if not isinstance(payload, Mapping):
        return False
    for key in ("error", "error_info"):
        if payload.get(key):
            return True
    for key in ("ok", "success"):
        if key in payload and payload[key] is False:
            return True
    status = str(payload.get("status") or payload.get("state") or "").lower()
    if status in FAILED_STATUSES:
        return True
    for key in ("result", "response"):
        if _payload_failed(payload.get(key)):
            return True
    return False


def _step_failed(step: Mapping[str, Any]) -> bool:
    if _payload_failed(step):
        return True
    for key in ("result", "response"):
        if _payload_failed(step.get(key)):
            return True
    return False


def _prune_reason(step: Mapping[str, Any]) -> str:
    if _is_read_only_step(step):
        return "read-only"
    if _step_failed(step):
        return "failed"
    return ""


def _topological_order(steps: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_id = {step["id"]: step for step in steps}
    indegree = {step["id"]: len(step["depends_on"]) for step in steps}
    children: Dict[str, List[str]] = {step["id"]: [] for step in steps}
    for step in steps:
        for dependency in step["depends_on"]:
            children.setdefault(dependency, []).append(step["id"])
    ready = [step["id"] for step in steps if indegree[step["id"]] == 0]
    ordered: List[Dict[str, Any]] = []
    while ready:
        current = ready.pop(0)
        ordered.append(by_id[current])
        for child in children.get(current, []):
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
    if len(ordered) != len(steps):
        remaining = [step_id for step_id, value in indegree.items() if value > 0]
        raise WorkflowCompileError(
            ["workflow contains a cycle among steps: %s" % ", ".join(remaining)]
        )
    return ordered


@dataclass
class WorkflowPlan:
    """A compiled workflow with symbolic input/output slots."""

    steps: List[Dict[str, Any]]
    inputs: Dict[str, Dict[str, str]] = field(default_factory=dict)
    artifacts: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    id: str = ""
    name: str = ""
    schema: str = WORKFLOW_SCHEMA
    created_at: str = field(default_factory=_now)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "id": self.id,
            "name": self.name,
            "created_at": self.created_at,
            "inputs": copy.deepcopy(self.inputs),
            "steps": copy.deepcopy(self.steps),
            "artifacts": copy.deepcopy(self.artifacts),
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "WorkflowPlan":
        if data.get("schema") != WORKFLOW_SCHEMA:
            raise WorkflowCompileError(
                ["unsupported workflow schema %r" % data.get("schema")]
            )
        plan = cls(
            steps=[dict(step) for step in data.get("steps", [])],
            inputs=copy.deepcopy(data.get("inputs", {})),
            artifacts=copy.deepcopy(data.get("artifacts", {})),
            warnings=list(data.get("warnings", [])),
            id=str(data.get("id", "")),
            name=str(data.get("name", "")),
            created_at=str(data.get("created_at", _now())),
        )
        plan.validate()
        return plan

    def validate(self) -> None:
        problems: List[str] = []
        if not self.steps:
            problems.append("workflow has no steps")
        seen_ids = set()
        for step in self.steps:
            step_id = step.get("id", "")
            if not step_id:
                problems.append("step has no id")
                continue
            if step_id in seen_ids:
                problems.append("duplicate step id %r" % step_id)
            seen_ids.add(step_id)
            if not step.get("op"):
                problems.append("step %s has no op" % step_id)
        known_placeholders = set(self.artifacts)
        known_placeholders.update(
            value.get("placeholder", "") for value in self.inputs.values()
        )
        for step in self.steps:
            for dependency in step.get("depends_on", []):
                if dependency not in seen_ids:
                    problems.append(
                        "step %s depends on unknown step %r"
                        % (step.get("id", ""), dependency)
                    )
            for placeholder in _placeholders_in(step.get("params", {})):
                if placeholder not in known_placeholders:
                    problems.append(
                        "step %s uses unknown placeholder %s"
                        % (step.get("id", ""), placeholder)
                    )
            for placeholder in _placeholders_in(step.get("outputs", {})):
                if placeholder not in known_placeholders:
                    problems.append(
                        "step %s declares unknown placeholder %s"
                        % (step.get("id", ""), placeholder)
                    )
        for placeholder, artifact in self.artifacts.items():
            if not PLACEHOLDER_RE.fullmatch(placeholder):
                problems.append("invalid placeholder %r" % placeholder)
                continue
            if artifact.get("step") not in seen_ids:
                problems.append(
                    "artifact %s names unknown step %r"
                    % (placeholder, artifact.get("step"))
                )
        if problems:
            raise WorkflowCompileError(problems)

    def bindings(self, output_base: Optional[str] = None) -> Dict[str, str]:
        """Map every input/artifact placeholder to a concrete path.

        ``output_base`` rewrites artifact paths under a run directory; inputs
        keep their declared original path.
        """
        bindings: Dict[str, str] = {}
        for value in self.inputs.values():
            bindings[value["placeholder"]] = value["original"]
        for placeholder, artifact in self.artifacts.items():
            if output_base:
                filename = artifact["role"] + artifact.get("ext", "")
                bindings[placeholder] = os.path.join(
                    os.path.abspath(output_base),
                    artifact["step"],
                    filename,
                )
            else:
                bindings[placeholder] = artifact["original"]
        return bindings

    def resolve(self, value: Any, bindings: Optional[Mapping[str, str]] = None) -> Any:
        """Substitute placeholders in a nested step parameter value."""
        bindings = self.bindings() if bindings is None else bindings

        def replace(text: str) -> str:
            def substitute(match: re.Match) -> str:
                key = match.group(0)
                if key not in bindings:
                    raise KeyError("unbound workflow placeholder %s" % key)
                return bindings[key]

            return PLACEHOLDER_RE.sub(substitute, text)

        if isinstance(value, dict):
            return {key: self.resolve(item, bindings) for key, item in value.items()}
        if isinstance(value, list):
            return [self.resolve(item, bindings) for item in value]
        if isinstance(value, tuple):
            return [self.resolve(item, bindings) for item in value]
        return replace(value) if isinstance(value, str) else value

    def resolved_steps(self, bindings: Optional[Mapping[str, str]] = None):
        """Yield steps in dependency order with placeholders substituted."""
        for step in _topological_order(self.steps):
            resolved = dict(step)
            resolved["params"] = self.resolve(step.get("params", {}), bindings)
            resolved["outputs"] = self.resolve(
                step.get("outputs", {}), bindings
            )
            yield resolved


def compile_workflow(
    spec: Mapping[str, Any] | Sequence[Mapping[str, Any]],
    *,
    name: str = "",
    allow_undeclared_inputs: bool = True,
    prune_non_executable: bool = True,
) -> WorkflowPlan:
    """Compile a workflow spec into a placeholder graph.

    ``spec`` is either ``{"inputs": ..., "steps": [...]}`` or a step sequence.
    Steps may be given out of order; dependencies are derived by matching
    declared output paths against later input values.
    """
    if isinstance(spec, Mapping):
        raw_steps = list(spec.get("steps", []))
        raw_inputs = spec.get("inputs", {})
        workflow_name = name or str(spec.get("name", ""))
        workflow_id = str(spec.get("id", ""))
    else:
        raw_steps = list(spec)
        raw_inputs = {}
        workflow_name = name
        workflow_id = ""

    warnings: List[str] = []
    if prune_non_executable:
        executable_steps = []
        for index, raw_step in enumerate(raw_steps):
            if not isinstance(raw_step, Mapping):
                raise WorkflowCompileError(
                    ["step %d must be an object" % (index + 1)]
                )
            reason = _prune_reason(raw_step)
            if reason:
                warnings.append(
                    "pruned %s step %s (%s)"
                    % (reason, raw_step.get("id", index + 1),
                       raw_step.get("op", "unknown"))
                )
            else:
                executable_steps.append(raw_step)
        raw_steps = executable_steps
    if not raw_steps:
        raise WorkflowCompileError([
            "workflow has no executable steps after pruning read-only or "
            "failed operations"
        ])
    if not workflow_id:
        workflow_id = "wf_%s_%s" % (
            datetime.now().strftime("%Y%m%d_%H%M%S"),
            uuid.uuid4().hex[:6],
        )

    input_index: Dict[str, str] = {}
    inputs: Dict[str, Dict[str, str]] = {}
    if isinstance(raw_inputs, Mapping):
        raw_input_items = list(raw_inputs.items())
    else:
        raw_input_items = [
            ("input_%d" % (index + 1), value)
            for index, value in enumerate(raw_inputs)
        ]
    for input_name, input_path in raw_input_items:
        if not isinstance(input_path, str) or not input_path.strip():
            continue
        role = _sanitize_token(str(input_name), "input")
        ext = _extension(input_path)
        if ext and not role.lower().endswith(ext.lower()):
            role += ext
        placeholder = "%%input.%s%%" % role
        input_index[_normalize_path(input_path)] = placeholder
        inputs[str(input_name)] = {
            "placeholder": placeholder,
            "original": input_path,
        }

    output_index: Dict[str, Dict[str, Any]] = {}
    artifacts: Dict[str, Dict[str, Any]] = {}
    prepared: List[Dict[str, Any]] = []
    used_ids = set()

    for index, raw_step in enumerate(raw_steps):
        step_id = _sanitize_token(
            str(raw_step.get("id", "")), "s%02d" % (index + 1)
        )
        if step_id in used_ids:
            raise WorkflowCompileError(["duplicate step id %r" % step_id])
        used_ids.add(step_id)
        params = copy.deepcopy(raw_step.get("params", {}))
        declared_outputs = raw_step.get("outputs")
        if declared_outputs is None:
            declared_outputs = _infer_outputs(params)
        if not isinstance(declared_outputs, Mapping):
            raise WorkflowCompileError(
                ["step %s outputs must be an object" % step_id]
            )
        step_artifacts: Dict[str, str] = {}
        for role, output_path in declared_outputs.items():
            if not isinstance(output_path, str) or not output_path.strip():
                continue
            normalized = _normalize_path(output_path)
            if not normalized:
                continue
            if normalized in output_index:
                raise WorkflowCompileError(
                    [
                        "duplicate output path %r (steps %s and %s)"
                        % (
                            output_path,
                            output_index[normalized]["step"],
                            step_id,
                        )
                    ]
                )
            role_token = _sanitize_token(str(role), "out")
            ext = _extension(output_path)
            placeholder = "%%%s.%s%s%%" % (
                step_id, role_token, ext
            )
            artifact = {
                "step": step_id,
                "role": role_token,
                "ext": ext,
                "placeholder": placeholder,
                "original": output_path,
            }
            output_index[normalized] = artifact
            artifacts[placeholder] = artifact
            step_artifacts[str(role)] = placeholder
        prepared.append({
            "id": step_id,
            "op": str(raw_step.get("op", "")),
            "params": params,
            "outputs": step_artifacts,
            "depends_on": [],
        })

    for step in prepared:
        dependencies = set()
        step["params"] = _rewrite_value(
            step["params"], output_index, input_index, dependencies
        )
        # A step may name its own output in ``params``; that is not a
        # dependency edge back to itself.
        dependencies.discard(step["id"])
        step["depends_on"] = sorted(dependencies)
        if not step["op"]:
            raise WorkflowCompileError(
                ["step %s has no op" % step["id"]]
            )
        for path in _unresolved_paths(step["params"]):
            if _normalize_path(path) in input_index:
                continue
            message = (
                "external path left literal in step %s: %s"
                % (step["id"], path)
            )
            if allow_undeclared_inputs:
                warnings.append(message)
            else:
                raise WorkflowCompileError([message])

    plan = WorkflowPlan(
        steps=prepared,
        inputs=inputs,
        artifacts=artifacts,
        warnings=warnings,
        id=workflow_id,
        name=workflow_name or "workflow",
    )
    plan.steps = _topological_order(plan.steps)
    plan.validate()
    return plan


def compile_jsonl(source: str | os.PathLike | Iterable[str], **kwargs) -> WorkflowPlan:
    """Compile one JSON command per line into a workflow plan.

    Each line may be a bridge-style command (``op``/``params``) with optional
    ``id`` and ``outputs``. Comments beginning with ``#`` and blank lines are
    ignored. ``source`` may be a path or an iterable of lines.
    """
    if isinstance(source, (str, os.PathLike)) and os.path.exists(source):
        with open(source, "r", encoding="utf-8") as handle:
            lines = handle.readlines()
    elif isinstance(source, (str, os.PathLike)):
        lines = str(source).splitlines()
    else:
        lines = list(source)
    steps = []
    header: Dict[str, Any] = {}
    for line_number, line in enumerate(lines, 1):
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:
            raise WorkflowCompileError(
                ["line %d is not JSON: %s" % (line_number, exc)]
            )
        if not isinstance(payload, dict):
            raise WorkflowCompileError(
                ["line %d must be a JSON object" % line_number]
            )
        if (payload.get("kind") == "workflow"
                or "steps" in payload
                or (payload.get("schema") == WORKFLOW_SCHEMA
                    and not payload.get("op"))):
            header = payload
            continue
        steps.append(payload)
    spec: Dict[str, Any] = {"steps": steps}
    if header:
        spec["name"] = header.get("name", "")
        spec["inputs"] = header.get("inputs", {})
    return compile_workflow(spec, **kwargs)


def _placeholders_in(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for item in value.values():
            yield from _placeholders_in(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _placeholders_in(item)
    elif isinstance(value, str):
        for match in PLACEHOLDER_RE.finditer(value):
            yield match.group(0)


def _unresolved_paths(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for item in value.values():
            yield from _unresolved_paths(item)
    elif isinstance(value, list):
        for item in value:
            yield from _unresolved_paths(item)
    elif isinstance(value, str):
        if PLACEHOLDER_RE.search(value):
            return
        if _looks_like_path(value):
            yield value


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Small CLI: compile a JSON or JSONL workflow spec to portable JSON."""
    import argparse
    import sys

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", help="workflow JSON/JSONL file, or '-' for stdin")
    parser.add_argument("--jsonl", action="store_true",
                        help="treat the input as one JSON command per line")
    parser.add_argument("--output", "-o", default="",
                        help="write compiled JSON here instead of stdout")
    parser.add_argument("--resolve", default="",
                        help="bind output placeholders under this base directory "
                             "and print resolved steps")
    args = parser.parse_args(argv)

    try:
        if args.spec == "-":
            text = sys.stdin.read()
        else:
            with open(args.spec, "r", encoding="utf-8") as handle:
                text = handle.read()
        is_jsonl = args.jsonl or (
            args.spec != "-" and args.spec.lower().endswith(".jsonl")
        )
        if is_jsonl:
            plan = compile_jsonl(text)
        else:
            plan = compile_workflow(json.loads(text))
        if args.resolve:
            payload = {
                "workflow": plan.to_dict(),
                "bindings": plan.bindings(args.resolve),
                "resolved_steps": list(
                    plan.resolved_steps(plan.bindings(args.resolve))
                ),
            }
        else:
            payload = plan.to_dict()
    except (OSError, json.JSONDecodeError, WorkflowCompileError) as exc:
        print("workflow compile failed: %s" % exc, file=sys.stderr)
        return 2

    rendered = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(rendered)
    else:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
