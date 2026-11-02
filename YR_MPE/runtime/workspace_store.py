"""Single source of truth for a YR-MPE workspace.

Layout::

    <workspace>/
      workspace.json
      imports/<import_id>/{import.json, files/...}
      runs/<run_id>/{manifest.json, console.jsonl, outputs/...}
      datasets/<dataset_id>/items/<item_id>/...
      recycled/

Nodes are data (items); edges are analysis methods (runs). A run can have
several inputs and several outputs. Imports are runs with no inputs and a
dedicated manifest, so every imported item can be traced to an exact source
path, hash, importer and parameter set.

The module is stdlib-only and has no PyQt dependency.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

WORKSPACE_SCHEMA = "yrmpe.workspace.v1"
RUN_SCHEMA = "yrmpe.run.v1"
IMPORT_SCHEMA = "yrmpe.import.v1"

_EXTENSIONS = {
    "fasta": ".fna",
    "phylip": ".phy",
    "nexus": ".nex",
    "newick": ".nwk",
    "json": ".json",
    "tsv": ".tsv",
}


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _short_id(prefix: str) -> str:
    return "%s_%s_%s" % (
        prefix,
        datetime.now().strftime("%Y%m%d_%H%M%S"),
        uuid.uuid4().hex[:8],
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _atomic_write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(_jsonable(payload), handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temporary, path)


@dataclass
class ItemRecord:
    """A data node."""

    id: str
    type: str
    name: str
    dataset_id: str
    format: str = ""
    file: str = ""
    length: int = 0
    sequence_count: int = 0
    is_aligned: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)
    data: Dict[str, Any] = field(default_factory=dict)
    import_id: Optional[str] = None
    run_id: Optional[str] = None
    parent_id: Optional[str] = None
    created_at: str = field(default_factory=_now)
    modified_at: str = field(default_factory=_now)

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(self.__dict__)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ItemRecord":
        return cls(**{key: data[key] for key in cls.__dataclass_fields__ if key in data})


@dataclass
class DatasetRecord:
    id: str
    name: str
    description: str = ""
    source: str = ""
    item_ids: List[str] = field(default_factory=list)
    stream_ids: List[str] = field(default_factory=list)
    settings: Dict[str, Any] = field(default_factory=dict)
    created_at: str = field(default_factory=_now)

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(self.__dict__)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DatasetRecord":
        return cls(**{key: data[key] for key in cls.__dataclass_fields__ if key in data})


@dataclass
class RunRecord:
    """An analysis method edge: inputs -> outputs."""

    id: str
    name: str
    method: str
    tool: str = ""
    tool_version: str = ""
    command: List[str] = field(default_factory=list)
    parameters: Dict[str, Any] = field(default_factory=dict)
    parameter_schema: Dict[str, Any] = field(default_factory=dict)
    input_ids: List[str] = field(default_factory=list)
    output_ids: List[str] = field(default_factory=list)
    status: str = "running"
    created_at: str = field(default_factory=_now)
    completed_at: str = ""
    cwd: str = ""
    log: str = ""
    artifacts: List[str] = field(default_factory=list)
    stream_id: Optional[str] = None
    manifest: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(self.__dict__)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RunRecord":
        return cls(**{key: data[key] for key in cls.__dataclass_fields__ if key in data})


@dataclass
class StreamRecord:
    """An ordered chain of runs inside one dataset."""

    id: str
    name: str
    dataset_id: str
    run_ids: List[str] = field(default_factory=list)
    created_at: str = field(default_factory=_now)

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(self.__dict__)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StreamRecord":
        return cls(**{key: data[key] for key in cls.__dataclass_fields__ if key in data})


@dataclass
class ImportRecord:
    """A precise import event (a run with no inputs)."""

    id: str
    source_paths: List[str] = field(default_factory=list)
    source_hashes: Dict[str, str] = field(default_factory=dict)
    importer: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)
    dataset_id: Optional[str] = None
    item_ids: List[str] = field(default_factory=list)
    created_at: str = field(default_factory=_now)
    manifest: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(self.__dict__)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ImportRecord":
        return cls(**{key: data[key] for key in cls.__dataclass_fields__ if key in data})


@dataclass
class WorkspaceIndex:
    schema: str = WORKSPACE_SCHEMA
    workspace_id: str = field(default_factory=lambda: uuid.uuid4().hex)
    name: str = ""
    created_at: str = field(default_factory=_now)
    updated_at: str = field(default_factory=_now)
    datasets: Dict[str, DatasetRecord] = field(default_factory=dict)
    items: Dict[str, ItemRecord] = field(default_factory=dict)
    runs: Dict[str, RunRecord] = field(default_factory=dict)
    streams: Dict[str, StreamRecord] = field(default_factory=dict)
    imports: Dict[str, ImportRecord] = field(default_factory=dict)
    settings: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "workspace_id": self.workspace_id,
            "name": self.name,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "datasets": {key: value.to_dict() for key, value in self.datasets.items()},
            "items": {key: value.to_dict() for key, value in self.items.items()},
            "runs": {key: value.to_dict() for key, value in self.runs.items()},
            "streams": {key: value.to_dict() for key, value in self.streams.items()},
            "imports": {key: value.to_dict() for key, value in self.imports.items()},
            "settings": _jsonable(self.settings),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkspaceIndex":
        index = cls(
            schema=data.get("schema", WORKSPACE_SCHEMA),
            workspace_id=data.get("workspace_id", uuid.uuid4().hex),
            name=data.get("name", ""),
            created_at=data.get("created_at", _now()),
            updated_at=data.get("updated_at", _now()),
            settings=data.get("settings", {}) or {},
        )
        index.datasets = {
            key: DatasetRecord.from_dict(value)
            for key, value in data.get("datasets", {}).items()
        }
        index.items = {
            key: ItemRecord.from_dict(value)
            for key, value in data.get("items", {}).items()
        }
        index.runs = {
            key: RunRecord.from_dict(value)
            for key, value in data.get("runs", {}).items()
        }
        index.streams = {
            key: StreamRecord.from_dict(value)
            for key, value in data.get("streams", {}).items()
        }
        index.imports = {
            key: ImportRecord.from_dict(value)
            for key, value in data.get("imports", {}).items()
        }
        return index


class WorkspaceStore:
    """Read/write a workspace directory."""

    def __init__(self, root: str, name: str = ""):
        self.root = Path(root).expanduser().resolve()
        self.index_path = self.root / "workspace.json"
        self.ensure_layout()
        existed = self.index_path.exists()
        self.index = self.load()
        if not existed:
            self.index.name = name or self.root.name

    # ------------------------------------------------------------------
    # paths / layout
    # ------------------------------------------------------------------
    def ensure_layout(self) -> None:
        for name in ("imports", "runs", "datasets", "recycled"):
            (self.root / name).mkdir(parents=True, exist_ok=True)

    def path(self, relative: str) -> Path:
        path = Path(relative)
        return path if path.is_absolute() else self.root / path

    def relative(self, path: Path) -> str:
        path = Path(path).resolve()
        try:
            return str(path.relative_to(self.root))
        except ValueError:
            return str(path)

    def _dataset_dir(self, dataset_id: str) -> Path:
        return self.root / "datasets" / dataset_id

    def _item_dir(self, dataset_id: str, item_id: str) -> Path:
        return self._dataset_dir(dataset_id) / "items" / item_id

    def item_path(self, dataset_id: str, item_id: str, filename: str = "data.fna") -> Path:
        path = self._item_dir(dataset_id, item_id) / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def _run_dir(self, run_id: str) -> Path:
        return self.root / "runs" / run_id

    def _import_dir(self, import_id: str) -> Path:
        return self.root / "imports" / import_id

    # ------------------------------------------------------------------
    # persistence
    # ------------------------------------------------------------------
    def load(self) -> WorkspaceIndex:
        if self.index_path.exists():
            with self.index_path.open("r", encoding="utf-8") as handle:
                return WorkspaceIndex.from_dict(json.load(handle))
        return WorkspaceIndex()

    def save(self) -> None:
        self.index.updated_at = _now()
        _atomic_write_json(self.index_path, self.index.to_dict())

    # ------------------------------------------------------------------
    # datasets / items
    # ------------------------------------------------------------------
    def add_dataset(
        self,
        name: str,
        description: str = "",
        source: str = "",
        settings: Optional[Dict[str, Any]] = None,
        dataset_id: Optional[str] = None,
    ) -> DatasetRecord:
        dataset_id = dataset_id or _short_id("ds")
        record = DatasetRecord(
            id=dataset_id,
            name=name,
            description=description,
            source=source,
            settings=settings or {},
        )
        self.index.datasets[dataset_id] = record
        self._dataset_dir(dataset_id).mkdir(parents=True, exist_ok=True)
        return record

    def add_item(
        self,
        dataset_id: str,
        item_type: str,
        name: str,
        source_path: Optional[str] = None,
        copy: bool = False,
        file_format: str = "",
        metadata: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        import_id: Optional[str] = None,
        run_id: Optional[str] = None,
        parent_id: Optional[str] = None,
        length: int = 0,
        sequence_count: int = 0,
        is_aligned: bool = False,
        item_id: Optional[str] = None,
    ) -> ItemRecord:
        if dataset_id not in self.index.datasets:
            raise KeyError("unknown dataset %r" % dataset_id)

        item_id = item_id or _short_id("item")
        item_dir = self._item_dir(dataset_id, item_id)
        item_dir.mkdir(parents=True, exist_ok=True)

        stored_file = ""
        item_metadata = dict(metadata or {})
        if source_path:
            source = Path(source_path)
            if source.exists():
                if source.is_dir():
                    target = item_dir / source.name
                    if copy:
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(source, target)
                        stored_file = self.relative(target)
                    else:
                        stored_file = str(source.resolve())
                else:
                    target = item_dir / source.name
                    if copy:
                        shutil.copy2(source, target)
                        stored_file = self.relative(target)
                    else:
                        stored_file = str(source.resolve())
                item_metadata.setdefault("source_path", str(source.resolve()))
                if source.is_file():
                    item_metadata.setdefault("sha256", _sha256(source))
            else:
                item_metadata["missing_source"] = str(source)

        record = ItemRecord(
            id=item_id,
            type=item_type,
            name=name,
            dataset_id=dataset_id,
            format=file_format,
            file=stored_file,
            length=length,
            sequence_count=sequence_count,
            is_aligned=is_aligned,
            metadata=item_metadata,
            data=data or {},
            import_id=import_id,
            run_id=run_id,
            parent_id=parent_id,
        )
        self.index.items[item_id] = record
        self.index.datasets[dataset_id].item_ids.append(item_id)
        if import_id and import_id in self.index.imports:
            import_record = self.index.imports[import_id]
            import_record.item_ids.append(item_id)
            self._write_import_manifest(import_record)
        return record

    def get_item(self, item_id: str) -> Optional[ItemRecord]:
        return self.index.items.get(item_id)

    def get_dataset(self, dataset_id: str) -> Optional[DatasetRecord]:
        return self.index.datasets.get(dataset_id)

    def items_in_dataset(self, dataset_id: str) -> List[ItemRecord]:
        dataset = self.index.datasets.get(dataset_id)
        if not dataset:
            return []
        return [self.index.items[item_id] for item_id in dataset.item_ids if item_id in self.index.items]

    def children_of(self, item_id: str) -> List[ItemRecord]:
        return [item for item in self.index.items.values() if item.parent_id == item_id]

    def parent_of(self, item_id: str) -> Optional[ItemRecord]:
        item = self.index.items.get(item_id)
        if not item or not item.parent_id:
            return None
        return self.index.items.get(item.parent_id)

    # ------------------------------------------------------------------
    # imports / runs / streams
    # ------------------------------------------------------------------
    def new_import(
        self,
        source_paths: Iterable[str],
        importer: str,
        parameters: Optional[Dict[str, Any]] = None,
        copy: bool = True,
        dataset_id: Optional[str] = None,
    ) -> ImportRecord:
        import_id = _short_id("import")
        import_dir = self._import_dir(import_id)
        files_dir = import_dir / "files"
        files_dir.mkdir(parents=True, exist_ok=True)

        record = ImportRecord(
            id=import_id,
            importer=importer,
            parameters=parameters or {},
            dataset_id=dataset_id,
            manifest=self.relative(import_dir / "import.json"),
        )

        for raw_path in source_paths:
            source = Path(raw_path)
            record.source_paths.append(str(source))
            if not source.exists():
                continue
            if source.is_file():
                record.source_hashes[str(source)] = _sha256(source)
                if copy:
                    shutil.copy2(source, files_dir / source.name)
            elif source.is_dir():
                if copy:
                    target = files_dir / source.name
                    if target.exists():
                        shutil.rmtree(target)
                    shutil.copytree(source, target)
                for child in sorted(source.rglob("*")):
                    if child.is_file():
                        record.source_hashes[str(child)] = _sha256(child)

        self.index.imports[import_id] = record
        self._write_import_manifest(record)
        return record

    def get_import(self, import_id: str) -> Optional[ImportRecord]:
        return self.index.imports.get(import_id)

    def new_run(
        self,
        method: str,
        tool: str,
        input_ids: Iterable[str],
        name: Optional[str] = None,
        tool_version: str = "",
        parameters: Optional[Dict[str, Any]] = None,
        parameter_schema: Optional[Dict[str, Any]] = None,
        command: Optional[Iterable[str]] = None,
        output_ids: Optional[Iterable[str]] = None,
        stream_id: Optional[str] = None,
        cwd: str = "",
    ) -> RunRecord:
        run_id = _short_id("run")
        run_dir = self._run_dir(run_id)
        (run_dir / "outputs").mkdir(parents=True, exist_ok=True)
        log_path = run_dir / "console.jsonl"
        log_path.touch(exist_ok=True)

        record = RunRecord(
            id=run_id,
            name=name or method,
            method=method,
            tool=tool,
            tool_version=tool_version,
            command=list(command or []),
            parameters=parameters or {},
            parameter_schema=parameter_schema or {},
            input_ids=list(input_ids),
            output_ids=list(output_ids or []),
            status="running",
            cwd=cwd,
            log=self.relative(log_path),
            stream_id=stream_id,
            manifest=self.relative(run_dir / "manifest.json"),
        )
        self.index.runs[run_id] = record
        if stream_id and stream_id in self.index.streams:
            self.index.streams[stream_id].run_ids.append(run_id)
        self._write_run_manifest(record)
        return record

    def finish_run(
        self,
        run_id: str,
        status: str = "completed",
        output_ids: Optional[Iterable[str]] = None,
        artifacts: Optional[Iterable[str]] = None,
    ) -> RunRecord:
        record = self.index.runs.get(run_id)
        if record is None:
            raise KeyError("unknown run %r" % run_id)
        record.status = status
        record.completed_at = _now()
        if output_ids is not None:
            record.output_ids = list(output_ids)
        if artifacts is not None:
            record.artifacts = [str(item) for item in artifacts]
        self._write_run_manifest(record)
        return record

    def get_run(self, run_id: str) -> Optional[RunRecord]:
        return self.index.runs.get(run_id)

    def rename_run(self, run_id: str, name: str) -> RunRecord:
        record = self.index.runs.get(run_id)
        if record is None:
            raise KeyError("unknown run %r" % run_id)
        record.name = name
        self._write_run_manifest(record)
        return record

    def add_stream(
        self,
        name: str,
        dataset_id: str,
        run_ids: Optional[Iterable[str]] = None,
    ) -> StreamRecord:
        if dataset_id not in self.index.datasets:
            raise KeyError("unknown dataset %r" % dataset_id)
        stream_id = _short_id("stream")
        record = StreamRecord(
            id=stream_id,
            name=name,
            dataset_id=dataset_id,
            run_ids=list(run_ids or []),
        )
        self.index.streams[stream_id] = record
        self.index.datasets[dataset_id].stream_ids.append(stream_id)
        for run_id in record.run_ids:
            if run_id in self.index.runs:
                self.index.runs[run_id].stream_id = stream_id
        return record

    def get_stream(self, stream_id: str) -> Optional[StreamRecord]:
        return self.index.streams.get(stream_id)

    def add_run_to_stream(self, stream_id: str, run_id: str) -> StreamRecord:
        stream = self.index.streams.get(stream_id)
        if stream is None:
            raise KeyError("unknown stream %r" % stream_id)
        if run_id not in stream.run_ids:
            stream.run_ids.append(run_id)
        if run_id in self.index.runs:
            self.index.runs[run_id].stream_id = stream_id
        return stream

    # ------------------------------------------------------------------
    # manifests
    # ------------------------------------------------------------------
    def _write_import_manifest(self, record: ImportRecord) -> None:
        payload = record.to_dict()
        payload["schema"] = IMPORT_SCHEMA
        payload["items"] = [
            {
                "id": item.id,
                "type": item.type,
                "name": item.name,
                "file": item.file,
            }
            for item_id in record.item_ids
            for item in [self.index.items.get(item_id)]
            if item is not None
        ]
        _atomic_write_json(self.path(record.manifest), payload)

    def _write_run_manifest(self, record: RunRecord) -> None:
        payload = record.to_dict()
        payload["schema"] = RUN_SCHEMA
        payload["inputs"] = record.input_ids
        payload["outputs"] = record.output_ids
        _atomic_write_json(self.path(record.manifest), payload)

    def refresh_manifests(self) -> None:
        """Rewrite manifests after item links change."""
        for record in self.index.imports.values():
            self._write_import_manifest(record)
        for record in self.index.runs.values():
            self._write_run_manifest(record)

    # ------------------------------------------------------------------
    # migration
    # ------------------------------------------------------------------
    def migrate_legacy_selection_state(self, path: str) -> int:
        """Import datasets/items from dataset_selection_state.json.

        Preserves UUIDs. Files are referenced, not copied, so migration is
        cheap; new imports should use copy=True.
        """
        source = Path(path)
        if not source.exists():
            return 0
        with source.open("r", encoding="utf-8") as handle:
            state = json.load(handle)

        datasets = state.get("datasets", {})
        items = state.get("items", {})
        migrated = 0

        for dataset_id, dataset_data in datasets.items():
            dataset = DatasetRecord.from_dict(dataset_data)
            dataset.id = dataset_id
            dataset.settings = dataset_data.get("settings", {}) or {}
            dataset.item_ids = []
            self.index.datasets[dataset_id] = dataset
            self._dataset_dir(dataset_id).mkdir(parents=True, exist_ok=True)

        def _import_item(item_id: str, item_data: Dict[str, Any], dataset_id: str):
            nonlocal migrated
            item_type = item_data.get("item_type", "alignment")
            record = ItemRecord(
                id=item_id,
                type=item_type,
                name=item_data.get("loci_name", "") or item_data.get("name", ""),
                dataset_id=dataset_id,
                format=item_data.get("format", ""),
                file=str(item_data.get("file_path", "") or ""),
                length=int(item_data.get("length", 0) or 0),
                sequence_count=int(item_data.get("sequence_count", 0) or 0),
                is_aligned=bool(item_data.get("is_aligned", False)),
                metadata=dict(item_data.get("metadata", {}) or {}),
                data=dict(item_data.get("data", {}) or {}),
                parent_id=item_data.get("parent_id"),
                created_at=item_data.get("created_at", _now()),
                modified_at=item_data.get("modified_at", _now()),
            )
            record.metadata.setdefault("migrated_from", source.name)
            self.index.items[item_id] = record
            if dataset_id in self.index.datasets:
                self.index.datasets[dataset_id].item_ids.append(item_id)
            migrated += 1

        for item_id, item_data in items.items():
            dataset_id = item_data.get("dataset_id", "")
            _import_item(item_id, item_data, dataset_id)

        # Older files only had dataset.settings['dataset_items'].
        for dataset_id, dataset in self.index.datasets.items():
            if dataset.item_ids:
                continue
            for item_data in dataset.settings.get("dataset_items", []):
                item_id = item_data.get("id") or _short_id("item")
                _import_item(item_id, item_data, dataset_id)

        self.save()
        backup = source.with_suffix(source.suffix + ".bak")
        if not backup.exists():
            os.replace(source, backup)
        return migrated
