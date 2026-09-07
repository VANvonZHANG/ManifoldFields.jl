"""Tests for the artifact-fetch helper (examples/python/_data.py).

Run: /home/zhangfan/anaconda3/envs/geo/bin/python -m pytest examples/python -q
"""
import tarfile
from pathlib import Path

import pytest

from _data import _read_entry, _sha256, ensure_data


def _make_tarball(tmp: Path, files: dict[str, bytes]) -> Path:
    root = tmp / "stage"
    root.mkdir(parents=True)
    for name, blob in files.items():
        (root / name).write_bytes(blob)
    tarball = tmp / "artifact.tar.gz"
    with tarfile.open(tarball, "w:gz") as tf:
        for name in files:
            tf.add(root / name, arcname=name)  # flat: files at archive root (see Task 5 layout note)
    return tarball


def _make_manifest(tmp: Path, tarball: Path) -> Path:
    manifest = tmp / "Artifacts.toml"
    manifest.write_text(
        "[example-data]\n"
        'git-tree-sha1 = "deadbeefdeadbeefdeadbeefdeadbeefdeadbeef"\n'
        "lazy = true\n\n"
        "[[example-data.download]]\n"
        f'url = "file://{tarball}"\n'
        f'sha256 = "{_sha256(tarball)}"\n'
    )
    return manifest


def test_read_entry(tmp_path):
    tarball = _make_tarball(tmp_path, {"a.nc": b"x"})
    url, sha = _read_entry(_make_manifest(tmp_path, tarball))
    assert url == f"file://{tarball}"
    assert len(sha) == 64


def test_ensure_data_downloads_verifies_extracts(tmp_path):
    tarball = _make_tarball(tmp_path, {"a.nc": b"hello", "b.nc": b"world"})
    root = ensure_data(_make_manifest(tmp_path, tarball), tmp_path / "dest")
    assert (root / "a.nc").read_bytes() == b"hello"
    assert (root / "b.nc").read_bytes() == b"world"


def test_ensure_data_idempotent_and_offline(tmp_path):
    tarball = _make_tarball(tmp_path, {"a.nc": b"hello"})
    manifest = _make_manifest(tmp_path, tarball)
    ensure_data(manifest, tmp_path / "dest")
    tarball.unlink()  # source gone: second run must not need it
    assert (ensure_data(manifest, tmp_path / "dest") / "a.nc").exists()


def test_bad_sha256_rejected(tmp_path):
    tarball = _make_tarball(tmp_path, {"a.nc": b"hello"})
    manifest = _make_manifest(tmp_path, tarball)
    manifest.write_text(manifest.read_text().replace(_sha256(tarball), "0" * 64))
    with pytest.raises(RuntimeError, match="sha256"):
        ensure_data(manifest, tmp_path / "dest")
