"""Fetch the shared example-data artifact for the marimo notebooks.

Reads examples/Artifacts.toml (the same manifest the Julia side resolves with
Pkg.ensure_artifact_installed), downloads the tarball, verifies its sha256,
extracts it under examples/data/ and returns the data root. Standard library
only — the Python track never depends on a Julia installation.
"""
from __future__ import annotations

import hashlib
import tarfile
import urllib.request
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

EXAMPLES_DIR = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = EXAMPLES_DIR / "Artifacts.toml"
DEFAULT_DEST = EXAMPLES_DIR / "data"


def _read_entry(manifest: Path, name: str = "example-data") -> tuple[str, str]:
    entry = tomllib.loads(manifest.read_text())[name]
    download = entry["download"][0]
    return download["url"], download["sha256"]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_data(
    manifest: Path | str = DEFAULT_MANIFEST,
    dest: Path | str = DEFAULT_DEST,
    name: str = "example-data",
) -> Path:
    """Return the directory holding the extracted example data files."""
    manifest, dest = Path(manifest), Path(dest)
    url, expected = _read_entry(manifest, name)
    marker = dest / f".{expected[:12]}.complete"
    if marker.exists():
        return dest  # already fetched
    dest.mkdir(parents=True, exist_ok=True)
    tarball = dest / "download.tar.gz"
    urllib.request.urlretrieve(url, tarball)  # supports file:// (tests)
    if _sha256(tarball) != expected:
        raise RuntimeError(f"sha256 mismatch for {url}")
    with tarfile.open(tarball) as tf:
        try:
            tf.extractall(dest, filter="data")
        except TypeError:  # Python < 3.12
            tf.extractall(dest)
    tarball.unlink()
    marker.touch()
    return dest
