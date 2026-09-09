"""
Unit tests for `sem.provenance.git_commit`, the dolfinx-free helper shared
by `sem.cli` (open-pore result header) and `sem.scripts.derive_geometry`
(`*_derivation.json`'s `provenance` object).
"""

import subprocess

from sem.provenance import git_commit


def test_git_commit_returns_short_sha_in_this_repo():
    sha = git_commit()
    assert sha != "unknown"
    assert 4 <= len(sha) <= 40
    int(sha, 16)  # must be valid hex


def test_git_commit_never_raises_on_subprocess_error(monkeypatch):
    def _boom(*args, **kwargs):
        raise OSError("git not found")

    monkeypatch.setattr(subprocess, "run", _boom)
    assert git_commit() == "unknown"


def test_git_commit_never_raises_on_timeout(monkeypatch):
    def _timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="git", timeout=2)

    monkeypatch.setattr(subprocess, "run", _timeout)
    assert git_commit() == "unknown"


def test_git_commit_returns_unknown_on_nonzero_exit(monkeypatch):
    class _Result:
        returncode = 128
        stdout = ""

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: _Result())
    assert git_commit() == "unknown"


def test_git_commit_returns_unknown_outside_a_repo(tmp_path):
    assert git_commit(cwd=tmp_path) == "unknown"
