#!/usr/bin/env python3
"""Build local version of spec2nii and wheel artifacts for this repository."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a local version of spec2nii and wheel artifacts."
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove previous build artifacts before building.",
    )
    parser.add_argument(
        "--out-dir",
        default="dist",
        help="Output directory for build artifacts (default: dist).",
    )
    parser.add_argument(
        "--venv-dir",
        default=".venv",
        help="Virtual environment directory (default: .venv).",
    )
    parser.add_argument(
        "--skip-venv",
        action="store_true",
        help="Use current interpreter instead of creating/reusing a local venv.",
    )
    parser.add_argument(
        "--skip-install-deps",
        action="store_true",
        help="Skip dependency installation in the selected Python environment.",
    )
    return parser.parse_args()


def run(cmd: list[str], cwd: Path) -> None:
    print(f"+ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def display_path(path: Path, repo_root: Path) -> str:
    try:
        return str(path.relative_to(repo_root))
    except ValueError:
        return str(path)


def venv_python_path(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def ensure_venv(venv_dir: Path, cwd: Path) -> Path:
    python_path = venv_python_path(venv_dir)
    if python_path.exists():
        print(f"Using virtual environment at {venv_dir}", flush=True)
        return python_path

    print(f"Creating virtual environment at {venv_dir}", flush=True)
    run([sys.executable, "-m", "venv", str(venv_dir)], cwd=cwd)
    if not python_path.exists():
        raise RuntimeError(f"Virtual environment creation failed: {python_path} missing")
    return python_path


def install_dependencies(python_exe: Path, cwd: Path) -> None:
    run(
        [
            str(python_exe),
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "setuptools",
            "wheel",
            "build",
        ],
        cwd=cwd,
    )
    run(
        [
            str(python_exe),
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "--no-build-isolation",
            "-e",
            ".",
        ],
        cwd=cwd,
    )


def can_import_module(python_exe: Path, module: str, cwd: Path) -> bool:
    check = subprocess.run(
        [str(python_exe), "-c", f"import {module}"],
        cwd=cwd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return check.returncode == 0


def clean_artifacts(repo_root: Path, out_dir: Path) -> None:
    paths_to_remove = [repo_root / "build", out_dir]
    paths_to_remove.extend(repo_root.glob("*.egg-info"))

    for path in paths_to_remove:
        if not path.exists():
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
        print(f"Removed {display_path(path, repo_root)}", flush=True)


def main() -> int:
    args = parse_args()
    script_path = Path(__file__).resolve()
    repo_root = script_path.parent.parent
    out_dir = (repo_root / args.out_dir).resolve()
    venv_dir = (repo_root / args.venv_dir).resolve()

    if args.clean:
        clean_artifacts(repo_root, out_dir)

    if args.skip_venv:
        python_exe = Path(sys.executable).resolve()
        print(f"Using current interpreter: {python_exe}", flush=True)
    else:
        python_exe = ensure_venv(venv_dir=venv_dir, cwd=repo_root)

    if args.skip_install_deps:
        print("Skipping dependency installation.", flush=True)
    else:
        install_dependencies(python_exe=python_exe, cwd=repo_root)

    out_dir.mkdir(parents=True, exist_ok=True)

    if can_import_module(python_exe=python_exe, module="build", cwd=repo_root):
        build_cmd = [
            str(python_exe),
            "-m",
            "build",
            "--sdist",
            "--wheel",
            "--no-isolation",
            "--outdir",
            str(out_dir),
        ]
        run(build_cmd, cwd=repo_root)
    else:
        print("Python package 'build' not found; using setup.py fallback.", flush=True)
        run([str(python_exe), "setup.py", "sdist", "bdist_wheel"], cwd=repo_root)

    artifacts = sorted(out_dir.glob("*"))
    if artifacts:
        print("\nBuild artifacts:")
        for item in artifacts:
            print(f"- {display_path(item, repo_root)}")
    else:
        print("\nNo artifacts produced.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
