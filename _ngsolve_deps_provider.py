"""
Custom scikit-build-core metadata provider for the ngsolve dependency.

This module is used instead of ``ngsolve._scikit_build_core_dependencies``
so that the ngsolve *native* extension module is never imported during the
build-metadata phase.  The native module requires shared libraries (e.g.
``libmkl_rt.so.2``) that are not present in standard manylinux build
containers or in isolated Windows build environments, causing the build to
fail.

By reading the installed-package metadata via :mod:`importlib.metadata` we
obtain the exact ngsolve version without loading any C extension.
"""

from __future__ import annotations

from importlib.metadata import version


def dynamic_metadata(
    field: str,
    settings: dict | None = None,
) -> list[str]:
    if field == "dependencies":
        ngsolve_version = version("ngsolve")
        return [f"ngsolve=={ngsolve_version}"]
    raise ValueError(f"Unknown metadata field: {field!r}")


def get_requires_for_dynamic_metadata(
    settings: dict | None = None,
) -> list[str]:
    return ["ngsolve"]
