#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

# import pytest
import sys
from pathlib import Path
import shutil
import os
import pytest

os.environ["PYPIPEGRAPH_NO_LOGGING"] = "1"
import pypipegraph as ppg

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))


@pytest.fixture
def new_pipeline(request):
    print("new_pipeline called")
    if request.cls is None:
        target_path = Path(__file__).parent / "run" / ("." + request.node.name)
    else:
        target_path = (
            Path(__file__).parent
            / "run"
            / (request.cls.__name__ + "." + request.node.name)
        )
    if target_path.exists():
        shutil.rmtree(target_path)
    Path(target_path).mkdir(parents=True, exist_ok=True)
    old_dir = Path(os.getcwd()).absolute()
    try:
        os.chdir(target_path)
        try:
            Path("logs").mkdir(parents=True, exist_ok=True)
            Path("cache").mkdir(parents=True, exist_ok=True)
            Path("results").mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        rc = ppg.resource_coordinators.LocalSystem(1)

        def np():
            ppg.new_pipegraph(rc, quiet=True)
            ppg.util.global_pipegraph.result_dir = Path("results")
            g = ppg.util.global_pipegraph
            g.new_pipeline = np
            return g

        yield np()
        try:
            # shutil.rmtree(Path(__file__).parent / "run")
            pass
        except OSError:
            pass
    finally:
        os.chdir(old_dir)


@pytest.fixture(scope="class")
def local_store():
    from mbf_externals import ExternalAlgorithmStore, change_global_store

    unpacked = Path(__file__).parent / "unpacked"
    if unpacked.exists():
        shutil.rmtree(unpacked)
    unpacked.mkdir()
    store = ExternalAlgorithmStore(Path(__file__).parent / "zipped", unpacked)
    change_global_store(store)
    yield store
    if unpacked.exists():
        shutil.rmtree(unpacked)


@pytest.fixture(scope="class")
def global_store():
    from mbf_externals import virtual_env_store

    store = virtual_env_store()
    yield store