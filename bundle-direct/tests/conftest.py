import importlib.util

import pytest
from pathlib import Path


@pytest.fixture
def script_module():
    """Import rdkit_ops.py as a module."""
    script_path = Path(__file__).parent.parent.joinpath("src", "rdkit_ops.py")
    spec = importlib.util.spec_from_file_location("rdkit_ops", str(script_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
