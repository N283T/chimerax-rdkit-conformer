import importlib.util

import pytest
from pathlib import Path


@pytest.fixture
def script_module():
    """Import input_to_sdf.py as a module (bypassing PEP 723 metadata)."""
    script_path = Path(__file__).parent.parent.joinpath("src", "input_to_sdf.py")
    spec = importlib.util.spec_from_file_location("input_to_sdf", str(script_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
