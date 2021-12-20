# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for mtcga.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import pytest
from mtcga.tcga import TCGA

from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    both_ppg_and_no_ppg,
    no_pipegraph,
    pytest_runtest_makereport,
)

@pytest.fixture
def small_tcga():
    disease_code = "ACC"
    request_name = "maf_and_htseq"
    data_types = ["htseq_counts", "maf", "htseq_FPKM"]
    tcga = TCGA(
        f"mini_{disease_code.lower()}_{request_name}",
        project_ids=[f"TCGA-{disease_code}"],
        data_types=data_types,
        sample_types=["primary tumor"],
        )
    return tcga
