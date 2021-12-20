# -*- coding: utf-8 -*-
import pytest


__author__ = "Marco Mernberger"
__copyright__ = "Marco Mernberger"
__license__ = "mit"


def test_init(small_tcga):
    tcga = small_tcga
    print(tcga.name)
    assert tcga.name == "mini_acc_maf_and_htseq"
    assert tcga.file_endpt == "https://api.gdc.cancer.gov/files/"
    assert tcga.cases_endpt == "https://api.gdc.cancer.gov/cases"
    assert tcga.projects_endpt == "https://api.gdc.cancer.gov/projects"
    assert tcga.gdc_client_command == "/project/cache/gdc_client/gdc-client"
    assert tcga.file_dir = Path("/machine/ffs/datasets/tcga")  #  "/project/incoming/tcga"
    assert tcga.cache_dir = Path("/project/cache/tcga") / self.name
    assert tcga.sample_types = sample_types
    assert tcga.data_types = data_types
    assert tcga.project_ids = project_ids
    assert tcga.client_file = client_file
