#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""tcga.py: Contains the TCGA class that offers some functionality for accessing
the GDC protal though the GDC API."""
from pathlib import Path
from typing import Callable, List, Optional, Dict
from pypipegraph import Job, FileGeneratingJob
import pandas as pd
import pypipegraph as ppg
import requests
import json
import os
import subprocess
import numpy as np
import tempfile
import zipfile

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class TCGA:
    def __init__(
        self,
        name: str,
        project_ids: List[str] = ["TCGA-BRCA"],
        data_types: List[str] = ["htseq_counts", "maf"],
        sample_types: Optional[List[str]] = None,
        client_path: Path = Path("/project/cache/gdc_client"),
        client_file: str = "gdc-client_v1.6.0_Ubuntu_x64-py3.7_0.zip",
    ) -> None:
        """
        Helper class to donwload and filter the TCGA files that we need.

        Parameters
        ----------
        name : str
            Unique name of the TCGA set.
        project_ids : List[str], optional
            List of poject ids, by default ["TCGA-BRCA"].
        datatypes : List[str], optional
            list of file types to select, by default ["htseq_counts", "maf"].
        sample_types : Optional[List[str]], optional, by default None.
            List of acceptable tissue types, by default ["primary tumor"].
        """
        self.name = name
        self.file_endpt = "https://api.gdc.cancer.gov/files/"
        self.cases_endpt = "https://api.gdc.cancer.gov/cases"
        self.projects_endpt = "https://api.gdc.cancer.gov/projects"
        self.client_path = client_path
        self.client_path.mkdir(parents=True, exist_ok=True)
        self.gdc_client_command = str(client_path / "gdc-client")
        self.file_dir = Path("/machine/ffs/datasets/tcga")  # "/project/incoming/tcga"
        self.cache_dir = Path("/project/cache/tcga") / self.name
        self.sample_types = sample_types
        if sample_types is None:
            self.sample_types = []
        self.data_types = data_types
        self.project_ids = project_ids
        self.client_file = client_file

    def download_client(self):
        """
        Downloads the GDC client to the client path.
        """
        url = f"https://gdc.cancer.gov/files/public/file/{self.client_file}"
        r = requests.get(url)
        f = tempfile.TemporaryFile()
        f.write(r.content)
        with zipfile.ZipFile(f, "r") as zip_ref:
            zip_ref.extractall(self.client_path)
        os.chmod(self.gdc_client_command, 0o770)

    def ensure_client(self) -> FileGeneratingJob:
        """
        Generates a FileGeneratingJob that downloads the gdc client.

        Dependency for the actual file download.

        Returns
        -------
        FileGeneratingJob
            Job to download GDC client.
        """
        self.client_path.mkdir(parents=True, exist_ok=True)
        outfile = self.client_path / "gdc-client"
        return ppg.FileGeneratingJob(outfile, self.download_client)

    def get_dependencies(self) -> List[Job]:
        """
        List of dependency jobs.

        Returns
        -------
        List[Job]
            List of dependency jobs.
        """
        return [
            self.ensure_client(),
            ppg.FunctionInvariant(self.name + "_ensure_client", self.ensure_client),
        ]

    def load(self) -> Job:
        """
        Returns a CachedDataLoadingJob that when done ensures that the TCGA class
        has some additional attributes.

        This adds several DataFrames to the class instance as attributes that
        need to be accessed by some methods in order to synchronize the files
        downloaded from the GDC data portal.

        Returns
        -------
        Job
            CachedDataLoadingJob that adds attributes to the class instance.
        """

        def __calc():
            """reads cached DataFrames after they are created."""
            dictionary_with_attributes = {}
            for data_type in self.data_types:
                for what in ["files", "cases"]:
                    filename = f"{self.name}_{data_type}"
                    fullpath = self.cache_dir / f"{filename}_{what}.tsv"
                    df = pd.read_csv(fullpath, sep="\t")
                    if "case_id" in df.columns:
                        df.index = df["case_id"]
                    dictionary_with_attributes[f"df_{filename}_{what}"] = df
            return dictionary_with_attributes

        def __load(dictionary_with_attributes):
            """loads pickled attributes form cache."""
            for attr_name in dictionary_with_attributes:
                print(attr_name)
                setattr(self, attr_name, dictionary_with_attributes[attr_name])

        jobs = []
        for data_type in self.data_types:
            # write DataFrame with requested cases
            job = self.write_df_from_request(
                f"{self.name}_{data_type}_cases",
                self.get_cases(self.project_ids, data_type, self.sample_types),
            )
            jobs.append(job)
            # write DataFrame with requested files
            job = self.write_df_from_request(
                f"{self.name}_{data_type}_files",
                self.get_files(self.project_ids, data_type, self.sample_types),
            )
            jobs.append(job)

        return (
            ppg.CachedDataLoadingJob(
                self.cache_dir / (self.name + "_load"), __calc, __load
            )
            .depends_on(self.get_dependencies())
            .depends_on(jobs)
            .depends_on(
                [
                    ppg.FunctionInvariant(self.name + "_get_filters", self._get_filter),
                    ppg.FunctionInvariant(self.name + "_get_files", self.get_files),
                    ppg.FunctionInvariant(self.name + "_get_cases", self.get_cases),
                    ppg.FunctionInvariant(
                        self.name + "_write_df_from_request", self.write_df_from_request
                    ),
                    ppg.ParameterInvariant(
                        self.name,
                        self.data_types + self.project_ids + self.sample_types,
                    ),
                ]
            )
        )

    def write_df_from_request(
        self, filename: str, request_function: Callable
    ) -> FileGeneratingJob:
        """
        Writes a df from a request to file.

        Returns a FileGeneratingJob that uses a callable for querying the GDC
        and a filename.

        Parameters
        ----------
        filename : str
            File name stem for result file.
        request_function : Callable
            callable that returns a DataFrame.

        Returns
        -------
        ppg.FileGeneratingJob
            Job that writes the file.
        """
        outfile = self.cache_dir / f"{filename}.tsv"
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __dump():
            result = request_function()
            if type(result) == tuple:
                df = result[0]
                with outfile.with_suffix(".log").open("w") as op:
                    op.write(result[1])
            else:
                df = result
            df.to_csv(outfile, sep="\t", index=False)

        return (
            ppg.FileGeneratingJob(outfile, __dump)
            .depends_on(
                [
                    ppg.FunctionInvariant(self.name + "_get_filters", self._get_filter),
                    ppg.FunctionInvariant(self.name + "_get_files", self.get_files),
                    ppg.FunctionInvariant(self.name + "_get_cases", self.get_cases),
                ]
            )
            .depends_on(self.get_dependencies())
        )

    def _get_filter(
        self, project_ids: List[str], datatype: str, sample_types: List[str] = []
    ) -> Dict[str, object]:
        """
        Returns a filter function for querying the gdc data protal using the GDC API.

        This returns a filter intended for filtering relevant cases and files
        based on the project id and the file type that we need.
        In addition, we can select for relevant tissue types.
        Currently that just supports htseq counts and maf.

        Parameters
        ----------
        project_ids : List[str]
            List of poject ids, e.g. ["TCGA-BRCA"].
        datatype : str
            list of file types to select, e.g. "htseq_counts".
        sample_types : List[str], by default [].
            List of acceptable tissue types, e.g. ["primary tumor"].

        Returns
        -------
        Dict[str, Object]
            A json dict.

        Raises
        ------
        ValueError
            "If the file type to filter for is not understood."
        """
        contents = []
        if len(sample_types) > 0:
            contents.append(
                {
                    "op": "=",
                    "content": {
                        "field": "cases.samples.sample_type",
                        "value": sample_types,
                    },
                }

            )
        if len(project_ids) > 0:
            contents.append(
                {
                    "op": "in",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": project_ids,
                    },
                }
            )
        if datatype == "htseq_counts":
            contents.extend(
                [
                    {
                        "op": "in",
                        "content": {
                            "field": "files.data_category",
                            "value": ["Transcriptome Profiling"],
                        },
                    },
                    {
                        "op": "=",
                        "content": {
                            "field": "files.analysis.workflow_type",
                            "value": "HTSeq - Counts",
                        },
                    },
                ]
            )
        elif datatype == "maf":
            contents.append(
                {"op": "in", "content": {"field": "files.data_format", "value": "maf"}}
            )
        elif datatype == "htseq_FPKM":
            contents.extend(
                [
                    {
                        "op": "in",
                        "content": {
                            "field": "files.data_category",
                            "value": ["Transcriptome Profiling"],
                        },
                    },
                    {
                        "op": "=",
                        "content": {
                            "field": "files.analysis.workflow_type",
                            "value": "HTSeq - FPKM",
                        },
                    },
                ]
            )
        else:
            raise ValueError(f"Don't know what to do with this datatype: {datatype}.")
        filters = {"op": "and", "content": contents}
        return filters

    def get_cases(
        self, project_ids: List[str], datatype: str, sample_types: List[str],
    ) -> Callable:
        """
        Retrieves case entities from the cases endpoint and returns them
        in a DataFrame.

        Depending on the filters to be specified, this queries the GDC data portal
        and converts the response to a DataFrame.

        Parameters
        ----------
        project_ids : List[str], optional
            List of poject ids, by default ["TCGA-BRCA"].
        datatype : str, optional
            list of file types to select, by default "htseq_counts".
        sample_types : List[str], optional
            List of acceptable tissue types, by default ["primary tumor"].

        Returns
        -------
        Callable
            Callable that returns a DataFrame with cases.
        """

        def call():
            _filter = self._get_filter(project_ids, datatype, sample_types)
            fields = [
                "case_id",
                "submitter_id",
                "primary_site",
                "project.project_id",
                "disease_type",
            ]
            params = {
                "filters": _filter,
                "fields": ",".join(fields),
                "format": "JSON",
            }
            df = self.fetch_df_from_endpoint(self.cases_endpt, params)
            return df, str(params)

        return call

    def fetch_df_from_endpoint(self, endpoint: str, params: dict):
        p = params.copy()
        params["filters"] = json.dumps(params["filters"])
        response = requests.get(endpoint, params=params)
        result_json = response.json()
        df = pd.read_json(json.dumps(result_json["data"]["hits"]))
        if df.empty:
            print("Paremeters:\n")
            print(json.dumps(p, indent=2))
            raise ValueError("The parameter provided did not reeeetieve any entities from the endpoint ({endpoint})")
        df = df.set_index("id")
        return df

    def get_files(
        self, project_ids: List[str], datatype: str, sample_types: List[str],
    ) -> Callable:
        """
        Retrieves file entities from the cases endpoint and returns them
        in a DataFrame.

        Depending on the filters to be specified, this queries the GDC data portal
        and converts the response to a DataFrame. In addition, it counts the
        number of cases associated with the files. If that is a 1-to-1
        relationship, the case_ids are used as indices of the DataFrame.

        Parameters
        ----------
        project_ids : List[str], optional
            List of poject ids, by default ["TCGA-BRCA"].
        datatype : str, optional
            list of file types to select, by default "htseq_counts".
        sample_types : List[str], optional
            List of acceptable tissue types, by default ["primary tumor"].

        Returns
        -------
        Callable
            Callable that returns a DataFrame with files.
        """

        def call():
            _filter = self._get_filter(project_ids, datatype, sample_types)
            retained = [
                "access",
                "data_category",
                "data_format",
                "data_type",
                "md5sum",
                "file_size",
                "state",
                "file_id",
                "file_name",
            ]
            fields = retained + [
                "cases.submitter_id",
                "cases.case_id",
                "cases.samples.tumor_descriptor",
                "cases.samples.tissue_type",
                "cases.samples.sample_type",
                "cases.samples.submitter_id",
                "cases.samples.sample_id",
                "cases.samples.portions.analytes.aliquots.aliquot_id",
                "cases.samples.portions.analytes.aliquots.submitter_id",
                "cases",
                "cases.project.project_id",
            ]
            params = {
                "filters": _filter,
                "fields": ",".join(fields),
                "format": "JSON",
            }
            df = self.fetch_df_from_endpoint(self.file_endpt, params)
            df = df[df["access"] == "open"]
            return df, str(params)

        return call

    def intersect_queries(self) -> Job:
        """
        Returns an AttributeLoadingJob that assures that a DataFrame with all
        usable cases depending on the queries specified is present as a class
        instance attribute.

        Intersects the cases for all data types obtained from queries to the
        GDC data portal and creates a DataFrame with all cases for which
        all data types are available.

        Returns
        -------
        Job
            Job that loads a DataFrame with all usable cases.
        """

        def get_cases_to_use():
            # intersect cases
            attr_name = f"df_{self.name}_{self.data_types[0]}_cases"
            df = getattr(self, attr_name)
            index = df.index
            for data_type in self.data_types[1:]:
                attr_name = f"df_{self.name}_{data_type}_cases"
                df2 = getattr(self, attr_name)
                index = index.intersection(df2.index)
            df_ret = df.loc[index]
            return df_ret

        return [
            ppg.CachedAttributeLoadingJob(
                self.cache_dir / f"{self.name}_df_cases",
                self,
                "df_cases",
                get_cases_to_use,
            ).depends_on(self.load()),
            self.write_df_from_request("df_cases_to_use", get_cases_to_use).depends_on(
                self.load()
            ),
        ]

    def write_manifest(self) -> None:
        """
        Writes a manifest file for the gdc download client.

        Filters the file DataFrames by case_id if applicable and removes
        files already present in the result dir. Remaining files are written to
        manifest.
        """
        output_files = []
        for project_id in self.project_ids:
            for data_type in self.data_types:
                manifest = (
                    self.cache_dir / project_id.lower() / data_type / "manifest.txt"
                )
                manifest_missing = (
                    self.cache_dir
                    / project_id.lower()
                    / data_type
                    / "manifest_missing.txt"
                )
                output_files.append(manifest)
                output_files.append(manifest_missing)
        invariant = []
        infolder = manifest.parent / "files"
        jobs = []
        if infolder.exists():
            for s in infolder.iterdir():
                invariant.append(s.name)
            jobs.append(
                ppg.ParameterInvariant(self.name + "_manifest_check", invariant)
            )

        def __write():
            df_cases = self.df_cases
            for project_id in self.project_ids:
                for data_type in self.data_types:
                    manifest = (
                        self.cache_dir / project_id.lower() / data_type / "manifest.txt"
                    )
                    manifest_missing = (
                        self.cache_dir
                        / project_id.lower()
                        / data_type
                        / "manifest_missing.txt"
                    )
                    manifest.parent.mkdir(exist_ok=True, parents=True)
                    output_folder = (
                        self.file_dir / project_id.lower() / data_type / "files"
                    )
                    output_folder.mkdir(parents=True, exist_ok=True)
                    df_files = getattr(self, f"df_{self.name}_{data_type}_files")
                    for c in df_files.columns:
                        print(c)
                        
                    df_files = df_files[
                        df_files["cases.0.project.project_id"] == project_id
                    ]
                    df = self.filter_by_case_id(df_files, df_cases.index, keep="first")
                    df = df[["file_id", "file_name", "md5sum", "file_size", "state"]]
                    df = df.rename(
                        columns={
                            "file_id": "id",
                            "file_name": "filename",
                            "file_size": "size",
                        }
                    )
                    df.to_csv(manifest, sep="\t", index=False)
                    dropme = set()

                    for subdir in output_folder.iterdir():
                        if subdir.name in df["id"].values:
                            dropme.update(df[df["id"] == subdir.name].index.values)
                    df = df.drop(dropme)
                    df.to_csv(manifest_missing, sep="\t", index=False)

        return (
            ppg.MultiFileGeneratingJob(output_files, __write)
            .depends_on(self.get_dependencies())
            .depends_on(self.load())
            .depends_on(self.intersect_queries())
            .depends_on(jobs)
        )

    def filter_by_case_id(
        self, df: pd.DataFrame, case_id_index: List, keep: str = "first"
    ) -> pd.DataFrame:
        """
        Filters a DataFrame by case_id and removes duplicates (e.g. two samples
        of the same case with same primary site).

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame to be filtered.
        case_id_index : Iterable
            The case_ids to filter by.
        keep : str, optional
            Keep duplicates, by default "first".

        Returns
        -------
        pd.DataFrame
            Filtered DataFrame.

        Raises
        ------
        ValueError
            If case_ids are not unique.
        """
        if len(df["no_cases"].unique()) == 1 and df["no_cases"].unique()[0] == 1:
            if keep is not all:
                df["filter_duplicate"] = df.duplicated("case_id", keep=keep)
                df = df[~df["filter_duplicate"]]
            df = df[df["case_id"].isin(case_id_index)]
            return df
        else:
            # ids are ambiguous, makes no sense to filter
            return df

    def fetch_files(self) -> ppg.job:
        """
        Download job for the remaining files.

        Returns
        -------
        ppg.job
            MultiFileGeneratingJob that starts the GDC downlöoad client.
        """
        outfiles_and_destinations = {}
        for project_id in self.project_ids:
            for data_type in self.data_types:
                stdout = (
                    self.cache_dir
                    / project_id.lower()
                    / data_type
                    / f"stdout.txt"
                )
                stdout.parent.mkdir(exist_ok=True, parents=True)
                output_folder = (
                    self.file_dir / project_id.lower() / data_type / "files"
                    )
                outfiles_and_destinations[stdout] = output_folder

        def __download():
            for outfile in outfiles_and_destinations:
                manifest = outfile.parent / "manifest_missing.txt"
                self.download(manifest, outfiles_and_destinations[outfile])

        return (
            ppg.MultiFileGeneratingJob(list(outfiles_and_destinations.keys()), __download)
            .depends_on(self.get_dependencies())
            .depends_on(self.write_manifest)
            .depends_on(self.load())
            .depends_on(self.intersect_queries())
        )

    def download(self, manifest: Path, outfolder: Path = None):
        """
        Starts the gdc client to download the files from a manifest.

        Parameters
        ----------
        manifest : Path
            Manifest file.
        """
        if outfolder is not None:
            output_folder = outfolder
        else:
            output_folder = manifest.parent
        stdout = manifest.parent / "stdout.txt"
        cmd = [self.gdc_client_command] + [
            "download",
            "-m",
            str(manifest.absolute()),
            "-d",
            str(output_folder.absolute()),
        ]
        try:
            with stdout.open("w") as outp:
                subprocess.check_call(cmd, stdout=outp)
        except subprocess.CalledProcessError:
            print(" ".join(cmd))
            raise ValueError()

    def write_htseq_meta(self) -> ppg.FileGeneratingJob:
        """
        Writes a tsv file containing file meta information for the cases to
        be used.

        Filters the files DataFrame by case_id.
        """
        outfiles = []
        for datatype in self.data_types:
            if "htseq" in datatype:
                outfile = (
                    Path("results")
                    / self.name
                    / f"df_{self.name}_{datatype}_files_used.tsv"
                )
                outfiles.append(outfile)
        if len(outfiles) == 0:
            raise ValueError("No htseq data requested.")
        outfiles[0].parent.mkdir(parents=True, exist_ok=True)

        def __write():
            df_cases = self.df_cases
            for data_type in self.data_types:
                if "htseq" in data_type:
                    outfile = (
                        Path("results")
                        / self.name
                        / f"df_{self.name}_{data_type}_files_used.tsv"
                    )
                    df_files = getattr(self, f"df_{self.name}_{data_type}_files")
                    columns_pretty = {
                        "cases.0.project.project_id": "project_id",
                    }
                    df_files = df_files.rename(columns=columns_pretty)
                    df = self.filter_by_case_id(df_files, df_cases.index, keep="first")
                    df.to_csv(str(outfile), sep="\t", index=False)

        return (
            ppg.MultiFileGeneratingJob(outfiles, __write)
            .depends_on(self.get_dependencies())
            .depends_on(self.load())
            .depends_on(self.intersect_queries())
            .depends_on(self.fetch_files())
        )

    def write_htseq_counts_df(self) -> ppg.Job:
        """
        Combines all htseq count data into a DataFrame.

        This depends on the job returned by self.write_htseq_meta().

        Returns
        -------
        ppg.Job
            Job that generates the dataframe and writes it to file.

        Raises
        ------
        ValueError
            If no htseq data was specified in the constructor.
        """
        outfiles = []
        for data_type in self.data_types:
            if "htseq" in data_type:
                outfiles.append(
                    Path("results") / self.name / f"df_{self.name}_{data_type}.tsv"
                )
        if len(outfiles) == 0:
            raise ValueError("No HTseq data requested.")
        outfiles[0].parent.mkdir(parents=True, exist_ok=True)

        def __write():
            for data_type in self.data_types:
                if "htseq" in data_type:
                    inpath = (
                        Path("results")
                        / self.name
                        / f"df_{self.name}_{data_type}_files_used.tsv"
                    )
                    df_meta = pd.read_csv(str(inpath), sep="\t")
                    dfs = []
                    for _, row in df_meta.iterrows():
                        case_id = row["case_id"]
                        project_id = row["project_id"]
                        file_id = row["file_id"]
                        file_name = row["file_name"]
                        path2file = (
                            self.file_dir
                            / project_id.lower()
                            / data_type
                            / "files"
                            / file_id
                            / file_name
                        )
                        df = pd.read_csv(
                            path2file,
                            compression="gzip",
                            sep="\t",
                            names=["ensembl", case_id],
                        )
                        df.index = [x.split(".")[0] for x in df["ensembl"]]
                        df = df[[case_id]]
                        dfs.append(df)
                    df_counts = dfs[0]
                    for df in dfs[1:]:
                        df_counts = df_counts.join(df)
                    df_counts["gene_stable_id"] = df_counts.index
                    outfile = (
                        Path("results") / self.name / f"df_{self.name}_{data_type}.tsv"
                    )
                    df_counts.to_csv(str(outfile), sep="\t", index=False)

        return (
            ppg.MultiFileGeneratingJob(outfiles, __write)
            .depends_on(self.get_dependencies())
            .depends_on(self.load())
            .depends_on(self.intersect_queries())
            .depends_on(self.write_htseq_meta())
            .depends_on(self.fetch_files())
        )

    def write_maf_df(self) -> ppg.Job:
        """
        Combines the maf mutations analysis data into a single dataframe and
        writes it to file.

        Returns
        -------
        ppg.Job
            The job that generates the file.
        """
        outfile = Path("results") / self.name / f"df_{self.name}_maf_used.tsv"
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __write():
            df_maf_files = getattr(self, f"df_{self.name}_maf_files")
            dfs = []
            print(len(self.df_cases["case_id"].unique()))
            for _, row in df_maf_files.iterrows():
                project_id = row["cases.0.project.project_id"]
                file_id = row["file_id"]
                file_name = row["file_name"]
                path2file = (
                    self.file_dir
                    / project_id.lower()
                    / "maf"
                    / "files"
                    / file_id
                    / file_name
                )
                print(path2file)
                caller = path2file.name.split(".")[2]
                df_caller = pd.read_csv(
                    path2file, compression="gzip", sep="\t", comment="#"
                )
                print(len(df_caller["case_id"].unique()))
                df_caller["Caller"] = [caller] * len(df_caller)
                #df_caller = df_caller[
                #    df_caller["case_id"].isin(self.df_cases["case_id"])
                #]
                dfs.append(df_caller)
            df = pd.concat(dfs)
            print(len(df["case_id"].unique()))
            
            raise ValueError()
            df.to_csv(str(outfile), sep="\t", index=False)

        return (
            ppg.FileGeneratingJob(outfile, __write)
            .depends_on(self.get_dependencies())
            .depends_on(self.load())
            .depends_on(self.intersect_queries())
            #            .depends_on(self.write_htseq_meta())
            .depends_on(self.fetch_files())
        )
