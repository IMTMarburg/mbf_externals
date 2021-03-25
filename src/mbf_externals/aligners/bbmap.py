from .base import Aligner
import pypipegraph as ppg
import mbf_align
from pathlib import Path
from ..util import download_file, Version, download_tar_bz2_and_turn_into_tar_gzip
import subprocess
import tempfile
import os


class BBMap(Aligner):
    def __init__(self, version="38.86", store=None):
        super().__init__(version, store)

    @property
    def name(self):
        return "BBMap"

    @property
    def multi_core(self):
        return False

    def get_latest_version(self):
        return "38.86"

    def fetch_version(self, version, target_filename):  # pragma: no cover
        url = f"https://sourceforge.net/projects/bbmap/files/BBMap_{version}.tar.gz/download"
#        cmd = ["curl", url, "-L", "--output", "bb.tar.gz"]
#        subprocess.check_call(cmd)
        with open(target_filename, "wb") as op:
            download_file(url, op)

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        raise NotImplementedError

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_basename,
        output_bam_filename,
        parameters,
    ):
        raise NotImplementedError

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        raise NotImplementedError


class ExtendCigarBBMap(mbf_align.post_process._PostProcessor):
    def __init__(self, samformat="1.4"):
        self.samformat = samformat
        self.bbmap = BBMap()
        self.name = "BBMap_reformat"
        self.result_folder_name = Path("results") / "aligned" / "ExtendCigarBBMap"

    def process(self, input_bam_name, output_bam_name, result_dir):
        self.bbmap.store.unpack_version(self.bbmap.name, self.bbmap.version)
        cmd = [
            str(self.bbmap.path / "bbmap" / "reformat.sh"),
            f"in={str(input_bam_name.absolute().resolve())}",
            f"out={str(output_bam_name.absolute().resolve())}",
            f"sam={self.samformat}",
        ]
        print(" ".join(cmd))
        subprocess.check_call(cmd)

    def register_qc(self, new_lane):
        pass  # pragma: no cover

    def get_version(self):
        return self.bbmap.version

    def get_parameters(self):
        return self.samformat

