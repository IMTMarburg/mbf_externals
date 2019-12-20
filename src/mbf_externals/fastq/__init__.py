#  import pypipegraph as ppg
#  import time
#  from pathlib import Path
#  from .. import find_code_path
from ..externals import ExternalAlgorithm, reproducible_tar
from ..util import download_file


class FASTQC(ExternalAlgorithm):
    @property
    def name(self):
        return "FASTQC"

    def build_cmd(self, output_directory, ncores, arguments):
        input_files = arguments
        return [
            str(self.path / "FastQC" / "fastqc"),
            "-t",
            str(ncores),
            "--noextract",
            "--quiet",
            "-o",
            str(output_directory),
        ] + [str(x) for x in input_files]

    @property
    def multi_core(self):
        return False  # fastqc has a threads option - and does not make use of it

    def get_latest_version(self):
        return "0.11.8"

    def fetch_version(self, version, target_filename):  # pragma: no cover
        import tempfile
        from pathlib import Path
        import subprocess

        v = version

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            url = f"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v{v}.zip"
            with (tmpdir / "fastqc.zip").open("wb") as zip_file:
                download_file(url, zip_file)
            import zipfile

            with zipfile.ZipFile(zip_file.name, "r") as zip_ref:
                zip_ref.extractall(tmpdir / "target")
            subprocess.check_call(
                ["chmod", "+x", str(tmpdir / "target" / "FastQC" / "fastqc")]
            )
            reproducible_tar(target_filename.absolute(), "./", cwd=tmpdir / "target")
            print(f"done downloading FASTQC version {v}")
