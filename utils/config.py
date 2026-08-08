from dataclasses import dataclass
from pathlib import Path
import yaml


@dataclass
class WorkflowConfig:
    output_dir: Path
    pdb_dir: Path
    ccd: Path
    valtrends_file: Path

    @property
    def density_dir(self) -> Path:
        return self.output_dir / "ccp4"

    @property
    def log_dir(self):
        return self.output_dir / "logs"
    
    @property
    def main_dir(self):
        return self.output_dir / "validation_data"
    
    @property
    def patterns_dir(self) -> Path:
        return self.main_dir / "result" / "RingsInHetResidues"

    @property
    def state_dir(self) -> Path:
        return self.main_dir / "last_state_data"

    @property
    def methods_info(self) -> Path:
        return self.state_dir / "methods_and_resolution.tsv"


def load_config(filename="config.yaml"):

    with open(filename) as f:
        data = yaml.safe_load(f)

    return WorkflowConfig(
        output_dir=Path(data["output_dir"]),
        pdb_dir=Path(data["pdb_dir"]),
        ccd=Path(data["ccd"]),
        valtrends_file=Path(data["valtrends_file"])
    )