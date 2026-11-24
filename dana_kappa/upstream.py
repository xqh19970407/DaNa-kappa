"""Load the separately obtained ALiEGNN baseline; its source is not vendored."""
import os
from pathlib import Path
import sys

UPSTREAM_REPOSITORY = "https://github.com/dannyzekunren/Dataset_thermoconductivity_pred"
UPSTREAM_COMMIT = "4b41cb2986402288c94d9b42a4565e59b48edef4"
FILES = ["aliegnn/__init__.py", "aliegnn/utils.py", "aliegnn/models/__init__.py",
         "aliegnn/models/utils.py", "aliegnn/models/alignn_egnn.py",
         "aliegnn/models/alignn_egnn_simple.py", "aliegnn/models/egnnconv.py",
         "aliegnn/models/egnnconv_GATv2.py", "aliegnn/models/egnn_e3nn.py"]


def configure(path=None):
    root = Path(path or os.environ.get("DANA_ALIEGNN_SOURCE", Path.home()/".cache/dana-kappa/aliegnn"))
    if not (root/"aliegnn/models/alignn_egnn.py").is_file():
        raise FileNotFoundError("Obtain the ALiEGNN source first: python scripts/get_aliegnn.py; "
                                "or pass --upstream /path/to/author/source.")
    sys.path.insert(0, str(root))
    return root
