"""Execute the original notebook code in an isolated directory and compare its tables.

Run with the existing analysis Python: python analysis_review/reproduce_notebook.py
No kernel server, downloads, or edits to original results are required. The scratch
directory is temporary. Full console output and comparison results stay here.
"""
from pathlib import Path
import contextlib
import json
import os
import shutil
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def compare_csv(original, reproduced):
    left, right = pd.read_csv(original), pd.read_csv(reproduced)
    if left.shape != right.shape or list(left) != list(right):
        return {"status": "different_schema", "original_shape": left.shape,
                "reproduced_shape": right.shape}
    errors = []
    max_absolute_error = 0.0
    for column in left:
        a, b = left[column], right[column]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            if not np.allclose(a, b, rtol=1e-6, atol=1e-8, equal_nan=True):
                errors.append(column)
            delta = (a.astype(float) - b.astype(float)).abs().dropna()
            if len(delta):
                max_absolute_error = max(max_absolute_error, float(delta.max()))
        elif not a.fillna("<NA>").equals(b.fillna("<NA>")):
            errors.append(column)
    return {"status": "match" if not errors else "different_values",
            "rows": len(left), "columns": len(left.columns),
            "different_columns": errors, "max_numeric_absolute_difference": max_absolute_error}


def main():
    notebook = json.loads((ROOT / "damp_secondary_analysis.ipynb").read_text())
    start = time.monotonic()
    result = {"executed_cells": [], "tables": {}, "errors": []}
    with tempfile.TemporaryDirectory(prefix="venomics-reproduction-") as temp:
        work = Path(temp)
        (work / "scratchpad").mkdir()
        (work / "analysis_outputs").mkdir()
        for name in ["df_damp_195625x57.tsv", "wet_lab_validated_amp_sources.xlsx"]:
            (work / name).symlink_to(ROOT / name)
        for file in (ROOT / "scratchpad").iterdir():
            target = work / "scratchpad" / file.name
            if file.suffix == ".py":
                shutil.copy2(file, target)
            else:
                target.symlink_to(file, target_is_directory=file.is_dir())
        namespace = {"__name__": "__main__"}
        # Notebook displays are irrelevant to numerical reproducibility and can be large.
        import IPython.display
        old_display = IPython.display.display
        IPython.display.display = lambda *args, **kwargs: None
        old_show = plt.show
        plt.show = lambda *args, **kwargs: plt.close("all")
        previous = Path.cwd()
        os.chdir(work)
        try:
            with (OUT / "notebook_execution.log").open("w") as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    for index, cell in enumerate(notebook["cells"]):
                        if cell["cell_type"] != "code":
                            continue
                        print(f"Executing original cell {index}", flush=True)
                        try:
                            exec(compile("".join(cell["source"]), f"notebook_cell_{index}", "exec"), namespace)
                            result["executed_cells"].append(index)
                        except Exception as error:
                            import traceback
                            traceback.print_exc()
                            result["errors"].append({"cell": index, "error": repr(error)})
                            break
            for file in sorted((work / "analysis_outputs").glob("*.csv")):
                original = ROOT / "analysis_outputs" / file.name
                if original.exists():
                    result["tables"][file.name] = compare_csv(original, file)
            originals = {p.name for p in (ROOT / "analysis_outputs").glob("*.csv")}
            result["tables_not_regenerated"] = sorted(originals - set(result["tables"]))
        finally:
            os.chdir(previous)
            IPython.display.display = old_display
            plt.show = old_show
    result["elapsed_seconds"] = time.monotonic() - start
    result["matching_tables"] = sum(r["status"] == "match" for r in result["tables"].values())
    (OUT / "notebook_reproduction.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "tables"}, indent=2))
    print("Table discrepancies:", {k: v for k, v in result["tables"].items() if v["status"] != "match"})


if __name__ == "__main__":
    main()
