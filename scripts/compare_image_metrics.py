"""Compare UHE image metrics between two ray-tracing output files."""

from pathlib import Path
import sys

import numpy as np


def load_image(path):
    return np.loadtxt(path, comments="#")


def metrics(data):
    tau = data[:, 4]
    psurv = data[:, 5]
    intensity = data[:, 6]
    return {
        "mean_tau": float(np.mean(tau)),
        "max_tau": float(np.max(tau)),
        "mean_survival_probability": float(np.mean(psurv)),
        "integrated_image_intensity": float(np.sum(intensity)),
    }


def rel_diff(new, ref):
    if ref == 0.0:
        return 0.0 if new == 0.0 else np.inf
    return (new - ref) / ref


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: python scripts/compare_image_metrics.py "
            "REFERENCE.dat NEW.dat REPORT.txt"
        )
        raise SystemExit(1)

    ref_path = Path(sys.argv[1])
    new_path = Path(sys.argv[2])
    report_path = Path(sys.argv[3])
    report_path.parent.mkdir(parents=True, exist_ok=True)

    ref = metrics(load_image(ref_path))
    new = metrics(load_image(new_path))

    lines = [
        "# Gaussian backward-compatibility image metric comparison",
        f"reference_file {ref_path}",
        f"new_file {new_path}",
        "# metric reference new relative_difference",
    ]

    for key in ref:
        lines.append(f"{key} {ref[key]:.16e} {new[key]:.16e} {rel_diff(new[key], ref[key]):.16e}")

    report_path.write_text("\n".join(lines) + "\n")
    print(report_path)


if __name__ == "__main__":
    main()
