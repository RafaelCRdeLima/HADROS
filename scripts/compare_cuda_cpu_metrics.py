"""Compare CPU and CUDA image metrics for validation runs."""

from pathlib import Path
import sys

import numpy as np


def load_image(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def metrics(data):
    tau = data[:, 4]
    psurv = data[:, 5]
    intensity = data[:, 6]
    captured = data[:, 7]
    finite = np.isfinite(tau) & np.isfinite(psurv) & np.isfinite(intensity)
    valid = finite & (captured == 0)

    return {
        "mean_tau": float(np.mean(tau[finite])),
        "max_tau": float(np.max(tau[finite])),
        "mean_P_surv": float(np.mean(psurv[finite])),
        "total_intensity": float(np.sum(intensity[finite])),
        "valid_rays": int(np.count_nonzero(valid)),
    }


def rel_diff(cuda_value, cpu_value):
    if isinstance(cpu_value, int):
        return cuda_value - cpu_value
    if cpu_value == 0.0:
        return 0.0 if cuda_value == 0.0 else np.inf
    return (cuda_value - cpu_value) / cpu_value


def main():
    if len(sys.argv) < 5 or (len(sys.argv) - 2) % 3 != 0:
        print(
            "Usage: python compare_cuda_cpu_metrics.py REPORT.txt "
            "LABEL CPU.dat CUDA.dat [LABEL CPU.dat CUDA.dat ...]"
        )
        raise SystemExit(1)

    report_path = Path(sys.argv[1])
    report_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# CPU/CUDA validation comparison",
        "# metric cpu cuda relative_difference",
        "# valid_rays reports cuda-cpu as an integer difference",
    ]

    for i in range(2, len(sys.argv), 3):
        label = sys.argv[i]
        cpu_path = Path(sys.argv[i + 1])
        cuda_path = Path(sys.argv[i + 2])

        cpu = metrics(load_image(cpu_path))
        cuda = metrics(load_image(cuda_path))

        lines.append("")
        lines.append(f"[{label}]")
        lines.append(f"cpu_file {cpu_path}")
        lines.append(f"cuda_file {cuda_path}")

        for key in ("mean_tau", "max_tau", "mean_P_surv", "total_intensity", "valid_rays"):
            diff = rel_diff(cuda[key], cpu[key])
            if key == "valid_rays":
                lines.append(f"{key} {cpu[key]} {cuda[key]} {diff}")
            else:
                lines.append(f"{key} {cpu[key]:.16e} {cuda[key]:.16e} {diff:.16e}")

    report_path.write_text("\n".join(lines) + "\n")
    print(report_path)


if __name__ == "__main__":
    main()
