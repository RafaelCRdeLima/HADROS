#!/usr/bin/env python3
"""Build a static dashboard for BH_Torus_RTX plots and output products."""

from __future__ import annotations

import html
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLOTS_DIR = ROOT / "plots"
OUTPUT_DIR = ROOT / "output"
DASHBOARD_DIR = ROOT / "dashboard"
INDEX_PATH = DASHBOARD_DIR / "index.html"
MANIFEST_PATH = DASHBOARD_DIR / "plot_manifest.json"

IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".svg", ".pdf"}
DATA_EXTENSIONS = {".csv", ".dat", ".txt", ".md", ".vtk", ".vti", ".vtu"}

CATEGORY_ORDER = [
    "Kerr images",
    "Density backgrounds",
    "UHE source morphologies",
    "Robustness scans",
    "Opacity surfaces",
    "DIS model validation",
    "Spectral models",
    "MeV neutrino physics",
    "Collapsar/NDAF-like preset",
    "ParaView exports",
    "Other",
]

EXPECTED_CATEGORIES = CATEGORY_ORDER[:-1]


def rel_to_root(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def rel_from_dashboard(path: Path) -> str:
    return path.relative_to(DASHBOARD_DIR).as_posix() if path.is_relative_to(DASHBOARD_DIR) else f"../{rel_to_root(path)}"


def title_from_path(path: Path) -> str:
    title = path.stem.replace("_", " ").replace("-", " ")
    replacements = {
        "uhe": "UHE",
        "mev": "MeV",
        "kerr": "Kerr",
        "tau": "tau",
        "rgb": "RGB",
        "gbw": "GBW",
        "iim": "IIM",
        "paraview": "ParaView",
    }
    words = []
    for word in title.split():
        lower = word.lower()
        words.append(replacements.get(lower, word.capitalize()))
    return " ".join(words)


def categorize(path: Path) -> tuple[str, str, str, str]:
    rel = rel_to_root(path).lower()
    name = path.name.lower()

    if "collapsar_ndaf_like" in rel or "collapsar_ndaf_like" in name:
        return (
            "Collapsar/NDAF-like preset",
            "Fiducial versus semi-analytic collapsar/NDAF-like diagnostic comparison",
            "Comparison product for the literature-guided collapsar/NDAF-like preset.",
            "make audit_collapsar_ndaf_like",
        )

    if "dis_validation" in rel or "sigma_nun_model" in name or "dis_model" in name:
        return (
            "DIS model validation",
            "GBW/IIM/PDF-reference cross-section validation",
            describe_dis_validation(name),
            "make validate_dis_models",
        )

    if "mev_physics" in rel or name.startswith("mev_"):
        return (
            "MeV neutrino physics",
            "Thermal MeV neutrino post-processing",
            describe_mev(name),
            "make validate_mev_physics",
        )

    if path.suffix.lower() in {".vtk", ".vti", ".vtu"} or "paraview" in rel:
        return (
            "ParaView exports",
            "3D visualization export products",
            describe_paraview(name),
            "make paraview_fields",
        )

    if "opacity_surfaces" in rel or "opacity" in name or "tau_surface" in name or "p_surv_energy_theta" in name:
        return (
            "Opacity surfaces",
            "UHE optical-depth structure",
            describe_opacity(name),
            "make opacity_surfaces",
        )

    if "plots/scans" in rel or "output/scans" in rel:
        return (
            "Robustness scans",
            "Point-3 robustness and sensitivity scans",
            describe_scan(name),
            "make robustness_scans",
        )

    if "plots/spectra" in rel or "output/spectra" in rel or "spectrum" in name or "spectral" in name or "multiband" in name:
        command = "make plot_multiband_image" if "multiband" in name else "make validate_spectra"
        if "emitted_vs_observed" in name:
            command = "make plot_spectrum_observed"
        return (
            "Spectral models",
            "UHE spectral emission and attenuation",
            describe_spectral(name),
            command,
        )

    if "validation_uhe_sources" in rel or "uhe_source" in name:
        return (
            "UHE source morphologies",
            "Phenomenological UHE source prescriptions",
            describe_source(name),
            "make validate_source_plots",
        )

    if (
        "validation_density" in rel
        or "density" in name
        or "semi_analytic_background" in name
        or ("torus" in name and "rho0_torus" not in rel)
    ):
        return (
            "Density backgrounds",
            "Semi-analytic density morphology validation",
            describe_density(name),
            "make validate_density_plots",
        )

    if (
        "kerr_" in name
        or "combined_" in name
        or "rgb" in name
        or "rho0_torus" in rel
        or "enu" in rel
        or "energy_diagnostics" in rel
        or "geometry_schematic" in name
    ):
        command = "make plot-energy-diagnostics" if "diagnostic_" in name else "make plot-kerr-image"
        return (
            "Kerr images",
            "Kerr ray-traced image and overlay products",
            describe_kerr(name),
            command,
        )

    return ("Other", "Unclassified project output", "No known description pattern.", "unknown")


def describe_kerr(name: str) -> str:
    if "combined_contours" in name:
        return "Side-by-side UHE/MeV image with contour overlays."
    if "combined_side_by_side" in name:
        return "Side-by-side UHE and MeV image comparison."
    if "combined_rgb_overlay" in name or "kerr_rgb" in name:
        return "False-color UHE plus MeV composite image."
    if "kerr_uhe_neutrino_image_sqrt" in name:
        return "Square-root scaled UHE neutrino image."
    if "kerr_uhe_neutrino_image" in name:
        return "UHE neutrino image from cached Kerr geodesics."
    if "kerr_mev_thermal" in name:
        return "Thermal MeV neutrino image."
    if "diagnostic_tau" in name:
        return "Energy or radial optical-depth diagnostic."
    if "diagnostic_survival" in name:
        return "Survival probability diagnostic versus energy."
    if "diagnostic_total_flux" in name:
        return "Total image flux diagnostic versus energy."
    if "diagnostic_shadow" in name:
        return "Shadow contrast diagnostic versus energy."
    return "Kerr ray-tracing plot or image diagnostic."


def describe_density(name: str) -> str:
    if "radial" in name:
        return "Radial density profile diagnostic."
    if "equatorial" in name:
        return "Equatorial density profile diagnostic."
    if "polar" in name:
        return "Polar density profile diagnostic."
    if "map" in name:
        return "Two-dimensional density morphology map."
    return "Semi-analytic density background diagnostic."


def describe_source(name: str) -> str:
    if "radial" in name:
        return "Radial profile of UHE source morphology."
    if "angular" in name:
        return "Angular profile of UHE source morphology."
    if "map" in name:
        return "Two-dimensional UHE source morphology map."
    return "UHE source morphology validation plot."


def describe_scan(name: str) -> str:
    if "source_model" in name:
        return "Comparison across UHE source prescriptions."
    if "density_background" in name:
        return "Comparison across density background morphologies."
    if "energy_dependence" in name:
        return "Energy-dependence robustness scan."
    if "inclination" in name:
        return "Observer-inclination robustness scan."
    if "spin" in name:
        return "Black-hole spin robustness scan."
    if "parameter_stability" in name:
        return "Numerical and physical stability under source-parameter changes."
    return "Robustness scan output."


def describe_opacity(name: str) -> str:
    if "tau_surface_comparison" in name:
        return "Overlay of tau=2/3, tau=1, and tau=3 opacity surfaces."
    if "energy_dependence" in name or "energy_tau_surface" in name:
        return "Energy dependence of extracted opacity surfaces."
    if "phase" in name:
        return "Opacity phase diagram."
    if "p_surv" in name or "survival" in name:
        return "Survival probability as a function of energy and angle."
    if "classification" in name:
        return "Transparent, transition, semi-opaque, and opaque regions."
    if "source_independence" in name:
        return "Validation that opacity surfaces are independent of source morphology."
    if "background" in name:
        return "Opacity-surface comparison across density backgrounds."
    return "Opacity-surface or optical-depth structure product."


def describe_spectral(name: str) -> str:
    if "emitted_vs_observed" in name:
        return "Comparison between emitted and attenuated observed UHE spectrum."
    if "spectral_models" in name:
        return "Spectral weights for monochromatic, power-law, and cutoff models."
    if "survival_probability" in name:
        return "Survival probability versus UHE neutrino energy."
    if "weighted_vs_monochromatic" in name:
        return "Spectral-weighted image diagnostics compared with monochromatic output."
    if "multiband" in name:
        return "Energy-band composite image; false colors encode UHE energy intervals."
    return "UHE spectral model diagnostic."


def describe_dis_validation(name: str) -> str:
    if "sigma_nun_model_comparison" in name:
        return "Direct comparison of GBW, IIM, and PDF-reference UHE nuN cross sections."
    if "sigma_nun_model_ratios" in name:
        return "Cross-section ratios between DIS models as a function of neutrino energy."
    if "observed_spectrum" in name:
        return "Observed UHE spectrum propagated through different DIS cross-section models."
    if "attenuation_ratio" in name:
        return "Observed/emitted attenuation ratio for each DIS model."
    if "uhe_image_dis_model_comparison" in name:
        return "Side-by-side UHE image comparison using identical astrophysics and different DIS tables."
    if "uhe_image_dis_model_ratio" in name:
        return "UHE image ratio relative to the PDF-reference DIS model."
    return "DIS model validation output."


def describe_mev(name: str) -> str:
    if "emissivity_channels" in name:
        return "URCA, pair, bremsstrahlung, and total MeV emissivity versus temperature."
    if "opacity_vs_energy" in name:
        return "MeV absorption and scattering opacity versus neutrino energy."
    if "channel_dominance" in name:
        return "Dominant MeV emissivity channel in the rho-T plane."
    if "physical_vs_toy" in name:
        return "Image comparison between previous toy and physical MeV prescriptions."
    return "Thermal MeV neutrino physics diagnostic."


def describe_paraview(name: str) -> str:
    if "bh_torus_fields" in name:
        return "3D Cartesian sampling of the axisymmetric density/source model."
    if "tau_surface" in name:
        return "Extracted opacity surface for ParaView visualization."
    if "classification" in name:
        return "Opacity classification regions for ParaView visualization."
    return "ParaView-readable VTK/VTI/VTU export."


def known_data_candidates(plot: Path, data_files: list[Path]) -> list[Path]:
    plot_rel = rel_to_root(plot).lower()
    stem = plot.stem.lower()
    candidates = []

    for data in data_files:
        data_rel = rel_to_root(data).lower()
        data_stem = data.stem.lower()
        if data_stem == stem:
            candidates.append(data)
            continue
        if "spectra" in plot_rel and "spectra" in data_rel:
            candidates.append(data)
        elif "scans" in plot_rel and "scans" in data_rel:
            candidates.append(data)
        elif "opacity_surfaces" in plot_rel and "opacity_surfaces" in data_rel:
            candidates.append(data)
        elif "mev_physics" in plot_rel and "mev_physics_validation" in data_stem:
            candidates.append(data)

    return sorted(candidates, key=lambda p: (len(rel_to_root(p)), rel_to_root(p)))[:3]


def discover_files() -> tuple[list[Path], list[Path]]:
    images = []
    data = []
    for base in [PLOTS_DIR, OUTPUT_DIR]:
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            suffix = path.suffix.lower()
            if suffix in IMAGE_EXTENSIONS:
                images.append(path)
            elif suffix in DATA_EXTENSIONS:
                data.append(path)
    return images, data


def file_created_at(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")


def build_manifest(images: list[Path], data_files: list[Path]) -> list[dict[str, str]]:
    entries = []
    for image in images:
        category, module, description, command = categorize(image)
        data_candidates = known_data_candidates(image, data_files)
        entries.append(
            {
                "title": title_from_path(image),
                "category": category,
                "plot_path": rel_to_root(image),
                "data_path": "; ".join(rel_to_root(p) for p in data_candidates),
                "description": description,
                "related_module": module,
                "command": command,
                "created_at": file_created_at(image),
            }
        )

    for data in data_files:
        if data.suffix.lower() in {".vtk", ".vti", ".vtu"}:
            category, module, description, command = categorize(data)
            entries.append(
                {
                    "title": title_from_path(data),
                    "category": category,
                    "plot_path": "",
                    "data_path": rel_to_root(data),
                    "description": description,
                    "related_module": module,
                    "command": command,
                    "created_at": file_created_at(data),
                }
            )

    return sorted(entries, key=lambda e: (CATEGORY_ORDER.index(e["category"]) if e["category"] in CATEGORY_ORDER else 999, e["title"]))


def group_entries(entries: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped = defaultdict(list)
    for entry in entries:
        grouped[entry["category"]].append(entry)
    return grouped


def render_card(entry: dict[str, str]) -> str:
    plot = entry.get("plot_path", "")
    data = entry.get("data_path", "")
    plot_link = f"../{html.escape(plot)}" if plot else ""
    data_links = []
    for item in [d.strip() for d in data.split(";") if d.strip()]:
        data_links.append(f'<a href="../{html.escape(item)}">{html.escape(item)}</a>')

    preview = ""
    if plot:
        suffix = Path(plot).suffix.lower()
        if suffix in {".png", ".jpg", ".jpeg", ".svg"}:
            preview = f'<a href="{plot_link}"><img src="{plot_link}" alt="{html.escape(entry["title"])}"></a>'
        elif suffix == ".pdf":
            preview = f'<a class="file-preview" href="{plot_link}">Open PDF preview</a>'

    if not preview and data:
        first_data = data.split(";")[0].strip()
        preview = f'<a class="file-preview" href="../{html.escape(first_data)}">Open data product</a>'

    data_html = ", ".join(data_links) if data_links else "unknown"

    return f"""
      <article class="card">
        <div class="preview">{preview}</div>
        <div class="card-body">
          <h3>{html.escape(entry["title"])}</h3>
          <p>{html.escape(entry["description"])}</p>
          <dl>
            <dt>Path</dt><dd>{html.escape(plot) if plot else "data product"}</dd>
            <dt>Data</dt><dd>{data_html}</dd>
            <dt>Module</dt><dd>{html.escape(entry.get("related_module", "unknown"))}</dd>
            <dt>Command</dt><dd><code>{html.escape(entry.get("command", "unknown"))}</code></dd>
            <dt>Updated</dt><dd>{html.escape(entry.get("created_at", "unknown"))}</dd>
          </dl>
        </div>
      </article>
    """


def render_html(entries: list[dict[str, str]], image_count: int, data_count: int) -> str:
    grouped = group_entries(entries)
    category_counts = Counter(entry["category"] for entry in entries)
    missing = [cat for cat in EXPECTED_CATEGORIES if category_counts.get(cat, 0) == 0]
    updated = datetime.now().isoformat(timespec="seconds")

    category_pills = "\n".join(
        f'<span class="pill">{html.escape(cat)}: {category_counts.get(cat, 0)}</span>'
        for cat in CATEGORY_ORDER
        if category_counts.get(cat, 0) > 0
    )

    warning = ""
    if missing:
        warning = (
            '<div class="warning"><strong>Missing expected categories:</strong> '
            + ", ".join(html.escape(cat) for cat in missing)
            + "</div>"
        )

    sections = []
    for category in CATEGORY_ORDER:
        items = grouped.get(category, [])
        if not items:
            continue
        cards = "\n".join(render_card(entry) for entry in items)
        sections.append(
            f"""
            <details open>
              <summary>
                <span>{html.escape(category)}</span>
                <span class="count">{len(items)}</span>
              </summary>
              <div class="grid">{cards}</div>
            </details>
            """
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BH_Torus_RTX Plot Dashboard</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f7f9;
      --panel: #ffffff;
      --ink: #1d2733;
      --muted: #657182;
      --line: #d9dee7;
      --accent: #245b85;
      --soft: #e9f2f8;
      --warn: #fff5d6;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--ink);
    }}
    header {{
      padding: 28px clamp(18px, 4vw, 48px);
      background: var(--panel);
      border-bottom: 1px solid var(--line);
    }}
    h1 {{ margin: 0 0 8px; font-size: clamp(1.8rem, 4vw, 3rem); letter-spacing: 0; }}
    header p {{ margin: 0; color: var(--muted); max-width: 900px; line-height: 1.5; }}
    main {{ padding: 24px clamp(18px, 4vw, 48px) 48px; }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin-bottom: 18px;
    }}
    .metric {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px 16px;
    }}
    .metric strong {{ display: block; font-size: 1.6rem; color: var(--accent); }}
    .metric span {{ color: var(--muted); font-size: 0.9rem; }}
    .pills {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 16px 0; }}
    .pill {{
      background: var(--soft);
      color: var(--accent);
      border: 1px solid #c8ddea;
      padding: 6px 10px;
      border-radius: 999px;
      font-size: 0.9rem;
    }}
    .warning {{
      background: var(--warn);
      border: 1px solid #e6c96e;
      border-radius: 8px;
      padding: 12px 14px;
      margin: 14px 0 20px;
    }}
    details {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      margin: 16px 0;
      overflow: hidden;
    }}
    summary {{
      cursor: pointer;
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 14px 18px;
      font-weight: 700;
      color: var(--accent);
      border-bottom: 1px solid var(--line);
    }}
    .count {{
      background: var(--soft);
      border-radius: 999px;
      padding: 3px 10px;
      color: var(--accent);
      font-weight: 600;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
      gap: 14px;
      padding: 16px;
    }}
    .card {{
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      background: #fff;
      display: flex;
      flex-direction: column;
      min-height: 100%;
    }}
    .preview {{
      min-height: 160px;
      background: #eef1f5;
      display: grid;
      place-items: center;
      border-bottom: 1px solid var(--line);
    }}
    .preview img {{
      width: 100%;
      max-height: 260px;
      object-fit: contain;
      display: block;
      background: #fff;
    }}
    .file-preview {{
      color: var(--accent);
      text-decoration: none;
      font-weight: 700;
    }}
    .card-body {{ padding: 14px; }}
    .card h3 {{ margin: 0 0 8px; font-size: 1rem; line-height: 1.3; }}
    .card p {{ margin: 0 0 12px; color: var(--muted); line-height: 1.45; }}
    dl {{ margin: 0; display: grid; grid-template-columns: 82px 1fr; gap: 6px 10px; font-size: 0.84rem; }}
    dt {{ color: var(--muted); }}
    dd {{ margin: 0; overflow-wrap: anywhere; }}
    a {{ color: var(--accent); }}
    code {{ background: #eef1f5; padding: 2px 4px; border-radius: 4px; }}
    footer {{ color: var(--muted); font-size: 0.86rem; margin-top: 24px; }}
  </style>
</head>
<body>
  <header>
    <h1>BH_Torus_RTX Plot Dashboard</h1>
    <p>Static index of existing plots and output products. This dashboard only scans files; it does not run simulations or modify the physics pipeline.</p>
  </header>
  <main>
    <section class="summary">
      <div class="metric"><strong>{image_count}</strong><span>plots found</span></div>
      <div class="metric"><strong>{data_count}</strong><span>data files found</span></div>
      <div class="metric"><strong>{len(category_counts)}</strong><span>categories found</span></div>
      <div class="metric"><strong>{html.escape(updated)}</strong><span>last update</span></div>
    </section>
    <div class="pills">{category_pills}</div>
    {warning}
    {''.join(sections)}
    <footer>
      Manifest: <a href="plot_manifest.json">plot_manifest.json</a>. Links are relative to this dashboard directory.
    </footer>
  </main>
</body>
</html>
"""


def validate_links(entries: list[dict[str, str]]) -> list[str]:
    missing = []
    for entry in entries:
        for key in ["plot_path", "data_path"]:
            value = entry.get(key, "")
            for item in [part.strip() for part in value.split(";") if part.strip()]:
                if not (ROOT / item).exists():
                    missing.append(item)
    return missing


def main() -> int:
    DASHBOARD_DIR.mkdir(parents=True, exist_ok=True)
    images, data_files = discover_files()
    entries = build_manifest(images, data_files)
    missing = validate_links(entries)
    if missing:
        raise SystemExit("Broken dashboard paths:\n" + "\n".join(missing[:50]))

    MANIFEST_PATH.write_text(json.dumps(entries, indent=2) + "\n", encoding="utf-8")
    INDEX_PATH.write_text(render_html(entries, len(images), len(data_files)), encoding="utf-8")

    categories = Counter(entry["category"] for entry in entries)
    print(f"Dashboard written: {INDEX_PATH.relative_to(ROOT)}")
    print(f"Manifest written: {MANIFEST_PATH.relative_to(ROOT)}")
    print(f"Plots indexed: {len(images)}")
    print(f"Data files found: {len(data_files)}")
    print("Categories found: " + ", ".join(f"{cat}={count}" for cat, count in sorted(categories.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
