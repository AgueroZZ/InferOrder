#!/usr/bin/env python3
"""Check the published surface and saved-result claims of the MPCurve site."""

import csv
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
PAGES = {
    "index.html": ("Method", "Simulation", "Analysis"),
    "method.html": ("Model", "Variational inference", "Several orderings"),
    "simulation_summary.html": ("Simulation design", "PCA initialization", "Isomap initialization"),
    "simulation_m1.html": ("Question and design", "Ordering recovery"),
    "simulation_m1_comparison.html": ("Question and design", "Ordering recovery", "Paired method differences", "Increasing feature count to 50", "Common Isomap initialization"),
    "simulation_m2.html": ("Feature assignments", "Sample orderings"),
    "estimate_intrinsic_m.html": ("Simulation design", "Methods", "Recovery of M", "Summary and reproduction"),
    "estimate_intrinsic_m_smooth.html": ("Simulation design", "Methods", "Recovery of M", "Summary and reproduction"),
    "fitness.html": ("One ordering", "Two orderings and environment groups"),
    "pancreas.html": ("One versus two orderings", "Feature assignments in the two-ordering fit"),
}
EXPECTED_IMAGES = {
    "index.html": 0,
    "method.html": 0,
    "simulation_summary.html": 3,
    "simulation_m1.html": 2,
    "simulation_m1_comparison.html": 7,
    "simulation_m2.html": 3,
    "estimate_intrinsic_m.html": 11,
    "estimate_intrinsic_m_smooth.html": 11,
    "fitness.html": 3,
    "pancreas.html": 3,
}
EXPECTED_TITLES = {
    "index.html": "MPCurve: Methods and Results",
    "method.html": "MPCurve and CAVI",
    "simulation_summary.html": "Simulation Summary: Single-Ordering Recovery",
    "simulation_m1.html": "Simulation: One Latent Ordering",
    "simulation_m1_comparison.html": "Simulation: Comparing Single-Ordering Recovery",
    "simulation_m2.html": "Simulation: Two Latent Orderings",
    "estimate_intrinsic_m.html": "Simulation: Estimating M with All-Monotone Trajectories",
    "estimate_intrinsic_m_smooth.html": "Simulation: Estimating M with One Monotone Anchor per Ordering",
    "fitness.html": "Analysis: Mutant Fitness Across Environments",
    "pancreas.html": "Analysis: Pancreatic Cell Loadings",
}


class PageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.headings = []
        self.links = []
        self.images = []
        self.captions = []
        self.text = []
        self._heading = None
        self._caption = None
        self._title = None
        self.title = ""
        self._skip_text = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "h2":
            self._heading = []
        elif tag == "caption":
            self._caption = []
        elif tag == "title":
            self._title = []
        elif tag in {"script", "style"}:
            self._skip_text = True
        elif tag == "a" and "href" in attrs:
            self.links.append(attrs["href"])
        elif tag == "img" and "src" in attrs:
            self.images.append(attrs["src"])

    def handle_endtag(self, tag):
        if tag == "h2" and self._heading is not None:
            self.headings.append("".join(self._heading).strip())
            self._heading = None
        elif tag == "caption" and self._caption is not None:
            self.captions.append("".join(self._caption).strip())
            self._caption = None
        elif tag == "title" and self._title is not None:
            self.title = "".join(self._title).strip()
            self._title = None
        elif tag in {"script", "style"}:
            self._skip_text = False

    def handle_data(self, data):
        if self._heading is not None:
            self._heading.append(data)
        if self._caption is not None:
            self._caption.append(data)
        if self._title is not None:
            self._title.append(data)
        if not self._skip_text:
            self.text.append(data)


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def page_text(page):
    return re.sub(r"\s+", " ", " ".join(page.text)).strip()


def check_auto_study(study_name, page, baseline_runs, show_comparators=True):
    auto_dir = ROOT / "experiments" / study_name / "auto_m_v034"
    validation = json.loads((auto_dir / "validation.json").read_text())
    auto_runs = rows(auto_dir / "summary" / "auto_runs.csv")
    auto_overall = rows(auto_dir / "summary" / "overall_summary.csv")
    auto_comparison = rows(auto_dir / "summary" / "comparison_summary.csv")[0]

    assert validation["complete"] and validation["datasets"] == 90
    assert validation["package_version"] == "0.3.4"
    assert validation["package_source_commit"] == (
        "15f2b0bbe5dfa61cd46da5160b2bc251e75a0475")
    assert validation["package_archive_sha256"] == (
        "58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6")
    assert validation["converged"] == 90 and validation["warnings"] == 0
    assert validation["independent_isomap_initialization_checks"] == 9
    assert len(auto_runs) == 90 and all(row["status"] == "success" for row in auto_runs)

    initial_exact = sum(int(row["selected_initial_M"]) == int(row["true_M"])
                        for row in auto_runs)
    effective_exact = sum(int(row["effective_M"]) == int(row["true_M"])
                          for row in auto_runs)
    assert validation["exact_initial_M"] == initial_exact
    assert validation["exact_effective_M"] == effective_exact
    assert min(int(row["minimum_initial_cluster_size"]) for row in auto_runs) >= 2

    baseline_exact = {
        method: sum(row["status"] == "success" and row["estimated_M"] == row["true_M"]
                    for row in baseline_runs if row["method"] == method)
        for method in ("adaptive", "forward")
    }
    reported_exact = {row["method"]: int(row["exact"]) for row in auto_overall}
    assert reported_exact == {
        "adaptive": baseline_exact["adaptive"],
        "auto_adaptive": effective_exact,
        "forward": baseline_exact["forward"],
    }

    baseline_adaptive = {
        row["id"]: row for row in baseline_runs if row["method"] == "adaptive"
    }
    gained_exact = 0
    lost_exact = 0
    for row in auto_runs:
        old = baseline_adaptive[row["id"]]
        old_exact = old["status"] == "success" and old["estimated_M"] == old["true_M"]
        new_exact = int(row["effective_M"]) == int(row["true_M"])
        gained_exact += new_exact and not old_exact
        lost_exact += old_exact and not new_exact
    assert int(auto_comparison["gained_exact"]) == gained_exact
    assert int(auto_comparison["lost_exact"]) == lost_exact

    text = page_text(page)
    assert "MPCurver 0.3.4" in text
    assert f"{initial_exact}/90" in text and f"{effective_exact}/90" in text
    if show_comparators:
        assert f'{baseline_exact["adaptive"]}/90' in text
        assert f'{baseline_exact["forward"]}/90' in text
    else:
        assert f'{baseline_exact["adaptive"]}/90' not in text
        assert "Original adaptive EB" not in text
        assert "Uniform + forward" not in text
    assert "minimum cluster size" in text and "two" in text
    assert "independently within" in text and "selected feature group" in text


def check_smooth_study(page, baseline_runs):
    study = ROOT / "experiments/estimate_intrinsic_m_smooth_v032"
    auto_dir = study / "auto_m_v034"
    auto_runs = rows(auto_dir / "summary" / "auto_runs.csv")
    current_conditions = rows(auto_dir / "summary" / "current_condition_summary.csv")
    assert len(auto_runs) == 90 and len(current_conditions) == 9
    assert all(row["status"] == "success" and row["method"] == "auto_adaptive"
               for row in auto_runs)
    assert all(int(row["datasets"]) == 10 for row in current_conditions)
    text = page_text(page)
    assert "MPCurver 0.3.4" in text and "effective M" in text
    assert "90 current-package fits" in text
    assert "MPCurver 0.3.2" not in text
    assert "Original adaptive EB" not in text
    assert "Uniform + forward" not in text
    check_auto_study("estimate_intrinsic_m_smooth_v032", page,
                     rows(study / "main_summary" / "runs.csv"),
                     show_comparators=False)


def main():
    actual = {path.name for path in DOCS.glob("*.html")}
    assert actual == set(PAGES), f"Unexpected public HTML pages: {actual ^ set(PAGES)}"

    parsed = {}
    for name, required_headings in PAGES.items():
        path = DOCS / name
        page = PageParser()
        html = path.read_text()
        page.feed(html)
        parsed[name] = page
        assert page.title == EXPECTED_TITLES[name], name
        assert all(heading in page.headings for heading in required_headings), name
        assert len(page.images) == EXPECTED_IMAGES[name], name
        assert {item for item in PAGES if item != name}.issubset(set(page.links)), name
        assert "id=\"workflowr-report\"" not in html, name
        assert "custom <code>fig.path</code>" not in html, name

        text = page_text(page)
        assert ("MPCurve" if name == "simulation_summary.html" else "CAVI") in text, name
        assert "smooth-em" not in text.lower(), name
        assert "smoothemr" not in text.lower(), name
        for item in page.links + page.images:
            if item.startswith("data:image/png;base64,"):
                assert item[22:34].startswith("iVBORw0KGgo"), name
                continue
            parsed_url = urlsplit(item)
            if parsed_url.scheme or parsed_url.netloc or not parsed_url.path:
                continue
            target = path.parent / unquote(parsed_url.path)
            assert target.exists(), f"Broken local reference in {name}: {item}"

    comparison_text = page_text(parsed["simulation_m1_comparison.html"])
    for study_name in ("m1_bspline_comparison", "m1_bspline_comparison_p50"):
        study = ROOT / "experiments" / study_name
        metrics = rows(study / "metrics.csv")
        summary = rows(study / "summary.csv")
        assert len(metrics) == 300 and len(summary) == 10
        assert len({(row["id"], row["method"]) for row in metrics}) == 300
        assert all(row["valid"] == "TRUE" for row in metrics)
        assert all(f'{float(row["median_rho"]):.3f}' in comparison_text for row in summary)
    isomap_study = ROOT / "experiments/m1_bspline_isomap"
    isomap_metrics = rows(isomap_study / "metrics.csv")
    isomap_summary = rows(isomap_study / "summary.csv")
    assert len(isomap_metrics) == 600 and len(isomap_summary) == 20
    assert len({(row["id"], row["method"]) for row in isomap_metrics}) == 600
    for row in isomap_summary:
        assert f'{float(row["median_rho"]):.4f}' in comparison_text
        if row["method"] != "Isomap":
            assert f'{float(row["median_total_seconds"]):.2f}' in comparison_text
    dimension_summary = rows(ROOT / "experiments/m1_bspline_comparison_p50/dimension_summary.csv")
    assert len(dimension_summary) == 10
    for row in dimension_summary:
        assert f'{float(row["mean_change"]):.3f}' in comparison_text
        if row["method"] != "PCA":
            assert f'{float(row["median_seconds_p50"]):.2f}' in comparison_text

    assert parsed["index.html"].headings == ["Method", "Simulation", "Analysis"]
    index_text = page_text(parsed["index.html"])
    assert "Each study records the package version" in index_text
    index_source = (ROOT / "analysis" / "index.Rmd").read_text()
    for label in ("- **Fixed $M=1$**", "- **Fixed $M=2$**",
                  "- **Estimate $M$ from the data**"):
        assert label in index_source
    for label in ("All trajectories monotone", "One monotone anchor per ordering"):
        assert label in index_text
    intrinsic = page_text(parsed["estimate_intrinsic_m.html"])
    assert "MPCurver 0.3.2" in intrinsic and "effective M" in intrinsic
    study_runs = rows(ROOT / "experiments/estimate_intrinsic_m_v032/main_summary/runs.csv")
    assert len(study_runs) == 180 and all(row["status"] == "success" for row in study_runs)
    for method, expected in (("adaptive", 87), ("forward", 90)):
        assert sum(row["estimated_M"] == row["true_M"] for row in study_runs if row["method"] == method) == expected
    check_auto_study("estimate_intrinsic_m_v032", parsed["estimate_intrinsic_m.html"], study_runs)
    check_smooth_study(parsed["estimate_intrinsic_m_smooth.html"], study_runs)
    assert len(parsed["fitness.html"].captions) == 1
    assert len(parsed["pancreas.html"].captions) == 2
    assert "Environment counts" in parsed["fitness.html"].captions[0]
    assert "Fixed model" in parsed["pancreas.html"].captions[0]
    assert "Posterior ordering probabilities" in parsed["pancreas.html"].captions[1]

    current_study = ROOT / "experiments/mpcurve_v034_site"
    source_hash = (current_study / "source/archive_sha256.txt").read_text().split()[0]
    assert source_hash == (
        "58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6"
    )

    sim_dir = current_study / "results/simulations"
    sim_provenance = rows(sim_dir / "package_provenance.csv")[0]
    assert sim_provenance["package_version"] == "0.3.4"
    assert sim_provenance["package_source_commit"] == (
        "15f2b0bbe5dfa61cd46da5160b2bc251e75a0475"
    )
    assert sim_provenance["source_archive_sha256"] == source_hash
    m1 = rows(sim_dir / "simulation_m1_summary.csv")[0]
    m2 = rows(sim_dir / "simulation_m2_summary.csv")[0]
    assert m1["converged"] == m2["converged"] == "TRUE"
    assert f"{float(m1['abs_spearman']):.4f}" in page_text(parsed["simulation_m1.html"])
    m2_text = page_text(parsed["simulation_m2.html"])
    for column in ("abs_spearman_A", "abs_spearman_B"):
        assert f"{float(m2[column]):.4f}" in m2_text
    assert int(m2["d"]) == 20 and float(m2["partition_accuracy"]) == 1
    assert "20 of 20" in m2_text
    for name in ("index.html", "method.html", "simulation_m1.html", "simulation_m2.html",
                 "fitness.html", "pancreas.html"):
        assert "0.3.4" in page_text(parsed[name]), name

    pancreas_dir = current_study / "results/pancreas"
    pancreas_provenance = rows(pancreas_dir / "package_provenance.csv")[0]
    assert pancreas_provenance == sim_provenance
    pancreas = rows(pancreas_dir / "fit_summary.csv")
    assert len(pancreas) == 2
    assert all(int(row["samples"]) == 865 and int(row["factors"]) == 8 for row in pancreas)
    assert all(row["converged"] == "TRUE" for row in pancreas)
    assignments = rows(pancreas_dir / "factor_assignments.csv")
    assert len(assignments) == 8
    assert {row["assignment"] for row in assignments} == {"A", "B"}
    assert sum(row["assignment"] == "A" for row in assignments) == 7
    assert "Seven factors follow ordering A" in page_text(parsed["pancreas.html"])

    fitness_dir = ROOT / "data/site_assets/fitness"
    counts = rows(fitness_dir / "assignment_counts.csv")
    totals = {ordering: sum(int(row["environments"]) for row in counts if row["ordering"] == ordering)
              for ordering in ("A", "B")}
    assert totals == {"A": 27, "B": 18}
    correlations = rows(fitness_dir / "ordering_correlation.csv")
    assert correlations[0]["B"] == "0.857"
    fitness_provenance = rows(fitness_dir / "source_provenance.csv")[0]
    assert fitness_provenance["package_version"] == "0.3.4"
    assert fitness_provenance["package_source_commit"] == (
        "15f2b0bbe5dfa61cd46da5160b2bc251e75a0475"
    )
    assert fitness_provenance["article_source_commit"] == (
        "eeaf65bcf4c8b021d590c88ad411506359202f1f"
    )
    comparison = rows(fitness_dir / "single_ordering_comparison.csv")
    assert len(comparison) == 1 and comparison[0]["metric"] == "Pearson correlation"
    fitness_text = page_text(parsed["fitness.html"])
    assert "27 environments" in fitness_text and "18 to" in fitness_text
    assert "0.857" in fitness_text
    assert f"{float(comparison[0]['value']):.3f}" in fitness_text

    print(f"Checked {len(PAGES)} public pages, local links, images, terminology, and saved-result claims.")


if __name__ == "__main__":
    main()
