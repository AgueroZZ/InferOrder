#!/usr/bin/env python3
"""Check the published surface and saved-result claims of the MPCurve site."""

import csv
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
PAGES = {
    "index.html": ("Method", "Simulation", "Analysis"),
    "method.html": ("Model", "Variational inference", "Several orderings"),
    "simulation_m1.html": ("Question and design", "Ordering recovery"),
    "simulation_m2.html": ("Feature assignments", "Sample orderings"),
    "fitness.html": ("One ordering", "Two orderings and environment groups"),
    "pancreas.html": ("One versus two orderings", "Feature assignments in the two-ordering fit"),
}
EXPECTED_IMAGES = {
    "index.html": 0,
    "method.html": 0,
    "simulation_m1.html": 2,
    "simulation_m2.html": 3,
    "fitness.html": 3,
    "pancreas.html": 3,
}
EXPECTED_TITLES = {
    "index.html": "MPCurve: Methods and Results",
    "method.html": "MPCurve and CAVI",
    "simulation_m1.html": "Simulation: One Latent Ordering",
    "simulation_m2.html": "Simulation: Two Latent Orderings",
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
        assert "CAVI" in text, name
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

    assert parsed["index.html"].headings == ["Method", "Simulation", "Analysis"]
    assert "MPCurver 0.3.0" in page_text(parsed["index.html"])
    assert len(parsed["fitness.html"].captions) == 1
    assert len(parsed["pancreas.html"].captions) == 2
    assert "Environment counts" in parsed["fitness.html"].captions[0]
    assert "Fixed model" in parsed["pancreas.html"].captions[0]
    assert "Posterior ordering probabilities" in parsed["pancreas.html"].captions[1]

    sim_dir = ROOT / "experiments/mpcurve_v030_site/results/simulations"
    m1 = rows(sim_dir / "simulation_m1_summary.csv")[0]
    m2 = rows(sim_dir / "simulation_m2_summary.csv")[0]
    assert f"{float(m1['abs_spearman']):.4f}" in page_text(parsed["simulation_m1.html"])
    m2_text = page_text(parsed["simulation_m2.html"])
    for column in ("abs_spearman_A", "abs_spearman_B"):
        assert f"{float(m2[column]):.4f}" in m2_text
    assert int(m2["d"]) == 20 and float(m2["partition_accuracy"]) == 1
    assert "20 of 20" in m2_text

    pancreas_dir = ROOT / "experiments/mpcurve_v030_site/results/pancreas"
    pancreas = rows(pancreas_dir / "fit_summary.csv")
    assert len(pancreas) == 2
    assert all(int(row["samples"]) == 865 and int(row["factors"]) == 8 for row in pancreas)
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
    assert correlations[0]["B"] == "0.866"
    comparison = rows(fitness_dir / "single_ordering_comparison.csv")
    assert len(comparison) == 1 and comparison[0]["metric"] == "Pearson correlation"
    fitness_text = page_text(parsed["fitness.html"])
    assert "27 environments" in fitness_text and "18 to" in fitness_text
    assert "0.866" in fitness_text
    assert f"{float(comparison[0]['value']):.3f}" in fitness_text

    print("Checked six public pages, local links, images, terminology, and saved-result claims.")


if __name__ == "__main__":
    main()
