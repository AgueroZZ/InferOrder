#!/usr/bin/env python3
"""Copy selected MPCurver 0.3.0 fitness results into InferOrder."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from html.parser import HTMLParser
from pathlib import Path


class TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self.table: list[list[str]] | None = None
        self.row: list[str] | None = None
        self.cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table":
            self.table = []
        elif self.table is not None and tag == "tr":
            self.row = []
        elif self.row is not None and tag in {"td", "th"}:
            self.cell = []

    def handle_data(self, data: str) -> None:
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag in {"td", "th"} and self.cell is not None and self.row is not None:
            self.row.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.row is not None and self.table is not None:
            if self.row:
                self.table.append(self.row)
            self.row = None
        elif tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mpcurver-root", type=Path, default=Path("../MPCurver"),
        help="Path to the MPCurver 0.3.0 repository",
    )
    args = parser.parse_args()

    source_root = args.mpcurver_root.resolve()
    description = (source_root / "DESCRIPTION").read_text()
    if "Version: 0.3.0" not in description:
        raise RuntimeError("The source repository must be MPCurver 0.3.0.")

    article = source_root / "docs/articles/fitness.html"
    figures = source_root / "docs/articles/fitness_files/figure-html"
    target = Path("data/site_assets/fitness")
    target.mkdir(parents=True, exist_ok=True)

    html = article.read_text()
    comparison_marker = "compare-single-orderings-1.png"
    if comparison_marker not in html:
        raise RuntimeError("Expected the single-ordering comparison figure.")
    comparison_prefix = html[:html.index(comparison_marker)]
    comparison_values = re.findall(r"#&gt; \[1\] ([0-9.]+)", comparison_prefix)
    if not comparison_values:
        raise RuntimeError("Expected the single-ordering comparison correlation.")
    single_ordering_correlation = float(comparison_values[-1])
    if not 0 <= single_ordering_correlation <= 1:
        raise RuntimeError("Invalid single-ordering comparison correlation.")
    with (target / "single_ordering_comparison.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["Pearson correlation", single_ordering_correlation])

    table_parser = TableParser()
    table_parser.feed(html)
    assignment_tables = [
        table for table in table_parser.tables
        if table and table[0] == ["environment", "category", "ordering"]
    ]
    if len(assignment_tables) != 1 or len(assignment_tables[0]) != 46:
        raise RuntimeError("Expected one environment-assignment table with 45 rows.")

    assignments = assignment_tables[0]
    with (target / "environment_assignments.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(assignments)

    counts: dict[tuple[str, str], int] = {}
    for _, category, ordering in assignments[1:]:
        key = (category, ordering)
        counts[key] = counts.get(key, 0) + 1
    with (target / "assignment_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["category", "ordering", "environments"])
        for (category, ordering), count in sorted(counts.items()):
            writer.writerow([category, ordering, count])

    correlation_tables = [
        table for table in table_parser.tables
        if table and table[0] == ["", "A", "B"] and len(table) == 3
    ]
    if len(correlation_tables) != 1:
        raise RuntimeError("Expected one two-ordering correlation table.")
    with (target / "ordering_correlation.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(correlation_tables[0])

    names = [
        "compare-single-orderings-1.png",
        "partition-plot-1.png",
        "known-sd-partition-trajectories-1.png",
    ]
    manifest = [(str(article.relative_to(source_root)), sha256(article))]
    for name in names:
        source = figures / name
        shutil.copy2(source, target / name)
        manifest.append((str(source.relative_to(source_root)), sha256(source)))

    with (target / "source_manifest.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["source_path_in_mpcurver", "sha256"])
        writer.writerows(manifest)

    print("Copied single-ordering correlation, 45 assignments, and three fitness figures.")


if __name__ == "__main__":
    main()
