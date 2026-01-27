#!/usr/bin/env python3

import sys
import yaml
from pathlib import Path


BEGIN = "<!-- BEGIN HIGHLIGHTS -->"
END   = "<!-- END HIGHLIGHTS -->"


CARD_TEMPLATE = """\
<div class="highlight-card" style="background-image:url('{image}');">
  <div class="highlight-overlay"></div>
  <div class="highlight-content">
    <h3>{title}</h3>
    <p>
      <a href="https://doi.org/{doi}" target="_blank">
        {authors}, {journal} ({year})
      </a>
    </p>
  </div>
</div>
"""


def render_cards(highlights):
    blocks = []
    for h in highlights:
        blocks.append(
            CARD_TEMPLATE.format(
                title=h["title"],
                authors=h["authors"],
                journal=h["journal"],
                year=h["year"],
                doi=h["doi"],
                image=h["image"],
            )
        )
    return "\n".join(blocks)


def inject(html, rendered):
    if BEGIN not in html or END not in html:
        raise RuntimeError("Highlight markers not found in HTML.")

    before, rest = html.split(BEGIN, 1)
    _, after = rest.split(END, 1)

    return (
        before
        + BEGIN
        + "\n\n<div class=\"highlight-grid\">\n"
        + rendered
        + "\n</div>\n\n"
        + END
        + after
    )


def main():
    if len(sys.argv) != 3:
        print("Usage: render_highlights.py highlights.yaml index.html")
        sys.exit(1)

    yaml_file = Path(sys.argv[1])
    html_file = Path(sys.argv[2])

    highlights = yaml.safe_load(yaml_file.read_text())
    html = html_file.read_text()

    rendered = render_cards(highlights)
    updated = inject(html, rendered)

    html_file.write_text(updated)
    print(f"Injected {len(highlights)} highlights into {html_file}")


if __name__ == "__main__":
    main()
