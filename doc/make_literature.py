#!/usr/bin/env python3
"""Generate the "literature using ngsxfem" pages from ``doc/literature.yaml``.

Outputs (both are committed to the repository):

* ``doc/literature.md`` -- plain Markdown (rendered on GitHub, compiled to a
  PDF with pandoc in the extras workflow)
* ``doc/sphinx/xfem_misc/literature_entries.rst`` -- reStructuredText with
  embedded HTML for the Sphinx documentation (styled by
  ``doc/sphinx/_static/literature.css``)

Usage::

    python3 doc/make_literature.py            # regenerate both files
    python3 doc/make_literature.py --check    # exit code 1 if files are outdated

The script only depends on PyYAML.
"""
import argparse
import html
import os
import sys
from collections import Counter

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(HERE, "literature.yaml")
MD_FILE = os.path.join(HERE, "literature.md")
RST_FILE = os.path.join(HERE, "sphinx", "xfem_misc", "literature_entries.rst")

GITHUB_DATA_URL = "https://github.com/ngsxfem/ngsxfem/blob/master/doc/literature.yaml"
GITHUB_ISSUES_URL = "https://github.com/ngsxfem/ngsxfem/issues"

# (label for the markdown output, emoji icon, label for the HTML type badge)
TYPES = {
    "article": ("Article", "📄", "journal article"),
    "preprint": ("Preprint", "📝", "preprint"),
    "inproceedings": ("Proceedings", "🗣️", "proceedings"),
    "incollection": ("Book chapter", "📖", "book chapter"),
    "phdthesis": ("PhD thesis", "🎓", "PhD thesis"),
    "mastersthesis": ("Master's thesis", "🎓", "Master's thesis"),
    "bachelorsthesis": ("Bachelor's thesis", "🎓", "Bachelor's thesis"),
    "software": ("Software", "🧩", "software"),
    "dataset": ("Dataset", "💾", "dataset"),
}
THESES = ("phdthesis", "mastersthesis", "bachelorsthesis")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def load_data():
    with open(DATA_FILE, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    validate(data)
    return data


def validate(data):
    cat_ids = [c["id"] for c in data["categories"]]
    if len(set(cat_ids)) != len(cat_ids):
        sys.exit("duplicate category ids in literature.yaml")
    keys = Counter(e["key"] for e in data["entries"])
    dups = [k for k, n in keys.items() if n > 1]
    if dups:
        sys.exit("duplicate entry keys in literature.yaml: " + ", ".join(dups))
    for e in data["entries"]:
        for field in ("key", "type", "authors", "title", "year", "categories"):
            if field not in e:
                sys.exit(f"entry {e.get('key', '?')} misses field '{field}'")
        if e["type"] not in TYPES:
            sys.exit(f"entry {e['key']}: unknown type '{e['type']}'")
        for c in e["categories"]:
            if c not in cat_ids:
                sys.exit(f"entry {e['key']}: unknown category '{c}'")
        if not isinstance(e["authors"], list):
            sys.exit(f"entry {e['key']}: 'authors' must be a list")
        if "arxiv" in e and not isinstance(e["arxiv"], str):
            sys.exit(f"entry {e['key']}: quote the arXiv id (\"{e['arxiv']}\"), "
                     "otherwise YAML reads it as a number and drops trailing zeros")


def author_string(authors):
    authors = [str(a) for a in authors]
    if len(authors) == 1:
        return authors[0]
    return ", ".join(authors[:-1]) + " and " + authors[-1]


def code_label(link):
    if "label" in link:
        return link["label"]
    url = link["url"]
    if "github.com" in url:
        return "GitHub"
    if "gitlab" in url:
        return "GitLab"
    if "zenodo" in url or str(link.get("doi", "")).startswith("10.5281"):
        return "Zenodo"
    if "bitbucket" in url:
        return "Bitbucket"
    return "Code"


def venue_parts(e):
    """Return (venue, details) as plain text, e.g. ('Numer. Math.', '135(2), 313–332, 2017')."""
    t = e["type"]
    venue = ""
    details = []
    if t in ("article", "inproceedings", "incollection"):
        venue = e.get("journal", "")
        if t == "incollection":
            if e.get("series"):
                details.append(e["series"] + (f" {e['volume']}" if e.get("volume") else ""))
            if e.get("publisher"):
                details.append(e["publisher"])
        else:
            vol = ""
            if e.get("volume") is not None:
                vol = str(e["volume"])
                if e.get("number") is not None:
                    vol += f"({e['number']})"
            if vol:
                details.append(vol)
        if e.get("pages") is not None:
            details.append(str(e["pages"]))
        details.append(str(e["year"]))
    elif t == "preprint":
        venue = ""
        details.append(f"preprint, {e['year']}")
    elif t in THESES:
        venue = TYPES[t][0]
        if e.get("school"):
            details.append(e["school"])
        details.append(str(e["year"]))
    else:  # software, dataset
        venue = e.get("journal", "")
        details.append(f"{TYPES[t][0].lower()}, {e['year']}")
    return venue, ", ".join(details)


def primary_url(e):
    if e.get("doi"):
        return "https://doi.org/" + e["doi"]
    if e.get("arxiv"):
        return "https://arxiv.org/abs/" + str(e["arxiv"])
    if e.get("hdl"):
        return "https://hdl.handle.net/" + e["hdl"]
    return e.get("url")


def sort_entries(entries):
    # newest first, then alphabetically by first author
    return sorted(entries, key=lambda e: (-int(e["year"]), str(e["authors"][0]).split()[-1].lower(), e["title"].lower()))


# --------------------------------------------------------------------------
# Markdown output
# --------------------------------------------------------------------------
def md_entry(e):
    venue, details = venue_parts(e)
    parts = [f"{author_string(e['authors'])}. “{e['title']}”."]
    if venue:
        parts.append(f"*{venue}* {details}.")
    else:
        parts.append(f"{details}.")
    links = []
    if e.get("doi"):
        links.append(f"doi: [{e['doi']}](https://doi.org/{e['doi']})")
    if e.get("arxiv"):
        links.append(f"arXiv: [{e['arxiv']}](https://arxiv.org/abs/{e['arxiv']})")
    if e.get("hdl"):
        links.append(f"hdl: [{e['hdl']}](https://hdl.handle.net/{e['hdl']})")
    if e.get("url") and not e.get("doi"):
        links.append(f"url: <{e['url']}>")
    if links:
        parts.append(", ".join(links) + ".")
    if e.get("code"):
        clinks = []
        for link in e["code"]:
            label = code_label(link)
            if link.get("doi"):
                clinks.append(f"[{label} (doi: {link['doi']})]({link['url']})")
            else:
                clinks.append(f"[{label}]({link['url']})")
        parts.append("Code: " + ", ".join(clinks) + ".")
    if e.get("note"):
        parts.append(f"({e['note']})")
    return "* " + " ".join(parts)


def write_markdown(data):
    cite = data["cite"]
    out = []
    out.append("---")
    out.append("Scientific literature using `ngsxfem`")
    out.append("---")
    out.append("")
    out.append("<!-- This file is generated from doc/literature.yaml by doc/make_literature.py. Do not edit by hand. -->")
    out.append("")
    out.append("This list collects scientific works (journal articles, preprints, theses) in which `ngsxfem` has been used. "
               "Entries are sorted by category and year (newest first); where available, links to reproduction code and data are given. "
               f"If you used `ngsxfem` in your work and it is missing here, please open an issue or a pull request on [GitHub]({GITHUB_ISSUES_URL}) "
               "(the list is generated from `doc/literature.yaml`).")
    out.append("")
    out.append("### Citing `ngsxfem`")
    out.append("")
    out.append(f"If you use `ngsxfem` for your research, please cite: "
               f"{author_string(cite['authors'])}. “{cite['title']}”. *{cite['journal']}* "
               f"{cite['volume']}({cite['number']}), {cite['pages']}, {cite['year']}. "
               f"doi: [{cite['doi']}](https://doi.org/{cite['doi']}).")
    if cite.get("software_doi"):
        out[-1] += f" Software releases are archived on Zenodo: doi: [{cite['software_doi']}](https://doi.org/{cite['software_doi']})."
    out.append("")
    out.append("```bibtex")
    out.append(cite["bibtex"].rstrip())
    out.append("```")
    out.append("")
    for cat in data["categories"]:
        entries = [e for e in data["entries"] if cat["id"] in e["categories"]]
        if not entries:
            continue
        out.append(f"### {cat['title']}")
        if cat.get("description"):
            out.append("")
            out.append(cat["description"].strip())
        out.append("")
        for e in sort_entries(entries):
            out.append(md_entry(e))
            out.append("")
    return "\n".join(out).rstrip() + "\n"


# --------------------------------------------------------------------------
# Sphinx (rst + html) output
# --------------------------------------------------------------------------
def h(s):
    return html.escape(str(s), quote=True)


def html_entry(e, anchor=None):
    t = e["type"]
    label, icon, badge = TYPES[t]
    venue, details = venue_parts(e)
    url = primary_url(e)
    title = h(e["title"])
    title_html = f'<a class="pub-title" href="{h(url)}">{title}</a>' if url else f'<span class="pub-title">{title}</span>'
    cls = f"pub pub-{t}" + (" pub-thesis" if t in THESES else "")
    lines = [f'<li class="{cls}" id="{h(anchor or e["key"])}">']
    lines.append(f'  <span class="pub-icon" title="{h(label)}" aria-hidden="true">{icon}</span>')
    lines.append('  <div class="pub-body">')
    lines.append(f'    <div class="pub-head">{title_html}<span class="pub-year">{h(e["year"])}</span></div>')
    lines.append(f'    <div class="pub-authors">{h(author_string(e["authors"]))}</div>')
    venue_html = f'<em>{h(venue)}</em> ' if venue else ""
    lines.append(f'    <div class="pub-venue">{venue_html}{h(details)}</div>')
    badges = [f'<span class="pub-badge pub-badge-type">{h(badge)}</span>']
    if e.get("doi"):
        badges.append(f'<a class="pub-badge pub-badge-doi" href="https://doi.org/{h(e["doi"])}">doi:{h(e["doi"])}</a>')
    if e.get("arxiv"):
        badges.append(f'<a class="pub-badge pub-badge-arxiv" href="https://arxiv.org/abs/{h(e["arxiv"])}">arXiv:{h(e["arxiv"])}</a>')
    if e.get("hdl"):
        badges.append(f'<a class="pub-badge pub-badge-doi" href="https://hdl.handle.net/{h(e["hdl"])}">hdl:{h(e["hdl"])}</a>')
    if e.get("url") and not e.get("doi"):
        badges.append(f'<a class="pub-badge pub-badge-url" href="{h(e["url"])}">{"PDF" if t in THESES else "Website"}</a>')
    for link in e.get("code") or []:
        lab = code_label(link)
        tip = f' title="{h(link["doi"])}"' if link.get("doi") else ""
        badges.append(f'<a class="pub-badge pub-badge-code" href="{h(link["url"])}"{tip}>💾 {h(lab)}</a>')
    lines.append('    <div class="pub-links">' + " ".join(badges) + "</div>")
    if e.get("note"):
        lines.append(f'    <div class="pub-note">{h(e["note"])}</div>')
    lines.append("  </div>")
    lines.append("</li>")
    return lines


def rst_raw_html(lines):
    # every line of the raw block must be indented (also lines inside <pre>)
    flat = []
    for l in lines:
        flat += l.split("\n")
    return [".. raw:: html", ""] + ["   " + l for l in flat] + [""]


def write_rst(data):
    cite = data["cite"]
    entries = data["entries"]
    n_art = sum(1 for e in entries if e["type"] in ("article", "inproceedings", "incollection"))
    n_pre = sum(1 for e in entries if e["type"] == "preprint")
    n_th = sum(1 for e in entries if e["type"] in THESES)
    n_code = sum(1 for e in entries if e.get("code"))

    out = [".. This file is generated from doc/literature.yaml by doc/make_literature.py.",
           "   Do not edit by hand.", ""]

    # intro + statistics + citation box
    intro = [
        '<div class="pub-intro">',
        '<p>This page collects scientific works in which <code>ngsxfem</code> has been used: journal articles, preprints and theses, '
        'sorted by topic and year (newest first). Where available, links to reproduction code and data are given.</p>',
        '<div class="pub-stats">',
        f'  <span class="pub-stat">📄 <b>{n_art}</b> articles</span>',
        f'  <span class="pub-stat">📝 <b>{n_pre}</b> preprints</span>',
        f'  <span class="pub-stat">🎓 <b>{n_th}</b> theses</span>',
        f'  <span class="pub-stat">💾 <b>{n_code}</b> with code / data</span>',
        '</div>',
        f'<p class="pub-missing">Missing a publication? Please open an <a href="{GITHUB_ISSUES_URL}">issue</a> '
        f'or a pull request: the list is generated from <a href="{GITHUB_DATA_URL}"><code>doc/literature.yaml</code></a>.</p>',
        '</div>',
    ]
    out += rst_raw_html(intro)

    out += ["Citing ngsxfem", "=" * len("Citing ngsxfem"), ""]
    cite_html = [
        '<div class="pub-cite">',
        '<p>If you use <code>ngsxfem</code> for your research, please cite the software paper:</p>',
        '<ul class="pub-list">',
    ]
    cite_entry = dict(cite)
    cite_entry.setdefault("type", "article")
    cite_entry.setdefault("key", "ngsxfem")
    cite_entry.setdefault("categories", [])
    cite_html += ["  " + l for l in html_entry(cite_entry)]
    cite_html += ['</ul>']
    if cite.get("software_doi"):
        cite_html.append(f'<p>Software releases are archived on Zenodo under the concept DOI '
                         f'<a href="https://doi.org/{h(cite["software_doi"])}">{h(cite["software_doi"])}</a>.</p>')
    cite_html.append('<details class="pub-bibtex"><summary>BibTeX</summary>')
    cite_html.append('<pre>' + h(cite["bibtex"].rstrip()) + '</pre>')
    cite_html.append('</details>')
    cite_html.append('</div>')
    out += rst_raw_html(cite_html)

    for cat in data["categories"]:
        cat_entries = [e for e in entries if cat["id"] in e["categories"]]
        if not cat_entries:
            continue
        out += [cat["title"], "=" * len(cat["title"]), ""]
        lines = []
        if cat.get("description"):
            lines.append(f'<p class="pub-cat-desc">{h(cat["description"].strip())}</p>')
        lines.append('<ul class="pub-list">')
        for e in sort_entries(cat_entries):
            # entries listed in several categories get a unique anchor per category
            anchor = e["key"] if e["categories"][0] == cat["id"] else f'{e["key"]}-{cat["id"]}'
            lines += ["  " + l for l in html_entry(e, anchor)]
        lines.append("</ul>")
        out += rst_raw_html(lines)
    return "\n".join(out).rstrip() + "\n"


# --------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--check", action="store_true", help="only check whether the generated files are up to date")
    args = p.parse_args(argv)
    data = load_data()
    outputs = {MD_FILE: write_markdown(data), RST_FILE: write_rst(data)}
    outdated = []
    for path, content in outputs.items():
        current = open(path, encoding="utf-8").read() if os.path.exists(path) else None
        if current != content:
            outdated.append(path)
            if not args.check:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(content)
    if args.check:
        if outdated:
            print("outdated (run doc/make_literature.py):", *outdated, sep="\n  ")
            return 1
        print("literature files are up to date")
        return 0
    print(f"{len(data['entries'])} entries written to", *outputs, sep="\n  ")
    return 0


if __name__ == "__main__":
    sys.exit(main())
