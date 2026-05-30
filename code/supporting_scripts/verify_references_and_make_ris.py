from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from difflib import SequenceMatcher
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT_MD = ROOT / "analysis" / "manuscript_assets" / "draft_manuscript" / "JTM_dedifferentiation_state_thyroid_cancer_revised_v2.md"
OUT_DIR = ROOT / "analysis" / "references"

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CROSSREF_WORKS = "https://api.crossref.org/works/"
USER_AGENT = "JTM-thyroid-reference-verifier/1.0 (mailto:example@example.com)"


def clean_text(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", value).strip()


def elem_text(elem: ET.Element | None) -> str:
    if elem is None:
        return ""
    return clean_text("".join(elem.itertext()))


def normalize_doi(value: str | None) -> str:
    value = clean_text(value or "").lower()
    value = re.sub(r"^\s*(https?://(dx\.)?doi\.org/)", "", value)
    return value.rstrip(".")


def normalize_title(value: str | None) -> str:
    value = clean_text(value or "").lower()
    value = re.sub(r"[^\w]+", " ", value)
    return clean_text(value)


def title_similarity(left: str, right: str) -> float:
    return SequenceMatcher(None, normalize_title(left), normalize_title(right)).ratio()


def fetch_url(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def parse_reference_lines() -> list[dict[str, str]]:
    text = MANUSCRIPT_MD.read_text(encoding="utf-8")
    ref_block = text.split("## References", 1)[1].split("## Tables", 1)[0]
    refs: list[dict[str, str]] = []
    for line in ref_block.splitlines():
        line = line.strip()
        if not re.match(r"^\d+\.\s+", line):
            continue
        number = re.match(r"^(\d+)\.", line).group(1)
        pmid_match = re.search(r"PMID:(\d+)", line)
        doi_match = re.search(r"doi:([^\.]+(?:\.[^\.]+)*?)\.\s+PMID:", line)
        before_pmid = re.sub(r"\s+PMID:\d+\.\s*$", "", line)
        before_pmid = re.sub(r"^\d+\.\s+", "", before_pmid)
        before_identifiers = re.sub(r"\s+doi:[^\.]+(?:\.[^\.]+)*\.$", "", before_pmid)
        parsed = re.match(r"^(?P<first_author>.+?)\. (?P<title>.+?)[\.\?\!] (?P<journal>.+?)\. (?P<year>\d{4})\.$", before_identifiers)
        refs.append(
            {
                "number": number,
                "line": line,
                "pmid": pmid_match.group(1) if pmid_match else "",
                "doi": normalize_doi(doi_match.group(1)) if doi_match else "",
                "first_author": parsed.group("first_author") if parsed else "",
                "title": parsed.group("title") if parsed else "",
                "journal": parsed.group("journal") if parsed else "",
                "year": parsed.group("year") if parsed else "",
            }
        )
    return refs


def fetch_pubmed(pmids: list[str]) -> dict[str, dict[str, object]]:
    params = urllib.parse.urlencode(
        {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "tool": "JTM-thyroid-reference-verifier",
        }
    )
    root = ET.fromstring(fetch_url(f"{EFETCH_URL}?{params}"))
    records: dict[str, dict[str, object]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = elem_text(article.find(".//PMID"))
        art = article.find(".//Article")
        journal = art.find("./Journal") if art is not None else None
        pub_date = journal.find("./JournalIssue/PubDate") if journal is not None else None
        year = elem_text(pub_date.find("./Year")) if pub_date is not None else ""
        if not year and pub_date is not None:
            medline_date = elem_text(pub_date.find("./MedlineDate"))
            year_match = re.search(r"\d{4}", medline_date)
            year = year_match.group(0) if year_match else ""
        authors = []
        for author in article.findall(".//AuthorList/Author"):
            collective = elem_text(author.find("./CollectiveName"))
            if collective:
                authors.append(collective)
                continue
            last = elem_text(author.find("./LastName"))
            fore = elem_text(author.find("./ForeName"))
            initials = elem_text(author.find("./Initials"))
            if last and fore:
                authors.append(f"{last}, {fore}")
            elif last and initials:
                authors.append(f"{last}, {initials}")
        doi = ""
        for article_id in article.findall(".//ArticleIdList/ArticleId"):
            if article_id.attrib.get("IdType") == "doi":
                doi = normalize_doi(elem_text(article_id))
                break
        if not doi:
            for elocation in article.findall(".//ELocationID"):
                if elocation.attrib.get("EIdType") == "doi":
                    doi = normalize_doi(elem_text(elocation))
                    break
        abstract_parts = [elem_text(x) for x in article.findall(".//Abstract/AbstractText")]
        records[pmid] = {
            "pmid": pmid,
            "title": elem_text(art.find("./ArticleTitle")) if art is not None else "",
            "journal": elem_text(journal.find("./Title")) if journal is not None else "",
            "journal_abbrev": elem_text(journal.find("./ISOAbbreviation")) if journal is not None else "",
            "year": year,
            "volume": elem_text(journal.find("./JournalIssue/Volume")) if journal is not None else "",
            "issue": elem_text(journal.find("./JournalIssue/Issue")) if journal is not None else "",
            "pages": elem_text(art.find("./Pagination/MedlinePgn")) if art is not None else "",
            "doi": doi,
            "authors": authors,
            "abstract": clean_text(" ".join([x for x in abstract_parts if x])),
        }
    return records


def fetch_crossref_by_doi(doi: str) -> dict[str, str]:
    if not doi:
        return {"status": "skipped_no_doi"}
    try:
        url = CROSSREF_WORKS + urllib.parse.quote(doi, safe="")
        payload = json.loads(fetch_url(url, timeout=30).decode("utf-8"))
        msg = payload.get("message", {})
        titles = msg.get("title") or []
        dates = msg.get("published-print") or msg.get("published-online") or msg.get("issued") or {}
        year = ""
        date_parts = dates.get("date-parts") or []
        if date_parts and date_parts[0]:
            year = str(date_parts[0][0])
        return {
            "status": "found",
            "title": clean_text(titles[0]) if titles else "",
            "journal": clean_text((msg.get("container-title") or [""])[0]),
            "year": year,
            "doi": normalize_doi(msg.get("DOI", "")),
        }
    except Exception as exc:  # noqa: BLE001
        return {"status": f"error:{type(exc).__name__}", "note": str(exc)[:180]}


def ris_escape(value: object) -> str:
    return clean_text(str(value or "")).replace("\n", " ")


def split_pages(pages: str) -> tuple[str, str]:
    pages = clean_text(pages)
    if "-" not in pages:
        return pages, ""
    first, last = pages.split("-", 1)
    return first, last


def format_ris(records: list[dict[str, object]]) -> str:
    blocks = []
    for rec in records:
        start_page, end_page = split_pages(str(rec.get("pages", "")))
        lines = ["TY  - JOUR"]
        for author in rec.get("authors", []):
            lines.append(f"AU  - {ris_escape(author)}")
        fields = [
            ("TI", rec.get("title")),
            ("JO", rec.get("journal")),
            ("JA", rec.get("journal_abbrev")),
            ("PY", rec.get("year")),
            ("VL", rec.get("volume")),
            ("IS", rec.get("issue")),
            ("SP", start_page),
            ("EP", end_page),
            ("DO", rec.get("doi")),
        ]
        for tag, value in fields:
            if value:
                lines.append(f"{tag}  - {ris_escape(value)}")
        doi = rec.get("doi")
        if doi:
            lines.append(f"UR  - https://doi.org/{ris_escape(doi)}")
        pmid = rec.get("pmid")
        if pmid:
            lines.append(f"UR  - https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            lines.append(f"AN  - PMID:{pmid}")
        if rec.get("abstract"):
            abstract = ris_escape(rec.get("abstract"))
            lines.append(f"N2  - {abstract[:1200]}")
        lines.append("DB  - PubMed")
        lines.append("ER  -")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    refs = parse_reference_lines()
    pubmed = fetch_pubmed([ref["pmid"] for ref in refs if ref["pmid"]])
    report_rows = []
    ordered_pubmed_records = []
    for ref in refs:
        rec = pubmed.get(ref["pmid"])
        if rec:
            ordered_pubmed_records.append(rec)
        crossref = fetch_crossref_by_doi(ref["doi"])
        time.sleep(0.15)
        manuscript_doi = normalize_doi(ref["doi"])
        pubmed_doi = normalize_doi(str(rec.get("doi", ""))) if rec else ""
        title_match = title_similarity(ref["title"], str(rec.get("title", ""))) if rec else 0.0
        doi_match = (not manuscript_doi) or (manuscript_doi == pubmed_doi)
        year_match = bool(rec and ref["year"] == str(rec.get("year", "")))
        crossref_similarity = title_similarity(ref["title"], crossref.get("title", "")) if crossref.get("title") else ""
        crossref_doi_match = (
            ""
            if not manuscript_doi or crossref.get("status") != "found"
            else str(manuscript_doi == normalize_doi(crossref.get("doi", "")))
        )
        manuscript_title_norm = normalize_title(ref["title"])
        pubmed_title_norm = normalize_title(str(rec.get("title", ""))) if rec else ""
        title_prefix_match = bool(manuscript_title_norm and pubmed_title_norm.startswith(manuscript_title_norm))
        if not rec:
            status = "not_found"
        elif (title_match >= 0.88 or title_prefix_match or (doi_match and title_match >= 0.65)) and year_match and doi_match:
            status = "verified"
        elif title_match >= 0.80 and year_match:
            status = "manual_check"
        else:
            status = "mismatch"
        report_rows.append(
            {
                "ref_no": ref["number"],
                "status": status,
                "pmid": ref["pmid"],
                "doi_in_manuscript": manuscript_doi,
                "doi_in_pubmed": pubmed_doi,
                "doi_match_pubmed": str(doi_match),
                "title_similarity_pubmed": f"{title_match:.3f}",
                "year_match_pubmed": str(year_match),
                "crossref_status": crossref.get("status", ""),
                "crossref_doi_match": crossref_doi_match,
                "title_similarity_crossref": f"{crossref_similarity:.3f}" if isinstance(crossref_similarity, float) else "",
                "manuscript_title": ref["title"],
                "pubmed_title": rec.get("title", "") if rec else "",
                "pubmed_journal": rec.get("journal", "") if rec else "",
                "pubmed_year": rec.get("year", "") if rec else "",
                "pubmed_authors_first_6": "; ".join(rec.get("authors", [])[:6]) if rec else "",
            }
        )
    report_csv = OUT_DIR / "JTM_references_verification_report.csv"
    with report_csv.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0].keys()))
        writer.writeheader()
        writer.writerows(report_rows)
    (OUT_DIR / "JTM_references_pubmed_metadata.json").write_text(
        json.dumps(ordered_pubmed_records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (OUT_DIR / "JTM_references_verified.ris").write_text(format_ris(ordered_pubmed_records), encoding="utf-8-sig")
    summary = {
        "total_references_in_manuscript": len(refs),
        "pubmed_records_retrieved": len(ordered_pubmed_records),
        "status_counts": {status: sum(row["status"] == status for row in report_rows) for status in sorted({row["status"] for row in report_rows})},
        "ris_records": len(ordered_pubmed_records),
        "report_csv": str(report_csv.relative_to(ROOT)),
        "ris_file": str((OUT_DIR / "JTM_references_verified.ris").relative_to(ROOT)),
    }
    (OUT_DIR / "JTM_references_verification_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
