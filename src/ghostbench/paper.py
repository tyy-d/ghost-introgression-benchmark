"""Evidence-driven course report and approval-proposal PDF export.

No result is invented to fill a report section. Diagnostic observations and
unresolved convergence remain labeled, and the November draft is not a claim
of project approval or final submission readiness.
"""

import json
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

from .analysis import bootstrap_summary, collect
from .figures import (
    export_study,
    figure_bytes,
    static_cost,
    static_performance,
    static_tracks,
    track_data,
)
from .registry import ROOT, validate_complete

INTRODUCTION = [
    [
        "Recombination makes a genome a mosaic of histories rather than a single family tree. "
        "When formerly separated populations interbreed, descendants can inherit segments whose "
        "ancestry differs from that of neighboring sequence. Introgression from an unsampled donor "
        "is particularly difficult to detect because no donor genome is available for direct matching. "
        "The word ghost describes that absence of a sampled reference. It does not establish that "
        "the donor was extinct, a separate species, or identifiable from a fossil.",
        "Two computational approaches can use different aspects of the same descendant sequences. "
        "ArchIE combines population-genetic summaries in a simulation-trained binomial classifier. "
        "Its features include allele frequencies, distances among haplotypes, divergence from a "
        "comparison population, and patterns of private derived variation. Training links those "
        "features to known simulated ancestry. Consequently, the fit between training and test "
        "histories is part of the scientific question, not merely an installation detail.",
        "TRACE instead examines genealogical information in ancestral recombination graphs. "
        "Deep branches that persist across linked sequence can support an introgression interpretation. "
        "The method combines information along the chromosome with a hidden Markov model. In this "
        "project, SINGER reconstructs graphs from the observed genomes before TRACE runs. These "
        "graphs are estimates with uncertainty, so a detector can lose information before ancestry "
        "classification begins. A true-graph diagnostic separates that reconstruction step from "
        "the detector itself.",
        "Durvasula and Sankararaman (2020) used ArchIE alongside a separate allele-frequency "
        "analysis to investigate archaic ancestry in African populations. Zhang and colleagues "
        "(2026) introduced TRACE using inferred genealogies. The present comparison concerns "
        "ancestry-segment detection with adapted software. It does not reproduce the 2020 paper's "
        "separate demographic-model analysis, and it does not presume either method will be superior."
    ],
    [
        "The precise analysis question is how sensitivity, achieved false discovery rate, and "
        "computational cost change when donor ancestry becomes older or rarer and when assumptions "
        "about a comparison population or demographic parameters fail. A secondary question asks "
        "how much TRACE changes when it receives inferred rather than true graphs. The primary "
        "comparison freezes thresholds using independent validation data, then evaluates paired "
        "methods on exactly the same held-out sequences.",
        "The main dataset is generated with msprime. Each independent region contains ten diploid "
        "target and ten diploid reference individuals across one million bases. The donor has no "
        "sampled genomes. Migration records identify ancestry along the simulated chromosomes, "
        "including intervals ancestral to multiple descendants. Phased biallelic sequence data and "
        "sample order are exported separately from the true graphs and ancestry labels. Training, "
        "validation and test seeds occupy disjoint ranges.",
        "The controlled baseline uses effective size 100,000, mutation and recombination rates of "
        "1e-8 per base per generation, donor divergence 200,000 generations ago, target/reference "
        "divergence 50,000 generations ago, and five percent admixture 2,000 generations ago. "
        "These are controlled software-study parameters, not estimates of a mosquito history. "
        "The data contain neutral simulated sequence and known migration histories; they do not "
        "contain selection, functional annotations, real sequencing errors, or direct evidence "
        "of a particular extinct organism.",
        "Public Ag1000G genomes motivate a conditional nonhuman demonstration, but empirical "
        "ancestry has no equivalent exact truth labels. Documented sample quality, relatedness, "
        "ancestral allele orientation and biological parameters are prerequisites. In particular, "
        "a reference allele is not automatically ancestral. The simulation study can quantify "
        "accuracy; empirical method agreement can only establish consistency. The assignment "
        "permits simulation projects, while approval of this specific proposal remains a separate "
        "course requirement."
    ]
]


def pdf_document(path, sections, title):
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        Image,
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="CourseBody", fontName="Times-Roman", fontSize=12, leading=18,
                             alignment=TA_LEFT, spaceAfter=9, splitLongWords=True))
    styles.add(ParagraphStyle(name="CourseTitle", fontName="Times-Bold", fontSize=16, leading=21, spaceAfter=14))
    styles.add(ParagraphStyle(name="CourseHeading", fontName="Times-Bold", fontSize=13, leading=18, spaceAfter=10, keepWithNext=True))
    styles.add(ParagraphStyle(name="CourseCaption", fontName="Times-Roman", fontSize=12, leading=18, spaceAfter=10))
    body = []
    for index, section in enumerate(sections):
        if index:
            body.append(PageBreak())
        body.append(Paragraph(escape(section["heading"]), styles["CourseTitle"] if index == 0 else styles["CourseHeading"]))
        for text in section.get("paragraphs", []):
            proposal_heading = text in {"Question and motivation", "Data and two approaches",
                                        "Tests and expected interpretation",
                                        "Optional mosquito demonstration and deliverables", "References"}
            style = styles["CourseHeading"] if proposal_heading else styles["CourseBody"]
            body.append(Paragraph(escape(text.replace("**", "")), style))
        if "table" in section:
            rows = [[Paragraph(escape(str(cell)), styles["CourseCaption"]) for cell in row] for row in section["table"]]
            table = Table(rows, colWidths=[1.4*inch, 1.2*inch, 1.3*inch, 1.3*inch, 1.3*inch][:len(rows[0])], repeatRows=1)
            table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
                                      ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                                      ("VALIGN", (0, 0), (-1, -1), "TOP"),
                                      ("BOTTOMPADDING", (0, 0), (-1, -1), 7)]))
            body.append(table)
            body.append(Spacer(1, 12))
        for filename, caption in section.get("figures", []):
            from PIL import Image as PILImage
            width, height = PILImage.open(filename).size
            body.append(Image(str(filename), width=6.5*inch, height=6.5*inch*height/width))
            body.append(Paragraph(escape(caption), styles["CourseCaption"]))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def footer(canvas, doc):
        canvas.setFont("Times-Roman", 9)
        canvas.drawString(inch, 0.55*inch, "BIOCB project | Draft evidence, not an approval record")
        canvas.drawRightString(7.5*inch, 0.55*inch, str(doc.page))

    SimpleDocTemplate(str(path), pagesize=(8.5*inch, 11*inch), leftMargin=inch, rightMargin=inch,
                      topMargin=0.8*inch, bottomMargin=0.8*inch, title=title).build(body, onFirstPage=footer, onLaterPages=footer)


def write_paper(root, destination):
    import matplotlib.pyplot as plt

    root, destination = Path(root), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    frame, summary, differences, runs = export_study(root, destination)
    _, all_rows = collect(root)
    diagnostic = all_rows[all_rows.role == "diagnostic"] if len(all_rows) else pd.DataFrame()
    representative = frame if len(frame) else diagnostic
    diagnostic_only = frame.empty
    if diagnostic_only and len(diagnostic):
        label = "DIAGNOSTIC ONLY: unfrozen threshold, single-seed functionality evidence"
        diagnostic_summary = bootstrap_summary(diagnostic)
        diagnostic_summary.to_csv(destination / "performance.tsv", sep="\t", index=False)
        diagnostic.to_csv(destination / "replicate_metrics.tsv", sep="\t", index=False)
        (destination / "selection.json").write_text(json.dumps({"threshold_status": "diagnostic-unfrozen", "role": "diagnostic",
            "preview": False, "run_fingerprints": sorted(diagnostic.fingerprint.unique().tolist())}, indent=2))
        fig = static_performance(diagnostic, diagnostic_summary)
        fig.suptitle(label, fontsize=9)
        for format in ["png", "svg", "pdf"]:
            (destination / f"performance.{format}").write_bytes(figure_bytes(fig, format))
        plt.close(fig)
        fig = static_cost(diagnostic)
        fig.suptitle(label, fontsize=9)
        for format in ["png", "svg", "pdf"]:
            (destination / f"cost_robustness.{format}").write_bytes(figure_bytes(fig, format))
        plt.close(fig)
        choice = diagnostic.sort_values(["scenario", "replicate", "method"]).iloc[0]
        data, _ = validate_complete(choice.path)
        metadata = pd.read_csv(data / "inputs/samples.tsv", sep="\t")
        haplotype = metadata[metadata.population == "target"].haplotype.iloc[0]
        args = track_data(choice.path, haplotype)
        args[0].to_csv(destination / "localization_data.tsv", sep="\t", index=False)
        (destination / "localization_selection.json").write_text(json.dumps(args[3], indent=2))
        fig = static_tracks(*args)
        fig.suptitle(label, fontsize=8)
        for format in ["png", "svg", "pdf"]:
            (destination / f"localization.{format}").write_bytes(figure_bytes(fig, format))
        plt.close(fig)
    if len(representative):
        table = [["Method", "Replicates", "Recall", "Achieved FDR", "Called kb"]]
        for method, group in representative.groupby("method", sort=True):
            from .analysis import pooled
            result = pooled(group)
            table.append([method, group[["scenario", "replicate"]].drop_duplicates().shape[0],
                          f"{result['recall']:.3f}" if result["recall"] is not None else "undefined",
                          f"{result['fdr']:.3f}" if result["fdr"] is not None else "undefined",
                          f"{result['called_bases']/1000:.1f}"])
    else:
        table = [["Method", "Status"], ["ArchIE", "No eligible result"], ["TRACE", "No eligible result"]]
    pilot = root / "pilot/pilot.json"
    pilot_record = json.loads(pilot.read_text()) if pilot.exists() else None
    pilot_status = pilot_record.get("status", "unresolved") if pilot_record else "not run"
    if pilot_record is not None and not pilot_record.get("finished"):
        pilot_status += " (running)"
    status = ("Only diagnostic results are available. All numerical values below use an unfrozen 0.5 threshold "
              "and do not estimate held-out calibrated performance." if diagnostic_only else
              f"There are {frame[['scenario', 'replicate']].drop_duplicates().shape[0]} usable paired test replicates. "
              "Thresholds were selected on independent validation regions and frozen before testing.")
    sections = [
        {"heading": "Detecting ancestry from unsampled lineages: introduction", "paragraphs": INTRODUCTION[0]},
        {"heading": "Introduction: question and data", "paragraphs": INTRODUCTION[1]},
        {"heading": "Results: evidence available and validity", "paragraphs": [
            status,
            f"The convergence pilot status is {pilot_status}. Two independent chains are screened using "
            "split R-hat and effective sample size for diversity fit, incompatible mutation sites, and "
            "TRACE ancestry amount on a predeclared focal haplotype. Constant statistics alone cannot "
            "establish mixing. Passing this pilot supports a choice of settings; it is not proof of "
            "convergence in every test history.",
            "The primary table pools exact callable-base counts. It must be read together with the "
            "per-condition tables and paired bootstrap results: pooling different conditions can "
            "conceal failure at low ancestry prevalence or old admixture. Undefined FDR or precision "
            "means there are no called bases, not perfect control of false discoveries. True-ARG "
            "TRACE is a diagnostic and is never a deployable empirical alternative."], "table": table},
        {"heading": "Results: comparisons, tests and resource cost", "paragraphs": [
            "The primary endpoint is recall at thresholds targeting estimated validation FDR no "
            "greater than five percent. Six independent one-megabase validation regions calibrate "
            "each method. If no nonempty call set qualifies, that method abstains at this operating "
            "target. The test result is achieved FDR alongside recall; a validation target cannot "
            "guarantee the same error rate under changed demographic histories.",
            "Uncertainty resamples independent simulation replicates and preserves method pairing. "
            "Five or more usable paired replicates are required for confidence intervals; smaller "
            "sets remain descriptive. Failed, incomplete and corrupted outputs are counted in the "
            "run-status table rather than imputed as zero ancestry. An alternative-history scenario "
            "contains ancestral migration, so its calls are not interpreted as ordinary false positives.",
            "The first principal figure summarizes recall and achieved FDR. The second shows an "
            "ancestry-localization example chosen by a deterministic largest-error rule for eligible "
            "test data; a diagnostic-only draft instead uses the first target haplotype explicitly. "
            "The third compares inference cost and available robustness conditions. ArchIE training "
            "is recorded once and excluded from its per-region inference time; SINGER reconstruction "
            "is included in TRACE's total. CPU and memory measurements include availability notes "
            "because the local sandbox may restrict process-tree sampling.",
            "Correctness checks cover exact mixed-base overlap, masks and adjacent intervals, null "
            "data, unattainable operating points, stale hashes, validation/test separation, paired "
            "resampling, and process-group cancellation. These checks verify scientific bookkeeping; "
            "they do not establish that the demographic model captures a natural population."]},
        {"heading": "Discussion: strengths, limitations and unanswered questions", "paragraphs": [
            "ArchIE turns interpretable population summaries into inexpensive repeated predictions, "
            "but its training distribution, comparison population and fixed window resolution can "
            "limit transfer. TRACE uses richer genealogical information and finer output intervals, "
            "but reconstructing those histories costs time and introduces uncertainty before detection. "
            "The true-graph diagnostic is useful precisely because it separates these two limitations.",
            "The controlled simulation supplies exact ancestry and excludes several natural processes. "
            "A high score can show success under specified assumptions, not species-wide validity. "
            "Reference admixture and demographic mismatch stress tests probe only two departures. "
            "A short chain, few independent regions, or no qualifying FDR threshold can make the "
            "answer inconclusive. Those are reportable limits rather than reasons to alter thresholds "
            "after seeing test labels or discard difficult regions.",
            "The present deliverable is a research draft for the November milestone. Written course "
            "approval has not been recorded. Final readiness is targeted for December 7, ahead of "
            "the assignment's December 8 deadline. Empirical mosquito inference remains conditional "
            "on QC, unrelated samples, polarization and documented parameters. Agreement would "
            "not identify an extinct donor or supply an accuracy estimate.",
            "References: Durvasula and Sankararaman (2020), Science Advances, DOI 10.1126/sciadv.aax5097. "
            "Zhang et al. (2026), Science, DOI 10.1126/science.aef8874. Software and full source/version "
            "records: github.com/sriramlab/ArchIE; github.com/YulinZhang9806/trace; "
            "github.com/popgenmethods/SINGER; tskit.dev/msprime. Public empirical resource: malariagen.net."]},
        {"heading": "Figures: evidence and reproducible exports", "figures": [
            (destination / "performance.png", "Figure 1. Recall and achieved FDR. Diagnostic-only displays are explicitly labeled; research intervals resample independent replicates."),
            (destination / "localization.png", "Figure 2. True and called ancestry at matching coordinates. Black outlines denote false positives and purple outlines false negatives; inaccessible sequence is excluded."),
            (destination / "cost_robustness.png", "Figure 3. Inference time and available robustness comparisons. Missing conditions are not implied to have been run.")]}]
    pdf_document(destination / "research-report-draft.pdf", sections, "Ghost introgression methods study - draft")
    markdown = ["# Research report draft", "", "Status: provisional; see run eligibility and convergence evidence. Not an approval record.", ""]
    for section in sections:
        markdown.extend(["## "+section["heading"], "", *[paragraph+"\n" for paragraph in section.get("paragraphs", [])]])
    (destination / "research-report-draft.md").write_text("\n".join(markdown))
    source = ROOT / "reports/approval-proposal.md"
    proposal_sections = []
    for part in source.read_text().split("\n## "):
        lines = part.splitlines()
        heading = lines[0].lstrip("# ")
        paragraphs = [" ".join(piece.splitlines()) for piece in "\n".join(lines[1:]).strip().split("\n\n") if piece.strip()]
        # Proposal headings flow naturally rather than forcing a page per heading.
        proposal_sections.append({"heading": heading, "paragraphs": paragraphs})
    proposal = {"heading": proposal_sections[0]["heading"], "paragraphs": []}
    for section in proposal_sections:
        if section != proposal_sections[0]:
            proposal["paragraphs"].append(section["heading"])
        proposal["paragraphs"].extend(section["paragraphs"])
    pdf_document(destination / "approval-proposal.pdf", [proposal], "BIOCB approval proposal")
    return destination
