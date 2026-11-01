# Project proposal: detecting ancestry from unsampled donor lineages

**Status:** proposal for written course approval; approval has not been recorded.
**November milestone:** minimal presentable version and written approval by November 14.
**Final target:** December 7 readiness; the assignment states December 8 submission.

## Question and motivation

How do ArchIE and TRACE differ in their ability to recover DNA inherited from an
unsampled donor population as admixture becomes older or rarer, and when their
demographic or reference-panel assumptions are violated? A secondary question is
how much TRACE's performance depends on reconstructing ancestral recombination
graphs rather than observing their true simulated histories.

Human introgression studies motivate these approaches. Durvasula and Sankararaman
(2020) used population-genetic summaries and ArchIE to map candidate archaic
segments without a donor genome. Zhang and colleagues (2026) introduced TRACE,
which uses features of inferred local genealogies. Their information requirements
differ: ArchIE uses a comparison population and simulation-trained parameters;
TRACE uses inferred genealogies. This project evaluates adapted software, not the
claim that unsampled ancestry automatically identifies an extinct species.

## Data and two approaches

The main dataset will consist of reproducible msprime simulations with donor
migration records. Ten diploid target and ten diploid reference individuals will
be simulated over independent 1-Mb regions. Neither detector receives donor
samples or test ancestry labels. Migration records provide exact truth for
evaluation, including mixed regions. Training, calibration and test seeds are
separate. The controlled model is not a fitted mosquito population history.

ArchIE will retain its independently trained binomial classifier. TRACE will use
SINGER-inferred ARGs; a separately labeled true-ARG analysis isolates reconstruction
losses. Existing source adaptations and upstream commit identifiers are documented.
An initial real software run has completed, but its short MCMC chain and one seed
establish functionality rather than comparative accuracy.

## Tests and expected interpretation

Primary scoring intersects predicted and true intervals across all callable bases.
Six independent baseline validation regions select thresholds targeting FDR at
most 5%, then thresholds are frozen. Test FDR may exceed that target under other
histories. If no nonempty validation call set qualifies, that absence is a result,
not a successful FDR guarantee.

Five to ten independent replicates will cover no introgression, different
introgression ages, lower ancestry proportion, and an alternative ancestral history.
The latter includes migration and is interpreted as an alternative-history call
rate rather than a pure false-positive control. Optional stress tests add reference
admixture and demographic mismatch. Confidence intervals resample independent
simulation replicates, retaining pairing between methods.

We will compare recall, achieved FDR, false-called sequence in the null model,
resource cost, and systematically selected localization errors. Null findings and
uncertain comparisons are valid outcomes. A read-only local explorer will link
performance plots to ancestry tracks; static figures and data remain reproducible.

## Optional mosquito demonstration and deliverables

Ag1000G provides public phased mosquito genomes and metadata. An empirical pilot
is conditional on defensible population parameters, QC/relatedness evidence, and
ancestral allele orientation. Reference alleles will not be treated as ancestral
alleles. Method agreement would be an empirical consistency check, not accuracy
or proof of a new ghost lineage. Biological observations are exploratory only.

Final deliverables follow the assignment: about two introduction pages, 2-3
results pages excluding figures, and 1-2 discussion pages, with 12-point font and
1.5 spacing; three principal figures; and annotated executed code and environment
records. Local experiments are capped at twelve elapsed hours with four workers.

## References

Durvasula A, Sankararaman S. Recovering signals of ghost archaic introgression in
African populations. Science Advances (2020). DOI: 10.1126/sciadv.aax5097.

Zhang Y et al. Recovering signatures of archaic hominin introgression using ancestral
recombination graphs. Science (2026). DOI: 10.1126/science.aef8874.

Software: https://github.com/sriramlab/ArchIE ; https://github.com/YulinZhang9806/trace ;
https://github.com/popgenmethods/SINGER . Data: https://www.malariagen.net/data_package/ag1000g-phase3-hap/ .
