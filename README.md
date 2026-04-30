# *Lineus longissimus* Functional Re-annotation

*Author: Anestis Gkanogiannis &lt;anestis@gkanogiannis.com&gt; - MNCN-CSIC, 2026*


Enhancement of the NCBI RefSeq annotation (`GCF_910592395.1-RS_2024_04`) of the bootlace worm
*Lineus longissimus* (Nemertea, taxID 88925), based on the Wellcome Sanger Institute
chromosome-level assembly **tnLinLong1.2** (RefSeq `GCF_910592395.1`).

The existing structural annotation covers 17,976 protein-coding genes (34,875 proteins), but
~25.9% of genes lack GO terms, ~35% are labelled "uncharacterized", and there is no KEGG
pathway mapping. This project layers additional functional evidence (GO, KEGG, Pfam/InterPro,
CAZymes, signal peptides, transmembrane topology, orthology) on top of the RefSeq proteins.

## Repository layout

```
assembly/        NCBI assembly report, stats, and genomic FASTA for tnLinLong1.2
annotation/      NCBI RefSeq GFF/GTF, CDS, RNA, protein FASTA, GO (GAF), expression counts
functional_annotation_osx.sh             Pipeline run on macOS (Apple Silicon)
functional_annotation_linux.sh           Full pipeline on Linux (conda)
functional_annotation_linux_singularity.sh   Same pipeline using Singularity containers
functional_annotation_osx_2026-03-31/    Outputs from the macOS run
REPORT_*.pdf                 Functional annotation report
annotation-improvement-plan.md           Tool selection and rationale
links.txt                                External data references
```

## Pipeline

**Input:** `annotation/GCF_910592395.1_tnLinLong1.2_protein.faa.gz` — 34,875 proteins from the
NCBI RefSeq release `GCF_910592395.1-RS_2024_04` (17,976 protein-coding genes).

**Final output:** a consolidated GFF3 layered on top of the RefSeq structural annotation,
adding GO terms, KEGG KO/pathway, Pfam/InterPro domains, CAZy families, signal-peptide and
transmembrane topology calls, and orthogroup IDs per gene.

### Tools and what they contribute

| # | Tool | Output | Notes |
|---|------|--------|-------|
| 1 | **eggNOG-mapper** v2.1 | GO, KEGG (KO/pathway/module), COG categories, preferred name, free-text description | tax_scope=Metazoa, evalue 1e-5, all GO evidence |
| 2 | **KofamScan** | KEGG KO assignment from HMM profiles (independent of eggNOG) | detail-tsv output, score-threshold filtering |
| 3 | **InterProScan** 5.77 | Pfam, SMART, PANTHER, CDD, SUPERFAMILY, MobiDB-lite, InterPro, GO | Linux only — JVM/native libs unavailable on ARM64 |
| 4 | **dbCAN** v4 | Carbohydrate-active enzyme families (GH, GT, PL, CE, AA, CBM) | Consensus from HMMER + DIAMOND + dbCAN-sub |
| 5 | **SignalP 6.0** | Secretory signal peptides (Sec/SPI, Sec/SPII, Tat/SPI) and cleavage site | `--organism eukarya` |
| 6 | **DeepTMHMM** v1.0 | Transmembrane topology (TM helices, signal peptides, inside/outside) | Linux + GPU recommended; ARM64 incompatible |
| 7 | **OrthoFinder** 3.1 + **DIAMOND** 2.1 | Orthogroups vs. 9 reference metazoans (Mollusca, Brachiopoda, Annelida, Echinodermata, Arthropoda, Chordata) | Used to transfer names to "uncharacterized" loci |

### Reference proteomes for orthology

Downloaded automatically via NCBI `datasets` CLI (see `functional_annotation_*.sh`):
*Crassostrea gigas*, *Lingula anatina*, *Homo sapiens*, *Drosophila melanogaster*,
*Strongylocentrotus purpuratus*, *Pecten maximus*, *Helobdella robusta*, *Lottia gigantea*,
*Octopus sinensis*.

### Pipeline variants

| Script | Host | Tool delivery | Coverage |
|--------|------|----------------|----------|
| `functional_annotation_osx.sh` | macOS, Apple Silicon (arm64), 16 GB / 10 CPU | conda + DeepTMHMM academic tarball | Skips InterProScan and DeepTMHMM (ARM64) |
| `functional_annotation_linux.sh` | Linux x86_64 | conda envs `funannotate`, `annotation`, `ortho` + native InterProScan | Full pipeline |
| `functional_annotation_linux_singularity.sh` | Linux x86_64 / HPC | Singularity/Apptainer images (BioContainers + custom defs) | Full pipeline, no Docker required |

Custom Singularity recipes for the two licensed tools are in this repo:

- `signalp6.def` — CPU-only SignalP 6.0 image (academic license tarball mounted at build time)
- `deeptmhmm.def` — GPU-enabled DeepTMHMM image (CUDA 11.2 / cuDNN 8, PyTorch 1.12+cu113)

### Required databases

Paths are configured at the top of each script. Approximate sizes:

| Database | Size | Used by |
|----------|------|---------|
| eggNOG v5.0 | ~44 GB | eggNOG-mapper |
| KEGG KOfam (profiles + ko_list) | ~7 GB | KofamScan |
| dbCAN v4 (HMM + DIAMOND + dbCAN-sub) | ~2 GB | dbCAN |
| InterProScan data (5.77-108.0) | ~45 GB | InterProScan (Linux) |

## Running

Edit the database paths and `CPUS` at the top of the chosen script, place the L. longissimus
protein FASTA at `annotation/GCF_910592395.1_tnLinLong1.2_protein.faa.gz` (already shipped),
then:

```bash
./functional_annotation_linux.sh              # conda: funannotate, annotation, ortho
./functional_annotation_linux_singularity.sh  # Singularity images (HPC-friendly)
./functional_annotation_osx.sh                # macOS arm64 subset (no IPS/DeepTMHMM)
```

Build the licensed-tool containers once before the Singularity run:

```bash
sudo -E singularity build signalp6.sif  signalp6.def
sudo -E singularity build deeptmhmm.sif deeptmhmm.def   # requires NVIDIA GPU at runtime
```

## Outputs

Per-tool results land under `functional_annotation_<host>_<date>/`:

```
eggnog/Llong.emapper.annotations          GO, KEGG, COG, descriptions (TSV)
kofamscan/Llong_kofam.txt                 KEGG KO calls (detail-tsv)
interproscan/Llong.tsv (+ .gff3, .xml)    Pfam/InterPro/PANTHER/MobiDB hits
dbcan/                                    overview.txt, CAZyme family calls
signalp/prediction_results.txt            signal peptide calls + cleavage sites
deeptmhmm/predicted_topologies.3line      TM topology (3-line format) + GFF3
orthofinder/Results_*/Orthogroups/        orthogroups TSV, single-copy orthologs
                                          and species-tree
Llong_functional.gff3                     consolidated annotation (final product)
```

See `REPORT_functional_annotation_L.longissimus.pdf` for the consolidated
summary, statistics (coverage gain in GO/KEGG/Pfam, "uncharacterized" rescue rate), and
figures.
