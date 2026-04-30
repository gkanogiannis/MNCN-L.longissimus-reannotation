#!/usr/bin/env bash
# =============================================================================
# Functional Annotation Pipeline for Lineus longissimus (OSX version with conda)
# =============================================================================
# Author: Anestis Gkanogiannis <anestis@gkanogiannis.com> - MNCN-CSIC, 2026
# =============================================================================
# This script was run on macOS (Apple Silicon arm64, 16 GB RAM, 10 CPUs)
# Date: 2026-03-31
#
# Prerequisites:
#   - conda/mamba with the following environments:
#     * funannotate: eggnog-mapper, kofamscan, dbcan, hmmer, diamond
#     * annotation:  signalp6, DeepTMHMM (academic license)
#     * ortho:       orthofinder, diamond
#   - Databases at /Volumes/SSD2TB/databases/:
#     * eggnog/ (eggNOG v5.0, ~44 GB)
#     * kofam/  (KEGG HMM profiles, ~7 GB)
#     * dbcan/  (dbCAN, ~2 GB)
#
# Note: InterProScan and DeepTMHMM failed on this machine due to ARM64
# compatibility issues. See functional_annotation_linux.sh for the full
# pipeline including those tools.
# =============================================================================

set -euo pipefail

# --- Configuration ---
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
PROTFILE="${WORKDIR}/annotation/GCF_910592395.1_tnLinLong1.2_protein.faa.gz"
OUTDIR="${WORKDIR}/functional_annotation_osx_2026-03-31"
CPUS=10

DB_EGGNOG="/Volumes/SSD2TB/databases/eggnog"
DB_KOFAM="/Volumes/SSD2TB/databases/kofam"
DB_DBCAN="/Volumes/SSD2TB/databases/dbcan"

# --- Setup ---
mkdir -p "${OUTDIR}"/{eggnog,kofamscan,signalp,deeptmhmm,dbcan,orthofinder/proteomes}

# Decompress proteins
if [ ! -f "${OUTDIR}/proteins.faa" ]; then
    echo "=== Decompressing protein sequences ==="
    gunzip -c "${PROTFILE}" > "${OUTDIR}/proteins.faa"
fi

PROTEINS="${OUTDIR}/proteins.faa"

# =============================================================================
# 1. eggNOG-mapper (GO, KEGG, COG, descriptions)
# =============================================================================
echo "=== Running eggNOG-mapper ==="
conda run -n funannotate emapper.py \
    -i "${PROTEINS}" \
    --output "${OUTDIR}/eggnog/Llong" \
    --data_dir "${DB_EGGNOG}" \
    --cpu ${CPUS} \
    --itype proteins \
    --tax_scope Metazoa \
    --go_evidence all \
    --target_orthologs all \
    --seed_ortholog_evalue 1e-5 \
    --override

# =============================================================================
# 2. KofamScan (KEGG KO assignment via HMM profiles)
# =============================================================================
echo "=== Running KofamScan ==="
conda run -n funannotate exec_annotation \
    -f detail-tsv \
    -o "${OUTDIR}/kofamscan/Llong_kofam.txt" \
    --profile "${DB_KOFAM}/profiles" \
    --ko-list "${DB_KOFAM}/ko_list" \
    --cpu ${CPUS} \
    "${PROTEINS}"

# =============================================================================
# 3. SignalP 6.0 (signal peptide prediction)
# =============================================================================
echo "=== Running SignalP 6.0 ==="
conda run -n annotation signalp6 \
    --fastafile "${PROTEINS}" \
    --output_dir "${OUTDIR}/signalp/" \
    --organism eukarya \
    --format txt

# =============================================================================
# 4. DeepTMHMM (transmembrane topology prediction)
# NOTE: Very slow on CPU (ARM64 Mac). May take 12+ hours.
#       Consider running on Linux with GPU instead.
# =============================================================================
echo "=== Running DeepTMHMM ==="
DEEPTMHMM_DIR="${HOME}/Software-Bioinfo/DeepTMHMM-Academic-License-v1.0"
ABSFASTA="$(cd "$(dirname "${PROTEINS}")" && pwd)/$(basename "${PROTEINS}")"
ABSOUT="$(cd "$(dirname "${OUTDIR}/deeptmhmm")" && pwd)/deeptmhmm"
(cd "${DEEPTMHMM_DIR}" && conda run -n annotation python3 predict.py \
    --fasta "${ABSFASTA}" \
    --output-dir "${ABSOUT}") || echo "WARNING: DeepTMHMM failed (ARM64 compatibility)"

# =============================================================================
# 5. dbCAN (CAZyme annotation)
# =============================================================================
echo "=== Running dbCAN ==="
conda run -n funannotate run_dbcan CAZyme_annotation \
    --mode protein \
    --input_raw_data "${PROTEINS}" \
    --output_dir "${OUTDIR}/dbcan/" \
    --db_dir "${DB_DBCAN}" \
    --threads ${CPUS}

# =============================================================================
# 6. OrthoFinder (orthology inference across 10 metazoan species)
# =============================================================================
echo "=== Running OrthoFinder ==="
ORTHO_DIR="${OUTDIR}/orthofinder/proteomes"

# Copy L. longissimus proteome
cp "${PROTEINS}" "${ORTHO_DIR}/Lineus_longissimus.faa"

# Download reference proteomes (NCBI datasets CLI required)
declare -A SPECIES_ACC=(
    ["Crassostrea_gigas"]="GCF_963853765.1"
    ["Lingula_anatina"]="GCF_001039355.2"
    ["Homo_sapiens"]="GCF_000001405.40"
    ["Drosophila_melanogaster"]="GCF_000001215.4"
    ["Strongylocentrotus_purpuratus"]="GCF_000002235.5"
    ["Pecten_maximus"]="GCF_902652985.1"
    ["Helobdella_robusta"]="GCF_000326865.1"
    ["Lottia_gigantea"]="GCF_000327385.1"
    ["Octopus_sinensis"]="GCF_006345805.1"
)

for species in "${!SPECIES_ACC[@]}"; do
    acc="${SPECIES_ACC[$species]}"
    if [ ! -f "${ORTHO_DIR}/${species}.faa" ]; then
        echo "Downloading ${species} (${acc})..."
        datasets download genome accession "${acc}" --include protein --filename "/tmp/${species}.zip" 2>/dev/null && \
        unzip -p "/tmp/${species}.zip" "*/protein.faa" > "${ORTHO_DIR}/${species}.faa" && \
        rm -f "/tmp/${species}.zip" || \
        echo "WARNING: Failed to download ${species}. Try manual download from NCBI FTP."
    fi
done

conda run -n ortho orthofinder \
    -f "${ORTHO_DIR}/" \
    -t ${CPUS} \
    -a 5 \
    -S diamond

# =============================================================================
# 7. Generate consolidated GFF3
# =============================================================================
echo "=== Generating consolidated functional annotation GFF3 ==="
python3 << 'PYEOF'
import gzip, re

protein_to_gene = {}
gene_to_mrna = {}
gene_coords = {}

with gzip.open("annotation/GCF_910592395.1_tnLinLong1.2_genomic.gff.gz", 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        feat_type, attrs = parts[2], parts[8]
        if feat_type == 'gene':
            m = re.search(r'ID=([^;]+)', attrs)
            if m: gene_coords[m.group(1)] = (parts[0], int(parts[3]), int(parts[4]), parts[6])
        elif feat_type == 'mRNA':
            m_id, m_par = re.search(r'ID=([^;]+)', attrs), re.search(r'Parent=([^;]+)', attrs)
            if m_id and m_par: gene_to_mrna[m_id.group(1)] = m_par.group(1)
        elif feat_type == 'CDS':
            m_pid, m_par = re.search(r'protein_id=([^;]+)', attrs), re.search(r'Parent=([^;]+)', attrs)
            if m_pid and m_par and m_pid.group(1) not in protein_to_gene:
                protein_to_gene[m_pid.group(1)] = gene_to_mrna.get(m_par.group(1), m_par.group(1))

import os
OUTDIR = os.environ.get('OUTDIR', 'functional_annotation_osx_2026-03-31')

eggnog = {}
with open(f"{OUTDIR}/eggnog/Llong.emapper.annotations") as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.strip().split('\t')
        if len(cols) < 21: continue
        eggnog[cols[0]] = {k: v for k, v in zip(['go','kegg_ko','kegg_pw','cog','desc'], [cols[9],cols[11],cols[12],cols[6],cols[7]]) if v != '-'}

kofam = {}
with open(f"{OUTDIR}/kofamscan/Llong_kofam.txt") as f:
    for line in f:
        if line.startswith('*'):
            p = line.strip().split()
            kofam.setdefault(p[1], []).append(p[2])

signalp = set()
with open(f"{OUTDIR}/signalp/output.gff3") as f:
    for line in f:
        if not line.startswith('#'):
            p = line.strip().split('\t')
            if len(p) >= 3 and p[2] == 'signal_peptide': signalp.add(p[0])

dbcan = {}
with open(f"{OUTDIR}/dbcan/overview.tsv") as f:
    next(f)
    for line in f:
        c = line.strip().split('\t')
        if len(c) >= 7 and int(c[5]) >= 2: dbcan[c[0]] = c[6]

all_pids = set(eggnog) | set(kofam) | signalp | set(dbcan)
with open(f"{OUTDIR}/Llong_functional_annotation.gff3", 'w') as out:
    out.write("##gff-version 3\n# Consolidated functional annotation for Lineus longissimus\n")
    for pid in sorted(all_pids):
        gid = protein_to_gene.get(pid)
        if not gid or gid not in gene_coords: continue
        ch, s, e, st = gene_coords[gid]
        a = [f"ID={pid}", f"gene_id={gid}"]
        eg = eggnog.get(pid, {})
        if 'go' in eg: a.append(f"Ontology_term={eg['go']}")
        if 'kegg_ko' in eg: a.append(f"kegg_ko={eg['kegg_ko']}")
        if 'kegg_pw' in eg: a.append(f"kegg_pathway={eg['kegg_pw']}")
        if 'cog' in eg: a.append(f"cog_category={eg['cog']}")
        if 'desc' in eg: a.append(f"eggnog_description={eg['desc'].replace(';','%3B').replace('=','%3D')}")
        if pid in kofam: a.append(f"kofamscan_ko={','.join(kofam[pid])}")
        if pid in signalp: a.append("signal_peptide=yes")
        if pid in dbcan: a.append(f"cazyme_family={dbcan[pid]}")
        out.write(f"{ch}\tfunctional_annotation\tprotein_match\t{s}\t{e}\t.\t{st}\t.\t{';'.join(a)}\n")
PYEOF

echo "=== Pipeline complete ==="
echo "Results in: ${OUTDIR}/"
