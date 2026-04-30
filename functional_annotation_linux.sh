#!/usr/bin/env bash
# =============================================================================
# Functional Annotation Pipeline for Lineus longissimus (Linux version with conda)
# =============================================================================
# Author: Anestis Gkanogiannis <anestis@gkanogiannis.com> - MNCN-CSIC, 2026
# =============================================================================
# Full pipeline including tools that failed on macOS ARM64:
#   - InterProScan (native, all analyses including PANTHER and MobiDB)
#   - DeepTMHMM (GPU-accelerated)
#
# Prerequisites:
#   - conda/mamba
#   - Java 11+ (for InterProScan)
#   - NVIDIA GPU + CUDA (for DeepTMHMM)
#   - NCBI datasets CLI (conda install -c conda-forge ncbi-datasets-cli)
#
# Installation (run once):
#   mamba create -n funannotate -y -c bioconda -c conda-forge python=3.11 \
#     eggnog-mapper hmmer diamond kofamscan dbcan interproscan
#   mamba create -n ortho -y -c bioconda -c conda-forge orthofinder diamond
#   mamba create -n annotation -y python=3.9 pytorch torchvision \
#     -c pytorch -c nvidia
#   # Then install SignalP 6.0 and DeepTMHMM from DTU academic licenses
#
# Database setup (run once):
#   # eggNOG (~44 GB)
#   conda run -n funannotate download_eggnog_data.py --data_dir ${DB_BASE}/eggnog -y
#   # KofamScan (~7 GB)
#   mkdir -p ${DB_BASE}/kofam && cd ${DB_BASE}/kofam
#   wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
#   wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
#   tar xzf profiles.tar.gz && gunzip ko_list.gz
#   # dbCAN (~2 GB)
#   mkdir -p ${DB_BASE}/dbcan && cd ${DB_BASE}/dbcan
#   conda run -n funannotate run_dbcan database --db_dir .
#   # InterProScan databases (auto-downloaded with conda install, ~15 GB)
# =============================================================================

set -euo pipefail

# --- Configuration ---
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
PROTFILE="${WORKDIR}/annotation/GCF_910592395.1_tnLinLong1.2_protein.faa.gz"
OUTDIR="${WORKDIR}/functional_annotation_linux"
CPUS=$(nproc)

# Adjust database paths for your Linux machine
DB_BASE="/data/databases"  # <-- CHANGE THIS
DB_EGGNOG="${DB_BASE}/eggnog"
DB_KOFAM="${DB_BASE}/kofam"
DB_DBCAN="${DB_BASE}/dbcan"

# DeepTMHMM path -- adjust for your installation
DEEPTMHMM_DIR="${HOME}/Software-Bioinfo/DeepTMHMM-Academic-License-v1.0"

# --- Setup ---
mkdir -p "${OUTDIR}"/{eggnog,interproscan,kofamscan,signalp,deeptmhmm,dbcan,orthofinder/proteomes}

# Decompress proteins
if [ ! -f "${OUTDIR}/proteins.faa" ]; then
    echo "=== Decompressing protein sequences ==="
    gunzip -c "${PROTFILE}" > "${OUTDIR}/proteins.faa"
fi

PROTEINS="${OUTDIR}/proteins.faa"

# =============================================================================
# 1. eggNOG-mapper (GO, KEGG, COG, descriptions)
# =============================================================================
echo "=== [1/7] Running eggNOG-mapper ==="
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
echo "=== [2/7] Running KofamScan ==="
conda run -n funannotate exec_annotation \
    -f detail-tsv \
    -o "${OUTDIR}/kofamscan/Llong_kofam.txt" \
    --profile "${DB_KOFAM}/profiles" \
    --ko-list "${DB_KOFAM}/ko_list" \
    --cpu ${CPUS} \
    "${PROTEINS}"

# =============================================================================
# 3. InterProScan (domains, GO, pathways -- FULL analysis)
#    This is the key addition vs. the OSX pipeline.
#    Includes PANTHER, MobiDB, Pfam, CDD, SMART, SUPERFAMILY, etc.
# =============================================================================
echo "=== [3/7] Running InterProScan ==="
conda run -n funannotate interproscan.sh \
    -i "${PROTEINS}" \
    -b "${OUTDIR}/interproscan/Llong_interpro" \
    -f TSV,GFF3,XML \
    -goterms -pa \
    -cpu ${CPUS} \
    -dp

# =============================================================================
# 4. SignalP 6.0 (signal peptide prediction)
# =============================================================================
echo "=== [4/7] Running SignalP 6.0 ==="
conda run -n annotation signalp6 \
    --fastafile "${PROTEINS}" \
    --output_dir "${OUTDIR}/signalp/" \
    --organism eukarya \
    --format txt

# =============================================================================
# 5. DeepTMHMM (transmembrane topology prediction -- GPU accelerated)
#    On a GPU machine this should complete in ~30-60 min vs 12+ hours on CPU.
# =============================================================================
echo "=== [5/7] Running DeepTMHMM (GPU) ==="
ABSFASTA="$(realpath "${PROTEINS}")"
ABSOUT="$(realpath "${OUTDIR}/deeptmhmm")"
rmdir "${ABSOUT}" 2>/dev/null || true  # DeepTMHMM requires non-existing output dir
(cd "${DEEPTMHMM_DIR}" && conda run -n annotation python3 predict.py \
    --fasta "${ABSFASTA}" \
    --output-dir "${ABSOUT}")

# =============================================================================
# 6. dbCAN (CAZyme annotation)
# =============================================================================
echo "=== [6/7] Running dbCAN ==="
conda run -n funannotate run_dbcan CAZyme_annotation \
    --mode protein \
    --input_raw_data "${PROTEINS}" \
    --output_dir "${OUTDIR}/dbcan/" \
    --db_dir "${DB_DBCAN}" \
    --threads ${CPUS}

# =============================================================================
# 7. OrthoFinder (orthology inference across 10 metazoan species)
# =============================================================================
echo "=== [7/7] Running OrthoFinder ==="
ORTHO_DIR="${OUTDIR}/orthofinder/proteomes"

# Copy L. longissimus proteome
cp "${PROTEINS}" "${ORTHO_DIR}/Lineus_longissimus.faa"

# Download reference proteomes
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
        datasets download genome accession "${acc}" --include protein --filename "/tmp/${species}.zip" && \
        unzip -p "/tmp/${species}.zip" "*/protein.faa" > "${ORTHO_DIR}/${species}.faa" && \
        rm -f "/tmp/${species}.zip" || \
        echo "WARNING: Failed to download ${species}. Try manual download."
    fi
done

conda run -n ortho orthofinder \
    -f "${ORTHO_DIR}/" \
    -t ${CPUS} \
    -a $(( CPUS / 2 )) \
    -S diamond

# =============================================================================
# 8. Generate consolidated GFF3
#    Includes InterProScan and DeepTMHMM results (not available in OSX version)
# =============================================================================
echo "=== Generating consolidated functional annotation GFF3 ==="
python3 << 'PYEOF'
import gzip, re, os, json

OUTDIR = os.environ.get('OUTDIR', 'functional_annotation_linux')

# --- Parse structural GFF for coordinate mapping ---
protein_to_gene, gene_to_mrna, gene_coords = {}, {}, {}
with gzip.open("annotation/GCF_910592395.1_tnLinLong1.2_genomic.gff.gz", 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        ft, attrs = parts[2], parts[8]
        if ft == 'gene':
            m = re.search(r'ID=([^;]+)', attrs)
            if m: gene_coords[m.group(1)] = (parts[0], int(parts[3]), int(parts[4]), parts[6])
        elif ft == 'mRNA':
            m_id, m_par = re.search(r'ID=([^;]+)', attrs), re.search(r'Parent=([^;]+)', attrs)
            if m_id and m_par: gene_to_mrna[m_id.group(1)] = m_par.group(1)
        elif ft == 'CDS':
            m_pid, m_par = re.search(r'protein_id=([^;]+)', attrs), re.search(r'Parent=([^;]+)', attrs)
            if m_pid and m_par and m_pid.group(1) not in protein_to_gene:
                protein_to_gene[m_pid.group(1)] = gene_to_mrna.get(m_par.group(1), m_par.group(1))

# --- Parse eggNOG ---
eggnog = {}
with open(f"{OUTDIR}/eggnog/Llong.emapper.annotations") as f:
    for line in f:
        if line.startswith('#'): continue
        c = line.strip().split('\t')
        if len(c) < 21: continue
        eggnog[c[0]] = {k: v for k, v in zip(['go','kegg_ko','kegg_pw','cog','desc'], [c[9],c[11],c[12],c[6],c[7]]) if v != '-'}

# --- Parse KofamScan ---
kofam = {}
with open(f"{OUTDIR}/kofamscan/Llong_kofam.txt") as f:
    for line in f:
        if line.startswith('*'):
            p = line.strip().split()
            kofam.setdefault(p[1], []).append(p[2])

# --- Parse SignalP ---
signalp = set()
with open(f"{OUTDIR}/signalp/output.gff3") as f:
    for line in f:
        if not line.startswith('#'):
            p = line.strip().split('\t')
            if len(p) >= 3 and p[2] == 'signal_peptide': signalp.add(p[0])

# --- Parse dbCAN (>=2 tools) ---
dbcan = {}
with open(f"{OUTDIR}/dbcan/overview.tsv") as f:
    next(f)
    for line in f:
        c = line.strip().split('\t')
        if len(c) >= 7 and int(c[5]) >= 2: dbcan[c[0]] = c[6]

# --- Parse InterProScan TSV ---
interpro = {}
ipro_tsv = f"{OUTDIR}/interproscan/Llong_interpro.tsv"
if os.path.exists(ipro_tsv):
    with open(ipro_tsv) as f:
        for line in f:
            c = line.strip().split('\t')
            pid = c[0]
            if pid not in interpro:
                interpro[pid] = {'domains': set(), 'go': set(), 'pathways': set()}
            if len(c) > 4:
                interpro[pid]['domains'].add(f"{c[3]}:{c[4]}")
            if len(c) > 13 and c[13]:
                for go in c[13].split('|'):
                    interpro[pid]['go'].add(go)
            if len(c) > 14 and c[14]:
                for pw in c[14].split('|'):
                    interpro[pid]['pathways'].add(pw)

# --- Parse DeepTMHMM ---
deeptmhmm = {}
tmhmm_file = f"{OUTDIR}/deeptmhmm/predicted_topologies.3line"
if os.path.exists(tmhmm_file):
    with open(tmhmm_file) as f:
        lines = f.readlines()
        for i in range(0, len(lines), 3):
            pid = lines[i].strip().lstrip('>')
            topology = lines[i+2].strip() if i+2 < len(lines) else ''
            if 'M' in topology:  # has TM region
                n_tm = topology.count('M') // 2  # rough count
                deeptmhmm[pid] = topology

# --- Write consolidated GFF3 ---
all_pids = set(eggnog) | set(kofam) | signalp | set(dbcan) | set(interpro) | set(deeptmhmm)
with open(f"{OUTDIR}/Llong_functional_annotation.gff3", 'w') as out:
    out.write("##gff-version 3\n")
    out.write("# Consolidated functional annotation for Lineus longissimus\n")
    out.write("# Tools: eggNOG-mapper, KofamScan, InterProScan, SignalP 6.0, DeepTMHMM, dbCAN\n")
    written = 0
    for pid in sorted(all_pids):
        gid = protein_to_gene.get(pid)
        if not gid or gid not in gene_coords: continue
        ch, s, e, st = gene_coords[gid]
        a = [f"ID={pid}", f"gene_id={gid}"]
        eg = eggnog.get(pid, {})
        # Merge GO from eggNOG + InterProScan
        all_go = set()
        if 'go' in eg: all_go.update(eg['go'].split(','))
        if pid in interpro: all_go.update(interpro[pid]['go'])
        if all_go: a.append(f"Ontology_term={','.join(sorted(all_go))}")
        if 'kegg_ko' in eg: a.append(f"kegg_ko={eg['kegg_ko']}")
        if 'kegg_pw' in eg: a.append(f"kegg_pathway={eg['kegg_pw']}")
        if 'cog' in eg: a.append(f"cog_category={eg['cog']}")
        if 'desc' in eg: a.append(f"eggnog_description={eg['desc'].replace(';','%3B').replace('=','%3D')}")
        if pid in kofam: a.append(f"kofamscan_ko={','.join(kofam[pid])}")
        if pid in interpro and interpro[pid]['domains']:
            a.append(f"interpro_domains={','.join(sorted(interpro[pid]['domains']))}")
        if pid in interpro and interpro[pid]['pathways']:
            a.append(f"interpro_pathways={','.join(sorted(interpro[pid]['pathways']))}")
        if pid in signalp: a.append("signal_peptide=yes")
        if pid in deeptmhmm: a.append("transmembrane=yes")
        if pid in dbcan: a.append(f"cazyme_family={dbcan[pid]}")
        out.write(f"{ch}\tfunctional_annotation\tprotein_match\t{s}\t{e}\t.\t{st}\t.\t{';'.join(a)}\n")
        written += 1
    print(f"Wrote {written} entries to {OUTDIR}/Llong_functional_annotation.gff3")
PYEOF

echo "=== Pipeline complete ==="
echo "Results in: ${OUTDIR}/"
