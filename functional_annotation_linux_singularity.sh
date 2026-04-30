#!/usr/bin/env bash
# =============================================================================
# Functional Annotation Pipeline for Lineus longissimus (Linux version with Singularity)
# =============================================================================
# Author: Anestis Gkanogiannis <anestis@gkanogiannis.com> - MNCN-CSIC, 2026
# =============================================================================
# Uses Singularity/Apptainer containers instead of conda environments.
#
# Container build (run once):
#   # --- BioContainers (pre-built, pulled from registries) ---
#   singularity pull ${SIF_DIR}/eggnog-mapper.sif docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0
#   singularity pull ${SIF_DIR}/kofamscan.sif      docker://quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2
#   singularity pull ${SIF_DIR}/interproscan.sif    docker://quay.io/biocontainers/interproscan:5.55_88.0--hec16e2b_0
#   singularity pull ${SIF_DIR}/dbcan.sif           docker://quay.io/biocontainers/dbcan:5.2.7--pyhdfd78af_0
#   singularity pull ${SIF_DIR}/orthofinder.sif     docker://quay.io/biocontainers/orthofinder:2.5.5--hdfd78af_2
#
#   # --- Licensed / custom-built containers ---
#   # SignalP 6.0 and DeepTMHMM require academic licenses from DTU.
#   # Build from your local installs:
#   #   sudo -E singularity build signalp6.sif signalp6.def
#   #   sudo -E singularity build deeptmhmm.sif deeptmhmm.def
#
# Database setup: same as conda version (eggNOG, KofamScan, dbCAN, InterProScan).
# =============================================================================

#SBATCH -J annotate-singularity
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH --gres=gpu:a100:1
#SBATCH -o logs/%x-%j.out
#SBATCH -e logs/%x-%j.err

set -euo pipefail

# --- Load Singularity module (CESGA) ---
module load cesga/2025 singularity/4.3.3

# --- Configuration ---
SUF="2026-04-17"
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
PROTFILE="${WORKDIR}/annotation/GCF_910592395.1_tnLinLong1.2_protein.faa.gz"
OUTDIR="${WORKDIR}/functional_annotation_linux_singularity_${SUF}"
CPUS=$(nproc)

DB_BASE="/${STORE}/databases"
DB_EGGNOG="${DB_BASE}/eggnog"
DB_KOFAM="${DB_BASE}/kofam"
DB_DBCAN="${DB_BASE}/dbcan"
# DeepTMHMM is bundled inside its container (see deeptmhmm.def in APPENDIX)

# --- Container images ---
SIF_DIR="${STORE}/containers"
SIF_EGGNOG="${SIF_DIR}/eggnog-mapper.sif"
SIF_KOFAM="${SIF_DIR}/kofamscan.sif"
SIF_INTERPRO="${SIF_DIR}/interproscan.sif"
SIF_SIGNALP="${SIF_DIR}/signalp6.sif"
SIF_DEEPTMHMM="${SIF_DIR}/deeptmhmm.sif"
SIF_DBCAN="${SIF_DIR}/dbcan.sif"
SIF_ORTHOFINDER="${SIF_DIR}/orthofinder.sif"

# --- Common Singularity options ---
# Bind the working directory, databases, and tmp so every container sees them.
SING_BIND="${WORKDIR}:${WORKDIR},${DB_BASE}:${DB_BASE},${OUTDIR}:${OUTDIR},/tmp:/tmp"
SING_OPTS="--cleanenv --bind ${SING_BIND}"

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
singularity exec ${SING_OPTS} "${SIF_EGGNOG}" \
    emapper.py \
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
singularity exec ${SING_OPTS} "${SIF_KOFAM}" \
    exec_annotation \
        -f detail-tsv \
        -o "${OUTDIR}/kofamscan/Llong_kofam.txt" \
        --profile "${DB_KOFAM}/profiles" \
        --ko-list "${DB_KOFAM}/ko_list" \
        --cpu ${CPUS} \
        "${PROTEINS}"

# =============================================================================
# 3. InterProScan (domains, GO, pathways -- FULL analysis)
# =============================================================================
echo "=== [3/7] Running InterProScan ==="
# InterProScan needs a writable temp directory; use a local one to avoid
# conflicts on shared /tmp.
IPRO_TMP="${OUTDIR}/interproscan/tmp"
mkdir -p "${IPRO_TMP}"

singularity exec ${SING_OPTS},${IPRO_TMP}:/tmp "${SIF_INTERPRO}" \
    interproscan.sh \
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
singularity exec ${SING_OPTS} "${SIF_SIGNALP}" \
    signalp6 \
        --fastafile "${PROTEINS}" \
        --output_dir "${OUTDIR}/signalp/" \
        --organism eukarya \
        --format txt

# =============================================================================
# 5. DeepTMHMM (transmembrane topology prediction -- GPU accelerated)
#    --nv enables NVIDIA GPU passthrough inside the container.
# =============================================================================
echo "=== [5/7] Running DeepTMHMM (GPU) ==="
ABSFASTA="$(realpath "${PROTEINS}")"
ABSOUT="$(realpath "${OUTDIR}/deeptmhmm")"
rmdir "${ABSOUT}" 2>/dev/null || true  # DeepTMHMM requires non-existing output dir

singularity run --nv ${SING_OPTS} "${SIF_DEEPTMHMM}" \
        --fasta "${ABSFASTA}" \
        --output-dir "${ABSOUT}"

# =============================================================================
# 6. dbCAN (CAZyme annotation)
# =============================================================================
echo "=== [6/7] Running dbCAN ==="
singularity exec ${SING_OPTS} "${SIF_DBCAN}" \
    run_dbcan CAZyme_annotation \
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

# Download reference proteomes (uses host NCBI datasets CLI or the container's)
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
        datasets download genome accession "${acc}" \
            --include protein --filename "/tmp/${species}.zip" && \
        unzip -p "/tmp/${species}.zip" "*/protein.faa" > "${ORTHO_DIR}/${species}.faa" && \
        rm -f "/tmp/${species}.zip" || \
        echo "WARNING: Failed to download ${species}. Try manual download."
    fi
done

singularity exec ${SING_OPTS} "${SIF_ORTHOFINDER}" \
    orthofinder \
        -f "${ORTHO_DIR}/" \
        -t ${CPUS} \
        -a $(( CPUS / 2 )) \
        -S diamond

# =============================================================================
# 8. Generate consolidated GFF3
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
            if 'M' in topology:
                n_tm = topology.count('M') // 2
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
