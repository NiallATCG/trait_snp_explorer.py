import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests

# Try to import cyvcf2 for real VCF parsing
use_real_vcf = False
vcf_ind, vcf_m, vcf_f = None, None, None
sample_ind, sample_mom, sample_dad = None, None, None
try:
    from cyvcf2 import VCF
except ImportError:
    VCF = None  # will fall back to mock

st.set_page_config(page_title="Phenome Query", layout="wide")

# ── 1. Trait definitions (Dimples removed) ──
traits_info = {
    "Freckles": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R encodes the melanocortin-1 receptor, which switches between "
            "eumelanin (brown/black) and pheomelanin (red/yellow). Variants "
            "rs1805007 and rs1805008 increase freckling."
        ),
        "inheritance": "dominant"
    },
    "Red-Green Colourblindness": {
        "gene": "OPN1LW",
        "snps": ["rs104894"],
        "description": (
            "Mutations in OPN1LW/OPN1MW on the X chromosome cause red-green colour defects. "
            "Males need one variant; females require two copies."
        ),
        "inheritance": "recessive"
    },
    "Hair Colour": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R variants reduce brown-black pigment (eumelanin), boosting "
            "red-yellow pigment (pheomelanin)."
        ),
        "inheritance": "dominant",
        "interpretation": (
            "Variants at rs1805007 and rs1805008 reduce the receptor's ability "
            "to produce eumelanin, leading to more pheomelanin. "
            "Heterozygotes often have auburn hair; homozygotes typically have true red hair."
        )
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 influences OCA2 expression. G/G at rs12913832 yields blue eyes; "
            "A/A or A/G yields brown."
        ),
        "inheritance": "recessive"
    },
    "Height": {
        "gene": None,
        "snps": [],
        "description": (
            "Height is polygenic. Mid-parental height adjusted by child sex gives an estimate."
        ),
        "inheritance": None
    },
    "Skin Tone": {
        "gene": "SLC24A5",
        "snps": ["rs1426654"],
        "description": (
            "rs1426654 A/A is associated with lighter skin tone."
        ),
        "inheritance": "recessive"
    },
    "Earwax Type": {
        "gene": "ABCC11",
        "snps": ["rs17822931"],
        "description": (
            "rs17822931 G allele → wet earwax; A/A → dry earwax."
        ),
        "inheritance": "dominant"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": (
            "rs4988235 T allele maintains lactase; C/C → lactose intolerance."
        ),
        "inheritance": "dominant"
    },
    "PTC Tasting": {
        "gene": "TAS2R38",
        "snps": ["rs713598", "rs1726866"],
        "description": "PAV haplotype (taster) is dominant over AVI (non-taster).",
        "inheritance": "dominant"
    },
    "Coriander Taste": {
        "gene": "OR6A2",
        "snps": ["rs72921001"],
        "description": (
            "rs72921001 C allele associates with soapy flavour perception."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "rs1815739 T allele introduces stop codon; CC/CT genotypes normal."
        ),
        "inheritance": "dominant"
    },
    "Alcohol Flush": {
        "gene": "ALDH2",
        "snps": ["rs671"],
        "description": (
            "rs671 A allele reduces enzyme activity, causing flush."
        ),
        "inheritance": "dominant"
    }
}

# ── 2. Mock genotype data ──
mock_vcf_data = {
    "rs1805007": {"ref":"C","alt":"T","mother":[0,1],"father":[1,1],"gt":[1,1]},
    "rs1805008": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,1]},
    "rs104894":  {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "rs12913832":{"ref":"A","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,1]},
    "rs1426654": {"ref":"G","alt":"A","mother":[1,0],"father":[0,0],"gt":[1,0]},
    "rs17822931":{"ref":"G","alt":"A","mother":[1,1],"father":[0,1],"gt":[1,1]},
    "rs4988235": {"ref":"C","alt":"T","mother":[1,1],"father":[0,1],"gt":[1,1]},
    "rs713598": {"ref":"C","alt":"G","mother":[0,1],"father":[1,1],"gt":[0,1]},
    "rs1726866": {"ref":"T","alt":"C","mother":[0,0],"father":[1,0],"gt":[0,0]},
    "rs72921001":{"ref":"T","alt":"C","mother":[1,0],"father":[1,1],"gt":[1,1]},
    "rs1815739":{"ref":"C","alt":"T","mother":[0,0],"father":[1,1],"gt":[1,1]},
    "rs671":     {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
}

# ── 3. ClinVar & gnomAD fetchers ──
@st.cache_data(ttl=24*3600)
def fetch_clinvar(rsid):
    rid = rsid.lstrip("rs")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rid}"
    r = requests.get(url)
    if not r.ok:
        return "Unavailable"
    data = r.json()
    cs = {
        term
        for anno in data.get("primary_snapshot_data", {}).get("allele_annotations", [])
        for term in anno.get("clinical_significances", [])
    }
    return ", ".join(cs) or "Not reported"

@st.cache_data(ttl=24*3600)
def fetch_gnomad(rsid):
    server = "https://gnomad.broadinstitute.org/api"
    query = """
    query($id: String!) {
      variant(variantId: $id) {
        genome { ac an populations { id ac an } }
      }
    }"""
    r = requests.post(server, json={"query": query, "variables": {"id": rsid}})
    if not r.ok:
        return None, {}
    v = r.json()["data"]["variant"]["genome"]
    if not v or v["an"] == 0:
        return None, {}
    global_af = v["ac"] / v["an"]
    pops = {
        pop["id"]: pop["ac"] / pop["an"]
        for pop in v["populations"] if pop["an"] > 0
    }
    return global_af, pops

# ── 4. Helper functions ──
def zygosity(gt):
    return "Homozygous" if gt[0] == gt[1] else "Heterozygous"

def alleles_from_gt(gt, ref, alt):
    return "/".join(ref if a == 0 else alt for a in gt)

def display_genotype(gt, ref, alt):
    return f"{gt[0]}/{gt[1]}", alleles_from_gt(gt, ref, alt)

def trait_present(gt, inheritance):
    if inheritance == "dominant":
        return any(a == 1 for a in gt)
    if inheritance == "recessive":
        return gt[0] == 1 and gt[1] == 1
    return None

def format_presence(gt, inheritance):
    pres = trait_present(gt, inheritance)
    if pres is None:
        return "", ""
    return ("Trait present", f"({inheritance})") if pres else ("Trait absent", "")

def child_genotype_probs(m_gt, f_gt):
    combos = [(m, f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    out = []
    for geno, cnt in sorted(counts.items()):
        pct = cnt / total * 100
        zyg = "Homozygous" if geno[0] == geno[1] else "Heterozygous"
        out.append({"geno": geno, "pct": pct, "zygosity": zyg})
    return out

def estimate_child_height(mom_cm, dad_cm, sex):
    return (mom_cm + dad_cm + (13 if sex == "Male" else -13)) / 2

def cm_to_ftin(cm):
    inches = cm / 2.54
    ft = int(inches // 12)
    inch = int(round(inches % 12))
    return ft, inch

def get_genotype(rsid, role):
    """
    role: "ind", "mom", or "dad"
    Returns: (gt_list, ref, alt)
    """
    # Real VCF branch
    if use_real_vcf and VCF:
        vobj, samp = {
            "ind": (vcf_ind, sample_ind),
            "mom": (vcf_m, sample_mom),
            "dad": (vcf_f, sample_dad),
        }[role]
        try:
            rec = next(vobj(f"{rsid}"), None)
            gt = rec.genotype(samp)["GT"]
            return gt, rec.REF, rec.ALT[0]
        except Exception:
            pass

    # Mock fallback branch
    data = mock_vcf_data[rsid]
    if role == "ind":
        return data["gt"], data["ref"], data["alt"]
    if role == "mom":
        return data["mother"], data["ref"], data["alt"]
    # role == "dad"
    return data["father"], data["ref"], data["alt"]

# ── 5. App UI ──
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate to:", ["Individual", "Child Phenome Predictor"])

# Data upload
st.sidebar.subheader("Data Upload")
if page == "Individual":
    vcf_file = st.sidebar.file_uploader("Upload Individual VCF", type=["vcf","vcf.gz"])
    if vcf_file and VCF:
        vcf_ind = VCF(vcf_file)
        sample_ind = st.sidebar.selectbox("Choose sample", vcf_ind.samples)
        use_real_vcf = True
else:
    vcf_mom = st.sidebar.file_uploader("Upload Mother VCF", type=["vcf","vcf.gz"])
    vcf_dad = st.sidebar.file_uploader("Upload Father VCF", type=["vcf","vcf.gz"])
    if vcf_mom and VCF:
        vcf_m = VCF(vcf_mom)
        sample_mom = st.sidebar.selectbox("Choose mother sample", vcf_m.samples)
        use_real_vcf = True
    if vcf_dad and VCF:
        vcf_f = VCF(vcf_dad)
        sample_dad = st.sidebar.selectbox("Choose father sample", vcf_f.samples)
        use_real_vcf = True

# Trait selection
selected = st.multiselect("Select traits:", list(traits_info.keys()))

# Height calculator on predictor page
if page == "Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    c1, c2 = st.columns(2)
    with c1:
        mom_cm = st.slider("Mother’s height (cm):", 140, 200, 165)
        dad_cm = st.slider("Father’s height (cm):", 140, 200, 180)
    with c2:
        ft_mom, in_mom = cm_to_ftin(mom_cm)
        ft_dad, in_dad = cm_to_ftin(dad_cm)
        st.write(f"Mother: {ft_mom} ft {in_mom} in")
        st.write(f"Father: {ft_dad} ft {in_dad} in")
    sex = st.selectbox("Child’s sex:", ["Male", "Female"])
    mean_h = estimate_child_height(mom_cm, dad_cm, sex)
    sigma = 4.7
    low, high = mean_h - 1.96*sigma, mean_h + 1.96*sigma
    ft_low, in_low = cm_to_ftin(low)
    ft_high, in_high = cm_to_ftin(high)
    st.markdown(f"**Predicted height**: {mean_h:.1f} cm")
    st.markdown(f"_95% interval_: {low:.1f}–{high:.1f} cm "
                f"(~{ft_low} ft {in_low} in to {ft_high} ft {in_high} in)")
    sims = np.random.normal(mean_h, sigma, 3000)
    df = pd.DataFrame({"Height (cm)": sims})
    chart = alt.Chart(df).mark_area(opacity=0.4).encode(
        alt.X("Height (cm):Q", bin=alt.Bin(maxbins=60)),
        alt.Y("count()", stack=None)
    ).properties(height=250, width=600)
    st.altair_chart(chart)
    st.markdown("---")

# Loop traits
for trait in selected:
    if page == "Individual" and trait == "Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        # Trait Gene Summary
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        # Hair Colour interpretation
        if trait == "Hair Colour":
            st.subheader("Trait Interpretation")
            st.write(info["interpretation"])

        # Genotypes & Inheritance
        if info["snps"]:
            st.subheader("Genotypes & Inheritance")
            for snp in info["snps"]:
                if page == "Individual":
                    gt, ref, alt = get_genotype(snp, "ind")
                    b, a = display_genotype(gt, ref, alt)
                    zg = zygosity(gt)
                    pres, mode = format_presence(gt, info["inheritance"])
                    st.markdown(f"**{snp}** (REF={ref}, ALT={alt})")
                    st.write(f"- Genotype: {b} → {a}, {zg}, {pres} {mode}")
                    if trait == "Hair Colour":
                        cnt_alt = gt.count(1)
                        result = ("True red hair" if cnt_alt == 2
                                  else "Auburn hair" if cnt_alt == 1
                                  else "Non-red hair")
                        st.write(f"**Interpretation:** {result}")

                else:  # Child Phenome Predictor
                    # Mother
                    m_gt, m_ref, m_alt = get_genotype(snp, "mom")
                    m_b, m_a = display_genotype(m_gt, m_ref, m_alt)
                    m_zg = zygosity(m_gt)
                    m_pres, m_mode = format_presence(m_gt, info["inheritance"])
                    # Father
                    f_gt, f_ref, f_alt = get_genotype(snp, "dad")
                    f_b, f_a = display_genotype(f_gt, f_ref, f_alt)
                    f_zg = zygosity(f_gt)
                    f_pres, f_mode = format_presence(f_gt, info["inheritance"])

                    st.markdown(f"**{snp}** (REF={m_ref}, ALT={m_alt})")
                    st.write(f"- Mother: {m_b} → {m_a}, {m_zg}, {m_pres} {m_mode}")
                    st.write(f"- Father: {f_b} → {f_a}, {f_zg}, {f_pres} {f_mode}")

                    # Annotations
                    with st.expander("Annotations", expanded=False):
                        clin = fetch_clinvar(snp)
                        st.write(f"- ClinVar: {clin}")
                        gaf, pops = fetch_gnomad(snp)
                        if gaf is not None:
                            st.write(f"- gnomAD AF: {gaf:.4f}")
                            for pop, af in pops.items():
                                st.write(f"  - {pop}: {af:.4f}")
                        else:
                            st.write("- gnomAD unavailable")

                    # Predicted child genotypes
                    st.subheader("Predicted Child Genotype Probabilities")
                    probs = child_genotype_probs(m_gt, f_gt)
                    for p in probs:
                        cb = f"{p['geno'][0]}/{p['geno'][1]}"
                        ca = alleles_from_gt(p["geno"], m_ref, m_alt)
                        pres, mode =	format_presence(p["geno"], info["inheritance"])
                        st.write(f"- {cb} ({ca}): {p['pct']:.0f}% → "
                                 f"{p['zygosity']}, {pres} {mode}")

                    # Hair colour prediction
                    if trait == "Hair Colour":
                        best = max(probs, key=lambda x: x["pct"])
                        cnt_alt = best["geno"].count(1)
                        result = ("True red hair" if cnt_alt == 2
                                  else "Auburn hair" if cnt_alt == 1
                                  else "Non-red hair")
                        st.write(f"**Child hair interpretation:** {result}")

                st.markdown("---")
        else:
            st.write("_No defined SNPs for this trait._")
            st.markdown("---")

# ── Delete VCF references after querying to free memory ──
# VCF files are deleted after use
if use_real_vcf:
    del vcf_ind, vcf_m, vcf_f
