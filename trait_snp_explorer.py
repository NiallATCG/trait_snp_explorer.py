import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests
from cyvcf2 import VCF

st.set_page_config(page_title="Phenome Query", layout="wide")

# 1. Trait definitions (Dimples removed)
traits_info = {
    "Freckles": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R encodes the melanocortin-1 receptor, which switches between "
            "eumelanin (brown/black) and pheomelanin (red/yellow). Variants "
            "rs1805007 (Arg151Cys) and rs1805008 (Arg160Trp) increase freckling. "
            "Heterozygotes often have mild freckling; homozygotes show pronounced freckling."
        ),
        "inheritance": "dominant"
    },
    "Red-Green Colourblindness": {
        "gene": "OPN1LW",
        "snps": ["rs104894"],
        "description": (
            "Mutations in OPN1LW/OPN1MW on the X chromosome cause red-green colour defects. "
            "Males (XY) need a single variant; females (XX) require two copies."
        ),
        "inheritance": "recessive"
    },
    "Hair Colour": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R variants at rs1805007/rs1805008 reduce eumelanin, boosting pheomelanin. "
            "Heterozygotes often have auburn hair; homozygotes typically have true red hair."
        ),
        "inheritance": "dominant",
        "interpretation": {
            "text": (
                "Variants in MC1R at rs1805007 and rs1805008 reduce the receptor's ability "
                "to stimulate brown-black pigment (eumelanin), leading to more red-yellow pigment (pheomelanin). "
                "Heterozygotes may have auburn hair; homozygotes are more likely to have true red hair."
            )
        }
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 regulates OCA2 expression, affecting iris melanin. "
            "G/G at rs12913832 yields blue eyes (recessive); A/A or A/G yields brown eyes (dominant)."
        ),
        "inheritance": "recessive"
    },
    "Height": {
        "gene": None,
        "snps": [],
        "description": (
            "Height is highly polygenic. A mid-parental estimate uses parental heights "
            "adjusted by child sex."
        ),
        "inheritance": None
    },
    "Skin Tone": {
        "gene": "SLC24A5",
        "snps": ["rs1426654"],
        "description": (
            "SLC24A5 variant rs1426654 A allele is associated with lighter skin tone. "
            "A/A yields lighter pigmentation (recessive)."
        ),
        "inheritance": "recessive"
    },
    "Earwax Type": {
        "gene": "ABCC11",
        "snps": ["rs17822931"],
        "description": (
            "ABCC11 rs17822931 G→A determines earwax type: G allele → wet earwax (dominant); "
            "A/A → dry earwax (recessive)."
        ),
        "inheritance": "dominant"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": (
            "MCM6 enhancer rs4988235 T allele maintains lactase into adulthood (dominant); "
            "C/C homozygotes are lactose intolerant."
        ),
        "inheritance": "dominant"
    },
    "PTC Tasting": {
        "gene": "TAS2R38",
        "snps": ["rs713598", "rs1726866"],
        "description": (
            "TAS2R38 haplotypes at rs713598/rs1726866 determine PTC tasting. "
            "PAV (taster) is dominant over AVI (non-taster)."
        ),
        "inheritance": "dominant"
    },
    "Coriander Taste": {
        "gene": "OR6A2",
        "snps": ["rs72921001"],
        "description": (
            "OR6A2 encodes a receptor responding to aldehydes in coriander. "
            "C allele at rs72921001 associates with soapy flavour (dominant)."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "ACTN3 encodes α-actinin-3 in fast-twitch fibres. "
            "T allele introduces a stop codon; CC/CT genotypes (normal) are dominant over TT."
        ),
        "inheritance": "dominant"
    },
    "Alcohol Flush": {
        "gene": "ALDH2",
        "snps": ["rs671"],
        "description": (
            "ALDH2 rs671 A allele reduces enzyme activity, causing alcohol flush. "
            "A is semi-dominant: heterozygotes flush moderately; A/A flush strongly."
        ),
        "inheritance": "dominant"
    }
}

# 2. Mock genotype data
mock_vcf_data = {
    "rs1805007": {"ref":"C","alt":"T","gt":[1,1]},  # mock male
    "rs1805008": {"ref":"G","alt":"A","gt":[0,1]},
    "rs104894":  {"ref":"A","alt":"G","gt":[0,1]},
    "rs12913832":{"ref":"A","alt":"G","gt":[0,1]},
    "rs1426654": {"ref":"G","alt":"A","gt":[1,0]},
    "rs17822931":{"ref":"G","alt":"A","gt":[1,1]},
    "rs4988235": {"ref":"C","alt":"T","gt":[1,1]},
    "rs713598": {"ref":"C","alt":"G","gt":[0,1]},
    "rs1726866": {"ref":"T","alt":"C","gt":[0,0]},
    "rs72921001":{"ref":"T","alt":"C","gt":[1,1]},
    "rs1815739":{"ref":"C","alt":"T","gt":[1,1]},
    "rs671":     {"ref":"G","alt":"A","gt":[0,1]},
}

# 3. Annotation fetchers
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

# 4. Helpers
def zygosity(gt): return "Homozygous" if gt[0]==gt[1] else "Heterozygous"
def alleles(gt, ref, alt): return "/".join(ref if a==0 else alt for a in gt)
def display(gt, ref, alt): return f"{gt[0]}/{gt[1]}", alleles(gt, ref, alt)
def present(gt, inh):
    if inh=="dominant": return any(a==1 for a in gt)
    if inh=="recessive": return gt[0]==1 and gt[1]==1
    return None
def fmt_presence(gt, inh):
    p = present(gt, inh)
    if p is None: return "", ""
    return ("Trait present","("+inh+")") if p else ("Trait absent","")
def child_probs(m, f):
    combos = [(mi, fi) for mi in m for fi in f]
    cnt = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    out = []
    for geno, c in sorted(cnt.items()):
        pct = c/total*100
        zog = "Homozygous" if geno[0]==geno[1] else "Heterozygous"
        out.append({"geno":geno,"pct":pct,"zygosity":zog})
    return out
def midparent(m, f, sex):
    return (m + f + (13 if sex=="Male" else -13)) / 2
def cm_to_ftin(cm):
    inches = cm/2.54
    ft = int(inches//12)
    inch = int(round(inches%12))
    return ft, inch

# 5. UI
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate:", ["Individual","Child Phenome Predictor"])

# Data upload
st.sidebar.subheader("Data Upload")
if page=="Individual":
    vcf_file = st.sidebar.file_uploader("Upload VCF (individual)", type=["vcf","vcf.gz"])
    sample_id = None
    if vcf_file:
        vcf = VCF(vcf_file)
        sample_id = st.sidebar.selectbox("Select sample", vcf.samples)
else:
    vcf_m = st.sidebar.file_uploader("Upload VCF (mother)", type=["vcf","vcf.gz"])
    vcf_f = st.sidebar.file_uploader("Upload VCF (father)", type=["vcf","vcf.gz"])
    mom_id = dad_id = None
    if vcf_m:
        vcfM = VCF(vcf_m)
        mom_id = st.sidebar.selectbox("Select mother sample", vcfM.samples)
    if vcf_f:
        vcfF = VCF(vcf_f)
        dad_id = st.sidebar.selectbox("Select father sample", vcfF.samples)

# Trait selection
selected = st.multiselect("Select traits:", list(traits_info.keys()))

# Height on predictor page
if page=="Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    c1, c2 = st.columns(2)
    with c1:
        mom_cm = st.slider("Mother’s height (cm):", 140, 200, 165)
        dad_cm = st.slider("Father’s height (cm):", 140, 200, 180)
    with c2:
        st.write(f"Mother: {cm_to_ftin(mom_cm)[0]} ft {cm_to_ftin(mom_cm)[1]} in")
        st.write(f"Father: {cm_to_ftin(dad_cm)[0]} ft {cm_to_ftin(dad_cm)[1]} in")
    sex = st.selectbox("Child’s sex:", ["Male","Female"])
    mean_h = midparent(mom_cm, dad_cm, sex)
    sigma = 4.7
    low, high = mean_h - 1.96*sigma, mean_h + 1.96*sigma
    lft, lin = cm_to_ftin(low); hft, hin = cm_to_ftin(high)
    st.markdown(f"**Predicted height**: {mean_h:.1f} cm ({lft} ft {lin} in)")
    st.markdown(f"_95% interval_: {low:.1f}–{high:.1f} cm (~{lft} ft {lin} in to {hft} ft {hin} in)")
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
    if page=="Individual" and trait=="Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        # Summary
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        # Interpret hair colour
        if trait=="Hair Colour":
            st.subheader("Trait Interpretation")
            st.write(info["interpretation"]["text"])

        # Genotype display
        if info["snps"]:
            st.subheader("Genotypes & Inheritance")
            for snp in info["snps"]:
                # load real or mock
                if page=="Individual":
                    if sample_id:
                        v = VCF(vcf_file)(f"{snp}")
                        rec = next(v, None)
                        gt = rec.genotype(sample_id)["GT"]
                        ref, alt = rec.REF, rec.ALT[0]
                    else:
                        data = mock_vcf_data[snp]
                        gt, ref, alt = data["gt"], data["ref"], data["alt"]
                    b, a = display(gt, ref, alt)
                    pres, mode = fmt_presence(gt, info["inheritance"])
                    st.markdown(f"**{snp}** (REF={ref}, ALT={alt})")
                    st.write(f"- Genotype: {b} → {a}, {zygosity(gt)}, {pres} {mode}")
                    # hair interpretation
                    if trait=="Hair Colour":
                        if any(gt.count(1) > 1 for _ in []): pass
                        cnt_alt = gt.count(1)
                        result = "True red hair" if cnt_alt==2 else ("Auburn hair" if cnt_alt==1 else "Non-red hair")
                        st.write(f"**Interpretation:** {result}")
                else:
                    # predictor page
                    # mother
                    if mom_id:
                        recM = next(VCF(vcf_m)(f"{snp}"), None)
                        m_gt = recM.genotype(mom_id)["GT"]; m_ref, m_alt = recM.REF, recM.ALT[0]
                    else:
                        d = mock_vcf_data[snp]; m_gt, m_ref, m_alt = d["mother"], d["ref"], d["alt"]
                    # father
                    if dad_id:
                        recF = next(VCF(vcf_f)(f"{snp}"), None)
                        f_gt = recF.genotype(dad_id)["GT"]; f_ref, f_alt = recF.REF, recF.ALT[0]
                    else:
                        d = mock_vcf_data[snp]; f_gt = d["father"]; f_ref, f_alt = d["ref"], d["alt"]
                    # display
                    st.markdown(f"**{snp}** (REF={m_ref}, ALT={m_alt})")
                    for label, gt, ref, alt in [
                        ("Mother", m_gt, m_ref, m_alt),
                        ("Father", f_gt, f_ref, f_alt),
                    ]:
                        b, a = display(gt, ref, alt)
                        pres, mode = fmt_presence(gt, info["inheritance"])
                        st.write(f"- {label}: {b} → {a}, {zygosity(gt)}, {pres} {mode}")
                    # annotations
                    with st.expander("Annotations", expanded=False):
                        cv = fetch_clinvar(snp)
                        st.write(f"- ClinVar: {cv}")
                        gaf, pops = fetch_gnomad(snp)
                        if gaf is not None:
                            st.write(f"- gnomAD AF: {gaf:.4f}")
                            for pop, af in pops.items():
                                st.write(f"  - {pop}: {af:.4f}")
                        else:
                            st.write("- gnomAD unavailable")
                    # child probabilities
                    st.subheader("Predicted Child Genotype Probabilities")
                    for p in child_probs(m_gt, f_gt):
                        cb = f"{p['geno'][0]}/{p['geno'][1]}"
                        ca = alleles(p["geno"], m_ref, m_alt)
                        pres, mode = fmt_presence(p["geno"], info["inheritance"])
                        st.write(f"- {cb} ({ca}): {p['pct']:.0f}% → {p['zygosity']}, {pres} {mode}")
                    # hair interpretation for child
                    if trait=="Hair Colour":
                        best = max(child_probs(m_gt, f_gt), key=lambda x: x["pct"])
                        cnt_alt = list(best["geno"]).count(1)
                        result = "True red hair" if cnt_alt==2 else ("Auburn hair" if cnt_alt==1 else "Non-red hair")
                        st.write(f"**Child hair interpretation:** {result}")
                st.markdown("---")
        else:
            st.write("_No defined SNPs for this trait._")
            st.markdown("---")
