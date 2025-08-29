import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests

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
            "Mutations in OPN1LW/OPN1MW on the X chromosome cause red-green colour vision defects. "
            "Males (XY) need only one mutated copy; females (XX) require two copies."
        ),
        "inheritance": "recessive"
    },
    "Hair Colour": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R variants influence hair pigmentation. ALT alleles at "
            "rs1805007/rs1805008 associate with red hair. Heterozygotes may have "
            "auburn shades; homozygotes often present true red hair."
        ),
        "inheritance": "dominant"
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 regulates OCA2 expression, affecting iris melanin. "
            "G/G at rs12913832 yields blue eyes (recessive); A/A or A/G yields "
            "brown eyes (dominant)."
        ),
        "inheritance": "recessive"
    },
    "Height": {
        "gene": None,
        "snps": [],
        "description": (
            "Height is highly polygenic. A rough estimate uses mid-parental height "
            "adjusted for the sex of the child."
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
            "ABCC11 rs17822931 G→A determines earwax: G allele → wet earwax (dominant); "
            "A/A → dry earwax (recessive)."
        ),
        "inheritance": "dominant"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": (
            "MCM6 enhancer variant rs4988235 T allele maintains lactase production "
            "into adulthood (dominant); C/C homozygotes are lactose intolerant."
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
            "C allele at rs72921001 associates with soapy flavour perception (dominant)."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "ACTN3 encodes α-actinin-3 in fast-twitch muscle fibres. "
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
    "rs1805007": {"ref":"C","alt":"T","mother":[0,1],"father":[1,1]},
    "rs1805008": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1]},
    "rs104894":  {"ref":"A","alt":"G","mother":[0,1],"father":[0,0]},
    "rs12913832":{"ref":"A","alt":"G","mother":[0,0],"father":[0,1]},
    "rs1426654": {"ref":"G","alt":"A","mother":[1,0],"father":[0,0]},
    "rs17822931":{"ref":"G","alt":"A","mother":[1,1],"father":[0,1]},
    "rs4988235": {"ref":"C","alt":"T","mother":[1,1],"father":[0,1]},
    "rs713598": {"ref":"C","alt":"G","mother":[0,1],"father":[1,1]},
    "rs1726866": {"ref":"T","alt":"C","mother":[0,0],"father":[1,0]},
    "rs72921001":{"ref":"T","alt":"C","mother":[1,0],"father":[1,1]},
    "rs1815739":{"ref":"C","alt":"T","mother":[0,0],"father":[1,1]},
    "rs671":     {"ref":"G","alt":"A","mother":[0,1],"father":[0,0]},
}

# 3. Annotation fetchers
@st.cache_data(ttl=24*3600)
def fetch_clinvar_annotation(rsid):
    rid = rsid.lstrip("rs")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rid}"
    r = requests.get(url, headers={"Accept":"application/json"})
    if not r.ok:
        return {"clinical_significance":"Unavailable"}
    data = r.json()
    cs = set()
    for anno in data.get("primary_snapshot_data", {}).get("allele_annotations", []):
        for term in anno.get("clinical_significances", []):
            cs.add(term)
    return {"clinical_significance": ", ".join(cs) if cs else "Not reported"}

@st.cache_data(ttl=24*3600)
def fetch_gnomad_freq(rsid):
    server = "https://gnomad.broadinstitute.org/api"
    query = """
    query($variantId: String!) {
      variant(variantId: $variantId) {
        genome {
          ac
          an
          populations {
            id
            ac
            an
          }
        }
      }
    }"""
    vars = {"variantId": rsid}
    r = requests.post(server, json={"query": query, "variables": vars})
    if not r.ok:
        return {"global_af": None, "populations": {}}
    v = r.json()["data"]["variant"]["genome"]
    if not v or v["an"] == 0:
        return {"global_af": None, "populations": {}}
    global_af = v["ac"] / v["an"]
    pop_af = {pop["id"]: pop["ac"]/pop["an"]
              for pop in v["populations"] if pop["an"] > 0}
    return {"global_af": global_af, "populations": pop_af}

# 4. Helper functions
def zygosity(gt):
    return "Homozygous" if gt[0]==gt[1] else "Heterozygous"

def alleles_from_gt(gt, ref, alt):
    return "/".join(ref if a==0 else alt for a in gt)

def display_genotype(gt, ref, alt):
    return f"{gt[0]}/{gt[1]}", alleles_from_gt(gt, ref, alt)

def trait_present(gt, inh):
    if inh=="dominant": return any(a==1 for a in gt)
    if inh=="recessive": return gt[0]==1 and gt[1]==1
    return None

def format_presence(gt, inh):
    pres = trait_present(gt, inh)
    if pres is None:
        return "", ""
    return ("Trait present", f"({inh})") if pres else ("Trait absent", "")

def child_genotype_probs(m_gt, f_gt):
    combos = [(m,f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    out = []
    for geno, cnt in sorted(counts.items()):
        pct = cnt/total*100
        zyg = "Homozygous" if geno[0]==geno[1] else "Heterozygous"
        out.append({"geno": geno, "pct": pct, "zygosity": zyg})
    return out

def estimate_child_height(mom_cm, dad_cm, sex):
    return (mom_cm + dad_cm + (13 if sex=="Male" else -13)) / 2

# 5. UI
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate to:", ["Individual", "Child Phenome Predictor"])
selected = st.multiselect("Select traits:", list(traits_info.keys()))

# Height calculator on predictor page
if page=="Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    col1, col2 = st.columns(2)
    with col1:
        mom_cm = st.slider("Mother’s height (cm):", 140, 200, 165)
        dad_cm = st.slider("Father’s height (cm):", 140, 200, 180)
    with col2:
        def cm_to_ftin(cm):
            total_inches = cm / 2.54
            ft = int(total_inches // 12)
            inch = int(round(total_inches % 12))
            return ft, inch
        mft, minch = cm_to_ftin(mom_cm)
        dft, dinch = cm_to_ftin(dad_cm)
        st.write(f"Mother: {mft} ft {minch} in")
        st.write(f"Father: {dft} ft {dinch} in")
    sex = st.selectbox("Child’s sex:", ["Male", "Female"])
    mean_h = estimate_child_height(mom_cm, dad_cm, sex)
    sigma = 4.7  # cm, real‐world SD around mid‐parental
    low, high = mean_h - 1.96*sigma, mean_h + 1.96*sigma
    lft, lin = cm_to_ftin(low)
    hft, hin = cm_to_ftin(high)
    st.markdown(f"**Predicted child height:** {mean_h:.1f} cm  ")
    st.markdown(f"_95% interval:_ {low:.1f}–{high:.1f} cm  ")
    st.markdown(f"_Which is ~ {lft} ft {lin} in  to  {hft} ft {hin} in_")

    sims = np.random.normal(loc=mean_h, scale=sigma, size=3000)
    df = pd.DataFrame({"Height (cm)": sims})
    chart = alt.Chart(df).mark_area(opacity=0.4).encode(
        alt.X("Height (cm):Q", bin=alt.Bin(maxbins=60)),
        alt.Y('count()', stack=None)
    ).properties(height=250, width=600)
    st.altair_chart(chart)
    st.markdown("---")

# Loop through traits
for trait in selected:
    # hide height on Individual page
    if page=="Individual" and trait=="Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        # Trait Gene Summary
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        # SNP Genotypes & Inheritance
        if info["snps"]:
            st.subheader("SNP Genotypes & Inheritance")
            for snp in info["snps"]:
                data = mock_vcf_data.get(snp)
                if not data:
                    st.write(f"- {snp}: no mock data")
                    continue

                ref, alt = data["ref"], data["alt"]
                m_gt, f_gt = data["mother"], data["father"]
                m_bin, m_alleles = display_genotype(m_gt, ref, alt)
                f_bin, f_alleles = display_genotype(f_gt, ref, alt)
                m_zyg, f_zyg = zygosity(m_gt), zygosity(f_gt)
                m_pres, m_mode = format_presence(m_gt, info["inheritance"])
                f_pres, f_mode = format_presence(f_gt, info["inheritance"])

                st.markdown(f"**SNP**: {snp} (REF={ref}, ALT={alt})")
                st.write(f"- Mother: {m_bin} → {m_alleles}, {m_zyg}, {m_pres} {m_mode}")
                st.write(f"- Father: {f_bin} → {f_alleles}, {f_zyg}, {f_pres} {f_mode}")

                # Annotations
                with st.expander("Annotations", expanded=False):
                    clin = fetch_clinvar_annotation(snp)
                    st.write(f"- ClinVar significance: {clin['clinical_significance']}")
                    gnomad = fetch_gnomad_freq(snp)
                    if gnomad["global_af"] is not None:
                        st.write(f"- gnomAD global AF: {gnomad['global_af']:.4f}")
                        for pop, af in gnomad["populations"].items():
                            st.write(f"  - {pop}: {af:.4f}")
                    else:
                        st.write("- gnomAD frequency unavailable")

                # Child probabilities
                if page=="Child Phenome Predictor":
                    st.write("**Predicted Child Genotype Probabilities**")
                    for p in child_genotype_probs(m_gt, f_gt):
                        c_bin = f"{p['geno'][0]}/{p['geno'][1]}"
                        c_alleles = alleles_from_gt(list(p["geno"]), ref, alt)
                        c_pres, c_mode = format_presence(p["geno"], info["inheritance"])
                        st.write(
                            f"- {c_bin} ({c_alleles}): {p['pct']:.0f}% → "
                            f"{p['zygosity']}, {c_pres} {c_mode}"
                        )
                st.markdown("---")
        else:
            st.write("_No defined SNPs for this trait._")
            st.markdown("---")
