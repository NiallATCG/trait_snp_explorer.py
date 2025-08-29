import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt

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
            "Red-green colour vision defects arise from mutations in OPN1LW/OPN1MW "
            "on the X chromosome. Males (XY) need only one mutated copy; females "
            "(XX) require two copies for the phenotype."
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
            "Height is highly polygenic. A rough estimate uses mid-parental height, "
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
    "rs1805007": {"ref": "C","alt": "T","mother":[0,1],"father":[1,1]},
    "rs1805008": {"ref": "G","alt": "A","mother":[0,0],"father":[0,1]},
    "rs104894":  {"ref": "A","alt": "G","mother":[0,1],"father":[0,0]},
    "rs12913832":{"ref": "A","alt": "G","mother":[0,0],"father":[0,1]},
    "rs1426654": {"ref": "G","alt": "A","mother":[1,0],"father":[0,0]},
    "rs17822931":{"ref": "G","alt": "A","mother":[1,1],"father":[0,1]},
    "rs4988235": {"ref": "C","alt": "T","mother":[1,1],"father":[0,1]},
    "rs713598": {"ref": "C","alt": "G","mother":[0,1],"father":[1,1]},
    "rs1726866": {"ref": "T","alt": "C","mother":[0,0],"father":[1,0]},
    "rs72921001":{"ref": "T","alt": "C","mother":[1,0],"father":[1,1]},
    "rs1815739":{"ref": "C","alt": "T","mother":[0,0],"father":[1,1]},
    "rs671":     {"ref": "G","alt": "A","mother":[0,1],"father":[0,0]},
}

# 3. Helper functions
def zygosity(gt): return "Homozygous" if gt[0]==gt[1] else "Heterozygous"
def alleles_from_gt(gt, ref, alt):
    return "/".join(ref if a==0 else alt for a in gt)
def display_genotype(gt, ref, alt):
    binary = f"{gt[0]}/{gt[1]}"
    return binary, alleles_from_gt(gt, ref, alt)
def trait_present(gt, inh):
    if inh=="dominant": return any(a==1 for a in gt)
    if inh=="recessive": return gt[0]==1 and gt[1]==1
    return None
def format_presence(gt, inh):
    pres = trait_present(gt, inh)
    if pres is None: return ("","")
    return ("Trait present","("+inh+")") if pres else ("Trait absent","")
def child_genotype_probs(m_gt, f_gt):
    combos = [(m,f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    out = []
    for geno,cnt in sorted(counts.items()):
        pct = cnt/total*100
        zyg = "Homozygous" if geno[0]==geno[1] else "Heterozygous"
        out.append({"geno":geno,"pct":pct,"zygosity":zyg})
    return out
def estimate_child_height(mom, dad, sex):
    return (mom + dad + (13 if sex=="Male" else -13)) / 2

# 4. App UI
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")

page = st.sidebar.radio("Navigate to:", ["Individual","Child Phenome Predictor"])
selected = st.multiselect("Select traits:", list(traits_info.keys()))

# Height calculator only on predictor page
if page=="Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    mom_h = st.slider("Mother’s height (cm):", 140, 200, 165)
    dad_h = st.slider("Father’s height (cm):", 140, 200, 180)
    sex = st.selectbox("Child’s sex:", ["Male","Female"])
    mean_h = estimate_child_height(mom_h, dad_h, sex)
    ci_low, ci_high = mean_h - 10, mean_h + 10  # ±10 cm ~95% interval
    st.markdown(f"**Predicted child height**: {mean_h:.1f} cm  ")
    st.markdown(f"_95% interval:_ {ci_low:.1f}–{ci_high:.1f} cm")

    # Generate distribution and plot
    sims = np.random.normal(loc=mean_h, scale=5, size=2000)
    df = pd.DataFrame({"Height": sims})
    chart = alt.Chart(df).mark_area(
        opacity=0.5,
        interpolate='step'
    ).encode(
        alt.X("Height:Q", bin=alt.Bin(maxbins=50)),
        alt.Y('count()', stack=None)
    ).properties(height=200, width=600)
    st.altair_chart(chart)
    st.markdown("---")

# Loop traits in collapsible panels
for trait in selected:
    # Skip height on Individual page
    if page=="Individual" and trait=="Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        # SNP details (if any)
        if info["snps"]:
            st.subheader("SNP Genotypes & Inheritance")
            for snp in info["snps"]:
                data = mock_vcf_data.get(snp)
                if not data:
                    st.write(f"- {snp}: no mock data")
                    continue

                ref, alt = data["ref"], data["alt"]
                m_gt, f_gt = data["mother"], data["father"]
                m_bin,m_alleles = display_genotype(m_gt, ref, alt)
                f_bin,f_alleles = display_genotype(f_gt, ref, alt)
                m_zyg,f_zyg = zygosity(m_gt), zygosity(f_gt)
                m_pres,m_mode = format_presence(m_gt, info["inheritance"])
                f_pres,f_mode = format_presence(f_gt, info["inheritance"])

                st.markdown(f"**SNP**: {snp}  (REF={ref}, ALT={alt})")
                st.write(f"- Mother: {m_bin} → {m_alleles}, {m_zyg}, {m_pres} {m_mode}")
                st.write(f"- Father: {f_bin} → {f_alleles}, {f_zyg}, {f_pres} {f_mode}")

                if page=="Child Phenome Predictor":
                    st.write("**Predicted Child Genotype Probabilities**")
                    for p in child_genotype_probs(m_gt, f_gt):
                        cb = f"{p['geno'][0]}/{p['geno'][1]}"
                        cal = alleles_from_gt(list(p["geno"]), ref, alt)
                        c_pres,c_mode = format_presence(list(p["geno"]), info["inheritance"])
                        st.write(
                            f"- {cb} ({cal}): {p['pct']:.0f}% → "
                            f"{p['zygosity']}, {c_pres} {c_mode}"
                        )
                st.markdown("")  # spacing
        else:
            st.write("_No defined SNPs for this trait._")
        st.markdown("---")
