import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests

# Attempt cyvcf2 import for real‐VCF parsing
use_real_vcf = False
vcf_ind, vcf_m, vcf_f = None, None, None
sample_ind, sample_mom, sample_dad = None, None, None
try:
    from cyvcf2 import VCF
except ImportError:
    VCF = None

st.set_page_config(page_title="Phenome Query", layout="wide")

# ── 1. Trait definitions ──
traits_info = {
    "Freckles": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R switches between brown/black eumelanin and red/yellow pheomelanin. "
            "Variants at rs1805007 & rs1805008 increase freckling."
        ),
        "inheritance": "dominant"
    },
    "Red-Green Colourblindness": {
        "gene": "OPN1LW",
        "snps": ["rs104894"],
        "description": (
            "X-linked defects in OPN1LW/OPN1MW cause red-green colourblindness. "
            "Males need one variant; females require two."
        ),
        "inheritance": "recessive"
    },
    "Hair Colour": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R variants reduce eumelanin, boosting pheomelanin. "
            "Associated with red-hued hair."
        ),
        "inheritance": "dominant",
        "interpretation": (
            "Variants at rs1805007 & rs1805008 reduce MC1R activity. "
            "Heterozygotes often have auburn hair; homozygotes typically have true red hair."
        )
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 regulates OCA2 expression. G/G at rs12913832 yields blue eyes; "
            "A/A or A/G yields brown."
        ),
        "inheritance": "recessive"
    },
    "Height": {
        "gene": None,
        "snps": [],
        "description": (
            "Height is polygenic. Mid-parental height adjusted by child’s sex gives an estimate."
        ),
        "inheritance": None
    },
    "Skin Tone": {
        "gene": "SLC24A5",
        "snps": ["rs1426654"],
        "description": "rs1426654 A/A associates with lighter skin tone.",
        "inheritance": "recessive"
    },
    "Earwax Type": {
        "gene": "ABCC11",
        "snps": ["rs17822931"],
        "description": "rs17822931 G→A: G allele → wet earwax; A/A → dry earwax.",
        "inheritance": "dominant"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": "rs4988235 T allele maintains lactase; C/C → lactose intolerance.",
        "inheritance": "dominant"
    },
    "PTC Tasting": {
        "gene": "TAS2R38",
        "snps": ["rs713598", "rs1726866"],
        "description": (
            "PAV (taster) is dominant over AVI (non-taster). "
            "PTC tasting refers to the inherited human trait of being able to taste the "
            "bitter compound phenylthiocarbamide (PTC), with some individuals perceiving "
            "it as strongly bitter while others taste little or nothing."
        ),
        "inheritance": "dominant"
    },
    "Coriander Taste": {
        "gene": "OR6A2",
        "snps": ["rs72921001"],
        "description": (
            "The perception of coriander (cilantro) as soapy is a genetic trait related "
            "to olfactory receptors that detect aldehydes, chemical compounds found in both "
            "cilantro and soap. A specific variation in an olfactory receptor gene, common "
            "in certain populations, makes these aldehydes taste and smell strongly of soap "
            "to affected individuals, overshadowing the herb's normal citrusy and herbal flavors."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "rs1815739 T allele introduces a stop codon; CC/CT normal performance. "
            "The \"sprint gene\" refers to a polymorphism in the ACTN3 gene, known as R577X "
            "(rs1815739), which influences a protein called alpha-actinin-3. Individuals "
            "with the \"R allele\" can produce this protein, which is important for fast-twitch "
            "muscle fibers involved in sprint and power performance. The \"XX genotype\" "
            "leads to a lack of functional alpha-actinin-3, which is associated with lower "
            "sprint performance and potentially higher injury risk in athletes."
        ),
        "inheritance": "dominant"
    },
    "Alcohol Flush": {
        "gene": "ALDH2",
        "snps": ["rs671"],
        "description": (
            "rs671 A allele reduces enzyme activity, causing flush. "
            "The alcohol flush reaction is an alcohol intolerance symptom where the face, "
            "neck, and chest turn red after consuming alcohol, caused by the buildup of the "
            "toxic compound acetaldehyde due to a genetic deficiency in the enzyme that breaks "
            "it down. Also known as Asian flush, this condition is common in people of East "
            "Asian descent and can lead to symptoms like nausea, rapid heartbeat, and, "
            "importantly, an increased risk of certain cancers, especially if alcohol and "
            "tobacco are consumed together."
        ),
        "inheritance": "dominant"
    }
}

# ── 2. Mock genotype data ──
mock_vcf_data = {
    snp: {"ref": d["ref"], "alt": d["alt"],
          "mother": d["mother"], "father": d["father"], "gt": d["gt"]}
    for snp, d in {
        "rs1805007":  {"ref":"C","alt":"T","mother":[0,1],"father":[1,1],"gt":[1,1]},
        "rs1805008":  {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,1]},
        "rs104894":   {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "rs12913832": {"ref":"A","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,1]},
        "rs1426654":  {"ref":"G","alt":"A","mother":[1,0],"father":[0,0],"gt":[1,0]},
        "rs17822931": {"ref":"G","alt":"A","mother":[1,1],"father":[0,1],"gt":[1,1]},
        "rs4988235":  {"ref":"C","alt":"T","mother":[1,1],"father":[0,1],"gt":[1,1]},
        "rs713598":   {"ref":"C","alt":"G","mother":[0,1],"father":[1,1],"gt":[0,1]},
        "rs1726866":  {"ref":"T","alt":"C","mother":[0,0],"father":[1,0],"gt":[0,0]},
        "rs72921001": {"ref":"T","alt":"C","mother":[1,0],"father":[1,1],"gt":[1,1]},
        "rs1815739":  {"ref":"C","alt":"T","mother":[0,0],"father":[1,1],"gt":[1,1]},
        "rs671":      {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    }.items()
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

# ── 4. Helpers ──
def zygosity(gt):
    return "Homozygous" if gt[0] == gt[1] else "Heterozygous"

def alleles_from_gt(gt, ref, alt):
    return "/".join(ref if a == 0 else alt for a in gt)

def display_genotype(gt, ref, alt):
    return f"{gt[0]}/{gt[1]}", alleles_from_gt(gt, ref, alt)

def trait_present(gt, inheritance):
    if gt is None or inheritance is None:
        return None
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
    if not m_gt or not f_gt or any(a is None for a in m_gt) or any(a is None for a in f_gt):
        return []
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
    d = mock_vcf_data[rsid]
    if role == "ind":
        return d["gt"], d["ref"], d["alt"]
    return (d["mother"], d["ref"], d["alt"]) if role == "mom" else (d["father"], d["ref"], d["alt"])

# ── 5. UI ──
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate to:", ["Individual", "Child Phenome Predictor"])
st.sidebar.subheader("Data Upload")
st.sidebar.markdown("_Disclaimer: all vcf data is not stored and is deleted after the query is run_")

if page == "Individual":
    vcf_file = st.sidebar.file_uploader("Upload Individual VCF", type=["vcf", "vcf.gz"])
    if vcf_file and VCF:
        vcf_ind = VCF(vcf_file)
        sample_ind = st.sidebar.selectbox("Select sample", vcf_ind.samples)
        use_real_vcf = True
else:
    vcf_mom = st.sidebar.file_uploader("Upload Mother VCF", type=["vcf", "vcf.gz"])
    vcf_dad = st.sidebar.file_uploader("Upload Father VCF", type=["vcf", "vcf.gz"])
    if vcf_mom and VCF:
        vcf_m = VCF(vcf_mom)
        sample_mom = st.sidebar.selectbox("Select mother sample", vcf_m.samples)
        use_real_vcf = True
    if vcf_dad and VCF:
        vcf_f = VCF(vcf_dad)
        sample_dad = st.sidebar.selectbox("Select father sample", vcf_f.samples)
        use_real_vcf = True

selected = st.multiselect("Select traits:", list(traits_info.keys()))

# Height on predictor page
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
    low, high = mean_h - 1.96 * sigma, mean_h + 1.96 * sigma
    ft_low, in_low = cm_to_ftin(low)
    ft_high, in_high = cm_to_ftin(high)
    st.markdown(f"**Predicted height**: {mean_h:.1f} cm")
    st.markdown(
        f"_95% interval_: {low:.1f}–{high:.1f} cm "
        f"(~{ft_low} ft {in_low} in to {ft_high} ft {in_high} in)"
    )
    sims = np.random.normal(mean_h, sigma, 3000)
    df = pd.DataFrame({"Height (cm)": sims})
    chart = alt.Chart(df).mark_area(opacity=0.4).encode(
        alt.X("Height (cm):Q", bin=alt.Bin(maxbins=60)),
        alt.Y("count()", stack=None),
    ).properties(height=250, width=600)
    st.altair_chart(chart)
    st.markdown("---")

# Loop traits
for trait in selected:
    if page == "Individual" and trait == "Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        if trait == "Hair Colour":
            st.subheader("Trait Interpretation")
            st.write(info["interpretation"])

        if info["snps"]:
            st.subheader("Genotypes & Inheritance")
            individual_present = False
            child_present_pcts = []

            for snp in info["snps"]:
                if page == "Individual":
                    gt, ref, alt = get_genotype(snp, "ind")
                    b, a = display_genotype(gt, ref, alt)
                    zg = zygosity(gt)
                    pres, mode = format_presence(gt, info["inheritance"])
                    if pres == "Trait present":
                        individual_present = True
                    st.markdown(f"**{snp}** (REF={ref}, ALT={alt})")
                    st.write(f"- Genotype: {b} → {a}, {zg}, {pres} {mode}")
                else:
                    m_gt, m_ref, m_alt = get_genotype(snp, "mom")
                    f_gt, f_ref, f_alt = get_genotype(snp, "dad")
                    m_b, m_a = display_genotype(m_gt, m_ref, m_alt)
                    f_b, f_a = display_genotype(f_gt, f_ref, f_alt)
                    m_zg = zygosity(m_gt); f_zg = zygosity(f_gt)
                    m_pres, m_mode = format_presence(m_gt, info["inheritance"])
                    f_pres, f_mode = format_presence(f_gt, info["inheritance"])

                    st.markdown(f"**{snp}** (REF={m_ref}, ALT={m_alt})")
                    st.write(f"- Mother: {m_b} → {m_a}, {m_zg}, {m_pres} {m_mode}")
                    st.write(f"- Father: {f_b} → {f_a}, {f_zg}, {f_pres} {f_mode}")

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

                    st.subheader("Predicted Child Genotype Probabilities")
                    probs = child_genotype_probs(m_gt, f_gt)
                    for p in probs:
                        cb = f"{p['geno'][0]}/{p['geno'][1]}"
                        ca = alleles_from_gt(p["geno"], m_ref, m_alt)
                        pres, mode = format_presence(p["geno"], info["inheritance"])
                        if pres == "Trait present":
                            child_present_pcts.append(p["pct"])
                        st.write(f"- {cb} ({ca}): {p['pct']:.0f}% → {p['zygosity']}, {pres} {mode}")

                st.markdown("")

            # Summary
            st.subheader("Summary")

            if trait == "Freckles":
                if page == "Individual":
                    total_alt = sum(get_genotype(s, "ind")[0].count(1) for s in info["snps"])
                    summary = (
                        "No freckles" if total_alt == 0
                        else "Mild freckling" if total_alt <= 2
                        else "Pronounced freckling"
                    )
                    st.write(summary)
                else:
                    m1,_,_ = get_genotype(info["snps"][0], "mom")
                    f1,_,_ = get_genotype(info["snps"][0], "dad")
                    p1 = child_genotype_probs(m1, f1)
                    m2,_,_ = get_genotype(info["snps"][1], "mom")
                    f2,_,_ = get_genotype(info["snps"][1], "dad")
                    p2 = child_genotype_probs(m2, f2)
                    P1_00 = next((x["pct"] for x in p1 if x["geno"]==(0,0)), 0)/100
                    P2_00 = next((x["pct"] for x in p2 if x["geno"]==(0,0)), 0)/100
                    P_no = P1_00 * P2_00 * 100
                    P1_hom = next((x["pct"] for x in p1 if x["geno"]==(1,1)), 0)/100
                    P2_hom = next((x["pct"] for x in p2 if x["geno"]==(1,1)), 0)/100
                    P_pron = (P1_hom + P2_hom - P1_hom * P2_hom) * 100
                    P_mild = 100 - P_no - P_pron
                    st.write(f"No freckles: {P_no:.1f}%")
                    st.write(f"Mild freckling: {P_mild:.1f}%")
                    st.write(f"Pronounced freckling: {P_pron:.1f}%")

            elif trait == "Hair Colour":
                if page == "Individual":
                    total_alt = sum(get_genotype(s, "ind")[0].count(1) for s in info["snps"])
                    hair_sum = (
                        "Non-red hair" if total_alt == 0
                        else "Auburn hair" if total_alt == 1
                        else "True red hair"
                    )
                    st.write(hair_sum)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    p00 = next((x["pct"] for x in p if x["geno"]==(0,0)),0)
                    p01 = next((x["pct"] for x in p if x["geno"]==(0,1)),0)
                    p11 = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    st.write(f"Non-red hair: {p00:.1f}%")
                    st.write(f"Auburn hair: {p01:.1f}%")
                    st.write(f"True red hair: {p11:.1f}%")

            elif trait == "Red-Green Colourblindness":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0], "ind")
                    phenotype = (
                        "Red green colour blind" if any(gt)
                        else "Not red green colour blind"
                    )
                    st.write(phenotype)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    present = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Red green colour blind: {present:.1f}%")
                    st.write(f"Not red green colour blind: {100-present:.1f}%")

            elif trait == "Eye Colour":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0], "ind")
                    phenotype = "Blue eyes" if gt == [0,0] else "Brown eyes"
                    st.write(phenotype)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    blue = next((x["pct"] for x in p if x["geno"]==(0,0)),0)
                    brown = 100 - blue
                    st.write(f"Blue eyes: {blue:.1f}%")
                    st.write(f"Brown eyes: {brown:.1f}%")

            elif trait == "Skin Tone":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0], "ind")
                    st.write(
                        "Lighter skin tone" if gt==[1,1]
                        else "Intermediate skin tone" if gt[0]!=gt[1]
                        else "Darker skin tone"
                    )
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    light = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    inter = next((x["pct"] for x in p if sorted(x["geno"])==[0,1]),0)
                    dark = 100 - light - inter
                    st.write(f"Lighter skin tone: {light:.1f}%")
                    st.write(f"Intermediate skin tone: {inter:.1f}%")
                    st.write(f"Darker skin tone: {dark:.1f}%")

            elif trait == "Earwax Type":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0], "ind")
                    ear = "Dry earwax" if gt.count(1)==2 else "Wet earwax"
                    st.write(ear)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    dry = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    wet = 100 - dry
                    st.write(f"Wet earwax: {wet:.1f}%")
                    st.write(f"Dry earwax: {dry:.1f}%")

            elif trait == "Lactose Intolerance":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0], "ind")
                    st.write("Lactose tolerant" if any(gt) else "Lactose intolerant")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0], "mom")
                    f_gt,_,_ = get_genotype(info["snps"][0], "dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    tol = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Lactose tolerant: {tol:.1f}%")
                    st.write(f"Lactose intolerant: {100-tol:.1f}%")

            elif trait == "PTC Tasting":
                if page == "Individual":
                    total_alt = sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                    st.write("Can taste PTC" if total_alt>0 else "Cannot taste PTC")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    taste = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Can taste PTC: {taste:.1f}%")
                    st.write(f"Cannot taste PTC: {100-taste:.1f}%")

            elif trait == "Coriander Taste":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Perceives coriander as soapy" if any(gt) else "Normal coriander taste")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    soapy = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Soapy perception: {soapy:.1f}%")
                    st.write(f"Normal taste: {100-soapy:.1f}%")

            elif trait == "Sprint Gene":
                if page == "Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    geno = "".join("R" if a==0 else "X" for a in sorted(gt))
                    if geno=="RR":
                        summary = "RR genotype: presence of R allele associated with better sprint and power performance."
                    elif geno in ("RX","XR"):
                        summary = "RX genotype: presence of R allele associated with better sprint and power performance."
                    else:
                        summary = (
                            "XX genotype: lack of functional alpha-actinin-3, "
                            "linked to reduced sprint performance and potentially higher injury risk."
                        )
                    st.write(summary)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt, f_gt)
                    rr = next((x["pct"] for x in p if x["geno"]==(0,0)),0)
                    rx = next((x["pct"] for x in p if sorted(x["geno"])==[0,1]),0)
                    xx = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    st.write(f"RR genotype: {rr:.1f}% → better sprint/power performance")
