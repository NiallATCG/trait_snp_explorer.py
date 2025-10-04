import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests

# Attempt cyvcf2 import for real VCF parsing
use_real_vcf = False
vcf_ind = vcf_m = vcf_f = None
sample_ind = sample_mom = sample_dad = None
try:
    from cyvcf2 import VCF
except ImportError:
    VCF = None

st.set_page_config(page_title="Phenome Query", layout="wide")

# 1. Trait definitions
traits_info = {
    "Freckles": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R switches between brown/black eumelanin and red/yellow pheomelanin. "
            "Variants at rs1805007 & rs1805008 increase freckling."
        ),
        "inheritance": "dominant",
        "overview": (
            "Freckles are small, pigmented spots on the skin that appear more prominently "
            "with sun exposure. They are strongly influenced by variants in the MC1R gene, "
            "which regulates the balance between eumelanin (brown/black pigment) and "
            "pheomelanin (red/yellow pigment). Certain MC1R variants, such as rs1805007 and "
            "rs1805008, reduce eumelanin production and increase pheomelanin, making the skin "
            "more prone to freckling. Inheritance is typically dominant, meaning one altered "
            "copy can increase the likelihood of freckles."
        )
    },
    "Red-Green Colourblindness": {
        "gene": "OPN1LW",
        "snps": ["rs104894"],
        "description": (
            "X-linked defects in OPN1LW/OPN1MW cause red-green colourblindness. "
            "Males need one variant; females require two."
        ),
        "inheritance": "recessive",
        "overview": (
            "Red‑green colourblindness is caused by defects in the OPN1LW and OPN1MW genes, "
            "which encode light‑sensitive opsins in the retina. The rs104894 variant disrupts "
            "normal function, impairing the ability to distinguish red from green. Because these "
            "genes are located on the X chromosome, inheritance is X‑linked recessive: males need "
            "only one altered copy to be affected, while females require two."
        )
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
        ),
        "overview": (
            "Hair colour is determined by the type and amount of melanin produced in hair follicles. "
            "Variants in the MC1R gene play a major role in red hair. The rs1805007 and rs1805008 "
            "variants reduce MC1R activity, shifting pigment production toward pheomelanin. "
            "Individuals with one copy often have auburn hair, while those with two copies usually "
            "have true red hair. This trait follows a dominant pattern, though expression can vary "
            "depending on other pigmentation genes."
        )
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 regulates OCA2 expression. G/G at rs12913832 yields blue eyes; "
            "A/A or A/G yields brown."
        ),
        "inheritance": "recessive",
        "overview": (
            "Eye colour is largely controlled by the HERC2 gene, which regulates expression of the "
            "neighbouring OCA2 gene involved in melanin production in the iris. A key SNP, rs12913832, "
            "determines whether eyes are blue or brown. Individuals with two G alleles typically have "
            "blue eyes due to reduced melanin, while those with at least one A allele usually have "
            "brown eyes. The inheritance is recessive, with blue eyes requiring two copies of the G allele."
        )
    },
    "Height": {
        "gene": None,
        "snps": [],
        "description": (
            "Height is polygenic. Mid-parental height adjusted by child’s sex gives an estimate."
        ),
        "inheritance": None,
        "overview": (
            "Height is a polygenic trait, meaning it is influenced by hundreds of genes as well as "
            "environmental factors like nutrition. While no single SNP determines height, a common "
            "clinical estimate uses mid‑parental height adjusted for the child’s sex. This provides a "
            "probabilistic range rather than a precise prediction, reflecting the complex inheritance of stature."
        )
    },
    "Skin Tone": {
        "gene": "SLC24A5",
        "snps": ["rs1426654"],
        "description": "rs1426654 A/A associates with lighter skin tone.",
        "inheritance": "recessive",
        "overview": (
            "Skin pigmentation is influenced by many genes, but SLC24A5 is one of the most important. "
            "The rs1426654 variant alters melanin production, with the A allele associated with lighter "
            "skin tones. Individuals with two A alleles tend to have lighter skin, heterozygotes show "
            "intermediate pigmentation, and those with two G alleles usually have darker skin. This SNP "
            "is a classic example of a recessive effect."
        )
    },
    "Earwax Type": {
        "gene": "ABCC11",
        "snps": ["rs17822931"],
        "description": "rs17822931 G→A: G allele → wet earwax; A/A → dry earwax.",
        "inheritance": "dominant",
        "overview": (
            "Earwax consistency is determined by the ABCC11 gene. The rs17822931 SNP distinguishes between "
            "wet and dry earwax. The G allele produces wet earwax, while individuals with two A alleles have "
            "dry earwax. This trait is inherited in a dominant fashion, with the presence of at least one G "
            "allele leading to wet earwax."
        )
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": "rs4988235 T allele maintains lactase; C/C → lactose intolerance.",
        "inheritance": "dominant",
        "overview": (
            "The ability to digest lactose in adulthood depends on regulatory variants near the MCM6 gene, "
            "which control expression of the lactase enzyme (LCT). The rs4988235 T allele maintains lactase "
            "activity, allowing lactose tolerance, while individuals with two C alleles typically lose lactase "
            "activity after childhood, leading to lactose intolerance. This trait is inherited in a dominant "
            "manner, with just one T allele sufficient for tolerance."
        )
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
        "inheritance": "dominant",
        "overview": (
            "The ability to taste the bitter compound phenylthiocarbamide (PTC) is controlled by the TAS2R38 gene. "
            "Different haplotypes, such as PAV (taster) and AVI (non‑taster), determine sensitivity. The PAV haplotype "
            "is dominant, so individuals with at least one copy can usually taste PTC as bitter, while AVI homozygotes "
            "cannot. This is a classic example of a simple Mendelian trait in taste perception."
        )
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
        "inheritance": "dominant",
        "overview": (
            "The perception of coriander (cilantro) as tasting soapy is influenced by the OR6A2 olfactory receptor gene. "
            "The rs72921001 variant alters sensitivity to aldehydes, compounds found in both coriander and soap. "
            "Individuals carrying the variant allele are more likely to perceive coriander as soapy, while others "
            "experience its normal citrusy flavour. The trait shows dominant inheritance, with one copy often enough "
            "to influence perception."
        )

    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "rs1815739 T allele introduces a stop codon; CC/CT normal performance. "
            "The “sprint gene” polymorphism (R577X) affects alpha-actinin-3. R allele → fast-twitch fibers; "
            "XX genotype → deficiency, reduced sprint performance."
        ),
        "inheritance": "dominant",
        "overview": (
            "The ACTN3 gene, often called the sprint gene, provides instructions for making the alpha-actinin-3 protein," 
            "which is crucial for fast-twitch muscle fibers that enable explosive,high-power movements like sprinting and jumping."
            "A specific variation, associated with the "R" allele and SNP rs1815739, is linked to greater muscle power and handgrip strength."
            "However, a non-functional variant, associated with the X allele, is linked to increased muscle damage and is"
            "common in certain populations, impacting athletic potential."
        )
    },
    "Alcohol Flush": {
        "gene": "ALDH2",
        "snps": ["rs671"],
        "description": (
            "rs671 A allele reduces enzyme activity, causing flush (Asian flush). "
            "Leads to facial redness, nausea, and higher cancer risk if alcohol+tobacco used."
        ),
        "inheritance": "dominant",
        "overview": (
            "Alcohol flush is an adverse reaction to alcohol caused by a genetic variant in the *ALDH2 gene" 
            "specifically the ALDH22 allele (rs671), which reduces the effectiveness of the acetaldehyde dehydrogenase enzyme."
            "This leads to a buildup of acetaldehyde in the body, producing symptoms like facial flushing, nausea, and headache, and is prevalent in East Asian populations."
        )
    }
}

# 2. Mock genotype data
mock_vcf_data = {
    snp: {
        "ref": d["ref"], "alt": d["alt"],
        "mother": d["mother"], "father": d["father"], "gt": d["gt"]
    }
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

# 3. ClinVar & gnomAD fetchers
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
    pops = {pop["id"]: pop["ac"] / pop["an"] for pop in v["populations"] if pop["an"]>0}
    return global_af, pops

# 4. Helpers
def zygosity(gt):
    return "Homozygous" if gt[0]==gt[1] else "Heterozygous"

def alleles_from_gt(gt, ref, alt):
    return "/".join(ref if a==0 else alt for a in gt)

def display_genotype(gt, ref, alt):
    return f"{gt[0]}/{gt[1]}", alleles_from_gt(gt, ref, alt)

def trait_present(gt, inheritance):
    if gt is None or inheritance is None:
        return None
    if inheritance=="dominant":
        return any(a==1 for a in gt)
    if inheritance=="recessive":
        return gt[0]==1 and gt[1]==1
    return None

def format_presence(gt, inheritance):
    pres = trait_present(gt, inheritance)
    if pres is None:
        return "", ""
    return ("Trait present", f"({inheritance})") if pres else ("Trait absent","")

def child_genotype_probs(m_gt, f_gt):
    if not m_gt or not f_gt or any(a is None for a in m_gt) or any(a is None for a in f_gt):
        return []
    combos = [(m, f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    out = []
    for geno, cnt in sorted(counts.items()):
        pct = cnt/total*100
        zyg = "Homozygous" if geno[0]==geno[1] else "Heterozygous"
        out.append({"geno":geno,"pct":pct,"zygosity":zyg})
    return out

def estimate_child_height(mom_cm, dad_cm, sex):
    return (mom_cm+dad_cm+(13 if sex=="Male" else -13))/2

def cm_to_ftin(cm):
    inches = cm/2.54
    ft = int(inches//12)
    inch = int(round(inches%12))
    return ft, inch

def get_genotype(rsid, role):
    if use_real_vcf and VCF:
        vobj, samp = {
            "ind": (vcf_ind, sample_ind),
            "mom": (vcf_m, sample_mom),
            "dad": (vcf_f, sample_dad)
        }[role]
        try:
            rec = next(vobj(f"{rsid}"), None)
            gt = rec.genotype(samp)["GT"]
            return gt, rec.REF, rec.ALT[0]
        except:
            pass
    d = mock_vcf_data[rsid]
    if role=="ind":
        return d["gt"], d["ref"], d["alt"]
    return (d["mother"], d["ref"], d["alt"]) if role=="mom" else (d["father"], d["ref"], d["alt"])

# 5. App UI
st.title("Genome Scan: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate to:", ["Individual","Child Phenome Predictor"])
st.sidebar.subheader("Data Upload")
st.sidebar.markdown("_Disclaimer: all vcf data is not stored and is deleted after the query is run_")

# VCF upload
if page=="Individual":
    vcf_file = st.sidebar.file_uploader("Upload Individual VCF", type=["vcf","vcf.gz"])
    if vcf_file and VCF:
        vcf_ind = VCF(vcf_file)
        sample_ind = st.sidebar.selectbox("Select sample", vcf_ind.samples)
        use_real_vcf = True
else:
    vcf_mom = st.sidebar.file_uploader("Upload Mother VCF", type=["vcf","vcf.gz"])
    vcf_dad = st.sidebar.file_uploader("Upload Father VCF", type=["vcf","vcf.gz"])
    if vcf_mom and VCF:
        vcf_m = VCF(vcf_mom)
        sample_mom = st.sidebar.selectbox("Mother sample", vcf_m.samples)
        use_real_vcf = True
    if vcf_dad and VCF:
        vcf_f = VCF(vcf_dad)
        sample_dad = st.sidebar.selectbox("Father sample", vcf_f.samples)
        use_real_vcf = True

# Define report groups
report_groups = {
    "SNP Associated diseases (Work in Progress)": [],
    "Externally Visible Characteristics": [
        "Freckles", "Hair Colour", "Eye Colour", "Skin Tone", "Earwax Type", "Height"
    ],
    "Nutrition and Fitness": [
        "Sprint Gene", "Lactose Intolerance", "PTC Tasting",
        "Coriander Taste", "Red-Green Colourblindness", "Alcohol Flush"
    ],
    "Genetic Reponse to Drugs (Work in Progress)": []
}

# Dropdown for report selection
report = st.selectbox("Select Report", list(report_groups.keys()), index=1)

# Larger multiselect for traits within chosen report
selected = st.multiselect(
    f"Select traits from {report}:", 
    report_groups[report],
    label_visibility="visible"
)

# ── Overall Summary Table ──
# ── Overall Summary Table ──
if selected:
    st.subheader("Overall Summary Table")
    rows = []
    for trait in selected:
        info = traits_info[trait]

        # Height is not SNP-based
        if trait == "Height":
            summary = "N/A"

        # Freckles: count ALT alleles across both SNPs
        elif trait == "Freckles":
            alt_count = sum(get_genotype(s, "ind")[0].count(1)
                            for s in info["snps"])
            if alt_count == 0:
                summary = "No freckles"
            elif alt_count <= 2:
                summary = "Mild freckling"
            else:
                summary = "Pronounced freckling"

        # Hair Colour: count ALT alleles
        elif trait == "Hair Colour":
            alt_count = sum(get_genotype(s, "ind")[0].count(1)
                            for s in info["snps"])
            if alt_count == 0:
                summary = "Non-red hair"
            elif alt_count == 1:
                summary = "Auburn hair"
            else:
                summary = "True red hair"

        # Eye Colour: G/G → blue, else brown
        elif trait == "Eye Colour":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = "Blue eyes" if gt == [0,0] else "Brown eyes"

        # Skin Tone: A/A → lighter, A/G → intermediate, G/G → darker
        elif trait == "Skin Tone":
            gt = get_genotype(info["snps"][0], "ind")[0]
            if gt == [1,1]:
                summary = "Lighter skin tone"
            elif gt[0] != gt[1]:
                summary = "Intermediate skin tone"
            else:
                summary = "Darker skin tone"

        # Earwax Type: A/A → dry, else wet
        elif trait == "Earwax Type":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = "Dry earwax" if gt.count(1) == 2 else "Wet earwax"

        # Lactose Intolerance: T present → tolerant, else intolerant
        elif trait == "Lactose Intolerance":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = "Lactose tolerant" if any(gt) else "Lactose intolerant"

        # PTC Tasting: any ALT → can taste, else cannot
        elif trait == "PTC Tasting":
            alt_count = sum(get_genotype(s, "ind")[0].count(1)
                            for s in info["snps"])
            summary = "Can taste PTC" if alt_count > 0 else "Cannot taste PTC"

        # Coriander Taste: any ALT → soapy perception, else normal
        elif trait == "Coriander Taste":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = ("Perceives coriander as soapy"
                       if any(gt) else "Normal coriander taste")

        # Red-Green Colourblindness: any ALT → colour blind, else not
        elif trait == "Red-Green Colourblindness":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = ("Red green colour blind" if any(gt)
                       else "Not red green colour blind")

        # Sprint Gene: any R allele (REF=0) → present, else absent
        elif trait == "Sprint Gene":
            gt = sorted(get_genotype(info["snps"][0], "ind")[0])
            # REF allele (0) is “R”
            sprint_present = 0 in gt
            summary = ("Sprint gene present: better sprint performance"
                       if sprint_present else "Sprint gene absent: reduced sprint")

        # Alcohol Flush: any ALT → flush present, else not present
        elif trait == "Alcohol Flush":
            gt = get_genotype(info["snps"][0], "ind")[0]
            summary = ("Alcohol flush present" if any(gt)
                       else "Alcohol flush not present")

        # Fallback
        else:
            summary = ""

        rows.append({"Trait": trait, "Summary": summary})

    st.table(pd.DataFrame(rows))


# Height on predictor page
if page=="Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    c1, c2 = st.columns(2)
    with c1:
        mom_cm = st.slider("Mother’s height (cm)", 140, 200, 165)
        dad_cm = st.slider("Father’s height (cm)", 140, 200, 180)
    with c2:
        ft_mom, in_mom = cm_to_ftin(mom_cm)
        ft_dad, in_dad = cm_to_ftin(dad_cm)
        st.write(f"Mother: {ft_mom} ft {in_mom} in")
        st.write(f"Father: {ft_dad} ft {in_dad} in")
    sex = st.selectbox("Child’s sex:", ["Male","Female"])
    mean_h = estimate_child_height(mom_cm, dad_cm, sex)
    sigma = 4.7
    low, high = mean_h-1.96*sigma, mean_h+1.96*sigma
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
    if page=="Individual" and trait=="Height":
        continue

    info = traits_info[trait]
    with st.expander(trait, expanded=True):
        # Gene summary
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])

        # Hair interpretation
        if trait=="Hair Colour":
            st.subheader("Trait Interpretation")
            st.write(info["interpretation"])

        # Genotypes & Inheritance
        if info["snps"]:
            st.subheader("Genotypes & Inheritance")
            individual_present = False
            child_present_pcts = []

            for snp in info["snps"]:
                if page=="Individual":
                    gt, ref, alt = get_genotype(snp, "ind")
                    b, a = display_genotype(gt, ref, alt)
                    zg = zygosity(gt)
                    pres, mode = format_presence(gt, info["inheritance"])
                    if pres=="Trait present":
                        individual_present = True
                    st.markdown(f"**{snp}** (REF={ref}, ALT={alt})")
                    st.write(f"- Genotype: {b} → {a}, {zg}, {pres} {mode}")
                else:
                    m_gt, m_ref, m_alt = get_genotype(snp,"mom")
                    f_gt, f_ref, f_alt = get_genotype(snp,"dad")
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
                        if pres=="Trait present":
                            child_present_pcts.append(p["pct"])
                        st.write(f"- {cb} ({ca}): {p['pct']:.0f}% → {p['zygosity']}, {pres} {mode}")

                st.markdown("")

            # Summary
            st.subheader("Summary")

            if trait == "Freckles":
                if page=="Individual":
                    total_alt = sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                    res = (
                        "No freckles" if total_alt==0
                        else "Mild freckling" if total_alt<=2
                        else "Pronounced freckling"
                    )
                    st.write(res)
                else:
                    m1,_,_ = get_genotype(info["snps"][0],"mom")
                    f1,_,_ = get_genotype(info["snps"][0],"dad")
                    p1 = child_genotype_probs(m1,f1)
                    m2,_,_ = get_genotype(info["snps"][1],"mom")
                    f2,_,_ = get_genotype(info["snps"][1],"dad")
                    p2 = child_genotype_probs(m2,f2)
                    P1_00 = next((x["pct"] for x in p1 if x["geno"]==(0,0)),0)/100
                    P2_00 = next((x["pct"] for x in p2 if x["geno"]==(0,0)),0)/100
                    P_no = P1_00*P2_00*100
                    P1_hom = next((x["pct"] for x in p1 if x["geno"]==(1,1)),0)/100
                    P2_hom = next((x["pct"] for x in p2 if x["geno"]==(1,1)),0)/100
                    P_pron = (P1_hom+P2_hom-P1_hom*P2_hom)*100
                    P_mild = 100-P_no-P_pron
                    st.write(f"No freckles: {P_no:.1f}%")
                    st.write(f"Mild freckling: {P_mild:.1f}%")
                    st.write(f"Pronounced freckling: {P_pron:.1f}%")

            elif trait == "Hair Colour":
                if page=="Individual":
                    total_alt = sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                    hair_sum = (
                        "Non-red hair" if total_alt==0
                        else "Auburn hair" if total_alt==1
                        else "True red hair"
                    )
                    st.write(hair_sum)
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    p00 = next((x["pct"] for x in p if x["geno"]==(0,0)),0)
                    p01 = next((x["pct"] for x in p if x["geno"]==(0,1)),0)
                    p11 = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    st.write(f"Non-red hair: {p00:.1f}%")
                    st.write(f"Auburn hair: {p01:.1f}%")
                    st.write(f"True red hair: {p11:.1f}%")

            elif trait == "Red-Green Colourblindness":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Red green colour blind" if any(gt) else "Not red green colour blind")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    pres = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Red green colour blind: {pres:.1f}%")
                    st.write(f"Not red green colour blind: {100-pres:.1f}%")

            elif trait == "Eye Colour":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Blue eyes" if gt==[0,0] else "Brown eyes")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    blue = next((x["pct"] for x in p if x["geno"]==(0,0)),0)
                    st.write(f"Blue eyes: {blue:.1f}%")
                    st.write(f"Brown eyes: {100-blue:.1f}%")

            elif trait == "Skin Tone":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    tone = (
                        "Lighter" if gt==[1,1]
                        else "Intermediate" if gt[0]!=gt[1]
                        else "Darker"
                    )
                    st.write(f"{tone} skin tone")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    light = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    inter = next((x["pct"] for x in p if sorted(x["geno"])==[0,1]),0)
                    dark = 100 - light - inter
                    st.write(f"Lighter skin tone: {light:.1f}%")
                    st.write(f"Intermediate skin tone: {inter:.1f}%")
                    st.write(f"Darker skin tone: {dark:.1f}%")

            elif trait == "Earwax Type":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Dry earwax" if gt.count(1)==2 else "Wet earwax")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    dry = next((x["pct"] for x in p if x["geno"]==(1,1)),0)
                    st.write(f"Dry earwax: {dry:.1f}%")
                    st.write(f"Wet earwax: {100-dry:.1f}%")

            elif trait == "Lactose Intolerance":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Lactose tolerant" if any(gt) else "Lactose intolerant")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    tol = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Lactose tolerant: {tol:.1f}%")
                    st.write(f"Lactose intolerant: {100-tol:.1f}%")

            elif trait == "PTC Tasting":
                if page=="Individual":
                    total_alt = sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                    st.write("Can taste PTC" if total_alt>0 else "Cannot taste PTC")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    taste = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Can taste PTC: {taste:.1f}%")
                    st.write(f"Cannot taste PTC: {100-taste:.1f}%")

            elif trait == "Coriander Taste":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Perceives coriander as soapy" if any(gt) else "Normal taste")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    soapy = sum(x["pct"] for x in p if any(x["geno"]))
                    st.write(f"Soapy perception: {soapy:.1f}%")
                    st.write(f"Normal taste: {100-soapy:.1f}%")

            elif trait == "Sprint Gene":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    geno = "".join("R" if a==0 else "X" for a in sorted(gt))
                    if "R" in geno:
                        st.write("Sprint gene present: better sprint/power performance")
                    else:
                        st.write("Sprint gene absent: reduced sprint performance")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    pres_pct = sum(x["pct"] for x in p if any(a==0 for a in x["geno"]))
                    st.write(f"Sprint gene present: {pres_pct:.1f}% chance")
                    st.write(f"Sprint gene absent: {100-pres_pct:.1f}% chance")

            elif trait == "Alcohol Flush":
                if page=="Individual":
                    gt,_,_ = get_genotype(info["snps"][0],"ind")
                    st.write("Alcohol flush present" if any(gt) else "Alcohol flush not present")
                else:
                    m_gt,_,_ = get_genotype(info["snps"][0],"mom")
                    f_gt,_,_ = get_genotype(info["snps"][0],"dad")
                    p = child_genotype_probs(m_gt,f_gt)
                    pres_pct = sum(x["pct"] for x in p if any(a==1 for a in x["geno"]))
                    st.write(f"Alcohol flush present: {pres_pct:.1f}% chance")
                    st.write(f"Alcohol flush not present: {100-pres_pct:.1f}% chance")

        else:
            st.write("_No defined SNPs for this trait._")

        st.markdown("---")

# Cleanup
if use_real_vcf:
    del vcf_ind, vcf_m, vcf_f
