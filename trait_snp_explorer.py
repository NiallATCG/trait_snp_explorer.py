import streamlit as st
from collections import Counter
import numpy as np
import pandas as pd
import altair as alt
import requests

# To raise upload limit to 500 MB, add a file at .streamlit/config.toml:
# [server]
# maxUploadSize = 500

# Attempt cyvcf2 import for real‐VCF parsing
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
        "snps": ["rs1805007","rs1805008"],
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
        "snps": ["rs1805007","rs1805008"],
        "description": (
            "MC1R variants reduce eumelanin and boost pheomelanin, associated with red hair."
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
            "HERC2 regulates OCA2. G/G at rs12913832 → blue eyes; A/A or A/G → brown eyes."
        ),
        "inheritance": "recessive"
    },
    "Height": {
        "gene": None, "snps": [],
        "description": "Height is polygenic; mid-parental estimate adjusted by sex.",
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
        "description": "rs17822931 G→A: G allele → wet; A/A → dry earwax.",
        "inheritance": "dominant"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": "rs4988235 T allele → lactase persistence; C/C → intolerance.",
        "inheritance": "dominant"
    },
    "PTC Tasting": {
        "gene": "TAS2R38",
        "snps": ["rs713598","rs1726866"],
        "description": (
            "PAV (taster) is dominant; AVI (non-taster). "
            "PTC tasting refers to tasting the bitter compound phenylthiocarbamide."
        ),
        "inheritance": "dominant"
    },
    "Coriander Taste": {
        "gene": "OR6A2",
        "snps": ["rs72921001"],
        "description": (
            "Variants in OR6A2 cause coriander to taste soapy by detecting aldehydes found in "
            "cilantro and soap, overpowering citrusy/herbal notes."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "R577X polymorphism: R allele → α-actinin-3 in fast-twitch fibers; "
            "XX genotype → deficiency, reduced sprint performance."
        ),
        "inheritance": "dominant"
    },
    "Alcohol Flush": {
        "gene": "ALDH2",
        "snps": ["rs671"],
        "description": (
            "rs671 A allele reduces ALDH2 activity, causing alcohol flush (face, neck red). "
            "Also known as Asian flush, linked to nausea, rapid heartbeat, cancer risk."
        ),
        "inheritance": "dominant"
    }
}

# 2. Mock genotype data
mock_vcf_data = {
    snp: dict(d)
    for snp,d in {
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

# 3. ClinVar & gnomAD
@st.cache_data(ttl=24*3600)
def fetch_clinvar(rsid):
    r = requests.get(f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid.lstrip('rs')}")
    if not r.ok:
        return "Unavailable"
    data = r.json()
    cs = {
        term
        for anno in data["primary_snapshot_data"]["allele_annotations"]
        for term in anno.get("clinical_significances", [])
    }
    return ", ".join(cs) or "Not reported"

@st.cache_data(ttl=24*3600)
def fetch_gnomad(rsid):
    query = """
    query($id:String!){variant(variantId:$id){genome{ac an populations{id ac an}}}}
    """
    r = requests.post("https://gnomad.broadinstitute.org/api", json={"query":query,"variables":{"id":rsid}})
    if not r.ok:
        return None, {}
    v = r.json()["data"]["variant"]["genome"]
    if not v or v["an"]==0:
        return None, {}
    af = v["ac"]/v["an"]
    pops = {pop["id"]:pop["ac"]/pop["an"] for pop in v["populations"] if pop["an"]>0}
    return af, pops

# 4. Helpers
def zygosity(gt): return "Homozygous" if gt[0]==gt[1] else "Heterozygous"
def alleles_from_gt(gt, ref, alt): return "/".join(ref if a==0 else alt for a in gt)
def display_genotype(gt, ref, alt): return f"{gt[0]}/{gt[1]}", alleles_from_gt(gt, ref, alt)

def trait_present(gt, inh):
    if gt is None or inh is None:
        return None
    if inh=="dominant": return any(a==1 for a in gt)
    if inh=="recessive": return gt[0]==1 and gt[1]==1
    return None

def format_presence(gt, inh):
    pres = trait_present(gt, inh)
    if pres is None: return "", ""
    return ("Trait present",f"({inh})") if pres else ("Trait absent","")

def child_genotype_probs(m_gt, f_gt):
    if not m_gt or not f_gt or any(a is None for a in m_gt) or any(a is None for a in f_gt):
        return []
    combos = [(m,f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    return [
        {"geno":geno, "pct":cnt/total*100, 
         "zygosity":("Homozygous" if geno[0]==geno[1] else "Heterozygous")}
        for geno,cnt in sorted(counts.items())
    ]

def estimate_child_height(mom,dad,sex): return (mom+dad+(13 if sex=="Male" else -13))/2
def cm_to_ftin(cm):
    inches=cm/2.54; ft=int(inches//12); inch=int(round(inches%12)); return ft, inch

def get_genotype(rsid, role):
    if use_real_vcf and VCF:
        vobj,samp=({"ind":(vcf_ind,sample_ind),
                    "mom":(vcf_m,sample_mom),
                    "dad":(vcf_f,sample_dad)})[role]
        try:
            rec=next(vobj(f"{rsid}"),None)
            gt=rec.genotype(samp)["GT"]
            return gt,rec.REF,rec.ALT[0]
        except:
            pass
    d=mock_vcf_data[rsid]
    if role=="ind": return d["gt"],d["ref"],d["alt"]
    return (d["mother"],d["ref"],d["alt"]) if role=="mom" else (d["father"],d["ref"],d["alt"])

# 5. UI
st.title("Phenome Query: Enhanced Trait-Based SNP Explorer")
page = st.sidebar.radio("Navigate to:",["Individual","Child Phenome Predictor"])
st.sidebar.subheader("Data Upload")
st.sidebar.markdown("_Disclaimer: all VCF data is not stored and is deleted after the query is run_")

if page=="Individual":
    vcf_file = st.sidebar.file_uploader("Upload Individual VCF", type=["vcf","vcf.gz"])
    if vcf_file and VCF:
        vcf_ind = VCF(vcf_file)
        sample_ind = st.sidebar.selectbox("Select sample",vcf_ind.samples)
        use_real_vcf = True
else:
    vcf_mom = st.sidebar.file_uploader("Upload Mother VCF",type=["vcf","vcf.gz"])
    vcf_dad = st.sidebar.file_uploader("Upload Father VCF",type=["vcf","vcf.gz"])
    if vcf_mom and VCF:
        vcf_m=VCF(vcf_mom); sample_mom=st.sidebar.selectbox("Mother sample",vcf_m.samples); use_real_vcf=True
    if vcf_dad and VCF:
        vcf_f=VCF(vcf_dad); sample_dad=st.sidebar.selectbox("Father sample",vcf_f.samples); use_real_vcf=True

selected = st.multiselect("Select traits:",list(traits_info.keys()))

# Summary table
if selected:
    if page=="Individual":
        st.subheader("Individual Summary")
        rows=[]
        for t in selected:
            info=traits_info[t]
            if not info["snps"]:
                summ="-"
            elif t=="Freckles":
                alt=sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                summ=("No freckles" if alt==0 else "Mild freckling" if alt<=2 else "Pronounced")
            elif t=="Hair Colour":
                alt=sum(get_genotype(s,"ind")[0].count(1) for s in info["snps"])
                summ=("Non-red" if alt==0 else "Auburn" if alt==1 else "True red")
            else:
                gt,_,_=get_genotype(info["snps"][0],"ind")
                pres=trait_present(gt,info["inheritance"])
                summ=("Present" if pres else "Absent") if pres is not None else "-"
            rows.append({"Trait":t,"Summary":summ})
        st.table(pd.DataFrame(rows))
    else:
        st.subheader("Child Predictor Summary")
        rows=[]
        for t in selected:
            info=traits_info[t]
            if not info["snps"]:
                rows.append({"Trait":t,"Present (%)":"-","Absent (%)":"-"})
            else:
                m_gt,_,_=get_genotype(info["snps"][0],"mom")
                f_gt,_,_=get_genotype(info["snps"][0],"dad")
                p=child_genotype_probs(m_gt,f_gt)
                pres_pct=sum(x["pct"] for x in p if trait_present(x["geno"],info["inheritance"]))
                rows.append({"Trait":t,"Present (%)":f"{pres_pct:.1f}","Absent (%)":f"{100-pres_pct:.1f}"})
        st.table(pd.DataFrame(rows))

# Height calculator on predictor page
if page=="Child Phenome Predictor" and "Height" in selected:
    st.subheader("Height Calculator")
    c1,c2=st.columns(2)
    with c1:
        mom_cm=st.slider("Mother’s height (cm)",140,200,165)
        dad_cm=st.slider("Father’s height (cm)",140,200,180)
    with c2:
        ft_m,in_m=cm_to_ftin(mom_cm)
        ft_d,in_d=cm_to_ftin(dad_cm)
        st.write(f"Mother: {ft_m} ft {in_m} in"); st.write(f"Father: {ft_d} ft {in_d} in")
    sex=st.selectbox("Child’s sex:",["Male","Female"])
    mean_h=estimate_child_height(mom_cm,dad_cm,sex)
    sigma=4.7; low,high=mean_h-1.96*sigma,mean_h+1.96*sigma
    ft_l,in_l=cm_to_ftin(low); ft_h,in_h=cm_to_ftin(high)
    st.markdown(f"**Predicted height**: {mean_h:.1f} cm")
    st.markdown(f"_95% interval_: {low:.1f}–{high:.1f} cm "
                f"(~{ft_l} ft {in_l} in to {ft_h} ft {in_h} in)")
    sims=np.random.normal(mean_h,sigma,3000)
    df=pd.DataFrame({"Height (cm)":sims})
    ch=alt.Chart(df).mark_area(opacity=0.4).encode(
        alt.X("Height (cm):Q",bin=alt.Bin(maxbins=60)),
        alt.Y("count()",stack=None)
    ).properties(height=250,width=600)
    st.altair_chart(ch)
    st.markdown("---")

# Detailed trait sections
for trait in selected:
    if page=="Individual" and trait=="Height":
        continue
    info=traits_info[trait]
    with st.expander(trait,expanded=True):
        st.subheader("Trait Gene Summary")
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")
        st.write(info["description"])
        if trait=="Hair Colour":
            st.subheader("Trait Interpretation"); st.write(info["interpretation"])
        if info["snps"]:
            st.subheader("Genotypes & Inheritance")
            for snp in info["snps"]:
                if page=="Individual":
                    gt,ref,alt=get_genotype(snp,"ind")
                    b,a=display_genotype(gt,ref,alt)
                    zg=zygosity(gt); pres,mode=format_presence(gt,info["inheritance"])
                    st.markdown(f"**{snp}** (REF={ref}, ALT={alt})")
                    st.write(f"- Genotype: {b} → {a}, {zg}, {pres} {mode}")
                else:
                    m_gt,m_ref,m_alt=get_genotype(snp,"mom")
                    f_gt,f_ref,f_alt=get_genotype(snp,"dad")
                    m_b,m_a=display_genotype(m_gt,m_ref,m_alt)
                    f_b,f_a=display_genotype(f_gt,f_ref,f_alt)
                    m_zg,f_zg=zygosity(m_gt),zygosity(f_gt)
                    m_pres,m_mode=format_presence(m_gt,info["inheritance"])
                    f_pres,f_mode=format_presence(f_gt,info["inheritance"])
                    st.markdown(f"**{snp}** (REF={m_ref}, ALT={m_alt})")
                    st.write(f"- Mother: {m_b} → {m_a}, {m_zg}, {m_pres} {m_mode}")
                    st.write(f"- Father: {f_b} → {f_a}, {f_zg}, {f_pres} {f_mode}")
                    with st.expander("Annotations",expanded=False):
                        st.write(f"- ClinVar: {fetch_clinvar(snp)}")
                        gaf,pops=fetch_gnomad(snp)
                        if gaf is not None:
                            st.write(f"- gnomAD AF: {gaf:.4f}")
                            for pop,af in pops.items(): st.write(f"  - {pop}: {af:.4f}")
                        else: st.write("- gnomAD unavailable")
                    st.subheader("Predicted Child Genotype Probabilities")
                    for p in child_genotype_probs(m_gt,f_gt):
                        cb=f"{p['geno'][0]}/{p['geno'][1]}"
                        ca=alleles_from_gt(p["geno"],m_ref,m_alt)
                        pr,md=format_presence(p["geno"],info["inheritance"])
                        st.write(f"- {cb} ({ca}): {p['pct']:.0f}% → {p['zygosity']}, {pr} {md}")
            # Summary below each trait expander
            st.subheader("Summary")
            # (Insert your summary logic here as before)
        else:
            st.write("_No defined SNPs for this trait._")
        st.markdown("---")

# Cleanup VCF objects
if use_real_vcf:
    del vcf_ind, vcf_m, vcf_f
