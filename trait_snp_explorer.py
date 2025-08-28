import streamlit as st
from collections import Counter

# 1. Trait definitions with summaries and inheritance mode
traits_info = {
    "Freckles": {
        "gene": "MC1R",
        "snps": ["rs1805007", "rs1805008"],
        "description": (
            "MC1R encodes the melanocortin-1 receptor, which regulates the switch "
            "between eumelanin (brown/black) and pheomelanin (red/yellow). Variants "
            "rs1805007 (Arg151Cys) and rs1805008 (Arg160Trp) are associated with "
            "increased freckling. Heterozygotes often have mild freckling, while "
            "homozygous carriers display pronounced freckling."
        ),
        "inheritance": "dominant"
    },
    "Eye Colour": {
        "gene": "HERC2",
        "snps": ["rs12913832"],
        "description": (
            "HERC2 regulates OCA2 expression, influencing melanin in the iris. "
            "At rs12913832, the G allele reduces OCA2 activity and is associated with blue eyes "
            "(recessive), whereas the A allele leads to brown eyes (dominant)."
        ),
        "inheritance": "recessive"
    },
    "Lactose Intolerance": {
        "gene": "MCM6",
        "snps": ["rs4988235"],
        "description": (
            "Variants in the MCM6 enhancer (rs4988235 C→T) control lactase gene persistence. "
            "The T allele maintains lactase into adulthood (dominant); homozygous CC "
            "individuals typically lose lactase activity and are lactose intolerant."
        ),
        "inheritance": "dominant"
    },
    "Sprint Gene": {
        "gene": "ACTN3",
        "snps": ["rs1815739"],
        "description": (
            "ACTN3 encodes α-actinin-3 in fast-twitch muscle fibres. "
            "The rs1815739 T allele introduces a stop codon, leading to deficiency. "
            "CC or CT genotypes have normal sprint performance (dominant), while "
            "TT homozygotes lack α-actinin-3."
        ),
        "inheritance": "dominant"
    },
    # Add other traits here…
}

# Mock SNP data for testing (0 = REF allele, 1 = ALT allele)
mock_vcf_data = {
    # Visual Traits
    "rs1805007": {"mother": [0, 1], "father": [1, 1]},   # MC1R – Freckles (Arg151Cys)
    "rs1805008": {"mother": [0, 0], "father": [0, 1]},   # MC1R – Freckles (Arg160Trp)
    "rs104894":  {"mother": [0, 1], "father": [0, 0]},   # OPN1LW – Red/Green colourblindness
    "rs12913832":{"mother": [0, 0], "father": [0, 1]},   # HERC2 – Eye colour
    # Hair colour – can reuse MC1R variants if desired
    "rs1426654": {"mother": [1, 0], "father": [0, 0]},   # SLC24A5 – Skin tone
    "rs17822931":{"mother": [1, 1], "father": [0, 1]},   # ABCC11 – Earwax type
    # Lifestyle Traits
    "rs4988235": {"mother": [1, 1], "father": [0, 1]},   # MCM6 – Lactose intolerance
    "rs713598": {"mother": [0, 1], "father": [1, 1]},    # TAS2R38 – PTC tasting
    "rs1726866": {"mother": [0, 0], "father": [1, 0]},   # TAS2R38 – PTC tasting
    "rs72921001":{"mother": [1, 0], "father": [1, 1]},   # OR6A2 – Coriander taste
    "rs1815739":{"mother": [0, 0], "father": [1, 1]},    # ACTN3 – Sprint gene
    "rs671":     {"mother": [0, 1], "father": [0, 0]}     # ALDH2 – Alcohol flush
}


# Helper functions
def zygosity(gt):
    return "Homozygous" if gt[0] == gt[1] else "Heterozygous"

def trait_present(gt, inheritance):
    has_alt = any(allele == 1 for allele in gt)
    if inheritance == "dominant":
        return has_alt
    else:  # recessive
        return gt[0] == 1 and gt[1] == 1

def format_presence(gt, inheritance):
    present = trait_present(gt, inheritance)
    return ("Trait present", f"({inheritance})") if present else ("Trait absent", "")

def child_genotype_probs(m_gt, f_gt):
    # Build all possible allele pairings
    combos = [(m, f) for m in m_gt for f in f_gt]
    counts = Counter(tuple(sorted(c)) for c in combos)
    total = len(combos)
    probs = []
    for geno, cnt in sorted(counts.items()):
        percent = cnt / total * 100
        zyg = "Homozygous" if geno[0] == geno[1] else "Heterozygous"
        probs.append({
            "geno_str": f"{geno[0]}/{geno[1]}",
            "percent": percent,
            "zygosity": zyg,
            "alleles": geno
        })
    return probs

# Streamlit UI
st.title("🧬 Enhanced Trait-Based SNP Explorer")

selected = st.multiselect("Choose traits:", list(traits_info.keys()))
st.write("Using mock genotypes for Mother & Father (0 = REF, 1 = ALT).")

for trait in selected:
    info = traits_info[trait]
    st.header(trait)
    
    # Trait Gene Summary
    st.subheader("Trait Gene Summary")
    st.write(f"**Gene**: {info['gene']}")
    st.write(info["description"])
    
    # SNP details
    st.subheader("SNP Genotypes & Inheritance")
    for snp in info["snps"]:
        data = mock_vcf_data.get(snp)
        if not data:
            st.write(f"- {snp}: no data")
            continue
        
        m_gt = data["mother"]
        f_gt = data["father"]
        m_zyg = zygosity(m_gt)
        f_zyg = zygosity(f_gt)
        m_pres, m_mode = format_presence(m_gt, info["inheritance"])
        f_pres, f_mode = format_presence(f_gt, info["inheritance"])
        
        st.markdown(f"**SNP**: {snp}")
        st.write(f"- Mother: {m_gt} → {m_zyg}, {m_pres} {m_mode}")
        st.write(f"- Father: {f_gt} → {f_zyg}, {f_pres} {f_mode}")
        
        # Child probabilities
        st.write("**Predicted Child Genotype Probabilities**")
        for p in child_genotype_probs(m_gt, f_gt):
            c_pres, c_mode = format_presence(p["alleles"], info["inheritance"])
            st.write(
                f"- {p['geno_str']}: {p['percent']:.0f}% → "
                f"{p['zygosity']}, {c_pres} {c_mode}"
            )
        st.markdown("---")
