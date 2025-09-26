# trait_snp_explorer.py
# Streamlit app: Trait SNP Explorer
import streamlit as st
import pandas as pd

# -----------------------
# Configuration / Mapping
# -----------------------

REPORT_TRAITS = {
    "SNP Associated Diseases (WIP)": [],
    "Nutrition": [
        "Lactose Intolerance",
        "PTC Tasting",
        "Coriander Taste",
        "Red-Green Colourblindness",
        "Alcohol Flush"
    ],
    "Health and Fitness": [
        "Freckles",
        "Hair Colour",
        "Eye Colour",
        "Skin Tone",
        "Earwax Type",
        "Height (non-SNP based)",
        "Sprint Gene"
    ],
    "Genetic Response to Drugs (WIP)": []
}

# A placeholder list of all supported traits (flat)
ALL_TRAITS = []
for v in REPORT_TRAITS.values():
    ALL_TRAITS.extend(v)
ALL_TRAITS = sorted(list(set(ALL_TRAITS)))

# -----------------------
# Helper functions
# -----------------------

def get_user_genotype_for_snp(snp_id, genotype_df):
    """
    Placeholder: lookup genotype_df for genotype string for given snp_id
    genotype_df assumed to have columns ['snp', 'genotype'].
    Returns genotype string or None.
    """
    if genotype_df is None:
        return None
    row = genotype_df[genotype_df['snp'] == snp_id]
    if row.empty:
        return None
    return str(row.iloc[0]['genotype'])

def simple_table_from_summary(summary_rows):
    """Make a simple pandas DataFrame from summary rows (list of dicts)."""
    if not summary_rows:
        return pd.DataFrame(columns=["Trait", "Result", "Notes"])
    return pd.DataFrame(summary_rows)

# -----------------------
# Trait interpretation logic
# -----------------------

def interpret_trait(trait, genotype_df=None):
    """
    Trait-level interpreter. Returns dict with keys: Trait, Result, Notes.
    For SNP-based traits the genotype_df is used. For non-SNP traits there may be
    static or derived text.
    """
    trait_key = trait.strip()
    res = {"Trait": trait_key, "Result": "Unknown", "Notes": ""}

    # Nutrition group
    if trait_key == "Lactose Intolerance":
        # Example SNPs: LCT -13910 C/T (rs4988235). This is placeholder logic.
        gt = get_user_genotype_for_snp("rs4988235", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs4988235 genotype not in uploaded data"
        else:
            # Common encoding: CC = lactose intolerant, CT = intermediate, TT = tolerant
            if gt in ["CC", "C/C"]:
                res["Result"] = "Predicted lactose intolerance"
                res["Notes"] = "Genotype indicates reduced lactase persistence"
            elif gt in ["CT", "C/T", "T/C"]:
                res["Result"] = "Intermediate risk"
                res["Notes"] = "Heterozygous genotype; phenotype may vary"
            elif gt in ["TT", "T/T"]:
                res["Result"] = "Predicted lactose tolerance"
                res["Notes"] = "Genotype indicates lactase persistence"
            else:
                res["Result"] = gt
                res["Notes"] = "Unrecognised genotype format"
        return res

    if trait_key == "PTC Tasting":
        # Example SNP: TAS2R38 haplotype inference (rs713598, rs1726866, rs10246939)
        # Placeholder: check presence of one sentinel SNP
        gt = get_user_genotype_for_snp("rs713598", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs713598 missing"
        else:
            # Using simplified assumption: CC = taster, GG = non-taster (placeholder)
            if gt in ["CC", "C/C"]:
                res["Result"] = "Predicted taster"
                res["Notes"] = "Genotype associated with PTC tasting"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Predicted non-taster"
                res["Notes"] = "Genotype associated with reduced PTC sensitivity"
            else:
                res["Result"] = gt
                res["Notes"] = "Complex haplotype; full haplotype needed for accuracy"
        return res

    if trait_key == "Coriander Taste":
        # Example: olfactory receptor variants; placeholder uses rs72921001
        gt = get_user_genotype_for_snp("rs72921001", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs72921001 missing"
        else:
            if gt in ["AA", "A/A"]:
                res["Result"] = "Predicted soapy perception"
                res["Notes"] = "Alleles associated with coriander tasting as soapy"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Predicted neutral/favourable perception"
                res["Notes"] = "Alleles associated with typical coriander perception"
            else:
                res["Result"] = gt
                res["Notes"] = "Uncertain; actual perception is multifactorial"
        return res

    if trait_key == "Red-Green Colourblindness":
        # Placeholder: Many variants on OPN1LW/OPN1MW; here we flag common absence
        # We'll look for a sentinel SNP as placeholder
        gt = get_user_genotype_for_snp("rs1800414", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "Sentinel SNP missing; colour vision genetics requires complex analysis"
        else:
            # This is illustrative only
            res["Result"] = "Complex analysis required"
            res["Notes"] = "Full gene copy number and sequence haplotype required for reliable call"
        return res

    if trait_key == "Alcohol Flush":
        # ALDH2 rs671 common in East Asian populations
        gt = get_user_genotype_for_snp("rs671", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs671 missing"
        else:
            if gt in ["AA", "A/A"]:
                res["Result"] = "Strong alcohol flush risk"
                res["Notes"] = "ALDH2 inactive variant; increased flush, alcohol sensitivity"
            elif gt in ["AG", "A/G", "GA", "G/A"]:
                res["Result"] = "Intermediate alcohol flush risk"
                res["Notes"] = "Heterozygous reduction in enzyme activity"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Typical response"
                res["Notes"] = "No common ALDH2 flushing variant detected"
            else:
                res["Result"] = gt
                res["Notes"] = "Unrecognised genotype"
        return res

    # Health and Fitness group
    if trait_key == "Freckles":
        # MC1R variants implicated. Placeholder: check rs1805007
        gt = get_user_genotype_for_snp("rs1805007", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs1805007 missing; MC1R haplotype needed for robust prediction"
        else:
            if gt in ["TT", "T/T", "T"]:
                res["Result"] = "Increased freckling likely"
                res["Notes"] = "MC1R variant associated with red hair / freckling phenotype"
            elif gt in ["CC", "C/C"]:
                res["Result"] = "Typical freckling"
                res["Notes"] = "No common MC1R variant at rs1805007"
            else:
                res["Result"] = gt
                res["Notes"] = "Haplotype-level phasing would improve accuracy"
        return res

    if trait_key == "Hair Colour":
        # Placeholder combining MC1R and other loci
        gt = get_user_genotype_for_snp("rs12913832", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "Important pigmentation SNPs missing"
        else:
            if gt in ["AA", "A/A"]:
                res["Result"] = "Higher probability of brown/black hair"
                res["Notes"] = "rs12913832 alleles associated with darker pigmentation"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Higher probability of lighter hair"
                res["Notes"] = "rs12913832 alleles associated with lighter pigmentation"
            else:
                res["Result"] = "Intermediate"
                res["Notes"] = "Complex polygenic trait; prediction is probabilistic"
        return res

    if trait_key == "Eye Colour":
        gt = get_user_genotype_for_snp("rs12913832", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "Key eye colour SNP missing"
        else:
            if gt in ["GG", "G/G"]:
                res["Result"] = "Blue eyes more likely"
                res["Notes"] = "Variant associated with lower pigmentation in iris"
            elif gt in ["AA", "A/A"]:
                res["Result"] = "Brown eyes more likely"
                res["Notes"] = "Variant associated with higher pigmentation in iris"
            else:
                res["Result"] = "Intermediate"
                res["Notes"] = "Eye colour prediction is polygenic"
        return res

    if trait_key == "Skin Tone":
        gt = get_user_genotype_for_snp("rs1426654", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "Key pigmentation SNP missing"
        else:
            if gt in ["AA", "A/A"]:
                res["Result"] = "Lighter skin tone alleles present"
                res["Notes"] = "rs1426654 influences skin pigmentation"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Darker skin tone alleles present"
                res["Notes"] = "Polygenic trait; single SNP gives limited info"
            else:
                res["Result"] = "Intermediate"
                res["Notes"] = "Prediction is probabilistic"
        return res

    if trait_key == "Earwax Type":
        # ABCC11 rs17822931; GG = wet, AA = dry (population-specific)
        gt = get_user_genotype_for_snp("rs17822931", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs17822931 missing"
        else:
            if gt in ["AA", "A/A"]:
                res["Result"] = "Dry earwax likely"
                res["Notes"] = "Allele associated with dry earwax (common in East Asians)"
            elif gt in ["GG", "G/G"]:
                res["Result"] = "Wet earwax likely"
                res["Notes"] = "Allele associated with wet earwax"
            else:
                res["Result"] = "Intermediate"
                res["Notes"] = "Heterozygote; phenotype may vary"
        return res

    if trait_key == "Height (non-SNP based)":
        # Non-SNP based: placeholder uses family history or simple average
        res["Result"] = "Non-SNP based estimate"
        res["Notes"] = "Height prediction not implemented from SNPs in this app"
        return res

    if trait_key == "Sprint Gene":
        # Example: ACTN3 rs1815739 R577X
        gt = get_user_genotype_for_snp("rs1815739", genotype_df)
        if gt is None:
            res["Result"] = "No genotype found"
            res["Notes"] = "rs1815739 missing"
        else:
            if gt in ["CC", "C/C", "R/R"]:
                res["Result"] = "Power/sprint-associated genotype"
                res["Notes"] = "Functional ACTN3 allele present"
            elif gt in ["TT", "T/T", "X/X"]:
                res["Result"] = "Reduced ACTN3 function"
                res["Notes"] = "May be associated with endurance performance traits"
            else:
                res["Result"] = gt
                res["Notes"] = "Heterozygous or ambiguous encoding"
        return res

    # Genetic Response to Drugs group (WIP)
    # Put placeholders here for future expansion
    if trait_key in ("SNP Associated Diseases (WIP)", "Genetic Response to Drugs (WIP)"):
        res["Result"] = "Report category selected"
        res["Notes"] = "Detailed analyses for this report are a work in progress"
        return res

    # Fallback
    res["Result"] = "Trait interpretation not implemented"
    res["Notes"] = "No specific interpretation code for this trait"
    return res

# -----------------------
# Streamlit UI
# -----------------------

def sidebar_upload_section():
    st.sidebar.header("Input data")
    uploaded = st.sidebar.file_uploader("Upload genotype CSV (columns: snp, genotype)", type=["csv", "txt"])
    genotype_df = None
    if uploaded is not None:
        try:
            genotype_df = pd.read_csv(uploaded)
            st.sidebar.success("Genotype file loaded")
        except Exception as e:
            st.sidebar.error(f"Failed to load file: {e}")
            genotype_df = None
    return genotype_df

def individual_tab(genotype_df):
    st.header("Individual report")

    # ----- CHANGED: Report selector replaces the old "select traits" label -----
    # Note: label uses exactly the text you requested: 'select report;:'
    report_choice = st.selectbox("select report;:", list(REPORT_TRAITS.keys()))
    available_traits = REPORT_TRAITS.get(report_choice, [])

    if not available_traits:
        st.info("No traits are currently assigned to this report.")
    # trait multiselect populated from selected report
    selected_traits = st.multiselect("select traits", available_traits, default=available_traits[:1])

    # Run interpretation only for selected traits
    summary_rows = []
    for trait in selected_traits:
        interpretation = interpret_trait(trait, genotype_df=genotype_df)
        summary_rows.append({
            "Trait": interpretation.get("Trait", trait),
            "Result": interpretation.get("Result", ""),
            "Notes": interpretation.get("Notes", "")
        })

    st.subheader("Report summary")
    df_summary = simple_table_from_summary(summary_rows)
    st.table(df_summary)

    # Expanders per trait for full detail
    for row in summary_rows:
        with st.expander(row["Trait"]):
            st.write("Result:", row["Result"])
            st.write("Notes:", row["Notes"])

def population_tab():
    st.header("Population / cohort analyses")
    st.write("This tab is a placeholder for cohort-level features.")
    # Implement cohort features as required

def main():
    st.set_page_config(page_title="Trait SNP Explorer", layout="wide")
    st.title("Trait SNP Explorer")

    genotype_df = sidebar_upload_section()

    tabs = st.tabs(["Individual", "Population"])
    with tabs[0]:
        individual_tab(genotype_df)
    with tabs[1]:
        population_tab()

if __name__ == "__main__":
    main()
