import streamlit as st
import pandas as pd

# Mock SNP data for testing
mock_vcf_data = {
    "rs1805007": {"mother": [0, 1], "father": [1, 1]},  # MC1R - freckles
    "rs12913832": {"mother": [0, 0], "father": [0, 1]},  # HERC2 - eye colour
    "rs4988235": {"mother": [1, 1], "father": [0, 1]},  # MCM6 - lactose intolerance
    "rs1815739": {"mother": [0, 0], "father": [1, 1]},  # ACTN3 - sprint gene
}

trait_snp_map = {
    "Freckles": {"gene": "MC1R", "snps": ["rs1805007"]},
    "Eye Colour": {"gene": "HERC2", "snps": ["rs12913832"]},
    "Lactose Intolerance": {"gene": "MCM6", "snps": ["rs4988235"]},
    "Sprint Gene": {"gene": "ACTN3", "snps": ["rs1815739"]}
}

def predict_child_alleles(mother_gt, father_gt):
    alleles = set(mother_gt + father_gt)
    return sorted(list(alleles))

# Streamlit UI
st.title("🧬 Trait-Based SNP Explorer")
selected_traits = st.multiselect("Select traits to explore:", list(trait_snp_map.keys()))
st.write("Mock parental genotypes loaded for testing.")

if selected_traits:
    st.subheader("Predicted Child Genotypes")
    for trait in selected_traits:
        snps = trait_snp_map[trait]["snps"]
        for snp in snps:
            data = mock_vcf_data.get(snp)
            if data:
                child_gt = predict_child_alleles(data["mother"], data["father"])
                st.markdown(f"**Trait:** {trait}")
                st.write(f"SNP: {snp}")
                st.write(f"Mother GT: {data['mother']}, Father GT: {data['father']}")
                st.write(f"Predicted Child GT: {child_gt}")
                st.markdown("---")

