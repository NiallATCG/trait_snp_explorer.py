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
