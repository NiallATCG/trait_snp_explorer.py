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
    "Tanning Response": {
        "gene": "MC1R; ASIP; IRF4; HERC2/OCA2; SLC45A2; TYR",
        "snps": ["rs1805007", "rs12203592", "rs12913832", "rs16891982", "rs1042602"],
        "description": (
            "Variants in pigmentation genes influence tanning ability. MC1R and IRF4 variants "
            "reduce tanning capacity; HERC2/OCA2 and SLC45A2 variants affect melanin type and amount."
        ),
        "inheritance": "polygenic",
        "overview": (
            "Tanning response reflects the skin’s ability to produce melanin after UV exposure. "
            "MC1R and IRF4 variants predispose to burning rather than tanning. HERC2/OCA2 rs12913832 "
            "and SLC45A2 rs16891982 influence melanin levels, while TYR rs1042602 affects pigment synthesis. "
            "Together, these variants determine tanning vs. burning tendency."
        )
    },
    "Lentigines (Sun Spots)": {
        "gene": "MC1R; HERC2/OCA2; ASIP",
        "snps": ["rs1805007", "rs12913832", "rs6058017"],
        "description": (
            "Pigmentation SNPs predispose to solar lentigines, persistent pigmented lesions caused by UV exposure."
        ),
        "inheritance": "polygenic",
        "overview": (
            "Lentigines, or sun spots, are darker patches that persist with age and UV exposure. "
            "MC1R variants increase susceptibility, while HERC2/OCA2 rs12913832 and ASIP rs6058017 "
            "modulate pigmentation. These variants interact with environmental UV exposure to influence risk."
        )
    },

    "Wrinkle & Collagen Degradation": {
        "gene": "MMP1; MMP16; COL17A1; SOD2",
        "snps": ["rs1799750", "rs6469206", "rs805698", "rs4880"],
        "description": (
            "Variants in collagen metabolism and oxidative stress genes influence wrinkle formation. "
            "MMP1 rs1799750 increases collagenase activity; SOD2 rs4880 alters oxidative stress handling."
        ),
        "inheritance": "polygenic",
        "overview": (
            "Wrinkle formation and collagen degradation are influenced by both UV exposure and genetics. "
            "MMP1 rs1799750 increases collagen breakdown, MMP16 and COL17A1 variants affect dermal structure, "
            "and SOD2 rs4880 alters oxidative stress resilience. These variants predispose to earlier or more "
            "pronounced wrinkling."
        )
    },
    "Stretch Marks (Striae Distensae)": {
        "gene": "ELN; FBN1; COL1A1; HMCN1",
        "snps": ["rs3757587", "rs2118181", "rs1800012", "rs7999168"],
        "description": (
            "Variants in connective tissue genes predispose to dermal tearing and striae formation."
        ),
        "inheritance": "polygenic",
        "overview": (
            "Stretch marks occur when dermal connective tissue is disrupted. ELN (elastin) and FBN1 (fibrillin‑1) "
            "variants reduce elasticity, COL1A1 rs1800012 affects collagen strength, and HMCN1 variants influence "
            "dermal resilience. These genetic factors increase susceptibility to striae under mechanical stress."
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
        
      # 🫀 Cardiology
    },
    "Warfarin response": {
        "gene": "VKORC1; CYP2C9; CYP4F2",
        "snps": ["rs9923231", "CYP2C9*2", "CYP2C9*3", "rs2108622"],
        "description": (
            "VKORC1 rs9923231 reduces VKORC1 expression (higher sensitivity); "
            "CYP2C9*2/*3 reduce warfarin clearance (higher exposure); "
            "CYP4F2 rs2108622 reduces vitamin K oxidation (slightly higher dose)."
        ),
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Warfarin dosing varies widely. VKORC1 variants increase sensitivity, CYP2C9 loss‑of‑function slows "
            "clearance and raises bleeding risk, and CYP4F2 variants can increase dose requirements. Combined "
            "genotyping guides initial dose and reduces adverse events."
        )
    },
    "Statin myopathy risk": {
        "gene": "SLCO1B1; ABCG2",
        "snps": ["rs4149056", "rs2231142"],
        "description": (
            "SLCO1B1 rs4149056 (Val174Ala) lowers hepatic uptake → higher plasma statin and myopathy risk; "
            "ABCG2 rs2231142 reduces efflux → increased exposure."
        ),
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Statin intolerance often stems from higher systemic exposure. SLCO1B1 and ABCG2 variants raise blood "
            "levels and myopathy risk, especially with simvastatin. Results can inform statin choice or lower dosing."
        )
    },
    "Clopidogrel response": {
        "gene": "CYP2C19",
        "snps": ["CYP2C19*2", "CYP2C19*3", "CYP2C19*17"],
        "description": (
            "CYP2C19 loss‑of‑function (*2, *3) reduces clopidogrel activation → reduced antiplatelet effect; "
            "*17 increases activity → higher bleeding risk."
        ),
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Clopidogrel is a prodrug requiring CYP2C19 activation. Poor metabolisers have higher thrombotic risk; "
            "ultrarapid metabolisers may bleed more. Genotyping supports switching to prasugrel/ticagrelor or dose changes."
        )
    },
    "Dabigatran activation": {
        "gene": "CES1",
        "snps": ["CES1 variants"],
        "description": (
            "CES1 encodes carboxylesterase‑1 that activates dabigatran etexilate; functional variants alter conversion and exposure."
        ),
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Dabigatran requires CES1 activation. Variants may change active drug levels and bleeding risk. Testing is emerging; "
            "clinical use is context‑dependent."
        )
    },

    # 🧠 Psychiatry & Neurology
    "Opioid analgesic response": {
        "gene": "CYP2D6",
        "snps": ["CYP2D6*3", "CYP2D6*4", "CYP2D6*5", "CYP2D6*6", "copy number"],
        "description": (
            "CYP2D6 status drives activation/clearance of codeine, tramadol, hydrocodone, oxycodone. Poor → lack of efficacy; "
            "ultrarapid → toxicity risk."
        ),
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Many opioids are prodrugs or CYP2D6‑dependent. Poor metabolisers get little pain relief; ultrarapid metabolisers "
            "can reach dangerous levels. Genotyping prevents treatment failure or overdose."
        )
    },
    "Atomoxetine response": {
        "gene": "CYP2D6",
        "snps": ["CYP2D6 variants"],
        "description": "CYP2D6 poor metabolisers have higher atomoxetine levels → more side effects; ultrarapid may need higher doses.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Atomoxetine clearance depends on CYP2D6. Dose reduction is advised in poor metabolisers; therapeutic failure can "
            "occur in ultrarapid metabolisers."
        )
    },
    "Tricyclic antidepressant response": {
        "gene": "CYP2D6; CYP2C19",
        "snps": ["CYP2D6 variants", "CYP2C19 variants"],
        "description": "CYP2D6 and CYP2C19 influence amitriptyline/nortriptyline/imipramine levels → adjust dose or switch.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "TCAs have narrow therapeutic windows. Poor metabolisers risk toxicity; ultrarapid may not respond. Combined genotypes "
            "guide dosing and drug choice."
        )
    },
    "SSRI response": {
        "gene": "CYP2D6; CYP2C19",
        "snps": ["CYP2D6 variants", "CYP2C19 variants"],
        "description": "Metaboliser status affects paroxetine, fluoxetine, sertraline, citalopram, escitalopram exposure and tolerability.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "SSRI efficacy and side effects track with CYP2D6/CYP2C19 activity. Genotyping supports dose adjustments or choosing "
            "alternatives with better metabolic fit."
        )
    },
    "Carbamazepine hypersensitivity": {
        "gene": "HLA-B; HLA-A",
        "snps": ["HLA-B*15:02", "HLA-A*31:01"],
        "description": "HLA-B*15:02 and HLA-A*31:01 associate with SJS/TEN risk on carbamazepine/oxcarbazepine.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Certain HLA alleles markedly increase risk of life‑threatening skin reactions. Screening before therapy prevents SJS/TEN, "
            "especially in at‑risk ancestries."
        )
    },
    "Phenytoin toxicity risk": {
        "gene": "CYP2C9; HLA-B",
        "snps": ["CYP2C9*2", "CYP2C9*3", "HLA-B*15:02"],
        "description": "CYP2C9 loss‑of‑function elevates phenytoin levels; HLA‑B*15:02 raises SJS/TEN risk.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Phenytoin has non‑linear kinetics; reduced CYP2C9 activity quickly leads to toxicity. Concurrent HLA risk mandates "
            "alternative therapy or careful dosing."
        )
    },
    "Valproic acid and POLG": {
        "gene": "POLG",
        "snps": ["POLG mutations"],
        "description": "POLG pathogenic variants predispose to valproate‑induced liver failure/encephalopathy.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "POLG mutations impair mitochondrial replication. Valproate can precipitate catastrophic liver failure; avoid in suspected "
            "mitochondrial disease."
        )
    },
    "Siponimod contraindication": {
        "gene": "CYP2C9",
        "snps": ["CYP2C9*2", "CYP2C9*3"],
        "description": "Siponimod is contraindicated in CYP2C9 poor metabolisers (*3/*3); dose adjustments for other genotypes.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "CYP2C9 genotype dictates siponimod exposure. Poor metabolisers are at higher risk; label recommends genotype‑guided use."
        )
    },

    # 🧬 Oncology
    "Fluoropyrimidine toxicity risk": {
        "gene": "DPYD",
        "snps": ["DPYD*2A", "DPYD*13", "rs67376798", "rs75017182"],
        "description": "Pathogenic DPYD variants reduce DPD activity → severe 5‑FU/capecitabine toxicity. Reduce dose or avoid.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "DPYD encodes dihydropyrimidine dehydrogenase. Deficiency causes life‑threatening toxicity with fluoropyrimidines. "
            "Pre‑treatment genotyping guides dose or alternative therapy."
        )
    },
    "Irinotecan toxicity risk": {
        "gene": "UGT1A1",
        "snps": ["UGT1A1*28"],
        "description": "UGT1A1*28 (TA repeat) reduces glucuronidation of SN‑38 → neutropenia risk; lower starting dose if homozygous.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Irinotecan’s active metabolite SN‑38 is cleared by UGT1A1. Reduced activity elevates toxicity. Testing identifies "
            "patients who benefit from dose reduction."
        )
    },
    "Tamoxifen efficacy": {
        "gene": "CYP2D6",
        "snps": ["CYP2D6 variants"],
        "description": "CYP2D6 poor metabolisers have reduced endoxifen formation → potentially reduced tamoxifen efficacy.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Tamoxifen requires CYP2D6 to form active metabolites. Low activity may compromise benefit; consider alternative endocrine "
            "therapy or careful monitoring."
        )
    },
    "Thiopurine toxicity risk": {
        "gene": "TPMT; NUDT15",
        "snps": ["TPMT activity alleles", "rs116855232"],
        "description": "Low TPMT/NUDT15 activity causes severe myelosuppression with thiopurines → reduce dose or avoid.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "TPMT and NUDT15 variants impair thiopurine inactivation. Even standard doses can be dangerous; phenotype/genotype "
            "guides safe dosing."
        )
    },
    "Anthracycline cardiotoxicity markers": {
        "gene": "RARG; SLC28A3",
        "snps": ["RARG variants", "SLC28A3 variants"],
        "description": "Genetic markers associated with higher cardiotoxicity risk with doxorubicin/daunorubicin.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Anthracyclines can damage the heart. Certain variants may increase risk; results support intensified monitoring or "
            "risk‑mitigation strategies."
        )
    },

    # 🦠 Infectious Disease
    "Abacavir hypersensitivity": {
        "gene": "HLA-B",
        "snps": ["HLA-B*57:01"],
        "description": "HLA‑B*57:01 confers high risk of abacavir hypersensitivity → contraindicated if positive.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Pre‑treatment HLA‑B*57:01 testing is standard for abacavir. Positive patients should not receive the drug."
        )
    },
    "Allopurinol severe skin reaction risk": {
        "gene": "HLA-B",
        "snps": ["HLA-B*58:01"],
        "description": "HLA‑B*58:01 increases risk of SCAR (SJS/TEN) with allopurinol → avoid if positive.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Screening for HLA‑B*58:01 helps prevent life‑threatening reactions, especially in high‑prevalence populations."
        )
    },
    "Flucloxacillin liver injury risk": {
        "gene": "HLA-B",
        "snps": ["HLA-B*57:01"],
        "description": "HLA‑B*57:01 associates with higher risk of flucloxacillin‑induced liver injury.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Genetic risk stratification can inform vigilance and alternative antibiotics in susceptible patients."
        )
    },
    "Efavirenz exposure": {
        "gene": "CYP2B6",
        "snps": ["rs3745274"],
        "description": "CYP2B6 rs3745274 reduces clearance → higher efavirenz levels and CNS side effects; dose reduction may help.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Efavirenz neuropsychiatric side effects correlate with CYP2B6 genotype. Lower doses can improve tolerability in poor metabolisers."
        )
    },
    "Atazanavir hyperbilirubinaemia": {
        "gene": "UGT1A1",
        "snps": ["UGT1A1*28"],
        "description": "UGT1A1*28 increases indirect bilirubin with atazanavir → consider monitoring or alternative agent.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Benign jaundice from atazanavir is common with UGT1A1*28; genotype informs expectations and switching decisions."
        )
    },
    "Voriconazole dosing": {
        "gene": "CYP2C19",
        "snps": ["CYP2C19*2", "CYP2C19*3", "CYP2C19*17"],
        "description": "CYP2C19 genotype strongly affects voriconazole exposure → adjust dose to avoid toxicity or failure.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Therapeutic drug monitoring plus genotype gives the best outcome; poor metabolisers need lower doses, ultrarapid may need higher."
        )
    },

    # 🧪 Immunosuppression & Transplant
    "Tacrolimus dosing": {
        "gene": "CYP3A5",
        "snps": ["CYP3A5*3 (rs776746)"],
        "description": "CYP3A5*3 non‑expressors have lower tacrolimus clearance → lower dose needed to reach target troughs.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Tacrolimus is highly variable. CYP3A5 expressors require higher doses; non‑expressors need less to achieve therapeutic levels."
        )
    },
    "Thiopurine dosing (transplant)": {
        "gene": "TPMT; NUDT15",
        "snps": ["TPMT activity alleles", "rs116855232"],
        "description": "Low TPMT/NUDT15 activity → severe myelosuppression with azathioprine/6‑MP/thioguanine; reduce dose or avoid.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Routine TPMT/NUDT15 testing prevents marrow toxicity. Dose individualisation is standard of care in many centres."
        )
    },
    "Mycophenolate response (research)": {
        "gene": "IMPDH1; IMPDH2",
        "snps": ["IMPDH variants"],
        "description": "IMPDH variants may influence mycophenolate efficacy/toxicity; evidence is emerging.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Mycophenolate targets IMPDH; genetic influences are under study and not yet used routinely in clinical decision‑making."
        )
    },

    # 🚬 Smoking Cessation
    "Smoking cessation pharmacogenetics": {
        "gene": "CYP2A6; CHRNA5",
        "snps": ["rs16969968", "CYP2A6 activity alleles"],
        "description": "CHRNA5 rs16969968 associates with nicotine dependence; CYP2A6 activity affects nicotine clearance and cessation outcomes.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "Genetics shape nicotine dependence and response to therapies. CHRNA5 risk alleles and slow CYP2A6 clearance may guide "
            "choice of varenicline, bupropion, or nicotine replacement."
        )
    },
    "Bupropion dosing": {
        "gene": "CYP2B6",
        "snps": ["rs3745274"],
        "description": "CYP2B6 poor metabolisers have higher bupropion exposure → adjust dose or monitor side effects.",
        "inheritance": "Pharmacogenetic",
        "overview": (
            "CYP2B6 genotype influences bupropion levels for depression and smoking cessation. Dose tailoring can improve tolerability."
 
        )

    }
}
# 2. Mock genotype data
# 2. Mock genotype data
mock_vcf_data = {
    snp: {
        "ref": d["ref"], "alt": d["alt"],
        "mother": d["mother"], "father": d["father"], "gt": d["gt"]
    }
    for snp, d in {
        # ── Original traits ──
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
        "rs12203592": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # IRF4
        # ── Placeholder entries for symbolic/star alleles ──
        "CYP2C9*2": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,1]},
        "CYP2C9*3": {"ref":"A","alt":"C","mother":[0,1],"father":[0,0],"gt":[0,1]},  # already have rs1057910, but keep alias
        "CYP2C19*2": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "CYP2C19*3": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
        "CYP2C19*17": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "CYP2D6*3": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "CYP2D6*4": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
        "CYP2D6*5": {"ref":"del","alt":"-","mother":[0,0],"father":[0,1],"gt":[0,0]},  # deletion allele
        "CYP2D6*6": {"ref":"T","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "CYP2A6 activity alleles": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
        "TPMT activity alleles": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "IMPDH variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
        "RARG variants": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "SLC28A3 variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
        "CES1 variants": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
        "POLG mutations": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},

        # ── Pharmacogenetic traits ──
        # Cardiology
        "rs9923231":  {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},  # VKORC1
        "rs2108622":  {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CYP4F2
        "rs4149056":  {"ref":"T","alt":"C","mother":[0,1],"father":[0,0],"gt":[0,1]},  # SLCO1B1
        "rs2231142":  {"ref":"G","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # ABCG2

        # Psychiatry & Neurology
        "rs3892097":  {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},  # CYP2D6*4 proxy
        "rs1065852":  {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CYP2D6*10 proxy
        "rs4244285":  {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},  # CYP2C19*2
        "rs4986893":  {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CYP2C19*3
        "rs12248560": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # CYP2C19*17
        "rs3909184":  {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # HLA-A*31:01 proxy
        "rs10484555": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # HLA-B*15:02 proxy
        "rs121918508":{"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # POLG proxy
        "rs1057910":  {"ref":"A","alt":"C","mother":[0,1],"father":[0,0],"gt":[0,1]},  # CYP2C9*3

        # Oncology
        "rs3918290":  {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # DPYD*2A
        "rs67376798": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # DPYD variant
        "rs75017182": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # DPYD variant
        "rs8175347":  {"ref":"TA6","alt":"TA7","mother":[0,1],"father":[0,0],"gt":[0,1]}, # UGT1A1*28
        "rs116855232":{"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # NUDT15
        "rs1800460":  {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},  # TPMT proxy

        # Infectious Disease
        "rs2395029":  {"ref":"T","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,0]},  # HLA-B*57:01 proxy
        "rs9263726":  {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # HLA-B*58:01 proxy
        "rs3745274":  {"ref":"G","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CYP2B6
        "rs8175347b": {"ref":"TA6","alt":"TA7","mother":[0,1],"father":[0,0],"gt":[0,1]}, # UGT1A1*28 (Atazanavir)

        # Immunosuppression & Transplant
        "rs776746":   {"ref":"A","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CYP3A5*3
        "rs2073838":  {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # IMPDH1 proxy

        # Smoking Cessation
        "rs16969968": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CHRNA5

        # ── Dermatology traits ──
        # Tanning response
        "rs16891982": {"ref":"C","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},  # SLC45A2
        "rs1042602":  {"ref":"C","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # TYR

        # Lentigines
        "rs6058017":  {"ref":"G","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # ASIP

         # Wrinkle & collagen degradation
        "rs1799750":  {"ref":"1","alt":"2","mother":[0,1],"father":[0,0],"gt":[0,1]},  # MMP1 promoter ins/del
        "rs6469206":  {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # MMP16
        "rs805698":   {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},  # COL17A1
        "rs4880":     {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # SOD2

        # Stretch marks
        "rs3757587":  {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # ELN
        "rs2118181":  {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},  # FBN1
        "rs1800012":  {"ref":"G","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},  # COL1A1
        "rs7999168":  {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},  # HMCN1
                                                     
    }.items()

}

# Add missing dermatology SNPs and pharmacogenetic symbolic variants
mock_vcf_data.update({
    "rs12203592": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "rs6058017":  {"ref":"G","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2C9*2": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,1]},
    "CYP2C9*3": {"ref":"A","alt":"C","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2C19*2": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2C19*3": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CYP2C19*17": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2D6*3": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2D6*4": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CYP2D6*5": {"ref":"del","alt":"-","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CYP2D6*6": {"ref":"T","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2D6 variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CYP2A6 activity alleles": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "TPMT activity alleles": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "IMPDH variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "RARG variants": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "SLC28A3 variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CES1 variants": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "POLG mutations": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
})

# Complete mock entries for Psychiatry/Neurology, Oncology, Infectious Disease, Immunosuppression
mock_vcf_data.update({
    # ── Psychiatry & Neurology ──
    "CYP2D6*3": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "CYP2D6*4": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "CYP2D6*5": {"ref":"del","alt":"-","mother":[0,0],"father":[0,1],"gt":[0,0]},  # deletion (null)
    "CYP2D6*6": {"ref":"T","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "copy number": {"ref":"2","alt":"3","mother":[0,0],"father":[0,1],"gt":[0,0]},  # CNV placeholder
    "CYP2D6 variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},

    "CYP2C19 variants": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},  # generic placeholder

    "HLA-B*15:02": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "HLA-A*31:01": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},

    "CYP2C9*2": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,1]},
    "CYP2C9*3": {"ref":"A","alt":"C","mother":[0,1],"father":[0,0],"gt":[0,1]},

    "POLG mutations": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},

    # ── Oncology ──
    "DPYD*2A": {"ref":"G","alt":"A","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "DPYD*13": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},

    "UGT1A1*28": {"ref":"TA6","alt":"TA7","mother":[0,1],"father":[0,0],"gt":[0,1]},

    "TPMT activity alleles": {"ref":"A","alt":"G","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "RARG variants": {"ref":"G","alt":"A","mother":[0,1],"father":[0,0],"gt":[0,1]},
    "SLC28A3 variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},

    # ── Infectious Disease ──
    "HLA-B*57:01": {"ref":"T","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,0]},
    "HLA-B*58:01": {"ref":"C","alt":"T","mother":[0,1],"father":[0,0],"gt":[0,1]},

    # Alias for Atazanavir use case
    "UGT1A1*28 (atazanavir)": {"ref":"TA6","alt":"TA7","mother":[0,1],"father":[0,0],"gt":[0,1]},

    # ── Immunosuppression & Transplant ──
    "CYP3A5*3": {"ref":"A","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,0]},  # alias for rs776746

    "IMPDH variants": {"ref":"C","alt":"T","mother":[0,0],"father":[0,1],"gt":[0,0]},
})

# Ensure the exact key used in traits_info exists in mock_vcf_data
mock_vcf_data.update({
    "CYP3A5*3 (rs776746)": {"ref":"A","alt":"G","mother":[0,0],"father":[0,1],"gt":[0,0]}
})

missing = []
for trait, info in traits_info.items():
    for s in info["snps"]:
        if s not in mock_vcf_data:
            missing.append((trait, s))
print("Missing SNPs:", missing)

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
        "Freckles", "Hair Colour", "Eye Colour", "Skin Tone", "Earwax Type", "Tanning Response",
        "Lentigines (Sun Spots)", "Wrinkle & Collagen Degradation", "Stretch Marks (Striae Distensae)", "Height"
    ],
    "Nutrition and Fitness": [
        "Sprint Gene", "Lactose Intolerance", "PTC Tasting",
        "Coriander Taste", "Red-Green Colourblindness", "Alcohol Flush"
    ],
    "Genetic Response to Drugs": {
        "Cardiology": [
            "Warfarin response",
            "Statin myopathy risk",
            "Clopidogrel response",
            "Dabigatran activation"
        ],
        "Psychiatry & Neurology": [
            "Opioid analgesic response",
            "Atomoxetine response",
            "Tricyclic antidepressant response",
            "SSRI response",
            "Carbamazepine hypersensitivity",
            "Phenytoin toxicity risk",
            "Valproic acid and POLG",
            "Siponimod contraindication"
        ],
        "Oncology": [
            "Fluoropyrimidine toxicity risk",
            "Irinotecan toxicity risk",
            "Tamoxifen efficacy",
            "Thiopurine toxicity risk",
            "Anthracycline cardiotoxicity markers"
        ],
        "Infectious Disease": [
            "Abacavir hypersensitivity",
            "Allopurinol severe skin reaction risk",
            "Flucloxacillin liver injury risk",
            "Efavirenz exposure",
            "Atazanavir hyperbilirubinaemia",
            "Voriconazole dosing"
        ],
        "Immunosuppression & Transplant": [
            "Tacrolimus dosing",
            "Thiopurine dosing (transplant)",
            "Mycophenolate response (research)"
        ],
        "Smoking Cessation": [
            "Smoking cessation pharmacogenetics",
            "Bupropion dosing"
        ]
    }
}
# 2. Update the UI logic

# Dropdown for report selection
report = st.selectbox("Select Report", list(report_groups.keys()), index=1)

if report == "Genetic Response to Drugs":
    subsection = st.selectbox("Select Indication", list(report_groups[report].keys()))
    selected = st.multiselect(
        f"Select traits from {subsection}:",
        report_groups[report][subsection]
    )
else:
    selected = st.multiselect(
        f"Select traits from {report}:",
        report_groups[report]
    )

# ── Overall Summary Table ──
if selected:
    st.subheader("Overall Summary Table")
    rows = []
    for trait in selected:
        info = traits_info[trait]

        # 🔎 Debug: show trait name and genotypes being used
        print("DEBUG:", trait, [get_genotype(s, "ind")[0] for s in info["snps"]])
        
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

        # ── New Dermatology Traits ──

        elif trait == "Tanning Response":
            has_alt = any(1 in get_genotype(s, "ind")[0] for s in info["snps"])
            summary = "Poor tanning / burns easily" if has_alt else "Good tanning ability"

        elif trait == "Lentigines (Sun Spots)":
            has_alt = any(1 in get_genotype(s, "ind")[0] for s in info["snps"])
            summary = "Higher likelihood of sun spots" if has_alt else "Lower likelihood of sun spots"

        elif trait == "Wrinkle & Collagen Degradation":
            has_alt = any(1 in get_genotype(s, "ind")[0] for s in info["snps"])
            summary = "Higher wrinkle susceptibility" if has_alt else "Lower wrinkle susceptibility"

        elif trait == "Stretch Marks (Striae Distensae)":
            has_alt = any(1 in get_genotype(s, "ind")[0] for s in info["snps"])
            summary = "Higher stretch mark susceptibility" if has_alt else "Lower stretch mark susceptibility"

                
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

        # Show plain-language overview first
        if "overview" in info:
            st.subheader("**Overview**")
            st.write(info["overview"])

        # Gene
        if info["gene"]:
            st.write(f"**Gene**: {info['gene']}")

        # Molecular mechanisms
        if info["description"]:    
            st.write(f"**Molecular Mechanisms**: {info['description']}")

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
