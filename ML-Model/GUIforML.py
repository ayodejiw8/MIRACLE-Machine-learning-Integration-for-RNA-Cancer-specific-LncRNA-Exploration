import streamlit as st
import pandas as pd
import joblib
import numpy as np
import os

# 1. PATH CONFIGURATION
base_dir = os.path.dirname(__file__)
model_path = os.path.join(base_dir, 'os_lncrna_detector.pkl')

# 2. LOAD THE MODEL
@st.cache_resource
def load_model():
    if os.path.exists(model_path):
        return joblib.load(model_path)
    return None

model = load_model()

# 3. GUI BRANDING
st.set_page_config(page_title="MIRACLE Multi-Omic Hub", page_icon="🧬", layout="wide")

st.title("🧬 MIRACLE: Multi-Omic Osteosarcoma Detector")
st.markdown("### Albany State University Undergraduate Research 2025-2026")
st.markdown("**Principal Investigator:** Ayodeji Williams | **Mentors:** Dr. Oojo & Dr. Hu")

st.write("---")

# 4. SIDEBAR WITH DYNAMIC FILE INFO
st.sidebar.header("Data Upload Center")
st.sidebar.markdown("""
**Supported GDC Formats:**
* 📂 **RNA-Seq:** STAR-Aligned Augmented Counts
* 📂 **Genomics:** ASCAT Copy Number Variation
* 📂 **Epigenomics:** DNA Methylation Beta Values
""")

uploaded_file = st.sidebar.file_uploader("Upload Profile (.tsv)", type="tsv")

# 5. UNIVERSAL PROCESSING ENGINE
if model is None:
    st.error("Missing model file (os_lncrna_detector.pkl).")

elif uploaded_file is not None:
    try:
        # Detect where the data table starts (Skip metadata headers)
        uploaded_file.seek(0)
        lines = uploaded_file.readlines()
        skip_count = 0
        for i, line in enumerate(lines):
            line_str = line.decode('utf-8')
            # Look for common GDC column headers
            if any(key in line_str for key in ['gene_id', 'gene_name', 'composite element ref', 'gene symbol']):
                skip_count = i
                break
        
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep='\t', skiprows=skip_count)
        df.columns = [c.lower().strip() for c in df.columns]

        # --- DYNAMIC DATA TYPE DETECTION ---
        file_category = "Unknown"
        val_col = None
        gene_col = None

        if 'tpm_unstranded' in df.columns:
            file_category = "Transcriptome (RNA-Seq)"
            val_col = 'tpm_unstranded'
            gene_col = 'gene_name'
        elif 'copy_number' in df.columns:
            file_category = "Genome (Copy Number Variation)"
            val_col = 'copy_number'
            gene_col = 'gene symbol'
        elif 'beta_value' in df.columns or 'value' in df.columns:
            file_category = "Epigenome (DNA Methylation)"
            val_col = 'beta_value' if 'beta_value' in df.columns else 'value'
            gene_col = df.columns[0] # Usually the probe ID or Gene Symbol

        # --- FEATURE EXTRACTION ---
        # Robust search for MALAT1 and NEAT1 across different naming conventions
        m1_row = df[df[gene_col].astype(str).str.contains('MALAT1', na=False, case=False)]
        n1_row = df[df[gene_col].astype(str).str.contains('NEAT1', na=False, case=False)]

        if not m1_row.empty and not n1_row.empty:
            m1_val = float(m1_row[val_col].values[0])
            n1_val = float(n1_row[val_col].values[0])

            # Display Analysis Dashboard
            st.info(f"**Detected Data Source:** {file_category}")
            
            main_col1, main_col2 = st.columns([2, 3])
            
            with main_col1:
                st.subheader("Key lncRNA Markers")
                st.metric(label="MALAT1 Signature", value=f"{m1_val:.4f}")
                st.metric(label="NEAT1 Signature", value=f"{n1_val:.4f}")

            with main_col2:
                # Run Prediction
                features = np.array([[m1_val, n1_val]])
                prob = model.predict_proba(features)[0][1]
                
                st.subheader("AI Diagnostic Intelligence")
                if prob > 0.75:
                    st.error(f"🚨 **HIGH RISK**: {prob*100:.1f}% OS Probability")
                    st.write("Current profile closely matches established Osteosarcoma genomic signatures.")
                elif prob > 0.4:
                    st.warning(f"⚠️ **MODERATE RISK**: {prob*100:.1f}% Probability")
                    st.write("Suggests early-stage genomic rearrangement or borderline expression.")
                else:
                    st.success(f"✅ **LOW RISK**: {prob*100:.1f}% Probability")
                    st.write("Profile aligns with healthy/control genomic background.")
                
                st.progress(prob)

        else:
            st.warning("⚠️ Target markers (MALAT1/NEAT1) not identified. Ensure the file contains standardized gene symbols.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
else:
    st.info("Ready for analysis. Please upload a genomic data file in the sidebar.")

# 6. RESEARCH FOOTER
st.markdown("---")
st.caption("MIRACLE Project v1.3 | Albany State University | Multi-Omic Integration for Osteosarcoma Exploration")
