import streamlit as st
import pandas as pd
import joblib
import numpy as np
import os

# 1. PATH CONFIGURATION
base_dir = os.path.dirname(__file__)
model_path = os.path.join(base_dir, 'os_lncrna_detector (1).pkl')

# 2. LOAD THE MULTI-OMIC MODEL (v1.3)
# This model now understands RNA Expression, CNV, and Methylation
if os.path.exists(model_path):
    model = joblib.load(model_path)
else:
    st.error("Model file (os_lncrna_detector.pkl) not found. Please upload it to the directory.")

# 3. GUI BRANDING & HEADER
st.set_page_config(page_title="MIRACLE Project | ASU", page_icon="🧬")

st.title("MIRACLE Project: Osteosarcoma lncRNA Detector")
st.markdown("### Albany State University Undergraduate Research 2025-2026 Cohort")
st.markdown("**Developed by:** Ayodeji Williams | Albany State University Bioinformatics")
st.markdown("**Mentors:** Dr. Olabisi Ojo & Dr. Wanjun Hu")
st.markdown("*College of Arts and Science | Dept. of Math, CS and Physics & Dept. of Natural Sciences*")

st.write("---")

# 4. SIDEBAR INPUT
st.sidebar.header("Data Input Control")
uploaded_file = st.sidebar.file_uploader("Upload Genomic Profile (.tsv)", type="tsv")
st.sidebar.info("Supported formats: STAR-Aligned Counts, ASCAT Copy Number, or Methylation Beta Values.")

# 5. PROCESSING LOGIC
if uploaded_file is not None:
    try:
        # Load data and standardize column names
        df = pd.read_csv(uploaded_file, sep='\t')
        df.columns = [c.lower() for c in df.columns]
        
        # Determine the Data Type (RNA vs DNA vs Epigenetic)
        # We check for standard GDC column headers
        if 'tpm_unstranded' in df.columns:
            data_type = "RNA-Seq (Expression)"
            val_col = 'tpm_unstranded'
            gene_col = 'gene_name'
        elif 'copy_number' in df.columns:
            data_type = "DNA (Copy Number)"
            val_col = 'copy_number'
            gene_col = 'gene symbol'
        elif 'value' in df.columns or any(df.iloc[:,0].str.startswith('cg', na=False)):
            data_type = "Epigenetic (Methylation)"
            # For methylation, we use a proxy mapping to MALAT1/NEAT1 probes
            val_col = df.columns[1] 
            gene_col = df.columns[0]
        else:
            data_type = "General Genomic"
            val_col = df.columns[1]
            gene_col = df.columns[0]

        # Extract MALAT1 and NEAT1 values
        m1_search = df[df[gene_col].str.contains('MALAT1', na=False, case=False)]
        n1_search = df[df[gene_col].str.contains('NEAT1', na=False, case=False)]

        if not m1_search.empty and not n1_search.empty:
            malat1 = m1_search[val_col].values[0]
            neat1 = n1_search[val_col].values[0]

            # Display Metadata
            st.success(f"Successfully identified {data_type} signature.")
            
            st.write(f"### Analysis for {uploaded_file.name}")
            col1, col2, col3 = st.columns(3)
            col1.metric("MALAT1 Value", f"{malat1:.2f}")
            col2.metric("NEAT1 Value", f"{neat1:.2f}")
            col3.metric("Data Type", data_type)

            # 6. PERFORM PREDICTION
            # Reshape for model input
            sample = np.array([[malat1, neat1]])
            prob = model.predict_proba(sample)[0][1]
            
            st.write("---")
            st.write("## Diagnostic AI Result")
            
            # Gauge-style result display
            if prob > 0.75:
                st.error(f"🚨 **HIGH RISK DETECTED**")
                st.write(f"The Multi-Omic signature indicates an **{prob*100:.1f}%** probability of Osteosarcoma.")
            elif prob > 0.40:
                st.warning(f"⚠️ **MODERATE RISK / BORDERLINE**")
                st.write(f"Signature suggests genomic instability: **{prob*100:.1f}%** probability.")
            else:
                st.success(f"✅ **NORMAL / HEALTHY SIGNATURE**")
                st.write(f"Confidence: **{(1-prob)*100:.1f}%** Normal.")

            # Progress bar for visual impact on poster
            st.progress(prob)
            
        else:
            st.error("Error: Could not find MALAT1 or NEAT1 markers in this file.")

    except Exception as e:
        st.error(f"File Processing Error: {e}")

# 7. FOOTER FOR ASU PRESENTATION
st.markdown("---")
st.caption("MIRACLE Project v1.3 | Powered by Random Forest Multi-Omic Integration | ASU 2026")

else:
    st.info("Please upload an RNA-Seq .tsv file in the sidebar to begin analysis.")
