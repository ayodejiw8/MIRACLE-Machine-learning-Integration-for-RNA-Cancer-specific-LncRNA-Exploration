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
st.set_page_config(page_title="MIRACLE Project | ASU", page_icon="🧬")

st.title("MIRACLE Project: Osteosarcoma lncRNA Detector")
st.markdown("### Albany State University Undergraduate Research 2025-2026 Cohort")
st.markdown("**Developed by:** Ayodeji Williams | Albany State University Bioinformatics")
st.markdown("**Mentors:** Dr. Olabisi Ojo & Dr. Wanjun Hu")
st.markdown("*College of Arts and Science | Dept. of Math, CS and Physics & Dept. of Natural Sciences*")

st.write("---")

# 4. SIDEBAR INSTRUCTIONS
st.sidebar.header("User Control Panel")
st.sidebar.markdown("""
**Instructions:**
1. Upload a **.tsv** STAR-aligned file.
2. The AI will scan for **MALAT1** & **NEAT1**.
3. Results will display in the main window.
""")

uploaded_file = st.sidebar.file_uploader("Upload Genomic Profile (.tsv)", type="tsv")

# 5. MAIN LOGIC
if model is None:
    st.error("Critical Error: Model file 'os_lncrna_detector.pkl' missing.")

elif uploaded_file is not None:
    try:
        # --- ROBUST FILE READING ---
        # We read the file line by line to find the actual header row
        # GDC files typically start the data after several metadata lines
        
        # Reset file pointer
        uploaded_file.seek(0)
        lines = uploaded_file.readlines()
        skip_count = 0
        for i, line in enumerate(lines):
            line_str = line.decode('utf-8')
            if 'gene_id' in line_str or 'gene_name' in line_str:
                skip_count = i
                break
        
        # Re-read the file with the correct skip-row count
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep='\t', skiprows=skip_count)
        df.columns = [c.lower().strip() for c in df.columns]

        # Identify Column Names
        gene_col = 'gene_name' if 'gene_name' in df.columns else 'gene symbol'
        val_col = 'tpm_unstranded' if 'tpm_unstranded' in df.columns else 'copy_number' if 'copy_number' in df.columns else df.columns[1]

        # Extract Values
        m1_data = df[df[gene_col].str.contains('MALAT1', na=False, case=False)]
        n1_data = df[df[gene_col].str.contains('NEAT1', na=False, case=False)]
        
        if not m1_data.empty and not n1_data.empty:
            malat1 = float(m1_data[val_col].values[0])
            neat1 = float(n1_data[val_col].values[0])
            
            st.write(f"### Analysis for {uploaded_file.name}")
            col1, col2 = st.columns(2)
            col1.metric("MALAT1 TPM", f"{malat1:.2f}")
            col2.metric("NEAT1 TPM", f"{neat1:.2f}")

            # Prediction
            sample = pd.DataFrame([[malat1, neat1]], columns=['MALAT1', 'NEAT1'])
            prob = model.predict_proba(sample)[0][1]
            
            st.write("---")
            st.write("## AI Diagnostic Result")
            
            if prob > 0.75:
                st.error(f"🚨 **High Risk Detected**: {prob*100:.1f}% probability")
            elif prob > 0.4:
                st.warning(f"⚠️ **Moderate Risk / Borderline**: {prob*100:.1f}% probability")
            else:
                st.success(f"✅ **Normal/Healthy Signature**: {prob*100:.1f}% probability")
            
            st.progress(prob)
            
        else:
            st.warning("⚠️ MALAT1 or NEAT1 not found in the parsed data. Please check file format.")

    except Exception as e:
        st.error(f"Error processing file: {e}")

else:
    st.info("Please upload a .tsv file in the sidebar to begin analysis.")

# 6. FOOTER
st.markdown("---")
st.caption("MIRACLE v1.3 | ASU 2026")
