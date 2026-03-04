import streamlit as st
import pandas as pd
import joblib
import numpy as np
import os

# 1. FIX THE PATH: Find the model file in the same folder as this script
base_dir = os.path.dirname(__file__)
model_path = os.path.join(base_dir, 'os_lncrna_detector.pkl')

# 2. LOAD THE MODEL ONCE (at the top)
model = joblib.load(model_path)

st.title("MIRACLE Project: Osteosarcoma lncRNA Detector")
st.markdown("Albany State University undergraduate research 2025-2026 cohort")
st.markdown("Developed by Ayodeji Williams | Albany State University Bioinformatics")
st.markdown("Mentors: Dr.Olabisi Ojo & Dr. Wanjun Hu")
st.markdown("College of Arts and Science| Department of Math, CS and Physics & Department of Natural Sciences")

st.sidebar.header("User Input")
uploaded_file = st.sidebar.file_uploader("Upload STAR-Aligned Gene Counts (.tsv)", type="tsv")

if uploaded_file is not None:
    # Read the dataset
    df = pd.read_csv(uploaded_file, sep='\t', skiprows=1)
    
    # Extract the necessary lncRNA values
    try:
        # We use tpm_unstranded as it's the standard for our model
        malat1 = df[df['gene_name'] == 'MALAT1']['tpm_unstranded'].values[0]
        neat1 = df[df['gene_name'] == 'NEAT1']['tpm_unstranded'].values[0]
        
        # Display extracted values
        st.write(f"### Analysis for {uploaded_file.name}")
        col1, col2 = st.columns(2)
        col1.metric("MALAT1 TPM", f"{malat1:.2f}")
        col2.metric("NEAT1 TPM", f"{neat1:.2f}")

        # Perform Prediction
        sample = pd.DataFrame([[malat1, neat1]], columns=['MALAT1', 'NEAT1'])
        prob = model.predict_proba(sample)[0][1]
        
        st.write("---")
        st.write("## Detection Result")
        
        if prob > 0.7:
            st.error(f"High Risk of Osteosarcoma Signature: {prob*100:.1f}% probability")
        elif prob > 0.4:
            st.warning(f"Moderate Risk / Borderline: {prob*100:.1f}% probability")
        else:
            st.success(f"Normal/Healthy Signature: {prob*100:.1f}% probability")
            
    except Exception as e:
        st.error(f"Error extracting data: Ensure columns 'gene_name' and 'tpm_unstranded' exist.")
            
    except IndexError:
        st.error("Error: Could not find MALAT1 or NEAT1 in the dataset. Please ensure the file follows GDC standards.")

else:
    st.info("Please upload an RNA-Seq .tsv file in the sidebar to begin analysis.")
