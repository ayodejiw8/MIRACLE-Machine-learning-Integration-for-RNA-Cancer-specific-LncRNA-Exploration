import joblib
import pandas as pd
import numpy as np

# Load the model you built with Gemini
model = joblib.load('os_lncrna_detector.pkl')

def detect_os(malat1_tpm, neat1_tpm):
    # Create a small dataframe for the model to read
    sample = pd.DataFrame([[malat1_tpm, neat1_tpm]], columns=['MALAT1', 'NEAT1'])
    
    # Run the prediction
    prediction = model.predict(sample)
    probability = model.predict_proba(sample)[0][1]
    
    print(f"--- Analysis Results ---")
    if prediction[0] == 1:
        print(f"Status: OSTEOSARCOMA SIGNATURE DETECTED")
        print(f"Probability: {probability*100:.2f}%")
    else:
        print(f"Status: NORMAL / BASELINE")
        print(f"Probability: {probability*100:.2f}%")

# TEST IT: Enter values from your recent STAR-aligned counts
# Example: Using Sample 3's values
detect_os(12684.09, 596.39)
