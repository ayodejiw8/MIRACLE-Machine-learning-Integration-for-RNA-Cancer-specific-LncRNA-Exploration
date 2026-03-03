import tarfile
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def analyze_gdc_rnaseq(tar_path):
    # 1. Extract Archive
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall()
        tsv_files = [f for f in tar.getnames() if f.endswith('.tsv')]
    
    if not tsv_files:
        print("No TSV files found in the archive.")
        return

    data_file = tsv_files[0]
    print(f"Analyzing file: {data_file}")

    # 2. Load and Clean Data
    # Skipping the first row (# gene-model) and handling the header
    df = pd.read_csv(data_file, sep='\t', skiprows=1)
    
    # Remove mapping statistics (rows starting with N_)
    df_clean = df[~df['gene_id'].str.startswith('N_')].copy()

    # 3. Gene Type Distribution
    type_counts = df_clean['gene_type'].value_counts()
    print("\nTop 5 Gene Types:\n", type_counts.head())

    # 4. Specific lncRNA Analysis
    lncrna_df = df_clean[df_clean['gene_type'] == 'lncRNA'].copy()
    top_lnc = lncrna_df.sort_values(by='tpm_unstranded', ascending=False).head(10)
    
    # 5. Visualizations
    # Top 10 lncRNAs
    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_lnc, x='tpm_unstranded', y='gene_name', palette='viridis')
    plt.title('Top 10 Expressed lncRNAs')
    plt.xlabel('TPM (Transcripts Per Million)')
    plt.ylabel('Gene Symbol')
    plt.savefig('top_lncrna_expression.png')
    
    # Expression Comparison (Boxplot)
    plt.figure(figsize=(12, 6))
    major_types = type_counts.head(4).index.tolist()
    plot_df = df_clean[df_clean['gene_type'].isin(major_types)].copy()
    plot_df['log_tpm'] = np.log10(plot_df['tpm_unstranded'] + 1)
    
    sns.boxplot(x='gene_type', y='log_tpm', data=plot_df)
    plt.title('Expression Level Distribution by Gene Type')
    plt.ylabel('Log10(TPM + 1)')
    plt.savefig('expression_by_type.png')

    # 6. Save results
    lncrna_df.to_csv('all_lncrna_expression.csv', index=False)
    print("\nResults saved to 'all_lncrna_expression.csv' and PNG images.")

# Execute the analysis
analyze_gdc_rnaseq('gdc_download_20260303_162919.355206.tar.gz')