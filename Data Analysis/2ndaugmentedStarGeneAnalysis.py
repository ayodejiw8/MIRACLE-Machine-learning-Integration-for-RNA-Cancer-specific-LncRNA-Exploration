import tarfile
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def perform_lncrna_analysis(tar_path):
    """
    Extracts GDC RNA-seq data, filters for lncRNAs, 
    and generates statistical summaries and visualizations.
    """
    
    # 1. Extraction
    print(f"--- Extracting {tar_path} ---")
    extract_dir = 'extracted_data_sample2'
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(path=extract_dir)
        # Identify the TSV data file
        tsv_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(extract_dir) 
                     for f in filenames if f.endswith('.tsv')]
    
    if not tsv_files:
        print("Error: No TSV data file found in the archive.")
        return

    data_file = tsv_files[0]
    print(f"Analyzing: {data_file}")

    # 2. Data Loading & Cleaning
    # GDC files typically have a comment on line 1, header on line 2
    df = pd.read_csv(data_file, sep='\t', skiprows=1)
    
    # Remove metadata rows (those starting with 'N_')
    df_clean = df[~df['gene_id'].str.startswith('N_')].copy()

    # 3. lncRNA Filtering
    lncrna_df = df_clean[df_clean['gene_type'] == 'lncRNA'].copy()
    protein_coding_df = df_clean[df_clean['gene_type'] == 'protein_coding'].copy()

    # 4. Statistical Summary
    print("\n--- Summary Statistics ---")
    print(f"Total Genes: {len(df_clean)}")
    print(f"lncRNA Count: {len(lncrna_df)}")
    print(f"Top 5 lncRNAs by TPM:")
    top_5 = lncrna_df.sort_values(by='tpm_unstranded', ascending=False).head(5)
    print(top_5[['gene_name', 'tpm_unstranded']])

    # 5. Visualizations
    sns.set_theme(style="whitegrid")

    # Plot A: Top 10 Expressed lncRNAs
    plt.figure(figsize=(12, 7))
    top_10 = lncrna_df.sort_values(by='tpm_unstranded', ascending=False).head(10)
    sns.barplot(data=top_10, x='tpm_unstranded', y='gene_name', palette='magma')
    plt.title('Top 10 Highly Expressed lncRNAs', fontsize=15)
    plt.xlabel('TPM (Transcripts Per Million)')
    plt.ylabel('Gene Symbol')
    plt.tight_layout()
    plt.savefig('top_lncrna_bar_chart.png')
    print("Saved: top_lncrna_bar_chart.png")

    # Plot B: Expression Distribution Comparison (lncRNA vs Protein Coding)
    plt.figure(figsize=(10, 6))
    compare_df = pd.concat([
        lncrna_df[['gene_type', 'tpm_unstranded']], 
        protein_coding_df[['gene_type', 'tpm_unstranded']]
    ])
    # Use log scale for better visibility of the spread
    compare_df['log_tpm'] = np.log10(compare_df['tpm_unstranded'] + 0.01)
    
    sns.violinplot(x='gene_type', y='log_tpm', data=compare_df, palette='Set2')
    plt.title('Expression Distribution: lncRNA vs Protein Coding', fontsize=14)
    plt.ylabel('log10(TPM + 0.01)')
    plt.savefig('expression_comparison_violin.png')
    print("Saved: expression_comparison_violin.png")

    # 6. Export Results
    lncrna_df.to_csv('lncRNA_expression_results.csv', index=False)
    print("Saved: lncRNA_expression_results.csv")

if __name__ == "__main__":
    # Path to your second dataset
    target_file = 'gdc_download_20260303_170024.760832.tar.gz'
    perform_lncrna_analysis(target_file)
