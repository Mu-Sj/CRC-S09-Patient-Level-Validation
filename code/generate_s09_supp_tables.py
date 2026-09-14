from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'results/final/s09_figures'
OUT = ROOT / 'results/final/s09_tables'
OUT.mkdir(parents=True, exist_ok=True)

# Table 1: frozen denominators and GSE294300 coverage.
den = pd.read_csv(FIG/'SupplementaryFigure1_denominators.tsv', sep='\t')
qc = pd.read_csv(FIG/'SupplementaryFigure1_GSE294300_cell_QC.tsv', sep='\t')
den.to_csv(OUT/'SupplementaryTable1_denominators.tsv', sep='\t', index=False)
qc.to_csv(OUT/'SupplementaryTable1_GSE294300_cell_QC.tsv', sep='\t', index=False)

# Tables 2-6 are frozen source-data exports used by Supplementary Figures 2-6.
mapping = {
    2: 'SupplementaryFigure2_input_inventory.tsv',
    3: 'SupplementaryFigure3_malignancy_CNV_evidence.tsv',
    4: 'SupplementaryFigure4_root_sensitivity_and_failures.tsv',
    5: 'SupplementaryFigure5_GSE178318_COL17LM_sensitivity.tsv',
    6: 'SupplementaryFigure6_progression_patient_deltas.tsv',
}
for num, src in mapping.items():
    pd.read_csv(FIG/src, sep='\t').to_csv(OUT/f'SupplementaryTable{num}.tsv', sep='\t', index=False)

manifest = pd.DataFrame([
    ['Supplementary Table 1a','SupplementaryTable1_denominators.tsv','Supplementary Figure 1','Frozen independent-unit denominators by analysis role'],
    ['Supplementary Table 1b','SupplementaryTable1_GSE294300_cell_QC.tsv','Supplementary Figure 1','GSE294300 per-patient cell and scoreability coverage'],
    ['Supplementary Table 2','SupplementaryTable2.tsv','Supplementary Figure 2','Candidate input matrix and structural QC inventory'],
    ['Supplementary Table 3','SupplementaryTable3.tsv','Supplementary Figure 3','Malignancy/CNV evidence tiers and cohort restrictions'],
    ['Supplementary Table 4','SupplementaryTable4.tsv','Supplementary Figure 4','Root sensitivity, bootstrap execution and failure records'],
    ['Supplementary Table 5','SupplementaryTable5.tsv','Supplementary Figure 5','GSE178318 COL17-LM low-count sensitivity layer'],
    ['Supplementary Table 6','SupplementaryTable6.tsv','Supplementary Figure 6','Patient/donor progression deltas and cohort identifiers'],
], columns=['table','file','linked_figure','description'])
manifest.to_csv(OUT/'S09_SUPPLEMENTARY_TABLE_MANIFEST.tsv', sep='\t', index=False)

with open(OUT/'README.md', 'w', encoding='utf-8') as f:
    f.write('S09 Supplementary Tables\n\n')
    f.write('Generated from frozen S09 figure source tables in WSL Ubuntu-24.04 using cast-py. No new analysis or parameter tuning was performed. Patients/donors are inferential units; cells are coverage/QC observations. Progression cohorts remain separate descriptive analyses.\n\n')
    f.write('See S09_SUPPLEMENTARY_TABLE_MANIFEST.tsv for the complete table-to-figure mapping.\n')
print(f'Generated Supplementary Tables 1-6 in {OUT}')
