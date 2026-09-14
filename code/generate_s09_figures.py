from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/final/s09_figures'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,
                     'axes.titlesize':10,'xtick.labelsize':7,'ytick.labelsize':7,
                     'figure.dpi':150,'savefig.dpi':600})
COLORS = {'normal':'#0072B2','primary':'#D55E00','metastasis':'#009E73',
          'GSE315534':'#E69F00','GSE178318':'#56B4E9','Moorman_2024':'#CC79A7'}

def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png', dpi=600, bbox_inches='tight')
    plt.close(fig)

# Figure 1: cohort design and inference units.
design = pd.DataFrame([
    ['GSE294300','normal-primary','16','patient','primary S09 endpoint'],
    ['GSE294300','normal-primary','2','patient','sensitivity layer'],
    ['GSE315534','primary-liver metastasis','3','patient','descriptive progression'],
    ['GSE178318','primary-liver metastasis','6','patient','descriptive progression'],
    ['Moorman 2024','primary-metastasis','6','donor','descriptive progression'],
], columns=['cohort','edge','n','inference_unit','role'])
design.to_csv(OUT/'figure1_study_design.tsv', sep='\t', index=False)
fig, ax = plt.subplots(figsize=(6.2,3.2))
y=np.arange(len(design)); vals=design.n.astype(int)
ax.barh(y, vals, color=[COLORS.get(c,'#666666') for c in design.cohort], edgecolor='none')
ax.set_yticks(y, [f"{r.cohort} | {r.edge}" for r in design.itertuples()])
ax.set_xlabel('Independent patients or donors (n)')
ax.set_title('S09 external-validation design')
for yi,v,r in zip(y,vals,design.itertuples()): ax.text(v+0.15, yi, f'n={v}; {r.inference_unit}', va='center', fontsize=7)
ax.set_xlim(0, max(vals)+4); ax.spines[['top','right']].set_visible(False)
save(fig,'Figure1_S09_study_design')

# Figure 2: paired GSE294300 endpoint.
pair = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_paired_state_concordance.tsv', sep='\t')
pair = pair[pair.cohort.eq('primary')].copy()
long=[]
for _,r in pair.iterrows():
    long += [[r.patient,'normal','ABSORPTIVE_BEST4_CA4',r.absorptive_normal], [r.patient,'primary','ABSORPTIVE_BEST4_CA4',r.absorptive_tumor],
             [r.patient,'normal','SECRETORY_MUC2_TFF3',r.secretory_normal], [r.patient,'primary','SECRETORY_MUC2_TFF3',r.secretory_tumor]]
pd.DataFrame(long,columns=['patient','arm','module','score']).to_csv(OUT/'figure2_gse294300_patient_scores.tsv',sep='\t',index=False)
fig,axs=plt.subplots(1,2,figsize=(6.2,3.0),sharey=True)
for ax,module in zip(axs,['ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3']):
    sub=pair
    a='absorptive' if 'ABSORPTIVE' in module else 'secretory'
    for _,r in sub.iterrows(): ax.plot([0,1],[r[f'{a}_normal'],r[f'{a}_tumor']],color='#999999',lw=.7,alpha=.7)
    ax.scatter(np.zeros(len(sub)),sub[f'{a}_normal'],color=COLORS['normal'],s=14,label='Normal')
    ax.scatter(np.ones(len(sub)),sub[f'{a}_tumor'],color=COLORS['primary'],s=14,label='Primary')
    ax.set_xticks([0,1],['Normal','Primary']); ax.set_title(module.replace('_',' ')); ax.set_xlabel('Arm')
    ax.spines[['top','right']].set_visible(False)
axs[0].set_ylabel('Module score (CP10K, log1p)'); axs[1].legend(frameon=False,fontsize=7,loc='best')
fig.suptitle('GSE294300 patient-level paired state endpoint',y=1.03)
save(fig,'Figure2_GSE294300_paired_endpoint')

# Figure 3: trajectory readiness.
ready=pd.read_csv(ROOT/'analysis/09_external_validation/S09_GSE294300_METHOD_READINESS.tsv',sep='\t')
ready.to_csv(OUT/'figure3_trajectory_readiness.tsv',sep='\t',index=False)
fig,ax=plt.subplots(figsize=(5.2,3.0)); x=np.arange(len(ready)); p=ready['pass']; f=ready['hard_failure']
ax.bar(x,p,color='#009E73',label='PASS'); ax.bar(x,f,bottom=p,color='#D55E00',label='HARD_FAILURE')
ax.set_xticks(x,[s.replace('_','\n') for s in ready.method_layer]); ax.set_ylabel('Replicates'); ax.set_title('Trajectory execution audit (11 replicates)')
ax.legend(frameon=False,fontsize=7); ax.spines[['top','right']].set_visible(False)
for i,r in ready.iterrows(): ax.text(i,(r['pass']+r['hard_failure'])+.2,f"{int(r['pass'])}/{int(r['replicates'])}",ha='center',fontsize=7)
save(fig,'Figure3_trajectory_execution_audit')

# Figure 4: separate progression cohorts with frozen intervals.
prog=pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260913_S09_direction_audit_v1/cross_dataset_primary_to_metastasis_summary.tsv',sep='\t')
prog.to_csv(OUT/'figure4_cross_cohort_progression.tsv',sep='\t',index=False)
fig,ax=plt.subplots(figsize=(6.2,3.5)); mods=['ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3']; ypos=[]; labels=[]
for j,m in enumerate(mods):
    sub=prog[prog.module.eq(m)].copy(); sub['y']=np.arange(len(sub))*1.0 + j*4
    ypos.extend(sub.y); labels.extend([f"{r.dataset}\n{r.n}" for r in sub.itertuples()])
    for _,r in sub.iterrows():
        ax.plot([r.bootstrap_ci95_low,r.bootstrap_ci95_high],[r.y,r.y],color=COLORS.get(r.dataset,'#444'),lw=2)
        ax.scatter(r.mean_delta,r.y,color=COLORS.get(r.dataset,'#444'),s=24,zorder=3)
ax.axvline(0,color='black',lw=.8); ax.set_yticks(ypos,labels); ax.set_xlabel('Primary-to-metastasis mean delta (95% bootstrap CI)'); ax.set_title('Progression estimates shown separately by cohort')
ax.text(0.02,0.97,'Absorptive',transform=ax.transAxes,va='top',fontsize=8,fontweight='bold'); ax.text(0.02,0.48,'Secretory',transform=ax.transAxes,va='top',fontsize=8,fontweight='bold')
ax.spines[['top','right']].set_visible(False); ax.grid(axis='x',alpha=.2)
save(fig,'Figure4_cross_cohort_progression')

manifest=pd.DataFrame([
 ['Figure1_S09_study_design','figure1_study_design.tsv'],['Figure2_GSE294300_paired_endpoint','figure2_gse294300_patient_scores.tsv'],
 ['Figure3_trajectory_execution_audit','figure3_trajectory_readiness.tsv'],['Figure4_cross_cohort_progression','figure4_cross_cohort_progression.tsv']],columns=['figure','data_table'])
manifest.to_csv(OUT/'S09_FIGURE_DATA_MANIFEST.tsv',sep='\t',index=False)

# Submission-level Figure 1: design, denominators, and evidence boundary.
cohort = pd.DataFrame([
    ['GSE294300','normal-primary',16,'patient','primary endpoint'],
    ['GSE294300','normal-primary',2,'patient','sensitivity only'],
    ['GSE315534','primary-liver metastasis',3,'patient','descriptive progression'],
    ['GSE178318','primary-liver metastasis',6,'patient','descriptive progression'],
    ['Moorman 2024','primary-metastasis',6,'donor','descriptive progression'],
], columns=['cohort','edge','n','unit','role'])
cohort.to_csv(OUT/'Figure1_v2_cohort_denominators.tsv',sep='\t',index=False)
evidence = pd.DataFrame([
    ['Frozen state program','S05 modules','reproducible expression state','primary/secondary'],
    ['Patient-level endpoint','GSE294300 paired scores','16 primary patients; CI reported','primary'],
    ['Trajectory execution','Monocle3/Slingshot/Python','technical stability and failures','secondary'],
    ['Malignancy evidence','CNV/marker/public labels','restricted or supportive proxy','limitation'],
    ['Progression direction','three separate cohorts','non-replicated; descriptive only','limitation'],
], columns=['evidence_layer','source','supported_interpretation','manuscript_role'])
evidence.to_csv(OUT/'Figure1_v2_evidence_boundary.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(7.2,5.2)); gs=fig.add_gridspec(2,2,hspace=.55,wspace=.38)
ax=fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('Study flow',loc='left',fontweight='bold')
boxes=[('Frozen S05\nstate programs',(0.05,.68)),('S09 contract\n& role freeze',(0.39,.68)),('Patient-level\nendpoints',(0.73,.68)),('Separate progression\ndescriptive analyses',(0.39,.18))]
for txt,(x0,y0) in boxes:
    ax.text(x0,y0,txt,ha='center',va='center',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.5',fc='#F3F4F6',ec='#444',lw=.8),fontsize=8)
for a,b in [((.18,.68),(.37,.68)),((.52,.68),(.71,.68)),((.52,.63),(.52,.28))]: ax.annotate('',xy=b,xytext=a,xycoords=ax.transAxes,arrowprops=dict(arrowstyle='->',lw=.9,color='#555'))
ax.text(.5,.02,'No retuning after validation unlock',ha='center',transform=ax.transAxes,fontsize=7,color='#555')
ax=fig.add_subplot(gs[0,1]); y=np.arange(len(cohort)); ax.barh(y,cohort.n,color=[COLORS.get(c,'#666') for c in cohort.cohort]); ax.set_yticks(y,[f'{r.cohort}\n{r.edge}' for r in cohort.itertuples()],fontsize=6); ax.set_xlabel('Patients or donors (n)'); ax.set_title('Frozen S09 denominators',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
for yi,r in zip(y,cohort.itertuples()): ax.text(r.n+.15,yi,f'{r.n}; {r.unit}',va='center',fontsize=7)
ax.set_xlim(0,8)
ax=fig.add_subplot(gs[1,:]); ax.axis('off'); ax.set_title('Evidence boundary',loc='left',fontweight='bold'); tbl=ax.table(cellText=evidence.values,colLabels=['Evidence layer','Source','Supported interpretation','Manuscript role'],loc='center',cellLoc='left',colWidths=[.18,.22,.42,.18]); tbl.auto_set_font_size(False); tbl.set_fontsize(6); tbl.scale(1,1.55)
for (i,j),cell in tbl.get_celld().items(): cell.set_edgecolor('#B0B0B0'); cell.set_linewidth(.4); cell.set_facecolor('#F3F4F6' if i==0 else 'white')
fig.suptitle('Figure 1. S09 study design and evidence boundaries',y=.99,fontsize=12)
save(fig,'Figure1_S09_study_design_v2')

# Submission-level Figure 3: endpoint, bootstrap, coverage, and sensitivity.
primary=pair.copy(); primary['delta_absorptive']=primary.absorptive_tumor-primary.absorptive_normal; primary['delta_secretory']=primary.secretory_tumor-primary.secretory_normal
primary[['patient','cells_normal','cells_tumor','scoreable_normal','scoreable_tumor','delta_absorptive','delta_secretory']].to_csv(OUT/'Figure3_v2_patient_endpoint.tsv',sep='\t',index=False)
boot=pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_bootstrap_ci.tsv',sep='\t'); boot.to_csv(OUT/'Figure3_v2_bootstrap.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(8.0,6.0)); gs=fig.add_gridspec(2,3,hspace=.62,wspace=.42)
ax=fig.add_subplot(gs[0,0]);
for _,r in primary.iterrows(): ax.plot([0,1],[r.absorptive_normal,r.absorptive_tumor],color='#BDBDBD',lw=.7); ax.scatter([0,1],[r.absorptive_normal,r.absorptive_tumor],c=[COLORS['normal'],COLORS['primary']],s=12)
ax.set_xticks([0,1],['Normal','Primary']); ax.set_ylabel('Absorptive score'); ax.set_title('A  Paired scores',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[0,1]); s=primary.sort_values('delta_absorptive'); ax.barh(np.arange(len(s)),s.delta_absorptive,color=np.where(s.delta_absorptive<0,COLORS['primary'],COLORS['normal'])); ax.axvline(0,color='black',lw=.7); ax.set_yticks(np.arange(len(s)),s.patient,fontsize=6); ax.set_xlabel('Tumor - normal'); ax.set_title('B  Patient delta',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[0,2]); ax.hist(boot.mean_delta_absorptive,bins=25,color=COLORS['primary'],alpha=.8,label='Absorptive'); ax.hist(boot.mean_delta_secretory,bins=25,color=COLORS['normal'],alpha=.55,label='Secretory'); ax.axvline(0,color='black',lw=.7); ax.set_xlabel('Bootstrap mean delta'); ax.set_ylabel('Count'); ax.set_title('C  Complete-patient bootstrap',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,0]); ax.scatter(primary.cells_tumor,primary.delta_absorptive,c=COLORS['primary'],s=18,label='Primary cells'); ax.scatter(primary.cells_normal,primary.delta_absorptive,c=COLORS['normal'],s=18,label='Normal cells'); ax.axhline(0,color='black',lw=.6); ax.set_xscale('log'); ax.set_xlabel('Cells per arm (log scale)'); ax.set_ylabel('Absorptive delta'); ax.set_title('D  Coverage audit',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6); ax.spines[['top','right']].set_visible(False); ax.tick_params(axis='x',labelsize=6)
ax=fig.add_subplot(gs[1,1:]); sens=pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_paired_state_concordance.tsv',sep='\t'); sens['layer']=np.where(sens.cohort.eq('primary'),'Primary (n=16)','Sensitivity (n=2)'); sm=sens.groupby('layer')[['delta_absorptive','delta_secretory']].mean().reindex(['Primary (n=16)','Sensitivity (n=2)']); xx=np.arange(len(sm)); w=.34; ax.bar(xx-w/2,sm.delta_absorptive,w,label='Absorptive',color=COLORS['primary']); ax.bar(xx+w/2,sm.delta_secretory,w,label='Secretory',color=COLORS['normal']); ax.axhline(0,color='black',lw=.7); ax.set_xticks(xx,sm.index); ax.set_ylabel('Mean paired delta'); ax.set_title('E  Prespecified sensitivity layer',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=7); ax.spines[['top','right']].set_visible(False)
fig.suptitle('Figure 3. GSE294300 patient-level independent validation',y=.995,fontsize=12)
save(fig,'Figure3_GSE294300_validation_v2')

with open(OUT/'Figure1_Figure3_v2_captions.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 1. S09 study design and evidence boundaries. (A) Frozen state programs, contract, patient-level endpoints, and separate descriptive progression analyses. (B) Frozen cohort denominators. (C) Evidence layers and their permitted interpretation.\n\nFigure 3. GSE294300 patient-level independent validation. (A) Paired module scores for 16 primary patients. (B) Patient-level absorptive deltas sorted by magnitude. (C) Complete-patient bootstrap distributions. (D) Cell-count coverage audit; cells are not independent replicates. (E) Prespecified sensitivity patients 1 and 17 shown separately. Bars and points are descriptive; no pooled-cell inference is used.\n')
# Submission-level Figure 4: technical trajectory stability and failures.
rtraj=pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_r_trajectory_formal_v9/r_trajectory.tsv',sep='\t')
rtraj.to_csv(OUT/'Figure4_v2_r_trajectory_replicates.tsv',sep='\t',index=False)
method=ready[['method_layer','replicates','pass','abstention','invalid_result','hard_failure','status']].copy(); method.to_csv(OUT/'Figure4_v2_method_readiness.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(8.0,5.6)); gs=fig.add_gridspec(2,2,hspace=.62,wspace=.38)
ax=fig.add_subplot(gs[0,0]); x=np.arange(len(method)); ax.bar(x,method['pass'],color='#009E73',label='PASS'); ax.bar(x,method['hard_failure'],bottom=method['pass'],color='#D55E00',label='HARD_FAILURE'); ax.set_xticks(x,[m.replace('_','\n') for m in method.method_layer],fontsize=6); ax.set_ylabel('Replicates (n)'); ax.set_title('A  Method readiness',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[0,1]); rtraj['slingshot_fraction']=rtraj.slingshot_finite_all/rtraj.cells; ax.plot(rtraj.replicate,rtraj.slingshot_fraction,'o-',color=COLORS['primary'],label='Slingshot all-cell finite'); ax.plot(rtraj.replicate,np.ones(len(rtraj)),':',color=COLORS['normal'],label='Monocle3 finite'); ax.set_ylim(0,1.08); ax.set_xlabel('Replicate'); ax.set_ylabel('Finite coverage fraction'); ax.set_title('B  Finite pseudotime coverage',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,0]); ax.scatter(rtraj.monocle3_graph_nodes,rtraj.monocle3_graph_leaves,c=rtraj.replicate,cmap='viridis',s=28); ax.set_xlabel('Monocle3 graph nodes'); ax.set_ylabel('Graph leaves'); ax.set_title('C  Graph structure across replicates',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,1]); mat=method[['pass','abstention','invalid_result','hard_failure']].to_numpy(); im=ax.imshow(mat,cmap='PuOr',aspect='auto',vmin=0,vmax=max(1,mat.max())); ax.set_xticks(range(4),['PASS','ABSTAIN','INVALID','HARD\nFAIL'],fontsize=6); ax.set_yticks(range(len(method)),[m.replace('_','\n') for m in method.method_layer],fontsize=6); ax.set_title('D  Failure-mode matrix',loc='left',fontweight='bold');
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]): ax.text(j,i,str(int(mat[i,j])),ha='center',va='center',fontsize=8)
fig.suptitle('Figure 4. Trajectory execution stability and failure modes',y=.995,fontsize=12); save(fig,'Figure4_trajectory_stability_v2')

# Submission-level Figure 5: separate progression cohorts.
prog.to_csv(OUT/'Figure5_v2_progression_estimates.tsv',sep='\t',index=False)
meta=pd.DataFrame([
 ['GSE315534',3,'patient','primary-liver metastasis','supportive CNV proxy'],
 ['GSE178318',6,'patient','primary-liver metastasis','restricted CNV evidence'],
 ['Moorman 2024',6,'donor','primary-metastasis','public malignant-cell labels'],
],columns=['dataset','n','inference_unit','edge','malignancy_evidence'])
meta.to_csv(OUT/'Figure5_v2_cohort_context.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(8.0,5.4)); gs=fig.add_gridspec(2,2,hspace=.6,wspace=.4)
for k,(m,title) in enumerate([('ABSORPTIVE_BEST4_CA4','A  Absorptive program'),('SECRETORY_MUC2_TFF3','B  Secretory program')]):
    ax=fig.add_subplot(gs[0,k]); sub=prog[prog.module.eq(m)].copy(); sub['y']=np.arange(len(sub));
    for _,r in sub.iterrows(): ax.plot([r.bootstrap_ci95_low,r.bootstrap_ci95_high],[r.y,r.y],lw=2.4,color=COLORS.get(r.dataset,'#555')); ax.scatter(r.mean_delta,r.y,s=26,color=COLORS.get(r.dataset,'#555'))
    ax.axvline(0,color='black',lw=.7); ax.set_yticks(sub.y,[f'{r.dataset} (n={r.n})' for r in sub.itertuples()],fontsize=6); ax.set_xlabel('Mean delta (95% bootstrap CI)'); ax.set_title(title,loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,0]); sign=prog.pivot(index='dataset',columns='module',values='mean_delta').reindex(['GSE315534','GSE178318','Moorman_2024']); im=ax.imshow(sign.values,cmap='PuOr',vmin=-.2,vmax=.2,aspect='auto'); ax.set_xticks(range(2),['Absorptive','Secretory'],fontsize=7); ax.set_yticks(range(3),sign.index,fontsize=7); ax.set_title('C  Direction sign by cohort',loc='left',fontweight='bold');
for i in range(sign.shape[0]):
    for j in range(sign.shape[1]): ax.text(j,i,f'{sign.iloc[i,j]:+.3f}',ha='center',va='center',fontsize=7)
ax=fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  Cohort context',loc='left',fontweight='bold'); short_meta=meta.copy(); short_meta['malignancy_evidence']=short_meta['malignancy_evidence'].replace({'supportive CNV proxy':'CNV proxy','restricted CNV evidence':'restricted CNV','public malignant-cell labels':'public labels'}); tbl=ax.table(cellText=short_meta.values,colLabels=['Dataset','n','Unit','Edge','Malignancy evidence'],loc='center',cellLoc='left',colWidths=[.23,.08,.14,.31,.24]); tbl.auto_set_font_size(False); tbl.set_fontsize(5.5); tbl.scale(1,1.45)
for (i,j),cell in tbl.get_celld().items(): cell.set_edgecolor('#B0B0B0'); cell.set_linewidth(.4); cell.set_facecolor('#F3F4F6' if i==0 else 'white')
fig.suptitle('Figure 5. Cross-cohort progression estimates without pooling',y=.995,fontsize=12); save(fig,'Figure5_cross_cohort_progression_v2')

with open(OUT/'Figure4_Figure5_v2_captions.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 4. Trajectory execution stability and failure modes. (A) Method readiness across 11 replicates. (B) Finite pseudotime coverage across formal replicates. (C) Monocle3 graph structure. (D) Failure-mode counts. These are technical execution results, not lineage evidence.\n\nFigure 5. Cross-cohort progression estimates without pooling. (A-B) Frozen patient/donor-level mean deltas and percentile bootstrap intervals shown separately by module and cohort. (C) Direction signs are heterogeneous. (D) Cohort context and malignancy evidence. No pooled estimate or direction vote is shown.\n')

# Submission-level Figure 1 v3: visual-first study design.
role_colors={'primary endpoint':'#0072B2','sensitivity only':'#56B4E9','descriptive progression':'#009E73'}
fig=plt.figure(figsize=(8.2,5.8)); gs=fig.add_gridspec(2,2,height_ratios=[1.0,1.18],hspace=.58,wspace=.34)
# Panel A: flow with concise labels.
ax=fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Frozen analysis path',loc='left',fontweight='bold',pad=8)
flow=[('S05','State\nprograms','#F0F4F8'),('S09','Contract\n& role freeze','#E6F2F8'),('Endpoint','Patient-level\nstate scores','#E8F5EF'),('Output','Descriptive\nexternal validation','#FFF3D6')]
for i,(tag,txt,fc) in enumerate(flow):
    x=.08+i*.29; ax.text(x,.58,tag,ha='center',va='center',fontsize=7,color='#555',transform=ax.transAxes)
    ax.text(x,.36,txt,ha='center',va='center',fontsize=8,fontweight='bold',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.55',fc=fc,ec='#667085',lw=.8))
    if i<3: ax.annotate('',xy=(x+.22,.36),xytext=(x+.12,.36),xycoords=ax.transAxes,arrowprops=dict(arrowstyle='-|>',lw=1.0,color='#667085'))
ax.text(.5,.08,'No retuning after validation unlock',ha='center',transform=ax.transAxes,fontsize=7,color='#667085')
# Panel B: cohort lanes, primary/sensitivity/progression.
ax=fig.add_subplot(gs[0,1]); ax.set_title('B  Cohort roles and denominators',loc='left',fontweight='bold',pad=8)
rows=[('GSE294300','Primary validation',16,'patient','primary endpoint'),('GSE294300','Sensitivity layer',2,'patient','sensitivity only'),('GSE315534','Progression',3,'patient','descriptive progression'),('GSE178318','Progression',6,'patient','descriptive progression'),('Moorman 2024','Progression',6,'donor','descriptive progression')]
yy=np.arange(len(rows)); ax.set_xlim(-1,8.7); ax.set_ylim(-.7,len(rows)-.3); ax.set_yticks(yy,[f'{r[0]}  |  {r[1]}' for r in rows],fontsize=6); ax.invert_yaxis(); ax.set_xlabel('Independent patients or donors (n)',fontsize=8); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
for i,(ds,role,n,unit,status) in enumerate(rows):
    ax.scatter(n, i, s=95, color=role_colors[status], zorder=3, edgecolor='white', linewidth=.8)
    ax.text(n+.28,i,f'{n}  {unit}',va='center',fontsize=7,color='#344054')
    ax.hlines(i,0,n,color=role_colors[status],lw=5,alpha=.22)
for label,x,c in [('primary',.02,'#0072B2'),('sensitivity',.02,'#56B4E9'),('progression',.02,'#009E73')]: pass
ax.text(.02,1.02,'Color encodes frozen role',transform=ax.transAxes,fontsize=6.5,color='#667085')
# Panel C: evidence boundary as aligned status cards.
ax=fig.add_subplot(gs[1,:]); ax.axis('off'); ax.set_title('C  Evidence boundary for interpretation',loc='left',fontweight='bold',pad=8)
cards=[('State programs','CAN ANSWER','reproducible expression states','#E8F5EF','#009E73'),('Patient endpoint','CAN ANSWER','paired state differences with bootstrap CI','#E6F2F8','#0072B2'),('Trajectory tools','TECHNICAL ONLY','execution stability, failure and abstention','#FFF3D6','#E69F00'),('Progression direction','DESCRIPTIVE ONLY','separate cohort estimates; no pooled direction','#FDECEC','#D55E00')]
for i,(head,status,desc,fc,edge) in enumerate(cards):
    x=.02+i*.245; ax.add_patch(plt.Rectangle((x,.18),.215,.62,transform=ax.transAxes,fc=fc,ec=edge,lw=1.0))
    ax.text(x+.015,.69,head,transform=ax.transAxes,fontsize=8,fontweight='bold',color='#1D2939')
    ax.text(x+.015,.52,status,transform=ax.transAxes,fontsize=7,fontweight='bold',color=edge)
    ax.text(x+.015,.34,desc,transform=ax.transAxes,fontsize=7,color='#344054',wrap=True,va='top')
ax.text(.02,.04,'Patients/donors are the inferential units; cells are not independent replicates. CNV/marker evidence is auxiliary or restricted.',transform=ax.transAxes,fontsize=7,color='#667085')
fig.suptitle('Figure 1. S09 study design, cohort roles and evidence boundaries',y=.995,fontsize=12)
save(fig,'Figure1_S09_study_design_v3')
pd.DataFrame(rows,columns=['dataset','role_label','n','inference_unit','frozen_role']).to_csv(OUT/'Figure1_v3_cohort_roles.tsv',sep='\t',index=False)
pd.DataFrame([(c[0],c[1],c[2]) for c in cards],columns=['evidence_layer','interpretation_class','supported_claim']).to_csv(OUT/'Figure1_v3_evidence_cards.tsv',sep='\t',index=False)
with open(OUT/'Figure1_v3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 1. S09 study design, cohort roles and evidence boundaries. (A) Frozen path from state-program definition through the S09 contract to patient-level descriptive validation. (B) Cohort roles and independent denominators; color encodes frozen role. (C) Evidence layers and the strongest interpretation permitted for each. Patients or donors are inferential units, and cells are not independent replicates. CNV and marker evidence remain auxiliary or restricted.\n')

# Figure 1 v4: single-column-safe layout with no text-heavy table.
fig=plt.figure(figsize=(8.2,7.1)); gs=fig.add_gridspec(3,1,height_ratios=[1.0,1.55,1.12],hspace=.62)
ax=fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Frozen analysis path',loc='left',fontweight='bold',pad=6)
flow=[('S05','State programs','#F0F4F8'),('S09','Contract +\nrole freeze','#E6F2F8'),('Endpoint','Patient-level\nstate scores','#E8F5EF'),('Output','Descriptive\nvalidation','#FFF3D6')]
for i,(tag,txt,fc) in enumerate(flow):
    x=.11+i*.26; ax.text(x,.73,tag,ha='center',va='center',fontsize=7,color='#667085',transform=ax.transAxes)
    ax.text(x,.43,txt,ha='center',va='center',fontsize=9,fontweight='bold',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.55',fc=fc,ec='#667085',lw=.9),linespacing=1.15)
    if i<3: ax.annotate('',xy=(x+.205,.43),xytext=(x+.145,.43),xycoords=ax.transAxes,arrowprops=dict(arrowstyle='-|>',lw=1.1,color='#667085'))
ax.text(.5,.08,'Frozen inputs → patient-level endpoint → descriptive handoff',ha='center',transform=ax.transAxes,fontsize=7,color='#667085')

ax=fig.add_subplot(gs[1,0]); ax.set_title('B  Cohort roles and independent denominators',loc='left',fontweight='bold',pad=6)
rows2=[('GSE294300','Primary validation',16,'patient','primary endpoint'),('GSE294300','Sensitivity layer',2,'patient','sensitivity only'),('GSE315534','Primary → liver metastasis',3,'patient','descriptive progression'),('GSE178318','Primary → liver metastasis',6,'patient','descriptive progression'),('Moorman 2024','Primary → metastasis',6,'donor','descriptive progression')]
yy=np.arange(len(rows2)); ax.set_xlim(0,18); ax.set_ylim(-.7,len(rows2)-.3); ax.invert_yaxis(); ax.set_yticks(yy,[f'{r[0]}  |  {r[1]}' for r in rows2],fontsize=7); ax.set_xlabel('Independent patients or donors (n)'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
for i,(ds,role,n,unit,status) in enumerate(rows2):
    c=role_colors[status]; unit_text = unit if n == 1 else unit + 's'; ax.hlines(i,0,n,color=c,lw=7,alpha=.20); ax.scatter(n,i,s=120,color=c,edgecolor='white',linewidth=.9,zorder=3); ax.text(n+.35,i,f'{n}  {unit_text}',va='center',fontsize=8,color='#344054')

ax=fig.add_subplot(gs[2,0]); ax.axis('off'); ax.set_title('C  Evidence boundary for interpretation',loc='left',fontweight='bold',pad=6)
cards2=[('State programs','SUPPORTED','Expression\nstates','#E8F5EF','#009E73'),('Patient endpoint','SUPPORTED','Paired delta\nwith 95% CI','#E6F2F8','#0072B2'),('Trajectory tools','TECHNICAL','Execution, failure\nand abstention','#FFF3D6','#E69F00'),('Progression direction','DESCRIPTIVE','Separate cohorts;\nno pooled direction','#FDECEC','#D55E00')]
for i,(head,status,desc,fc,edge) in enumerate(cards2):
    x=.015+i*.25; ax.add_patch(plt.Rectangle((x,.15),.22,.72,transform=ax.transAxes,fc=fc,ec=edge,lw=1.0)); ax.text(x+.014,.72,head,transform=ax.transAxes,fontsize=8,fontweight='bold',color='#1D2939'); ax.text(x+.014,.52,status,transform=ax.transAxes,fontsize=7,fontweight='bold',color=edge); ax.text(x+.014,.36,desc,transform=ax.transAxes,fontsize=7.5,color='#344054',va='top',linespacing=1.2)
ax.text(.015,.03,'Patients/donors are inferential units. CNV/marker evidence is auxiliary or restricted.',transform=ax.transAxes,fontsize=7,color='#667085')
fig.suptitle('Figure 1. S09 study design, cohort roles and evidence boundaries',y=.995,fontsize=12)
save(fig,'Figure1_S09_study_design_v4')
pd.DataFrame(rows2,columns=['dataset','role_label','n','inference_unit','frozen_role']).to_csv(OUT/'Figure1_v4_cohort_roles.tsv',sep='\t',index=False)
pd.DataFrame([(c[0],c[1],c[2]) for c in cards2],columns=['evidence_layer','interpretation_class','supported_claim']).to_csv(OUT/'Figure1_v4_evidence_cards.tsv',sep='\t',index=False)
with open(OUT/'Figure1_v4_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 1. S09 study design, cohort roles and evidence boundaries. (A) Frozen analysis path from state programs through the S09 contract to patient-level descriptive validation. (B) Independent patient/donor denominators by frozen role. (C) The strongest interpretation permitted for each evidence layer. Patients/donors are inferential units; cells are not independent replicates. CNV/marker evidence is auxiliary or restricted.\n')

# Figure 4 v3: visual-first technical stability figure.
fig=plt.figure(figsize=(8.2,6.2)); gs=fig.add_gridspec(2,2,hspace=.72,wspace=.42)
method_v=ready[['method_layer','replicates','pass','abstention','invalid_result','hard_failure']].copy()
ax=fig.add_subplot(gs[0,0]); labels_v=['R Slingshot','R Monocle3','Python PAGA +\nPalantir']; y=np.arange(len(method_v)); total=method_v.replicates.to_numpy(); pa=method_v['pass'].to_numpy(); hf=method_v.hard_failure.to_numpy(); ax.barh(y,pa,color='#009E73',height=.46,label='PASS'); ax.barh(y,hf,left=pa,color='#D55E00',height=.46,label='HARD FAILURE'); ax.set_yticks(y,labels_v,fontsize=7); ax.invert_yaxis(); ax.set_xlim(0,11); ax.set_xlabel('Formal replicates (n)'); ax.set_title('A  Execution readiness',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,handlelength=1.4,columnspacing=1.0); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
for i,(p,h,t) in enumerate(zip(pa,hf,total)): ax.text(min(t+.2,10.7),i,f'{p}/{t} pass',va='center',fontsize=7,color='#344054')
ax=fig.add_subplot(gs[0,1]); rtraj2=rtraj.copy(); rtraj2['slingshot_fraction']=rtraj2.slingshot_finite_all/rtraj2.cells; ax.plot(rtraj2.replicate,rtraj2.slingshot_fraction,'o-',lw=2,color='#D55E00',label='Slingshot'); ax.axhline(1,color='#0072B2',lw=2,ls=(0,(2,2)),label='Monocle3'); ax.fill_between(rtraj2.replicate,rtraj2.slingshot_fraction,1,color='#D55E00',alpha=.10); ax.set_ylim(.25,1.05); ax.set_xlabel('Bootstrap replicate'); ax.set_ylabel('Finite coverage fraction'); ax.set_title('B  Cell-level pseudotime coverage',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=6); ax.spines[['top','right']].set_visible(False); ax.set_xticks(range(0,11,2))
ax=fig.add_subplot(gs[1,0]); sc=ax.scatter(rtraj.monocle3_graph_nodes,rtraj.monocle3_graph_leaves,c=rtraj.replicate,cmap='viridis',s=46,edgecolor='white',linewidth=.6); ax.set_xlabel('Monocle3 graph nodes'); ax.set_ylabel('Graph leaves'); ax.set_title('C  Graph structure across replicates',loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False); cb=fig.colorbar(sc,ax=ax,pad=.02); cb.set_label('Replicate',fontsize=7); cb.ax.tick_params(labelsize=6)
for _,rr in rtraj.iterrows(): ax.text(rr.monocle3_graph_nodes+8,rr.monocle3_graph_leaves,str(int(rr.replicate)),fontsize=6,color='#344054')
ax=fig.add_subplot(gs[1,1]); fm=method_v[['pass','abstention','invalid_result','hard_failure']].to_numpy(); ax.set_xlim(-.5,3.5); ax.set_ylim(-.7,2.7); ax.set_xticks(range(4),['PASS','ABSTAIN','INVALID','HARD\nFAIL'],fontsize=6); ax.set_yticks(range(3),labels_v,fontsize=7); ax.invert_yaxis(); ax.set_title('D  Failure-mode record',loc='left',fontweight='bold'); ax.spines[:].set_visible(False); ax.grid(False)
for i in range(fm.shape[0]):
    for j in range(fm.shape[1]):
        val=int(fm[i,j]); color='#009E73' if j==0 else ('#D55E00' if j==3 and val else '#D0D5DD'); ax.scatter(j,i,s=380,color=color,alpha=.92,edgecolor='white',linewidth=.8); ax.text(j,i,str(val),ha='center',va='center',fontsize=8,color='white' if (j==0 or (j==3 and val)) else '#344054',fontweight='bold')
fig.suptitle('Figure 4. Trajectory execution stability and failure modes',y=.995,fontsize=12)
save(fig,'Figure4_trajectory_stability_v3')
with open(OUT/'Figure4_v3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 4. Trajectory execution stability and failure modes. (A) Formal replicate readiness for the three tested method layers. (B) Finite cell-level pseudotime coverage across replicates; the Monocle3 line denotes complete finite coverage in each replicate. (C) Monocle3 graph nodes and leaves across replicates, with replicate identifiers. (D) Preserved failure-mode counts. These panels describe technical execution only and do not establish lineage or biological direction.\n')
# Figure 5 v3: polished cross-cohort non-replication figure.
fig=plt.figure(figsize=(8.2,6.0)); gs=fig.add_gridspec(2,2,height_ratios=[1.28,1.0],hspace=.62,wspace=.42)
cohort_order=['GSE315534','GSE178318','Moorman_2024']; cohort_palette={'GSE315534':'#E69F00','GSE178318':'#56B4E9','Moorman_2024':'#CC79A7'}
ax=fig.add_subplot(gs[0,:]); ax.set_title('A  Independent cohort estimates',loc='left',fontweight='bold',pad=8)
ys=[]; ylabels=[]
for j,m in enumerate(['ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3']):
    sub=prog[prog.module.eq(m)].set_index('dataset').loc[cohort_order].reset_index(); base=j*4; 
    for i,rr in sub.iterrows():
        y=base+i; ys.append(y); ylabels.append(f'{str(rr.dataset).replace("Moorman_2024","Moorman 2024")}  (n={int(rr.n)})')
        c=cohort_palette[rr.dataset]; ax.plot([rr.bootstrap_ci95_low,rr.bootstrap_ci95_high],[y,y],color=c,lw=3,solid_capstyle='round'); ax.scatter(rr.mean_delta,y,s=58,color=c,edgecolor='white',linewidth=.9,zorder=3)
    ax.text(-.29,base+1,'Absorptive' if j==0 else 'Secretory',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8,fontweight='bold',color='#344054')
ax.axvline(0,color='#344054',lw=.9); ax.set_yticks(ys,ylabels,fontsize=6.5); ax.set_xlabel('Primary-to-metastasis mean delta (95% bootstrap CI)'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0,pad=4); ax.grid(axis='x',alpha=.18); ax.set_xlim(-.24,.09)
ax=fig.add_subplot(gs[1,0]); sign=prog.pivot(index='dataset',columns='module',values='mean_delta').reindex(cohort_order); im=ax.imshow(sign.values,cmap='PuOr',vmin=-.2,vmax=.2,aspect='auto'); ax.set_title('B  Direction by cohort',loc='left',fontweight='bold',pad=8); ax.set_xticks(range(2),['Absorptive','Secretory'],fontsize=7); ax.set_yticks(range(3),['GSE315534','GSE178318','Moorman 2024'],fontsize=6.5); ax.set_xlabel('Frozen module');
for i in range(sign.shape[0]):
    for j in range(sign.shape[1]): ax.text(j,i,f'{sign.iloc[i,j]:+.3f}',ha='center',va='center',fontsize=8,color='#1D2939',fontweight='bold')
ax.spines[:].set_visible(False)
ax=fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('C  Cohort context',loc='left',fontweight='bold',pad=8)
ctx=[('GSE315534','n=3 • patient','supportive CNV proxy','#FFF3D6','#E69F00'),('GSE178318','n=6 • patient','restricted CNV evidence','#E6F2F8','#56B4E9'),('Moorman 2024','n=6 • donor','public malignant labels','#F3EAF2','#CC79A7')]
for i,(name,note,evi,fc,edge) in enumerate(ctx):
    y=.78-i*.27; ax.add_patch(plt.Rectangle((.03,y-.17),.92,.21,transform=ax.transAxes,fc=fc,ec=edge,lw=.9)); ax.text(.07,y-.03,name,transform=ax.transAxes,fontsize=8,fontweight='bold',color='#1D2939'); ax.text(.07,y-.12,note+'  |  '+evi,transform=ax.transAxes,fontsize=6.5,color='#475467')
ax.text(.03,.015,'Cohorts shown separately; no pooled estimate or direction vote.',transform=ax.transAxes,fontsize=6.0,color='#667085')
fig.suptitle('Figure 5. Cross-cohort progression estimates without pooling',y=.995,fontsize=12)
save(fig,'Figure5_cross_cohort_progression_v3')
with open(OUT/'Figure5_v3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 5. Cross-cohort progression estimates without pooling. (A) Frozen patient/donor-level mean deltas and percentile bootstrap 95% intervals for each cohort and module. (B) Signed estimates demonstrate heterogeneous direction across cohorts. (C) Cohort context and malignancy-evidence tier. Cohorts were not pooled and no direction vote was performed; all estimates are descriptive.\n')

# Figure 2 v2: frozen state-program foundation and patient-level reproducibility.
s05 = ROOT/'analysis/05_state_discovery/runs/20260902_S05_normalization_corrective_audit_v1'
primary_genes = pd.read_csv(s05/'genesets/module_genes_primary.tsv', sep='\t')
primary_genes = primary_genes[primary_genes.module.isin(['ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3'])].copy()
primary_genes['module_label'] = primary_genes.module.map({'ABSORPTIVE_BEST4_CA4':'Absorptive','SECRETORY_MUC2_TFF3':'Secretory'})
primary_genes.to_csv(OUT/'Figure2_v2_frozen_primary_genes.tsv', sep='\t', index=False)

stab = pd.read_csv(s05/'stability/patient_bootstrap_marker_stability.tsv', sep='\t')
loo = pd.read_csv(s05/'stability/leave_one_patient_marker_stability.tsv', sep='\t')
# Summaries are restricted to discovery/support datasets and retain patient-level units.
stab_s = stab.groupby('dataset').agg(clusters=('cluster','nunique'), patients_evaluable=('patients_evaluable','sum'),
                                     bootstrap_jaccard_median=('marker_jaccard_median','median')).reset_index()
loo_s = loo.groupby('dataset').agg(loo_jaccard_median=('loo_jaccard_median','median'),
                                   loo_supported=('loo_status',lambda x: int((x=='SUPPORTED').sum()))).reset_index()
stability = stab_s.merge(loo_s,on='dataset',how='left')
stability.to_csv(OUT/'Figure2_v2_stability_summary.tsv', sep='\t', index=False)

ps = pd.read_csv(s05/'consolidation/patient_pseudobulk_module_scores.tsv', sep='\t')
ps = ps[ps.module.isin(['ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3'])].copy()
ps['module_label'] = ps.module.map({'ABSORPTIVE_BEST4_CA4':'Absorptive','SECRETORY_MUC2_TFF3':'Secretory'})
ps.to_csv(OUT/'Figure2_v2_patient_pseudobulk_scores.tsv', sep='\t', index=False)

fig = plt.figure(figsize=(8.2, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.22, 1.0], hspace=.62, wspace=.38)
module_colors = {'Absorptive':'#D55E00','Secretory':'#0072B2'}

# A. Gene-set architecture: compact rank-ordered lollipops with support encoded by size.
ax = fig.add_subplot(gs[0,0])
ax.set_title('A  Frozen primary gene programs', loc='left', fontweight='bold', pad=8)
gene_ticks, gene_labels = [], []
for i, lab in enumerate(['Absorptive','Secretory']):
    sub = primary_genes[primary_genes.module_label.eq(lab)].sort_values(['dataset_support','median_rank'], ascending=[False,True]).head(10).copy()
    y = np.arange(len(sub)) + (0 if i==0 else 13)
    gene_ticks.extend(y); gene_labels.extend(sub.gene)
    ax.hlines(y, 0, sub.dataset_support, color=module_colors[lab], alpha=.22, lw=5)
    ax.scatter(sub.dataset_support, y, s=30+sub.source_clusters*10, color=module_colors[lab], edgecolor='white', linewidth=.6, zorder=3)
    # Module labels sit in the inter-program whitespace, away from gene marks.
    label_y = 23.3 if lab == 'Secretory' else 10.7
    ax.text(.98, label_y, lab, ha='right', va='center', fontsize=8, fontweight='bold', color=module_colors[lab])
ax.set_yticks(gene_ticks, gene_labels, fontsize=6)
ax.set_xlim(0,4.05); ax.set_ylim(-1,23); ax.set_xlabel('Independent dataset support (0–3)'); ax.set_ylabel('Representative genes (ranked)')
ax.text(.98,.03,'Dot area = source-cluster support', transform=ax.transAxes, ha='right', fontsize=6.3, color='#667085')
ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y', length=0)

# B. Cross-dataset marker support, shown as a bounded distribution rather than a pooled claim.
ax = fig.add_subplot(gs[0,1])
ax.set_title('B  Cross-dataset marker stability', loc='left', fontweight='bold', pad=8)
datasets = ['GSE132465','GSE144735','GSE200997','GSE132257']
for j, ds in enumerate(datasets):
    ss = stab[stab.dataset.eq(ds)]
    vals = ss.marker_jaccard_median.dropna().to_numpy()
    if len(vals):
        vp = ax.violinplot(vals, positions=[j], widths=.72, showextrema=False)
        for body in vp['bodies']:
            body.set_facecolor('#D0D5DD'); body.set_edgecolor('#667085'); body.set_alpha(.8)
        ax.scatter(np.full(len(vals),j), vals, s=18, color='#344054', alpha=.72, edgecolor='white', linewidth=.35, zorder=3)
        ax.scatter(j, np.median(vals), s=46, color='#009E73', edgecolor='white', linewidth=.7, zorder=4)
ax.set_xticks(range(len(datasets)), ['GSE132465','GSE144735','GSE200997','GSE132257\ntechnical only'], rotation=28, ha='right', fontsize=6)
ax.set_ylabel('Patient-bootstrap marker Jaccard')
ax.set_ylim(.25,1.12); ax.axhline(.5, color='#98A2B3', lw=.7, ls=(0,(2,2)))
ax.text(.02,.04,'Dot = cluster; green = dataset median', transform=ax.transAxes, fontsize=6.5, color='#667085')
ax.spines[['top','right']].set_visible(False)

# C. Patient-level pseudobulk score distributions, with cell count encoded by point size.
ax = fig.add_subplot(gs[1,0])
ax.set_title('C  Patient-level pseudobulk scores', loc='left', fontweight='bold', pad=8)
plot_ps = ps.copy(); plot_ps['x'] = plot_ps.module_label.map({'Absorptive':0,'Secretory':1})
rng = np.random.default_rng(19)
for lab, x0 in [('Absorptive',0),('Secretory',1)]:
    sub = plot_ps[plot_ps.module_label.eq(lab)]
    jitter = rng.uniform(-.16,.16,len(sub))
    ax.scatter(x0+jitter, sub.module_score, s=18+np.sqrt(sub.cells)*2, color=module_colors[lab], alpha=.72, edgecolor='white', linewidth=.45, label=lab)
    if len(sub): ax.plot([x0-.22,x0+.22],[sub.module_score.median()]*2, color='#1D2939', lw=1.4)
ax.set_xticks([0,1],['Absorptive','Secretory']); ax.set_ylabel('Module score (patient pseudobulk)')
ax.text(.02,.04,'Point area scales with cells; horizontal line = median', transform=ax.transAxes, fontsize=6.5, color='#667085')
ax.spines[['top','right']].set_visible(False)

# D. Interpretation guardrail: explicit evidence contract.
ax = fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  What this figure supports', loc='left', fontweight='bold', pad=8)
cards = [('STATE','Reproducible expression programs','#E8F5EF','#009E73'),
         ('PATIENT','Patient-level score summaries','#E6F2F8','#0072B2'),
         ('BOUNDARY','Not lineage or transition probability','#FDECEC','#D55E00')]
for i,(head,desc,fc,edge) in enumerate(cards):
    y=.84-i*.245; ax.add_patch(plt.Rectangle((.04,y-.17),.90,.20,transform=ax.transAxes,fc=fc,ec=edge,lw=.9))
    ax.text(.08,y-.03,head,transform=ax.transAxes,fontsize=8,fontweight='bold',color=edge)
    ax.text(.08,y-.12,desc,transform=ax.transAxes,fontsize=7,color='#344054')
ax.text(.04,.045,'S05 modules frozen before S09; this figure summarizes expression-state evidence only.', transform=ax.transAxes, fontsize=6.3, color='#667085', wrap=True)
fig.suptitle('Figure 2. Frozen state programs and patient-level reproducibility', y=.995, fontsize=12)
save(fig,'Figure2_state_program_foundation_v3')
primary_genes.to_csv(OUT/'Figure2_v3_frozen_primary_genes.tsv', sep='\t', index=False)
stability.to_csv(OUT/'Figure2_v3_stability_summary.tsv', sep='\t', index=False)
ps.to_csv(OUT/'Figure2_v3_patient_pseudobulk_scores.tsv', sep='\t', index=False)
with open(OUT/'Figure2_v3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 2. Frozen state programs and patient-level reproducibility. (A) Representative genes from the two frozen primary programs, ordered by cross-dataset support; point size reflects source-cluster support. (B) Patient-bootstrap marker stability by dataset and cluster; green points denote dataset medians. (C) Patient-level pseudobulk module scores, with point area proportional to contributing cells. (D) Interpretation boundary: these data support reproducible expression states and patient-level score summaries, not lineage, transition probability or disease progression direction. S05 state programs were frozen before S09; this figure summarizes expression-state evidence only.\n')

# Figure 3 v3: polished patient-level endpoint figure.
primary3 = pair[pair.cohort.eq('primary')].copy()
primary3['delta_absorptive'] = primary3.absorptive_tumor - primary3.absorptive_normal
primary3['delta_secretory'] = primary3.secretory_tumor - primary3.secretory_normal
boot3 = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_bootstrap_ci.tsv', sep='\t')
sens3 = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_paired_state_concordance.tsv', sep='\t')
sens3['layer'] = np.where(sens3.cohort.eq('primary'), 'Primary (n=16)', 'Sensitivity (n=2)')
fig = plt.figure(figsize=(8.2, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=.64, wspace=.40)

# A. Paired slopes for both frozen programs, with shared patient ordering.
ax = fig.add_subplot(gs[0,0])
ax.set_title('A  Paired patient-level scores', loc='left', fontweight='bold', pad=8)
xvals = [0, 1, 3, 4]
for _, rr in primary3.iterrows():
    ax.plot([0,1], [rr.absorptive_normal, rr.absorptive_tumor], color='#B8C0CC', lw=.7, alpha=.48)
    ax.plot([3,4], [rr.secretory_normal, rr.secretory_tumor], color='#B8C0CC', lw=.7, alpha=.48)
ax.scatter(np.zeros(len(primary3)), primary3.absorptive_normal, s=18, color=COLORS['normal'], edgecolor='white', linewidth=.4)
ax.scatter(np.ones(len(primary3)), primary3.absorptive_tumor, s=18, color=COLORS['primary'], edgecolor='white', linewidth=.4)
ax.scatter(np.full(len(primary3),3), primary3.secretory_normal, s=18, color=COLORS['normal'], edgecolor='white', linewidth=.4)
ax.scatter(np.full(len(primary3),4), primary3.secretory_tumor, s=18, color=COLORS['primary'], edgecolor='white', linewidth=.4)
ax.set_xticks(xvals, ['N','P','N','P']); ax.set_xlabel('Normal (N) → primary (P)'); ax.set_ylabel('Module score (CP10K, log1p)')
ax.axvline(2, color='#D0D5DD', lw=.8); ax.text(.5,.98,'Absorptive',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7,color=COLORS['primary'],fontweight='bold'); ax.text(3.5,.98,'Secretory',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7,color=COLORS['normal'],fontweight='bold')
ax.spines[['top','right']].set_visible(False); ax.tick_params(axis='x',length=0)

# B. Main effect estimates with bootstrap intervals.
ax = fig.add_subplot(gs[0,1])
ax.set_title('B  Frozen endpoint estimates', loc='left', fontweight='bold', pad=8)
means = {'Absorptive': primary3.delta_absorptive.mean(), 'Secretory': primary3.delta_secretory.mean()}
ci = {'Absorptive': (boot3.mean_delta_absorptive.quantile(.025), boot3.mean_delta_absorptive.quantile(.975)),
      'Secretory': (boot3.mean_delta_secretory.quantile(.025), boot3.mean_delta_secretory.quantile(.975))}
yy = np.array([1,0]); labs = ['Absorptive','Secretory']; cols = [COLORS['primary'],COLORS['normal']]
for y0, lab, col in zip(yy,labs,cols):
    lo, hi = ci[lab]; ax.plot([lo,hi],[y0,y0], color=col, lw=4, solid_capstyle='round'); ax.scatter(means[lab], y0, s=72, color=col, edgecolor='white', linewidth=1.0, zorder=3); ax.text(hi+.012,y0,f'{means[lab]:+.3f}',va='center',fontsize=7,color='#344054')
ax.axvline(0,color='#344054',lw=.9); ax.set_yticks(yy,labs); ax.set_xlabel('Tumor − normal mean delta (95% bootstrap CI)'); ax.set_xlim(-.32,.18); ax.grid(axis='x',alpha=.18); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
ax.text(.02,.05,'Unit: complete patients; n=16',transform=ax.transAxes,fontsize=6.5,color='#667085')

# C. Bootstrap distributions with interval markers.
ax = fig.add_subplot(gs[1,0])
ax.set_title('C  Complete-patient bootstrap', loc='left', fontweight='bold', pad=8)
for vals, y0, col, lab in [(boot3.mean_delta_absorptive,1,COLORS['primary'],'Absorptive'),(boot3.mean_delta_secretory,0,COLORS['normal'],'Secretory')]:
    counts, edges = np.histogram(vals, bins=28); centers=(edges[:-1]+edges[1:])/2; heights=counts/counts.max()*.34
    ax.bar(centers, heights, bottom=y0-.17, width=np.diff(edges), color=col, alpha=.78, edgecolor='white', linewidth=.25)
    lo, hi = np.quantile(vals,[.025,.975]); ax.plot([lo,hi],[y0+.22,y0+.22],color='#344054',lw=1.4); ax.scatter(np.mean(vals),y0+.22,s=30,color='#344054',zorder=3)
ax.axvline(0,color='#344054',lw=.8); ax.set_yticks([0,1],['Secretory','Absorptive']); ax.set_xlabel('Bootstrap mean delta'); ax.set_ylabel('Program'); ax.set_ylim(-.32,1.48); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0); ax.grid(axis='x',alpha=.15)

# D. Coverage and prespecified sensitivity in one compact audit panel.
subgs = gs[1,1].subgridspec(1,2, width_ratios=[1.25,.75], wspace=.42)
ax = fig.add_subplot(subgs[0,0])
ax.set_title('D  Coverage QC', loc='left', fontweight='bold', pad=8)
for _, rr in primary3.iterrows():
    ax.plot([0,1],[rr.cells_normal,rr.cells_tumor],color='#B8C0CC',lw=.65,alpha=.45)
ax.scatter(np.zeros(len(primary3)),primary3.cells_normal,s=22,color=COLORS['normal'],edgecolor='white',linewidth=.4,label='Normal cells')
ax.scatter(np.ones(len(primary3)),primary3.cells_tumor,s=22,color=COLORS['primary'],edgecolor='white',linewidth=.4,label='Primary cells')
ax.set_yscale('log'); ax.set_xticks([0,1],['Normal','Primary']); ax.set_ylabel('Cells per patient (log scale)'); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.15,which='both')
ax.legend(frameon=False,fontsize=6,loc='upper left')
ax2 = fig.add_subplot(subgs[0,1]); sm=sens3.groupby('layer')[['delta_absorptive','delta_secretory']].mean().reindex(['Primary (n=16)','Sensitivity (n=2)']); xx2=np.arange(2)
for colname, lab, col in [('delta_absorptive','Absorptive',COLORS['primary']),('delta_secretory','Secretory',COLORS['normal'])]:
    ax2.plot(xx2, sm[colname].to_numpy(), 'o-', color=col, lw=2, ms=5, label=lab)
ax2.axhline(0,color='#344054',lw=.7); ax2.set_xticks(xx2,['16','2'],fontsize=6); ax2.set_xlabel('Patients',fontsize=6); ax2.set_ylabel('Mean Δ',fontsize=6); ax2.set_ylim(-.24,.08); ax2.set_title('Sensitivity layer',fontsize=7,loc='left'); ax2.spines[['top','right']].set_visible(False); ax2.tick_params(labelsize=6); ax2.grid(axis='y',alpha=.15); ax2.legend(frameon=False,fontsize=5.5,loc='lower left')
fig.suptitle('Figure 3. GSE294300 patient-level independent validation', y=.995, fontsize=12)
save(fig,'Figure3_GSE294300_validation_v3')
primary3[['patient','cells_normal','cells_tumor','scoreable_normal','scoreable_tumor','delta_absorptive','delta_secretory']].to_csv(OUT/'Figure3_v3_patient_endpoint.tsv',sep='\t',index=False)
boot3.to_csv(OUT/'Figure3_v3_bootstrap.tsv',sep='\t',index=False)
with open(OUT/'Figure3_v3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Figure 3. GSE294300 patient-level independent validation. (A) Paired normal-to-primary module scores for 16 primary patients; each line is a patient. (B) Frozen patient-level mean deltas with complete-patient bootstrap 95% intervals. (C) Bootstrap distributions and percentile intervals for the two programs. (D) Cell-count coverage is shown as a QC audit only, with a prespecified n=2 sensitivity layer inset; cells are not independent replicates. All estimates are descriptive and no pooled-cell inference is used.\n')

# Final manifest for the submission-level figure package.
manifest_final = pd.DataFrame([
    ['Figure1_S09_study_design_v4','Figure1_v4_cohort_roles.tsv; Figure1_v4_evidence_cards.tsv'],
    ['Figure2_state_program_foundation_v3','Figure2_v3_frozen_primary_genes.tsv; Figure2_v3_stability_summary.tsv; Figure2_v3_patient_pseudobulk_scores.tsv'],
    ['Figure3_GSE294300_validation_v3','Figure3_v3_patient_endpoint.tsv; Figure3_v3_bootstrap.tsv'],
    ['Figure4_trajectory_stability_v3','Figure4_v2_method_readiness.tsv; Figure4_v2_r_trajectory_replicates.tsv'],
    ['Figure5_cross_cohort_progression_v3','Figure5_v2_progression_estimates.tsv; Figure5_v2_cohort_context.tsv'],
    ['SupplementaryFigure1_disposition_and_denominators','SupplementaryFigure1_denominators.tsv; SupplementaryFigure1_GSE294300_cell_QC.tsv'],
    ['SupplementaryFigure2_candidate_input_and_QC','SupplementaryFigure2_input_inventory.tsv'],
    ['SupplementaryFigure3_malignancy_CNV_evidence','SupplementaryFigure3_malignancy_CNV_evidence.tsv'],
    ['SupplementaryFigure4_root_sensitivity_and_failures','SupplementaryFigure4_root_sensitivity_and_failures.tsv'],
    ['SupplementaryFigure5_GSE178318_COL17LM_sensitivity','SupplementaryFigure5_GSE178318_COL17LM_sensitivity.tsv'],
    ['SupplementaryFigure6_progression_patient_deltas','SupplementaryFigure6_progression_patient_deltas.tsv'],
], columns=['figure','data_tables'])
manifest_final.to_csv(OUT/'S09_SUBMISSION_FIGURE_MANIFEST.tsv', sep='\t', index=False)

# Supplementary Figure 1: sample/patient/cell disposition and denominator audit.
paired_all = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260912_S09_gse294300_state_concordance_v1/patient_paired_state_concordance.tsv', sep='\t')
paired_all['layer'] = np.where(paired_all.cohort.eq('primary'), 'Primary endpoint', 'Sensitivity layer')
denom_rows = pd.DataFrame([
    ['GSE294300', 'Primary endpoint', 'patient', 16, 'paired normal + primary'],
    ['GSE294300', 'Sensitivity layer', 'patient', 2, 'retained outside primary denominator'],
    ['GSE315534', 'Progression descriptive', 'patient', 3, 'primary + liver metastasis'],
    ['GSE178318', 'Progression descriptive', 'patient', 6, 'primary + liver metastasis'],
    ['Moorman 2024', 'Progression descriptive', 'donor', 6, 'primary + metastasis'],
], columns=['dataset','role','inference_unit','n_independent','design_note'])
denom_rows.to_csv(OUT/'SupplementaryFigure1_denominators.tsv', sep='\t', index=False)
cell_qc = paired_all[['patient','cohort','layer','cells_normal','cells_tumor','scoreable_normal','scoreable_tumor']].copy()
cell_qc.to_csv(OUT/'SupplementaryFigure1_GSE294300_cell_QC.tsv', sep='\t', index=False)

fig = plt.figure(figsize=(8.2, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=.64, wspace=.38)

# A. Disposition flow, using only frozen independent-unit counts.
ax = fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Independent-unit disposition', loc='left', fontweight='bold', pad=8)
flow = [('GSE294300\n18 patients assessed', .50, .77, '#E6F2F8', '#0072B2'),
        ('16 primary\nendpoint', .25, .42, '#E8F5EF', '#009E73'),
        ('2 sensitivity\nlayer', .75, .42, '#F0F4F8', '#56B4E9'),
        ('Progression cohorts\n3 + 6 + 6 units', .50, .12, '#FFF3D6', '#E69F00')]
for txt,x,y,fc,ec in flow:
    ax.text(x,y,txt,ha='center',va='center',transform=ax.transAxes,fontsize=8,fontweight='bold',bbox=dict(boxstyle='round,pad=.55',fc=fc,ec=ec,lw=1.0))
for a,b in [((.50,.68),(.31,.50)),((.50,.68),(.69,.50))]:
    ax.annotate('',xy=b,xytext=a,xycoords=ax.transAxes,arrowprops=dict(arrowstyle='-|>',lw=1.0,color='#667085'))
ax.text(.50,.27,'separate descriptive analysis',ha='center',transform=ax.transAxes,fontsize=6.5,color='#667085')
ax.text(.50,.93,'Patient/donor is the inferential unit',ha='center',transform=ax.transAxes,fontsize=7,color='#667085')

# B. Denominator bars, with role color and no cell-level aggregation.
ax = fig.add_subplot(gs[0,1]); ax.set_title('B  Frozen denominators', loc='left', fontweight='bold', pad=8)
roles = denom_rows.copy(); yy=np.arange(len(roles)); role_c={'Primary endpoint':'#0072B2','Sensitivity layer':'#56B4E9','Progression descriptive':'#009E73'}
ax.barh(yy, roles.n_independent, color=[role_c[r] for r in roles.role], alpha=.80, height=.56)
ax.set_yticks(yy,[f"{r.dataset}\n{r.role}" for r in roles.itertuples()],fontsize=6); ax.invert_yaxis(); ax.set_xlabel('Independent patients or donors (n)'); ax.set_xlim(0,18); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0)
for yi,r in zip(yy,roles.itertuples()): ax.text(r.n_independent+.28,yi,f'{r.n_independent} {r.inference_unit}{"s" if r.n_independent!=1 else ""}',va='center',fontsize=7,color='#344054')

# C. Per-patient cells and scoreable cells are an audit, not independent replicates.
ax = fig.add_subplot(gs[1,0]); ax.set_title('C  GSE294300 cell-count audit', loc='left', fontweight='bold', pad=8)
for x0, col, lab in [(0,'#0072B2','Normal'),(1,'#D55E00','Primary')]:
    vals = cell_qc['cells_normal' if x0==0 else 'cells_tumor'].to_numpy(); ax.scatter(np.full(len(vals),x0),vals,s=24,color=col,alpha=.72,edgecolor='white',linewidth=.4,label=lab)
for _,r in cell_qc.iterrows(): ax.plot([0,1],[r.cells_normal,r.cells_tumor],color='#B8C0CC',lw=.6,alpha=.45)
ax.set_yscale('log'); ax.set_xticks([0,1],['Normal','Primary']); ax.set_ylabel('Cells per patient (log scale)'); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.16,which='both'); ax.legend(frameon=False,fontsize=6,loc='upper left')
ax.text(.03,.04,'Lines pair the same patient; cells are observations.',transform=ax.transAxes,fontsize=6.5,color='#667085')

# D. Scoreability and role boundary.
ax = fig.add_subplot(gs[1,1]); ax.set_title('D  Scoreability and analysis role', loc='left', fontweight='bold', pad=8)
for i,(label,key, col) in enumerate([('Normal','scoreable_normal','#0072B2'),('Primary','scoreable_tumor','#D55E00')]):
    vals=cell_qc[key].to_numpy(); ax.scatter(np.full(len(vals),i),vals,s=24,color=col,alpha=.72,edgecolor='white',linewidth=.4,label=label)
ax.set_yscale('log'); ax.set_xticks([0,1],['Normal','Primary']); ax.set_ylabel('Scoreable cells per patient (log scale)'); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.16,which='both'); ax.legend(frameon=False,fontsize=6,loc='upper left')
ax.text(.03,.04,'Primary n=16 plus sensitivity n=2; sensitivity is not pooled.',transform=ax.transAxes,fontsize=6.5,color='#667085')
fig.suptitle('Supplementary Figure 1. Sample, patient and cell disposition', y=.995, fontsize=12)
save(fig,'SupplementaryFigure1_disposition_and_denominators')
with open(OUT/'SupplementaryFigure1_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 1. Sample, patient and cell disposition. (A) Frozen independent-unit disposition for GSE294300 and the three separate progression cohorts. (B) Independent patient/donor denominators by frozen analysis role. (C) Per-patient cell-count coverage in the GSE294300 paired cohort; lines connect the same patient across arms. (D) Scoreable-cell counts are shown as an input/QC audit. Patients or donors are the inferential units; cells are observations and are not treated as independent replicates. The two sensitivity patients remain outside the primary endpoint denominator and are not pooled with it.\n')

# Supplementary Figure 2: GSE294300 candidate input and QC limitations.
inv = pd.read_csv(ROOT/'data/processed/GSE294300_lenient_candidates_v1/filtered_input_inventory.tsv', sep='\t')
aud = pd.read_csv(ROOT/'data/processed/GSE294300_lenient_candidates_v1/filtered_input_audit.tsv', sep='\t')
inv['retention_fraction'] = inv.retained_barcodes / inv.input_barcodes
inv['arm'] = np.where(inv['sample'].str.contains('nor', case=False), 'Normal', 'Primary')
inv.to_csv(OUT/'SupplementaryFigure2_input_inventory.tsv', sep='\t', index=False)
fig = plt.figure(figsize=(8.2, 7.0)); gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=.64, wspace=.40)
ax = fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Candidate matrix contract', loc='left', fontweight='bold', pad=8)
steps = [('1','Barcode-rank input','1.6–2.1M barcodes','#667085'),('2','Retained barcodes','447–13,010 per sample','#0072B2'),('3','Feature space','62,700 rows','#009E73'),('4','State-score input','Paired patient summary','#E69F00')]
for i,(num,label,value,ec) in enumerate(steps):
    y=.72-i*.16; ax.add_patch(plt.Rectangle((.08,y-.057),.84,.114,transform=ax.transAxes,fc='#F8FAFC',ec='#D0D5DD',lw=.6)); ax.add_patch(plt.Rectangle((.08,y-.057),.035,.114,transform=ax.transAxes,fc=ec,ec=ec,lw=0)); ax.text(.098,y,num,ha='center',va='center',transform=ax.transAxes,fontsize=7,fontweight='bold',color='white'); ax.text(.16,y+.018,label,ha='left',va='center',transform=ax.transAxes,fontsize=7.0,fontweight='bold',color='#1D2939'); ax.text(.16,y-.027,value,ha='left',va='center',transform=ax.transAxes,fontsize=6.5,color='#475467')
ax.text(.5,.88,'Provisional barcode-rank-derived candidate input',ha='center',transform=ax.transAxes,fontsize=7.5,color='#667085')
ax.text(.5,.06,'Raw counts preserved separately; candidate layer shown.',ha='center',transform=ax.transAxes,fontsize=6.2,color='#667085')
ax = fig.add_subplot(gs[0,1]); ax.set_title('B  Barcode retention by sample', loc='left', fontweight='bold', pad=8)
order = inv.sort_values('retained_barcodes').reset_index(drop=True); yy=np.arange(len(order)); colors=np.where(order.arm.eq('Normal'),'#56B4E9','#D55E00')
ax.barh(yy, order['input_barcodes'], color='#D0D5DD', height=.72, label='Input barcodes'); ax.barh(yy, order['retained_barcodes'], color=colors, height=.46, label='Retained barcodes'); ax.set_xscale('log');
tick_idx = np.arange(0, len(order), 4); tick_labels = [str(i+1) for i in tick_idx]
ax.set_yticks(tick_idx, tick_labels, fontsize=6); ax.set_ylabel('Sample rank', fontsize=7); ax.set_xlabel('Barcodes (log scale)'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0); ax.legend(frameon=False,fontsize=6,loc='lower right')
ax.text(.02,.03,'Color: normal vs primary',transform=ax.transAxes,fontsize=6.2,color='#667085')
ax = fig.add_subplot(gs[1,0]); ax.set_title('C  Retention fraction', loc='left', fontweight='bold', pad=8)
for x0, arm, col in [(0,'Normal','#0072B2'),(1,'Primary','#D55E00')]:
    vals=inv.loc[inv.arm.eq(arm),'retention_fraction']; jitter=np.random.default_rng(7+x0).uniform(-.12,.12,len(vals)); ax.scatter(np.full(len(vals),x0)+jitter, vals, s=24, color=col, alpha=.75, edgecolor='white', linewidth=.4); ax.plot([x0-.2,x0+.2],[vals.median()]*2,color='#344054',lw=1.5)
ax.set_xticks([0,1],['Normal','Primary']); ax.set_ylabel('Retained / input barcodes'); ax.set_ylim(0,.012); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.15); ax.text(.03,.04,'Line = arm median',transform=ax.transAxes,fontsize=6.3,color='#667085')
ax = fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  Input limitations to carry forward', loc='left', fontweight='bold', pad=8)
cards=[('PROVISIONAL INPUT','Barcode-rank-derived candidate layer','#E6F2F8','#0072B2'),('NOT AUTHOR MATRIX','Do not claim reproduction of a filtered matrix','#FDECEC','#D55E00'),('QC STATUS','36 sample matrices passed structural audit','#E8F5EF','#009E73'),('INFERENCE','Patients are units; cells support coverage QC','#FFF3D6','#E69F00')]
for i,(head,desc,fc,edge) in enumerate(cards):
    y=.83-i*.22; ax.add_patch(plt.Rectangle((.04,y-.14),.92,.17,transform=ax.transAxes,fc=fc,ec=edge,lw=.9)); ax.text(.08,y-.02,head,transform=ax.transAxes,fontsize=7.5,fontweight='bold',color=edge); ax.text(.08,y-.095,desc,transform=ax.transAxes,fontsize=6.6,color='#344054')
fig.suptitle('Supplementary Figure 2. GSE294300 candidate input and QC boundaries', y=.995, fontsize=12); save(fig,'SupplementaryFigure2_candidate_input_and_QC')
with open(OUT/'SupplementaryFigure2_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 2. GSE294300 candidate input and QC boundaries. (A) The analyzed layer is a provisional barcode-rank-derived candidate matrix with 62,700 feature rows; raw counts remain preserved separately. (B) Input and retained barcode counts are shown for all 36 sample matrices. (C) Retention fractions are shown by arm, with horizontal lines denoting arm medians. (D) Structural matrix checks passed for all 36 samples, but the candidate layer is not asserted to reproduce an author-provided filtered matrix. Patients remain the inferential units; cells are used for coverage and scoreability audits only.\n')

# Supplementary Figure 3: malignancy/CNV evidence tiers and availability boundaries.
evidence_rows = pd.DataFrame([
    ['GSE294300','Auxiliary genomic consistency','CNV audit available; not malignant-cell truth','Normal/primary paired','Primary endpoint; conditional descriptive'],
    ['GSE315534','Supportive proxy','Author-reported CNV summary; no barcode-level CNV truth','Primary/liver metastasis','Descriptive progression only'],
    ['GSE178318','Restricted / uncertain','Paper-level InferCNV; no public recomputable cell-level label; COL17-LM low-count','Primary/liver metastasis','Descriptive; sensitivity limitation'],
    ['Moorman 2024','Stronger public label','Exact public malignant-cell label and non-malignant epithelial reference','Primary/metastasis','Descriptive progression only'],
], columns=['dataset','evidence_tier','evidence_basis','arm_context','permitted_role'])
evidence_rows.to_csv(OUT/'SupplementaryFigure3_malignancy_CNV_evidence.tsv', sep='\t', index=False)

fig = plt.figure(figsize=(8.2, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=.68, wspace=.40)
# A. Evidence ladder.
ax = fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Evidence tiers are not interchangeable', loc='left', fontweight='bold', pad=8)
ladder=[('Stronger public label','Moorman 2024',['Malignant-cell label + normal','epithelial reference'],'#E8F5EF','#009E73'),('Supportive proxy','GSE315534',['CNV summary; no barcode-level','truth'],'#FFF3D6','#E69F00'),('Restricted / uncertain','GSE178318',['Paper-level InferCNV; no public','cell-level label'],'#FDECEC','#D55E00'),('Auxiliary consistency','GSE294300',['CNV consistency audit; not','malignant truth'],'#E6F2F8','#0072B2')]
for i,(tier,ds,desc,fc,ec) in enumerate(ladder):
    y=.84-i*.205
    ax.add_patch(plt.Rectangle((.035,y-.155),.93,.20,transform=ax.transAxes,fc=fc,ec=ec,lw=1.1))
    ax.add_patch(plt.Rectangle((.035,y-.155),.055,.20,transform=ax.transAxes,fc=ec,ec=ec,lw=0))
    ax.text(.0625,y-.055,f'{i+1}',transform=ax.transAxes,ha='center',va='center',fontsize=9,fontweight='bold',color='white')
    ax.text(.115,y-.018,tier,transform=ax.transAxes,fontsize=8.0,fontweight='bold',color=ec)
    ax.text(.115,y-.078,ds,transform=ax.transAxes,fontsize=6.9,fontweight='bold',color='#344054')
    ax.text(.115,y-.123,desc[0]+' '+desc[1],transform=ax.transAxes,fontsize=6.35,color='#475467')
ax.text(.035,.025,'Tier describes evidence availability, not biological certainty.',transform=ax.transAxes,fontsize=6.6,color='#667085')

# B. Cohort-by-evidence matrix.
ax = fig.add_subplot(gs[0,1]); ax.set_title('B  Cohort evidence matrix', loc='left', fontweight='bold', pad=8)
cols=['Public\nmalignant','CNV /\ngenomic proxy','Normal\nepithelial ref.','Barcode-level\ntruth']
mat=np.array([[0,1,1,0],[0,1,0,0],[0,1,0,0],[1,0,1,1]])
im=ax.imshow(mat,cmap=plt.matplotlib.colors.ListedColormap(['#EAECF0','#009E73']),vmin=0,vmax=1,aspect='auto')
ax.set_xticks(range(4),cols,fontsize=6); ax.set_yticks(range(4),['GSE294300','GSE315534','GSE178318','Moorman 2024'],fontsize=6); ax.tick_params(length=0); ax.spines[:].set_visible(False)
for i in range(4):
    for j in range(4): ax.text(j,i,'available' if mat[i,j] else 'limited',ha='center',va='center',fontsize=5.5,color='white' if mat[i,j] else '#667085')
ax.text(.02,-.18,'Availability matrix; “limited” is not negative evidence.',transform=ax.transAxes,fontsize=6.2,color='#667085')

# C. Dataset-specific restrictions.
ax = fig.add_subplot(gs[1,0]); ax.set_title('C  Restrictions carried into analysis', loc='left', fontweight='bold', pad=8)
restr=[('GSE294300','Auxiliary CNV only'),('GSE315534','No barcode-level malignancy truth'),('GSE178318','COL17-LM low-count sensitivity'),('Moorman 2024','Small donor cohort; descriptive')]
for i,(ds,desc) in enumerate(restr):
    y=3-i; ax.barh(y,1,color=['#E6F2F8','#FFF3D6','#FDECEC','#F3EAF2'][i],edgecolor=['#0072B2','#E69F00','#D55E00','#CC79A7'][i],height=.58); ax.text(.04,y,ds,va='center',fontsize=7,fontweight='bold',color='#1D2939'); ax.text(.96,y,desc,ha='right',va='center',fontsize=6.5,color='#344054')
ax.set_xlim(0,1); ax.set_ylim(-.7,3.7); ax.axis('off')

# D. Boundary statement.
ax = fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  Interpretation boundary', loc='left', fontweight='bold', pad=8)
ax.add_patch(plt.Rectangle((.05,.23),.90,.50,transform=ax.transAxes,fc='#F8FAFC',ec='#667085',lw=1.0)); ax.text(.10,.61,'MALIGNANCY / CNV EVIDENCE',transform=ax.transAxes,fontsize=8,fontweight='bold',color='#344054'); ax.text(.10,.49,'Supports evidence grading and',transform=ax.transAxes,fontsize=7,color='#344054'); ax.text(.10,.42,'cohort context.',transform=ax.transAxes,fontsize=7,color='#344054'); ax.text(.10,.33,'It does not define state, lineage',transform=ax.transAxes,fontsize=7,color='#D55E00',fontweight='bold'); ax.text(.10,.26,'or transition.',transform=ax.transAxes,fontsize=7,color='#D55E00',fontweight='bold')
ax.text(.05,.08,'All progression results remain separate descriptive analyses.',transform=ax.transAxes,fontsize=6.5,color='#667085')
fig.suptitle('Supplementary Figure 3. Malignancy and CNV evidence tiers by cohort', y=.995, fontsize=12); save(fig,'SupplementaryFigure3_malignancy_CNV_evidence')
with open(OUT/'SupplementaryFigure3_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 3. Malignancy and CNV evidence tiers by cohort. (A) Evidence tiers reflect the availability and reproducibility of malignancy/CNV support and are not interchangeable. (B) Cohort-by-evidence availability matrix; limited evidence is not treated as negative evidence. (C) Restrictions carried into the descriptive analyses. (D) CNV and malignancy evidence support evidence grading and cohort context, but do not define expression state, lineage or transition. Progression cohorts remain separate descriptive analyses.\n')

# Supplementary Figure 4: S08 root sensitivity, resource layer, failures and abstention.
bootstrap4 = pd.DataFrame([
    ['PAGA', 10, 10, 0, 0, 'execution only'],
    ['Palantir', 10, 8, 2, 0, '2 ARPACK abstentions'],
    ['Slingshot', 10, 6, 4, 0, '4 graph/linear-system abstentions'],
    ['Monocle3', 10, 10, 0, 0, 'cap=2,000; all finite'],
], columns=['method','replicates','pass','abstention','invalid','note'])
root4 = pd.DataFrame([
    ['Palantir','full input',6,6,0,0,'PASS','0.437–0.999'],
    ['Monocle3','cap=500 cells/patient',6,6,0,0,'PASS','0.065–0.847'],
    ['Slingshot','cap=500 cells/patient',9,0,0,9,'INVALID_RESULT','0.768–0.948 finite cells only'],
], columns=['method','resource_layer','attempts','pass','abstention','invalid','status','absolute_spearman'])
audit4 = pd.DataFrame([
    ['S08 patient bootstrap','PAGA',10,10,0,'PASS','patient-level resampling'],
    ['S08 patient bootstrap','Palantir',10,8,2,'LIMITED','ARPACK starting vector zero'],
    ['S08 patient bootstrap','Slingshot',10,6,4,'LIMITED','singular/NA graph failures'],
    ['S08 patient bootstrap','Monocle3',10,10,0,'PASS','finite pseudotime in all replicates'],
    ['S08 root sensitivity','Palantir',6,6,0,'PASS_RESOURCE_LIMITED','full input'],
    ['S08 root sensitivity','Monocle3',6,6,0,'PASS_RESOURCE_LIMITED','cap=500 cells/patient'],
    ['S08 root sensitivity','Slingshot',9,0,0,'INVALID_RESULT','cap=500; incomplete finite pseudotime'],
], columns=['audit','method','attempts','pass','abstention','status','interpretation'])
pd.concat([bootstrap4.assign(audit='patient_bootstrap'), root4.assign(audit='root_sensitivity')], ignore_index=True, sort=False).to_csv(OUT/'SupplementaryFigure4_root_sensitivity_and_failures.tsv', sep='\t', index=False)

fig = plt.figure(figsize=(8.2, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=.68, wspace=.42)
# A. Resource layers and root sensitivity coverage.
ax = fig.add_subplot(gs[0,0]); ax.set_title('A  Root-sensitivity resource layers', loc='left', fontweight='bold', pad=8)
methods=['Palantir','Monocle3','Slingshot']; y=np.arange(3); attempts=[6,6,9]; passed=[6,6,0]; invalid=[0,0,9]
ax.barh(y, attempts, color='#EAECF0', height=.62, label='Attempts'); ax.barh(y, passed, color='#009E73', height=.42, label='PASS'); ax.barh(y, invalid, left=passed, color='#D55E00', height=.42, label='INVALID_RESULT')
ax.set_yticks(y,methods); ax.set_xlabel('Candidate root runs'); ax.set_xlim(0,10); ax.invert_yaxis(); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0); ax.grid(axis='x',alpha=.15); ax.legend(frameon=False,fontsize=6,loc='upper right')
ax.text(.02,.015,'Palantir: full input; Monocle3/Slingshot: cap=500 cells/patient',transform=ax.transAxes,fontsize=5.8,color='#667085')
# B. Patient bootstrap outcome counts.
ax = fig.add_subplot(gs[0,1]); ax.set_title('B  Patient-bootstrap execution status', loc='left', fontweight='bold', pad=8)
y=np.arange(4); ax.barh(y, bootstrap4['pass'], color='#009E73', height=.62, label='PASS'); ax.barh(y, bootstrap4['abstention'], left=bootstrap4['pass'], color='#E69F00', height=.62, label='ABSTENTION'); ax.set_yticks(y,bootstrap4['method']); ax.set_xlim(0,10); ax.set_xlabel('Replicates (n=10)'); ax.invert_yaxis(); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0); ax.grid(axis='x',alpha=.15); ax.legend(frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,handlelength=1.4,columnspacing=1.0)
for yi,row in bootstrap4.iterrows(): ax.text(10.15,yi,row['note'],va='center',fontsize=5.8,color='#475467',clip_on=False)
# C. Root sensitivity result and correlation range.
ax = fig.add_subplot(gs[1,0]); ax.set_title('C  Root-sensitivity result', loc='left', fontweight='bold', pad=8)
labels=['Palantir','Monocle3','Slingshot']; vals=[1,1,0]; cols=['#009E73','#009E73','#D55E00']; ax.bar(np.arange(3),vals,color=cols,width=.62); ax.set_ylim(0,1.25); ax.set_yticks([0,1],['not valid','all valid']); ax.set_xticks(np.arange(3),labels,rotation=20,ha='right'); ax.set_ylabel('Finite pseudotime / valid output'); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.15); ax.text(.03,.91,'Slingshot: 9/9 invalid under frozen nonfinite rule',transform=ax.transAxes,fontsize=6.0,color='#667085')
# D. Interpretation boundary.
ax = fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  What this audit supports', loc='left', fontweight='bold', pad=8)
cards=[('EXECUTION STABILITY','Reports whether a frozen run completed','#E8F5EF','#009E73'),('FAILURE / ABSTENTION','Failures remain visible; no repair after inspection','#FFF3D6','#E69F00'),('ROOT DEPENDENCE','Root changes are a technical sensitivity layer','#E6F2F8','#0072B2'),('NOT A BIOLOGICAL CLAIM','No software ranking, lineage or direction inference','#FDECEC','#D55E00')]
for i,(head,desc,fc,edge) in enumerate(cards):
    y=.86-i*.22; ax.add_patch(plt.Rectangle((.05,y-.14),.90,.17,transform=ax.transAxes,fc=fc,ec=edge,lw=.9)); ax.text(.09,y-.02,head,transform=ax.transAxes,fontsize=7.1,fontweight='bold',color=edge); ax.text(.09,y-.095,desc,transform=ax.transAxes,fontsize=6.2,color='#344054')
fig.suptitle('Supplementary Figure 4. S08 root sensitivity and execution failure audit', y=.995, fontsize=12); save(fig,'SupplementaryFigure4_root_sensitivity_and_failures')
with open(OUT/'SupplementaryFigure4_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 4. S08 root sensitivity and execution failure audit. (A) Candidate-root runs are shown by method and frozen resource layer; green denotes PASS and orange denotes INVALID_RESULT. (B) Patient-level bootstrap execution status for the 10 formal resamples; abstentions are retained rather than repaired or removed. (C) Root-sensitivity validity summary under the frozen nonfinite-output rule. (D) These outputs document technical execution stability and failure handling only; they do not rank software, establish lineage, or infer biological direction.\n')

# Supplementary Figure 5: GSE178318 COL17-LM low-cell sensitivity layer.
g178 = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260913_0005_S09_gse178318_progression_concordance_v1/sample_state_program_scores.tsv', sep='\t')
p178 = pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260913_0005_S09_gse178318_progression_concordance_v1/patient_paired_state_deltas.tsv', sep='\t')
summary178 = pd.DataFrame([
    ['all pairs', 6, p178['ABSORPTIVE_BEST4_CA4_primary_to_lm'].mean(), p178['SECRETORY_MUC2_TFF3_primary_to_lm'].mean()],
    ['LM candidates >=100 cells', int((~p178.lm_low_count).sum()), p178.loc[~p178.lm_low_count, 'ABSORPTIVE_BEST4_CA4_primary_to_lm'].mean(), p178.loc[~p178.lm_low_count, 'SECRETORY_MUC2_TFF3_primary_to_lm'].mean()],
], columns=['stratum','patients','absorptive_delta','secretory_delta'])
g178_out = g178[['patient','tissue','candidate_cells','low_count_flag','ABSORPTIVE_BEST4_CA4','SECRETORY_MUC2_TFF3']].merge(p178[['patient','lm_low_count','ABSORPTIVE_BEST4_CA4_primary_to_lm','SECRETORY_MUC2_TFF3_primary_to_lm']], on='patient', how='left')
g178_out.to_csv(OUT/'SupplementaryFigure5_GSE178318_COL17LM_sensitivity.tsv', sep='\t', index=False)
fig = plt.figure(figsize=(8.2, 7.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=.82, wspace=.52)
# A. Per-patient cell coverage.
ax = fig.add_subplot(gs[0,0]); ax.set_title('A  Paired cell coverage by patient', loc='left', fontweight='bold', pad=8)
patients=sorted(p178.patient.tolist()); x=np.arange(len(patients)); piv=g178.pivot(index='patient',columns='tissue',values='candidate_cells').reindex(patients); width=.34
ax.bar(x-width/2,piv['primary_CRC'],width,color='#D55E00',label='Primary CRC'); ax.bar(x+width/2,piv['liver_metastasis'],width,color='#56B4E9',label='Liver metastasis'); ax.axhline(100,color='#667085',lw=.8,ls=(0,(2,2))); ax.set_yscale('log'); ax.set_xticks(x,patients,rotation=45,ha='right',fontsize=6); ax.set_ylabel('Marker-gated cells (log scale)'); ax.legend(frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2,handlelength=1.4,columnspacing=1.0); ax.text(.02,.03,'Dashed line = 100-cell sensitivity floor',transform=ax.transAxes,fontsize=6.0,color='#667085')
# B. Low-count flag.
ax = fig.add_subplot(gs[0,1]); ax.set_title('B  Prespecified low-count flag', loc='left', fontweight='bold', pad=8)
flag_counts=g178.groupby(['tissue','low_count_flag']).size().unstack(fill_value=0).reindex(['primary_CRC','liver_metastasis']); flag_counts.plot(kind='bar',stacked=True,ax=ax,color=['#56B4E9','#D55E00'],width=.58); ax.set_xticklabels(['Primary CRC','Liver metastasis'],rotation=0); ax.set_ylabel('Samples (n)'); ax.set_xlabel(''); ax.legend(['>=100 cells','<100 cells'],frameon=False,fontsize=6,loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2,handlelength=1.4,columnspacing=1.0); ax.spines[['top','right']].set_visible(False); ax.text(.03,.03,'Only COL17 liver metastasis is <100 cells (60)',transform=ax.transAxes,fontsize=6.0,color='#667085')
# C. Patient-level delta with COL17 highlighted.
ax = fig.add_subplot(gs[1,0]); ax.set_title('C  Patient-level primary to LM deltas', loc='left', fontweight='bold', pad=8)
x=np.arange(len(p178)); c=np.where(p178.lm_low_count,'#D55E00','#56B4E9'); ax.axhline(0,color='#344054',lw=.8); ax.scatter(x-.12,p178.ABSORPTIVE_BEST4_CA4_primary_to_lm,s=36,color=c,edgecolor='white',linewidth=.5,label='Absorptive'); ax.scatter(x+.12,p178.SECRETORY_MUC2_TFF3_primary_to_lm,s=36,color=c,marker='s',edgecolor='white',linewidth=.5,label='Secretory'); ax.set_xticks(x,p178.patient,rotation=45,ha='right',fontsize=6); ax.set_ylabel('Primary -> liver metastasis delta'); ax.legend(frameon=False,fontsize=6,loc='lower right',handlelength=1.4); ax.text(.03,.93,'Orange = COL17 low-count LM pair',transform=ax.transAxes,fontsize=6.0,color='#667085')
# D. Stratum comparison and boundary.
ax = fig.add_subplot(gs[1,1]); ax.set_title('D  Sensitivity stratum comparison', loc='left', fontweight='bold', pad=8)
xx=np.arange(2); ax.axhline(0,color='#344054',lw=.8); ax.bar(xx-.18,summary178.absorptive_delta,.34,color='#D55E00',label='Absorptive'); ax.bar(xx+.18,summary178.secretory_delta,.34,color='#0072B2',label='Secretory'); ax.set_xticks(xx,['All pairs\n(n=6)','LM >=100\n(n=5)']); ax.set_ylabel('Mean delta'); ax.legend(frameon=False,fontsize=6,loc='center left',bbox_to_anchor=(1.02,.55),handlelength=1.4); ax.spines[['top','right']].set_visible(False)
fig.suptitle('Supplementary Figure 5. GSE178318 COL17-LM low-count sensitivity', y=.995, fontsize=12); save(fig,'SupplementaryFigure5_GSE178318_COL17LM_sensitivity')
with open(OUT/'SupplementaryFigure5_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 5. GSE178318 COL17-LM low-count sensitivity. (A) Marker-gated candidate-cell counts for six paired patients; the dashed line marks the prespecified 100-cell sensitivity floor. (B) The only low-count sample is the COL17 liver-metastasis block (60 cells), retained in the all-pair descriptive analysis. (C) Patient-level primary-to-liver-metastasis deltas are shown with COL17 highlighted. (D) Mean deltas are compared for all six pairs and the five-patient stratum excluding COL17; bootstrap intervals cross zero in both strata. This is a sensitivity and data-coverage audit, not evidence for absence of biological change or a stable progression direction.\n')

# Supplementary Figure 6: all progression patient/donor deltas and frozen bootstrap intervals.
progression_sources = {
    'GSE315534': ROOT/'analysis/09_external_validation/runs/20260912_2320_S09_gse315534_progression_concordance_v1/patient_paired_state_deltas.tsv',
    'GSE178318': ROOT/'analysis/09_external_validation/runs/20260913_0005_S09_gse178318_progression_concordance_v1/patient_paired_state_deltas.tsv',
    'Moorman 2024': ROOT/'analysis/09_external_validation/runs/20260913_0130_S09_moorman_progression_concordance_v1/patient_paired_state_deltas.tsv',
}
delta_rows=[]
for cohort, path in progression_sources.items():
    d=pd.read_csv(path,sep='\t')
    for _,r in d.iterrows():
        unit = r.get('donor_id', r.get('patient'))
        low = bool(r.get('lm_low_count', False))
        for module,col in [('Absorptive','ABSORPTIVE_BEST4_CA4_primary_to_lm'),('Secretory','SECRETORY_MUC2_TFF3_primary_to_lm'),('Absorptive','ABSORPTIVE_BEST4_CA4_primary_to_metastasis'),('Secretory','SECRETORY_MUC2_TFF3_primary_to_metastasis')]:
            if col in d.columns:
                delta_rows.append([cohort,str(unit),module,float(r[col]),low])
delta6=pd.DataFrame(delta_rows,columns=['cohort','unit','module','delta','low_count'])
delta6=delta6.drop_duplicates(['cohort','unit','module'])
delta6.to_csv(OUT/'SupplementaryFigure6_progression_patient_deltas.tsv',sep='\t',index=False)
prog6=pd.read_csv(ROOT/'analysis/09_external_validation/runs/20260913_S09_direction_audit_v1/cross_dataset_primary_to_metastasis_summary.tsv',sep='\t')
fig=plt.figure(figsize=(8.2,7.2)); gs=fig.add_gridspec(2,2,height_ratios=[1.1,1.0],hspace=.72,wspace=.48)
cohorts=['GSE315534','GSE178318','Moorman 2024']; palette={'GSE315534':'#E69F00','GSE178318':'#56B4E9','Moorman 2024':'#CC79A7'}
for ax,mod,title in [(fig.add_subplot(gs[0,0]),'Absorptive','A  Absorptive patient/donor deltas'),(fig.add_subplot(gs[0,1]),'Secretory','B  Secretory patient/donor deltas')]:
    sub=delta6[delta6.module.eq(mod)]
    for i,co in enumerate(cohorts):
        vals=sub[sub.cohort.eq(co)]; ax.scatter(np.full(len(vals),i),vals.delta,s=42,color=palette[co],edgecolor='white',linewidth=.6,zorder=3)
        sm=prog6[(prog6.dataset.eq(co.replace(' ','_')) if co=='Moorman 2024' else prog6.dataset.eq(co)) & prog6.module.str.startswith(mod.upper())]
        if len(sm): ax.plot([i-.16,i+.16],[sm.mean_delta.iloc[0]]*2,color=palette[co],lw=3); ax.vlines(i,sm.bootstrap_ci95_low.iloc[0],sm.bootstrap_ci95_high.iloc[0],color=palette[co],lw=2)
    ax.axhline(0,color='#344054',lw=.8); ax.set_xticks(range(3),['GSE315534','GSE178318','Moorman 2024'],rotation=25,ha='right',fontsize=6); ax.set_ylabel('Primary-to-metastasis delta'); ax.set_title(title,loc='left',fontweight='bold'); ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.15)
ax=fig.add_subplot(gs[1,0]); ax.set_title('C  Frozen bootstrap intervals',loc='left',fontweight='bold'); y=np.arange(6); labels=[]
for i,co in enumerate(cohorts):
    for j,mod in enumerate(['Absorptive','Secretory']):
        row=prog6[(prog6.dataset.eq(co.replace(' ','_')) if co=='Moorman 2024' else prog6.dataset.eq(co)) & prog6.module.str.startswith(mod.upper())].iloc[0]; yy=i*2+j; labels.append(f'{co} | {mod}'); ax.plot([row.bootstrap_ci95_low,row.bootstrap_ci95_high],[yy,yy],color=palette[co],lw=3); ax.scatter(row.mean_delta,yy,color=palette[co],s=38,edgecolor='white',zorder=3)
ax.axvline(0,color='#344054',lw=.8); ax.set_yticks(y,labels,fontsize=6); ax.set_xlabel('Mean delta (95% bootstrap CI)'); ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='y',length=0); ax.grid(axis='x',alpha=.15); ax.invert_yaxis()
ax=fig.add_subplot(gs[1,1]); ax.axis('off'); ax.set_title('D  Interpretation boundary',loc='left',fontweight='bold',pad=8)
cards=[('ALL PATIENTS / DONORS','Every available independent unit is shown','#E8F5EF','#009E73'),('INTERVALS','Bootstrap intervals are cohort-specific','#E6F2F8','#0072B2'),('NO POOLING','No pooled estimate or direction vote','#FFF3D6','#E69F00'),('DESCRIPTIVE ONLY','Heterogeneity is retained, not resolved','#FDECEC','#D55E00')]
for i,(head,desc,fc,edge) in enumerate(cards):
    yy=.84-i*.22; ax.add_patch(plt.Rectangle((.05,yy-.14),.90,.17,transform=ax.transAxes,fc=fc,ec=edge,lw=.9)); ax.text(.09,yy-.02,head,transform=ax.transAxes,fontsize=7.0,fontweight='bold',color=edge); ax.text(.09,yy-.095,desc,transform=ax.transAxes,fontsize=6.2,color='#344054')
fig.suptitle('Supplementary Figure 6. Patient-level progression deltas and uncertainty',y=.995,fontsize=12); save(fig,'SupplementaryFigure6_progression_patient_deltas')
with open(OUT/'SupplementaryFigure6_caption.txt','w',encoding='utf-8') as fh:
    fh.write('Supplementary Figure 6. Patient-level progression deltas and uncertainty. (A,B) All available patient/donor-level primary-to-metastasis deltas are shown for the absorptive and secretory programs; colored horizontal marks and vertical lines denote frozen cohort means and bootstrap 95% intervals. (C) Cohort-specific intervals are displayed without pooling. (D) These are separate descriptive progression analyses; patients/donors are independent units and no direction vote or stable progression claim is made.\n')
print(f'Generated submission figures 1-5 and data tables in {OUT}')
