#!/usr/bin/env python3
"""
Codon Usage Pipeline
--------------------
Calculates and visualizes codon usage metrics across multiple genomic FASTA files.
Metrics include RSCU, CAI, ENC (Effective Number of Codons), and GC3.
Outputs include raw data matrices, heatmaps, PCA, and regression analyses.
"""

from Bio import SeqIO
from Bio.Data import CodonTable
import os, sys, math
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- Setup codons ---
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
ALL_CODONS = [a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"]
SENSE_CODONS = [c for c in ALL_CODONS if c not in standard_table.stop_codons]
CODON2AA = {codon: standard_table.forward_table.get(codon, '*') for codon in ALL_CODONS}
AA2CODONS = defaultdict(list)
for cod in SENSE_CODONS:
    aa = CODON2AA.get(cod, None)
    if aa:
        AA2CODONS[aa].append(cod)

DEGENERACY_GROUPS = defaultdict(list)
for aa, cods in AA2CODONS.items():
    deg = len(cods)
    DEGENERACY_GROUPS[deg].append(aa)

# --- Utility functions ---
def gc_and_gc3(seq):
    seq = seq.upper()
    L = len(seq)
    if L == 0: return np.nan, np.nan
    gc = (seq.count("G") + seq.count("C")) / L
    thirds = seq[2::3]
    gc3 = (thirds.count("G") + thirds.count("C")) / len(thirds) if len(thirds)>0 else np.nan
    return gc, gc3

def count_codons(seq):
    seq = seq.upper()
    counts = Counter()
    L = len(seq) - (len(seq)%3)
    for i in range(0,L,3):
        cod = seq[i:i+3]
        if cod in SENSE_CODONS: counts[cod]+=1
    return counts

def compute_rscu_from_counts(codon_counts):
    rscu = {}
    for aa, cods in AA2CODONS.items():
        total = sum(codon_counts.get(c,0) for c in cods)
        m = len(cods)
        for c in cods:
            obs = codon_counts.get(c,0)
            expected = total/m if total>0 else np.nan
            rscu[c] = (obs/expected) if (expected and expected>0) else 0.0
    return rscu

def compute_codon_freqs_from_counts(codon_counts):
    total = sum(codon_counts.get(c,0) for c in SENSE_CODONS)
    return {c:(codon_counts.get(c,0)/total) if total>0 else 0.0 for c in SENSE_CODONS}

def cai_weights_from_rscu(rscu):
    w = {}
    eps = 1e-6
    for aa, cods in AA2CODONS.items():
        max_r = max([rscu.get(c,0.0) for c in cods]) if len(cods)>0 else 0.0
        for c in cods:
            w[c] = (rscu.get(c,0.0)/max_r) if max_r>0 else eps
    for c in SENSE_CODONS:
        if w.get(c,0)<=0: w[c]=eps
    return w

def compute_cai_for_codon_sequence(codon_list, w_table):
    logs=[]
    for cod in codon_list:
        if cod in w_table:
            val=w_table[cod]
            if val<=0: val=1e-6
            logs.append(math.log(val))
    return math.exp(sum(logs)/len(logs)) if len(logs)>0 else np.nan

def compute_enc_from_codon_counts(gene_counts):
    Fk_lists={2:[],3:[],4:[],6:[]}
    for k in [2,3,4,6]:
        for aa in DEGENERACY_GROUPS.get(k,[]):
            cods = AA2CODONS.get(aa,[])
            counts=[gene_counts.get(c,0) for c in cods]
            n=sum(counts)
            if n<=1: continue
            sum_sq=sum(x**2 for x in counts)
            denom=n*(n-1)
            if denom==0: continue
            F=(sum_sq-n)/denom
            F=max(0.0,min(1.0,F))
            Fk_lists[k].append(F)
    F_means={}
    for k, vals in Fk_lists.items():
        if len(vals)>0:
            avg=sum(vals)/len(vals)
            F_means[k]=avg if avg>0 else 1e-6
    if len(F_means)==0: return np.nan
    ENC=2.0
    ENC += 9.0/F_means[2] if 2 in F_means else 9.0
    ENC += 1.0/F_means[3] if 3 in F_means else 0.0
    ENC += 5.0/F_means[4] if 4 in F_means else 0.0
    ENC += 3.0/F_means[6] if 6 in F_means else 0.0
    ENC=max(20.0,min(61.0,ENC))
    return ENC

def plot_heatmap(df, outfile, title, cmap='vlag', annot_fmt=".3f"):
    plt.figure(figsize=(8, max(10, 0.35*df.shape[1]+6)))
    sns.set(font_scale=0.55)
    
    cg = sns.clustermap(df.T, cmap=cmap, metric='euclidean', method='average',
                        annot=True, fmt=annot_fmt, linewidths=0.3,
                        annot_kws={"size":9}, cbar_kws={"label": "Value"},
                        dendrogram_ratio=(0.2, 0.2))
    
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_yticklabels(), fontsize=9)
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xticklabels(), fontsize=11, rotation=45, ha='right')
    
    cbar = cg.ax_heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label("Value", fontsize=9)
    
    cg.fig.suptitle(title, y=1.02, fontsize=14)
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

def plot_pca_with_legend(df, outfile, title):
    X=df.values
    X=np.nan_to_num(X)
    pca=PCA(n_components=2)
    pcs=pca.fit_transform(X)
    pc_df=pd.DataFrame(pcs, index=df.index, columns=["PC1","PC2"])
    
    plt.figure(figsize=(8,6))
    genomes=df.index.tolist()
    palette=sns.color_palette("tab10", n_colors=max(10,len(genomes)))
    color_map={genomes[i]:palette[i%len(palette)] for i in range(len(genomes))}
    
    for g in genomes:
        plt.scatter(pc_df.loc[g,"PC1"], pc_df.loc[g,"PC2"], color=color_map[g], s=100)
        
    handles=[plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=color_map[g],
                        markersize=10, label=g) for g in genomes]
                        
    plt.legend(handles=handles, title="Genome", bbox_to_anchor=(1.02,1), loc='upper left')
    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    return pc_df

# --- Main Execution ---
if __name__ == "__main__":
    if len(sys.argv)<2:
        print("Usage: python codon_pipeline.py /path/to/cds_folder (use . for current dir)")
        sys.exit(1)
        
    folder=sys.argv[1] if sys.argv[1]!="." else os.getcwd()
    files=[f for f in os.listdir(folder) if f.lower().endswith((".fa",".fasta",".fna"))]
    if len(files)==0:
        print("No fasta files found in folder:", folder)
        sys.exit(1)

    codon_freq_list=[]
    rscu_list=[]
    meta_list=[]
    per_gene_records=[]
    data_quality_lines=[]

    for f in files:
        genome_name=os.path.splitext(f)[0]
        path=os.path.join(folder,f)
        total_counts=Counter()
        gc_list=[]
        gc3_list=[]
        total_codons=0
        per_gene_for_genome=[]
        
        for rec in SeqIO.parse(path,"fasta"):
            seq=str(rec.seq).upper()
            if len(seq)<9:
                gc,gc3=gc_and_gc3(seq)
                gc_list.append(gc)
                gc3_list.append(gc3)
                per_gene_for_genome.append({"genome":genome_name,"gene":rec.id,"GC":gc,"GC3":gc3,"codon_counts":Counter(),"codon_list":[]})
                data_quality_lines.append(f"Short sequence (<9bp): genome {genome_name}, gene {rec.id}")
                continue
                
            gc,gc3=gc_and_gc3(seq)
            gc_list.append(gc)
            gc3_list.append(gc3)
            cc=count_codons(seq)
            total_counts.update(cc)
            total_codons+=sum(cc.values())
            
            if sum(cc.values())==0:
                data_quality_lines.append(f"No codons found: genome {genome_name}, gene {rec.id}")
                
            codon_list=[seq[i:i+3] for i in range(0,len(seq)-(len(seq)%3),3) if seq[i:i+3] in SENSE_CODONS]
            per_gene_for_genome.append({"genome":genome_name,"gene":rec.id,"GC":gc,"GC3":gc3,"codon_counts":cc,"codon_list":codon_list})
            
        if len(per_gene_for_genome)<10:
            data_quality_lines.append(f"Few genes (<10) for genome {genome_name}")

        freq=compute_codon_freqs_from_counts(total_counts)
        codon_freq_list.append(pd.Series(freq,name=genome_name))
        rscu=compute_rscu_from_counts(total_counts)
        rscu_list.append(pd.Series(rscu,name=genome_name))
        meta_list.append({"genome":genome_name,"mean_GC":np.nanmean(gc_list),"mean_GC3":np.nanmean(gc3_list),
                          "total_codons":int(total_codons)})

        cai_w=cai_weights_from_rscu(rscu)
        for g in per_gene_for_genome:
            gene_counts=g["codon_counts"]
            codon_list=g["codon_list"]
            gc=g["GC"]
            gc3=g["GC3"]
            mean_rscu=float(np.nanmean([rscu.get(c,0.0) for c in codon_list])) if len(codon_list)>0 else np.nan
            enc=compute_enc_from_codon_counts(gene_counts)
            cai=compute_cai_for_codon_sequence(codon_list,cai_w)
            per_gene_records.append({"genome":genome_name,"gene":g["gene"],"GC":gc,"GC3":gc3,
                                     "mean_RSCU":mean_rscu,"ENC":enc,"CAI":cai,"codon_count_total":sum(gene_counts.values())})
                                     
        genome_mean_rscu=float(np.nanmean([rscu.get(c,0.0) for c in SENSE_CODONS]))
        if genome_mean_rscu<0.5 or genome_mean_rscu>2.0:
            data_quality_lines.append(f"Abnormal genome mean RSCU: {genome_name} = {genome_mean_rscu:.3f}")

    # --- Build DataFrames ---
    codon_freq_df=pd.DataFrame(codon_freq_list).fillna(0)
    rscu_df=pd.DataFrame(rscu_list).fillna(0)
    meta_df=pd.DataFrame(meta_list).set_index("genome")
    per_gene_df=pd.DataFrame(per_gene_records)

    # --- Save CSVs/XLSX ---
    codon_freq_df.to_csv("codon_frequencies.csv", float_format="%.4f")
    rscu_df.to_csv("rscu_matrix.csv", float_format="%.4f")
    meta_df.to_csv("genome_meta.csv", float_format="%.4f")
    per_gene_df.to_csv("per_gene_metrics.csv", index=False, float_format="%.4f")
    try:
        codon_freq_df.to_excel("codon_frequencies.xlsx", float_format="%.4f", engine='openpyxl')
        rscu_df.to_excel("rscu_matrix.xlsx", float_format="%.4f", engine='openpyxl')
        meta_df.to_excel("genome_meta.xlsx", float_format="%.4f", engine='openpyxl')
    except Exception as e:
        print("Warning: could not write .xlsx files. Error:", e)

    # --- Heatmaps ---
    plot_heatmap(codon_freq_df,"codon_freq_heatmap.png","Codon Frequencies (genome-wise)")
    plot_heatmap(rscu_df,"rscu_heatmap.png","RSCU Values (genome-wise)")

    # --- PCA plots ---
    plot_pca_with_legend(codon_freq_df,"codon_freq_pca.png","PCA of Codon Frequencies")
    plot_pca_with_legend(rscu_df,"rscu_pca.png","PCA of RSCU Values")

    # --- GC3 vs mean-RSCU regression ---
    plt.figure(figsize=(10,8))
    sns.set(style="whitegrid", font_scale=0.9)
    genomes_sorted=sorted(per_gene_df['genome'].unique())
    palette=sns.color_palette("tab10", n_colors=max(10,len(genomes_sorted)))
    color_map={genomes_sorted[i]: palette[i%len(palette)] for i in range(len(genomes_sorted))}
    legend_handles=[]
    
    for g in genomes_sorted:
        sub=per_gene_df[per_gene_df['genome']==g].dropna(subset=['GC3','mean_RSCU'])
        if len(sub)==0: continue
        plt.scatter(sub['GC3'],sub['mean_RSCU'],s=10,alpha=0.6,color=color_map[g])
        if len(sub)>=5:
            slope,intercept,r_val,p_val,stderr=linregress(sub['GC3'],sub['mean_RSCU'])
            xs=np.linspace(sub['GC3'].min(),sub['GC3'].max(),100)
            ys=intercept+slope*xs
            plt.plot(xs,ys,color=color_map[g],linewidth=1.5)
        legend_handles.append(plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=color_map[g],
                                         markersize=6,linestyle='-',label=g))
                                         
    all_valid=per_gene_df.dropna(subset=['GC3','mean_RSCU'])
    if len(all_valid)>=5:
        s,i,r,p,se=linregress(all_valid['GC3'],all_valid['mean_RSCU'])
        xs=np.linspace(all_valid['GC3'].min(),all_valid['GC3'].max(),100)
        plt.plot(xs,i+s*xs,color='k',linewidth=2.0,label='Overall fit')
        
    plt.xlabel("GC3 (fraction)")
    plt.ylabel("Mean RSCU (per gene)")
    plt.title("GC3 vs Mean RSCU per gene (colored by genome)")
    plt.legend(handles=legend_handles,title="Genome",bbox_to_anchor=(1.02,1),loc='upper left')
    plt.tight_layout()
    plt.savefig("gc3_vs_meanrscu_regression.png",dpi=300,bbox_inches='tight')
    plt.close()

    # --- ENC vs GC3 Regression ---
    plt.figure(figsize=(10, 8))
    sns.set(style="ticks", font_scale=0.9)
    legend_handles_enc=[]
    enc_gc3_report_lines=[]
    
    for g in genomes_sorted:
        sub = per_gene_df[per_gene_df['genome']==g].dropna(subset=['GC3', 'ENC'])
        if len(sub) == 0: continue
        plt.scatter(sub['GC3'], sub['ENC'], s=10, alpha=0.6, color=color_map[g])
        if len(sub) >= 5:
            slope, intercept, r_val, p_val, stderr = linregress(sub['GC3'], sub['ENC'])
            xs = np.linspace(sub['GC3'].min(), sub['GC3'].max(), 100)
            ys = intercept + slope * xs
            plt.plot(xs, ys, color=color_map[g], linewidth=1.5)
            enc_gc3_report_lines.append(f"{g}\tSlope={slope:.4f}\tIntercept={intercept:.4f}\tR²={r_val**2:.4f}")
        legend_handles_enc.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[g],
                                             markersize=6, linestyle='-', label=g))
                                             
    plt.xlabel("GC3 (fraction)")
    plt.ylabel("ENC")
    plt.title("Wright's ENC vs GC3 plot (per gene, colored by genome)")
    plt.legend(handles=legend_handles_enc, title="Genome", bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("enc_vs_gc3_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

    with open("enc_gc3_regression_report.txt", "w") as fh:
        if len(enc_gc3_report_lines)==0:
            fh.write("No valid regressions computed.\n")
        else:
            fh.write("\n".join(enc_gc3_report_lines)+"\n")

    # --- Summaries ---
    summary_rows=[]
    for g in genomes_sorted:
        sub=per_gene_df[per_gene_df['genome']==g]
        if len(sub)==0: continue
        sub_valid=sub.dropna(subset=['GC3','mean_RSCU'])
        r2_gc3_rscu=linregress(sub_valid['GC3'],sub_valid['mean_RSCU']).rvalue**2 if len(sub_valid)>=5 else np.nan
        sub_valid_enc=sub.dropna(subset=['GC3','ENC'])
        r2_enc=linregress(sub_valid_enc['GC3'],sub_valid_enc['ENC']).rvalue**2 if len(sub_valid_enc)>=5 else np.nan
        
        summary_rows.append({
            "genome":g,
            "num_genes":len(sub),
            "total_codons":sub['codon_count_total'].sum(),
            "mean_GC":sub['GC'].mean(),
            "mean_GC3":sub['GC3'].mean(),
            "mean_RSCU":sub['mean_RSCU'].mean(),
            "ENC_mean":sub['ENC'].mean(),
            "CAI_mean":sub['CAI'].mean(),
            "R2_GC3_vs_meanRSCU":r2_gc3_rscu,
            "R2_ENC_vs_GC3": r2_enc,
            "var_meanRSCU":sub['mean_RSCU'].var()
        })
        
    summary_df=pd.DataFrame(summary_rows)
    summary_df.to_csv("codon_summary.csv",index=False,float_format="%.4f")

    # --- QC log ---
    if len(data_quality_lines)>0:
        with open("data_quality_warnings.txt","w") as fh:
            fh.write("\n".join(data_quality_lines)+"\n")