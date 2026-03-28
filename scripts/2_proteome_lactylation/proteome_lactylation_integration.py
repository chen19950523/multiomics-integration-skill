#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework 2: Proteome + Lactylation Integration Analysis
蛋白质组 + 蛋白乳酸化 联合分析流程

分析模块:
1. 乳酸化位点鉴定 + 全局蛋白定量关联
2. 乳酸化与蛋白表达的关系分类 (Type A-G)
3. 乳酸化位点motif分析
4. 功能富集 (代谢通路、线粒体功能)

Author: Bio-Design Team
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact, chi2_contingency
import re
import warnings
warnings.filterwarnings('ignore')

# ============== 配置参数 ==============
CONFIG = {
    'fc_threshold': 1.5,           # 差异倍数阈值
    'pvalue_threshold': 0.05,      # P值阈值
    'lactylation_site_fc': 1.5,    # 乳酸化位点差异阈值
    'motif_window': 7,            # Motif分析窗口大小
    'species': 'human',
    'min_peptides': 2,             # 最少肽段数
    ' localization_threshold': 0.75 # 定位可信度阈值
}

# ============== 数据加载 ==============
def load_lactylation_data():
    """
    加载乳酸化组学数据
    输入数据格式:
    - lactylation_sites: 乳酸化位点数据 (蛋白ID, 位点, 修饰强度)
    - proteome_data: 全局蛋白定量数据
    - sample_info: 样本分组信息
    """
    lactylation_sites = pd.read_csv('input/lactylation_sites.csv', index_col=0)
    proteome_data = pd.read_csv('input/proteome_quantification.csv', index_col=0)
    sample_info = pd.read_csv('input/sample_info.csv')
    
    return lactylation_sites, proteome_data, sample_info


# ============== 1. 乳酸化位点鉴定 + 蛋白定量关联 ==============
def lactylation_proteome_association(lactylation_sites, proteome_data):
    """
    乳酸化位点与全局蛋白定量关联分析
    """
    # 提取蛋白ID (从乳酸化位点数据)
    lactylation_sites['ProteinID'] = lactylation_sites.index  # 或从Site列提取
    
    # 找到共同的蛋白
    common_proteins = list(set(lactylation_sites.index) & set(proteome_data.index))
    print(f"共同蛋白数量: {len(common_proteins)}")
    
    # 创建关联表
    association_results = []
    
    for protein in common_proteins:
        # 获取该蛋白的所有乳酸化位点
        protein_sites = lactylation_sites[lactylation_sites.index == protein]
        n_sites = len(protein_sites)
        
        # 获取蛋白表达量
        protein_expr = proteome_data.loc[protein]
        
        # 获取乳酸化强度 (汇总所有位点)
        site_intensities = protein_sites.filter(like='Intensity').mean(axis=0) if n_sites > 0 else 0
        
        association_results.append({
            'ProteinID': protein,
            'ProteinExpression': protein_expr.mean(),
            'N_Lactylation_Sites': n_sites,
            'Total_Lactylation_Intensity': site_intensities.mean() if isinstance(site_intensities, pd.Series) else site_intensities,
            'Mean_Site_Intensity': site_intensities.mean() if isinstance(site_intensities, pd.Series) and n_sites > 0 else 0
        })
    
    assoc_df = pd.DataFrame(association_results)
    assoc_df['Lactylation_Per_Protein'] = assoc_df['Total_Lactylation_Intensity'] / (assoc_df['ProteinExpression'] + 1)
    
    assoc_df.to_csv('results/1_lactylation_proteome_association.csv', index=False)
    
    print(f"=== 乳酸化-蛋白关联 ===")
    print(f"有乳酸化位点的蛋白: {len(assoc_df[assoc_df['N_Lactylation_Sites'] > 0])}")
    print(f"平均位点数: {assoc_df['N_Lactylation_Sites'].mean():.2f}")
    
    return assoc_df


# ============== 2. 乳酸化与蛋白表达关系分类 (Type A-G) ==============
def lactylation_expression_typing(lactylation_sites, proteome_data, deg_data=None):
    """
    乳酸化与蛋白表达的关系分类
    
    Type A: 蛋白上调 + 乳酸化上调 (一致激活)
    Type B: 蛋白下调 + 乳酸化下调 (一致抑制)
    Type C: 蛋白上调 + 乳酸化下调 (蛋白激活但修饰抑制)
    Type D: 蛋白下调 + 乳酸化上调 (蛋白抑制但修饰激活)
    Type E: 仅蛋白变化 (无乳酸化变化)
    Type F: 仅乳酸化变化 (无蛋白变化)
    Type G: 蛋白和乳酸化都无显著变化
    
    修改版 - 基于位点水平分析:
    - Site-level changes vs Protein-level changes
    """
    # 计算蛋白和乳酸化位点的变化
    results = []
    
    # 获取所有蛋白
    all_proteins = set(lactylation_sites.index) | set(proteome_data.index)
    
    for protein in all_proteins:
        # 蛋白表达信息
        if protein in proteome_data.index:
            prot_expr = proteome_data.loc[protein]
            prot_change = prot_expr.get('log2FC', 0)
            prot_significant = prot_expr.get('PValue', 1) < 0.05
        else:
            prot_expr = None
            prot_change = 0
            prot_significant = False
        
        # 乳酸化位点信息
        if protein in lactylation_sites.index:
            sites = lactylation_sites.loc[protein]
            
            # 支持多行(多站点)或单行(汇总)
            if isinstance(sites, pd.DataFrame):
                site_changes = sites['log2FC'].values if 'log2FC' in sites.columns else [0]
                site_ps = sites['PValue'].values if 'PValue' in sites.columns else [1]
                n_sites_total = len(sites)
                n_sites_sig = sum(1 for p in site_ps if p < 0.05)
                mean_site_change = np.mean(site_changes)
            else:
                site_changes = [sites.get('log2FC', 0)]
                site_ps = [sites.get('PValue', 1)]
                n_sites_total = 1
                n_sites_sig = 1 if site_ps[0] < 0.05 else 0
                mean_site_change = site_changes[0]
        else:
            n_sites_total = 0
            n_sites_sig = 0
            mean_site_change = 0
        
        # 分类判断
        if prot_significant and n_sites_sig > 0:
            if prot_change > 0 and mean_site_change > 0:
                lactype = 'A'  # 蛋白Up + Lactylation Up
            elif prot_change < 0 and mean_site_change < 0:
                lactype = 'B'  # 蛋白Down + Lactylation Down
            elif prot_change > 0 and mean_site_change < 0:
                lactype = 'C'  # 蛋白Up + Lactylation Down
            elif prot_change < 0 and mean_site_change > 0:
                lactype = 'D'  # 蛋白Down + Lactylation Up
            else:
                lactype = 'G'
        elif prot_significant:
            lactype = 'E'  # 仅蛋白变化
        elif n_sites_sig > 0:
            lactype = 'F'  # 仅乳酸化变化
        else:
            lactype = 'G'  # 都无显著变化
        
        results.append({
            'ProteinID': protein,
            'Protein_log2FC': prot_change,
            'Protein_Significant': prot_significant,
            'N_Lactylation_Sites': n_sites_total,
            'N_Significant_Sites': n_sites_sig,
            'Lactylation_log2FC': mean_site_change,
            'Lactylation_Type': lactype
        })
    
    type_df = pd.DataFrame(results)
    type_df.to_csv('results/2_lactylation_expression_typing.csv', index=False)
    
    # 统计各类型数量
    print(f"\n=== 乳酸化-表达关系分类 ===")
    type_counts = type_df['Lactylation_Type'].value_counts()
    for t, count in type_counts.items():
        pct = count / len(type_df) * 100
        print(f"Type {t}: {count} ({pct:.1f}%)")
    
    return type_df


def plot_lactylation_typing(type_df):
    """可视化分类结果"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 1. Type分布柱状图
    type_counts = type_df['Lactylation_Type'].value_counts().sort_index()
    colors = {
        'A': '#e74c3c', 'B': '#3498db', 'C': '#f39c12',
        'D': '#9b59b6', 'E': '#1abc9c', 'F': '#34495e', 'G': '#95a5a6'
    }
    
    bars = axes[0].bar(type_counts.index, type_counts.values,
                       color=[colors.get(t, '#95a5a6') for t in type_counts.index])
    axes[0].set_xlabel('Lactylation Type')
    axes[0].set_ylabel('Number of Proteins')
    axes[0].set_title('Distribution of Lactylation-Expression Relationship Types')
    
    # 添加数值标签
    for bar, count in zip(bars, type_counts.values):
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                    str(count), ha='center', fontsize=10)
    
    # 2. 蛋白变化 vs 乳酸化变化散点图
    valid_df = type_df[(type_df['Protein_log2FC'] != 0) | (type_df['Lactylation_log2FC'] != 0)]
    
    type_colors = valid_df['Lactylation_Type'].map(colors)
    axes[1].scatter(valid_df['Protein_log2FC'], valid_df['Lactylation_log2FC'],
                   c=type_colors, alpha=0.6, s=30)
    axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    axes[1].axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('Protein log2FC')
    axes[1].set_ylabel('Lactylation log2FC')
    axes[1].set_title('Protein vs Lactylation Changes')
    
    plt.tight_layout()
    plt.savefig('results/2_lactylation_typing.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/2_lactylation_typing.png'


# ============== 3. 乳酸化位点Motif分析 ==============
def motif_analysis(lactylation_sites, window=7):
    """
    乳酸化位点motif分析
    提取位点周围氨基酸序列，分析共有motif
    """
    from collections import Counter
    
    # 氨基酸单字母编码
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    # 提取序列窗口
    sequences = []
    for _, row in lactylation_sites.iterrows():
        if 'Sequence' in row and 'Position' in row:
            seq = row['Sequence']
            pos = int(row['Position']) - 1  # 转为0-index
            
            if isinstance(seq, str) and 0 <= pos < len(seq):
                # 提取窗口序列
                start = max(0, pos - window)
                end = min(len(seq), pos + window + 1)
                window_seq = seq[start:end]
                
                # 标记中心K
                center_k = f"[{window_seq[pos-start]}]"
                sequences.append({
                    'ProteinID': row.name,
                    'Site': f"{row.get('GeneSymbol', row.name)}-{pos+1}",
                    'Window_Sequence': window_seq,
                    'Center_Marked': center_k
                })
    
    seq_df = pd.DataFrame(sequences)
    seq_df.to_csv('results/3_motif_sequences.csv', index=False)
    
    # Position frequency matrix (PFM)
    if len(seq_df) > 0:
        pfm = {aa: [0] * (window * 2) for aa in amino_acids}
        
        for seq in seq_df['Window_Sequence']:
            for i, aa in enumerate(seq):
                if aa in amino_acids:
                    pfm[aa][i] += 1
        
        pfm_df = pd.DataFrame(pfm).T
        pfm_df.columns = [f'pos_{i-window}' for i in range(window * 2)]
        pfm_df.to_csv('results/3_position_frequency_matrix.csv')
        
        # 序列logos可视化
        plot_sequence_logo(pfm_df, window)
    
    print(f"=== Motif分析 ===")
    print(f"分析位点数: {len(seq_df)}")
    
    return seq_df, pfm_df


def plot_sequence_logo(pfm_df, window):
    """绘制序列logos (简化版)"""
    fig, ax = plt.subplots(figsize=(14, 4))
    
    # 转换为位置频率
    height = 10
    x_positions = np.arange(pfm_df.shape[1])
    
    # 颜色分类
    aa_colors = {
        'K': '#e74c3c', 'R': '#e74c3c',  # 碱性 - 红色
        'D': '#3498db', 'E': '#3498db',  # 酸性 - 蓝色
        'S': '#2ecc71', 'T': '#2ecc71',  # 极性 - 绿色
        'A': '#95a5a6', 'V': '#95a5a6', 'L': '#95a5a6', 'I': '#95a5a6',  # 疏水 - 灰色
    }
    
    bottom = np.zeros(pfm_df.shape[1])
    
    for aa in pfm_df.index:
        heights = pfm_df.loc[aa].values / pfm_df.sum(axis=0).values * height
        color = aa_colors.get(aa, '#95a5a6')
        ax.bar(x_positions, heights, bottom=bottom, label=aa if aa == 'K' else '',
              color=color, width=0.8, alpha=0.8)
        bottom += heights
    
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'pos_{i}' for i in range(-window, window)])
    ax.set_ylabel('Frequency')
    ax.set_title('Lactylation Site Motif Logo')
    ax.axvline(x=window, color='red', linestyle='--', alpha=0.5, label='Center K')
    
    plt.tight_layout()
    plt.savefig('results/3_motif_logo.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/3_motif_logo.png'


# ============== 4. 功能富集分析 (代谢通路、线粒体功能) ==============
def lactylation_functional_enrichment(gene_list, output_prefix='lactylation'):
    """
    乳酸化蛋白的功能富集分析
    特别关注:
    - 代谢通路
    - 线粒体功能
    - 表观遗传调控
    """
    try:
        from clusterProfiler import enrichGO, enrichKEGG, enricher
        import org.Hs.eg.db as orgdb
    except ImportError:
        print("请安装: pip install clusterProfiler org.Hs.eg.db")
        return None
    
    # GO富集
    go_results = enrichGO(
        gene=gene_list,
        OrgDb=orgdb,
        keyType='SYMBOL',
        ont='ALL',
        pvalueCutoff=0.05,
        qvalueCutoff=0.05
    )
    
    go_df = pd.DataFrame(go_results)
    go_df.to_csv(f'results/4_{output_prefix}_GO_enrichment.csv', index=False)
    
    # KEGG富集
    kegg_results = enrichKEGG(
        gene=gene_list,
        organism='hsa',
        pvalueCutoff=0.05,
        qvalueCutoff=0.05
    )
    
    kegg_df = pd.DataFrame(kegg_results)
    kegg_df.to_csv(f'results/4_{output_prefix}_KEGG_enrichment.csv', index=False)
    
    # 代谢通路和线粒体相关关键词筛选
    metabolic_keywords = ['glycolysis', ' TCA ', 'citrate', 'acetyl', 'metabol', 
                          'oxidative', 'pyruvate', 'fatty acid', 'lipid']
    mitochondrial_keywords = ['mitochondri', 'mitochondrial', 'electron transport',
                              'respiratory', 'ATP synthesis', 'NADH', 'NAD+']
    
    metabolic_genes = []
    mitochondrial_genes = []
    
    for gene in gene_list:
        # 简化的分类方法 - 实际应使用更精确的分类
        metabolic_genes.append(gene)
        mitochondrial_genes.append(gene)
    
    return {
        'go': go_df,
        'kegg': kegg_df,
        'metabolic_genes': metabolic_genes,
        'mitochondrial_genes': mitochondrial_genes
    }


def plot_lactylation_enrichment(enrich_result, type_df):
    """可视化乳酸化功能富集结果"""
    if enrich_result is None:
        return None
    
    go_df = enrich_result['go']
    kegg_df = enrich_result['kegg']
    
    if len(go_df) == 0:
        return None
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. GO BP 富集 (top 20)
    go_bp = go_df[go_df['ONTOLOGY'] == 'BP'].head(20)
    if len(go_bp) > 0:
        axes[0, 0].barh(range(len(go_bp)), -np.log10(go_bp['pvalue']))
        axes[0, 0].set_yticks(range(len(go_bp)))
        axes[0, 0].set_yticklabels([t[:40] for t in go_bp['Description']], fontsize=7)
        axes[0, 0].set_xlabel('-log10(P-value)')
        axes[0, 0].set_title('GO Biological Process Enrichment')
        axes[0, 0].invert_yaxis()
    
    # 2. KEGG通路富集
    kegg_top = kegg_df.head(20)
    if len(kegg_top) > 0:
        axes[0, 1].barh(range(len(kegg_top)), -np.log10(kegg_top['pvalue']))
        axes[0, 1].set_yticks(range(len(kegg_top)))
        axes[0, 1].set_yticklabels([t[:40] for t in kegg_top['Description']], fontsize=7)
        axes[0, 1].set_xlabel('-log10(P-value)')
        axes[0, 1].set_title('KEGG Pathway Enrichment')
        axes[0, 1].invert_yaxis()
    
    # 3. Type分布 (堆叠柱状图)
    type_counts = type_df['Lactylation_Type'].value_counts()
    axes[1, 0].pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%')
    axes[1, 0].set_title('Lactylation-Expression Type Distribution')
    
    # 4. 蛋白-位点关系
    # 统计每个蛋白的乳酸化位点数分布
    site_counts = type_df[type_df['N_Lactylation_Sites'] > 0]['N_Lactylation_Sites']
    axes[1, 1].hist(site_counts, bins=range(1, site_counts.max()+2), 
                    color='#3498db', alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('Number of Lactylation Sites per Protein')
    axes[1, 1].set_ylabel('Number of Proteins')
    axes[1, 1].set_title('Distribution of Lactylation Sites per Protein')
    
    plt.tight_layout()
    plt.savefig('results/4_lactylation_functional_enrichment.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/4_lactylation_functional_enrichment.png'


# ============== 主流程 ==============
def main():
    """主分析流程"""
    print("=" * 60)
    print("蛋白质组 + 乳酸化组 联合分析流程")
    print("=" * 60)
    
    # 1. 加载数据
    print("\n[1/4] 加载数据...")
    lactylation_sites, proteome_data, sample_info = load_lactylation_data()
    
    # 2. 乳酸化-蛋白关联分析
    print("\n[2/4] 乳酸化位点与蛋白表达关联...")
    assoc_df = lactylation_proteome_association(lactylation_sites, proteome_data)
    
    # 3. 乳酸化表达关系分类
    print("\n[3/4] 乳酸化-表达关系分类 (Type A-G)...")
    type_df = lactylation_expression_typing(lactylation_sites, proteome_data)
    plot_lactylation_typing(type_df)
    
    # 4. Motif分析
    print("\n[4/4] 乳酸化位点Motif分析...")
    seq_df, pfm_df = motif_analysis(lactylation_sites, window=CONFIG['motif_window'])
    
    # 5. 功能富集分析
    print("\n[Extra] 功能富集分析...")
    lactylated_genes = type_df[type_df['N_Lactylation_Sites'] > 0]['ProteinID'].tolist()
    enrich_result = lactylation_functional_enrichment(lactylated_genes)
    if enrich_result:
        plot_lactylation_enrichment(enrich_result, type_df)
    
    print("\n" + "=" * 60)
    print("分析完成! 结果保存在 results/ 目录")
    print("=" * 60)


if __name__ == '__main__':
    main()
