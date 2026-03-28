#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework 1: Transcriptome + Proteome Integration Analysis
转录组 + 蛋白质组 联合分析流程

分析模块:
1. 差异重叠分析 (DEGs ∩ DEPs)
2. mRNA-蛋白表达相关性
3. GO/KEGG联合富集
4. PPI网络构建
5. 翻译效率分析

Author: Bio-Design Team
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# ============== 配置参数 ==============
CONFIG = {
    'fc_threshold': 1.5,           # 差异倍数阈值
    'pvalue_threshold': 0.05,     # P值阈值
    'correlation_threshold': 0.6,  # 相关性阈值
    'ppi_score_threshold': 0.7,    # PPI相互作用置信度
    'top_n': 50,                   # 富集分析展示top N
    'species': 'human',            # human/mouse/rat
    'species_taxid': 9606,          # NCBI Taxonomy ID
}

# ============== 数据输入 ==============
def load_data():
    """
    加载转录组和蛋白质组数据
    输入格式要求:
    - rowname: gene symbol (转录组) / protein ID (蛋白组)
    - 列: 样本名
    - 必须有 group 列标识分组
    """
    # 示例数据结构，实际使用时替换为真实数据路径
    deg_data = pd.read_csv('input/DEGs.csv', index_col=0)
    dep_data = pd.read_csv('input/DEPs.csv', index_col=0)
    tpm_data = pd.read_csv('input/TPM.csv', index_col=0)  # 转录本表达量
    intensity_data = pd.read_csv('input/protein_intensity.csv', index_col=0)  # 蛋白表达量
    
    return deg_data, dep_data, tpm_data, intensity_data


# ============== 1. 差异重叠分析 ==============
def diff_overlap_analysis(deg_data, dep_data):
    """
    差异基因/蛋白重叠分析
    输出: Venn图、重叠基因列表、分类统计
    """
    degs = set(deg_data.index)
    deps = set(dep_data.index)
    
    # 计算重叠
    overlap = degs & deps
    deg_only = degs - deps
    dep_only = deps - degs
    
    print(f"=== 差异重叠分析 ===")
    print(f"DEGs数量: {len(degs)}")
    print(f"DEPs数量: {len(deps)}")
    print(f"重叠基因/蛋白: {len(overlap)}")
    print(f"仅DEGs: {len(deg_only)}")
    print(f"仅DEPs: {len(dep_only)}")
    
    # 创建重叠分类表
    overlap_df = pd.DataFrame({
        'GeneSymbol': list(overlap),
        'DEG_FC': [deg_data.loc[g, 'log2FC'] if 'log2FC' in deg_data.columns else np.nan for g in overlap],
        'DEP_FC': [dep_data.loc[g, 'log2FC'] if 'log2FC' in dep_data.columns else np.nan for g in overlap],
        'DEG_Direction': ['Up' if deg_data.loc[g, 'log2FC'] > 0 else 'Down' for g in overlap],
        'DEP_Direction': ['Up' if dep_data.loc[g, 'log2FC'] > 0 else 'Down' for g in overlap],
    })
    
    # 分类: 一致性方向
    overlap_df['Consistency'] = overlap_df.apply(
        lambda x: 'Consistent' if (x['DEG_Direction'] == x['DEP_Direction']) else 'Inconsistent', axis=1
    )
    
    # 保存结果
    overlap_df.to_csv('results/1_overlap_analysis.csv', index=False)
    
    return {
        'overlap_genes': overlap,
        'overlap_df': overlap_df,
        'deg_only': deg_only,
        'dep_only': dep_only
    }


def plot_overlap_venn(deg_data, dep_data, overlap_result):
    """绘制Venn图"""
    from matplotlib_venn import venn2, venn3
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 2-way Venn
    venn2([set(deg_data.index), set(dep_data.index)], 
          set_labels=['DEGs', 'DEPs'],
          ax=axes[0])
    axes[0].set_title('DEGs ∩ DEPs Overlap', fontsize=14)
    
    # 重叠基因方向分类饼图
    overlap_df = overlap_result['overlap_df']
    direction_counts = overlap_df.groupby(['DEG_Direction', 'DEP_Direction']).size()
    
    categories = ['Consistent Up', 'Consistent Down', 'Inconsistent']
    values = [
        len(overlap_df[(overlap_df['DEG_Direction']=='Up') & (overlap_df['DEP_Direction']=='Up')]),
        len(overlap_df[(overlap_df['DEG_Direction']=='Down') & (overlap_df['DEP_Direction']=='Down')]),
        len(overlap_df[overlap_df['Consistency']=='Inconsistent'])
    ]
    
    colors = ['#e74c3c', '#3498db', '#95a5a6']
    axes[1].pie(values, labels=categories, autopct='%1.1f%%', colors=colors)
    axes[1].set_title('Overlap Direction Consistency', fontsize=14)
    
    plt.tight_layout()
    plt.savefig('results/1_overlap_venn.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/1_overlap_venn.png'


# ============== 2. mRNA-蛋白表达相关性 ==============
def mrna_protein_correlation(tpm_data, intensity_data, min_samples=3):
    """
    mRNA-蛋白表达相关性分析
    支持Pearson、Spearman两种相关性计算
    """
    # 找到共同的基因/蛋白
    common_genes = list(set(tpm_data.index) & set(intensity_data.index))
    print(f"共同基因数量: {len(common_genes)}")
    
    tpm_common = tpm_data.loc[common_genes]
    intensity_common = intensity_data.loc[common_genes]
    
    # 计算相关性
    correlation_results = []
    
    for gene in common_genes:
        tpm_vals = tpm_common.loc[gene].values
        intensity_vals = intensity_common.loc[gene].values
        
        # 过滤低表达
        if np.mean(tpm_vals) < 1 or np.mean(intensity_vals) < np.median(intensity_common.values):
            continue
            
        # Pearson相关性
        try:
            r_p, p_p = pearsonr(tpm_vals, intensity_vals)
        except:
            r_p, p_p = np.nan, np.nan
            
        # Spearman相关性
        try:
            r_s, p_s = spearmanr(tpm_vals, intensity_vals)
        except:
            r_s, p_s = np.nan, np.nan
        
        correlation_results.append({
            'GeneSymbol': gene,
            'Pearson_R': r_p,
            'Pearson_P': p_p,
            'Spearman_R': r_s,
            'Spearman_P': p_s,
            'TPM_Mean': np.mean(tpm_vals),
            'Intensity_Mean': np.mean(intensity_vals)
        })
    
    corr_df = pd.DataFrame(correlation_results)
    corr_df['Pearson_FDR'] = stats.false_discovery_control(corr_df['Pearson_P'].fillna(1))
    corr_df['Spearman_FDR'] = stats.false_discovery_control(corr_df['Spearman_P'].fillna(1))
    
    # 分类
    corr_df['Correlation_Type'] = 'None'
    corr_df.loc[(corr_df['Pearson_R'] > 0.6) & (corr_df['Pearson_P'] < 0.05), 'Correlation_Type'] = 'Strong_Positive'
    corr_df.loc[(corr_df['Pearson_R'] > 0.3) & (corr_df['Pearson_P'] < 0.05), 'Correlation_Type'] = 'Moderate_Positive'
    corr_df.loc[(corr_df['Pearson_R'] < -0.6) & (corr_df['Pearson_P'] < 0.05), 'Correlation_Type'] = 'Strong_Negative'
    corr_df.loc[(corr_df['Pearson_R'] < -0.3) & (corr_df['Pearson_P'] < 0.05), 'Correlation_Type'] = 'Moderate_Negative'
    
    corr_df.to_csv('results/2_mRNA_protein_correlation.csv', index=False)
    
    print(f"=== mRNA-蛋白相关性 ===")
    print(corr_df['Correlation_Type'].value_counts())
    
    return corr_df


def plot_correlation_heatmap(corr_df, top_n=100):
    """绘制相关性热图"""
    # 按Pearson R值排序，取top N
    top_corr = corr_df.nlargest(top_n, 'Pearson_R')
    
    fig, ax = plt.subplots(figsize=(10, max(8, top_n/10)))
    
    # 颜色映射
    colors = ['#3498db' if r > 0 else '#e74c3c' for r in top_corr['Pearson_R']]
    
    ax.barh(range(len(top_corr)), top_corr['Pearson_R'], color=colors)
    ax.set_yticks(range(len(top_corr)))
    ax.set_yticklabels(top_corr['GeneSymbol'], fontsize=6)
    ax.set_xlabel('Pearson Correlation Coefficient')
    ax.set_title(f'Top {top_n} mRNA-Protein Correlations')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.axvline(x=0.6, color='green', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.axvline(x=-0.6, color='green', linestyle='--', linewidth=0.5, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('results/2_correlation_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/2_correlation_barplot.png'


# ============== 3. GO/KEGG联合富集分析 ==============
def enrichment_analysis(gene_list, gene_type='symbol', output_prefix='enrichment'):
    """
    GO/KEGG富集分析
    使用clusterProfiler进行富集分析
    """
    try:
        from clusterProfiler import enrichGO, enrichKEGG, enricher
        import org.Hs.eg.db as orgdb
    except ImportError:
        print("请安装: pip install clusterProfiler org.Hs.eg.db")
        return None
    
    # 设置物种数据库
    if CONFIG['species'] == 'human':
        organism_db = orgdb
        organism = 'hsa'
    elif CONFIG['species'] == 'mouse':
        organism_db = orgdb  # 需要org.Mm.eg.db
        organism = 'mmu'
    else:
        organism_db = orgdb
        organism = 'hsa'
    
    # GO富集分析
    go_results = enrichGO(
        gene=gene_list,
        OrgDb=organism_db,
        keyType='SYMBOL',
        ont='ALL',  # BP, MF, CC, ALL
        pvalueCutoff=0.05,
        qvalueCutoff=0.05
    )
    
    go_df = pd.DataFrame(go_results)
    go_df.to_csv(f'results/3_{output_prefix}_GO_enrichment.csv', index=False)
    
    # KEGG富集分析
    kegg_results = enrichKEGG(
        gene=gene_list,
        organism=organism,
        pvalueCutoff=0.05,
        qvalueCutoff=0.05
    )
    
    kegg_df = pd.DataFrame(kegg_results)
    kegg_df.to_csv(f'results/3_{output_prefix}_KEGG_enrichment.csv', index=False)
    
    return {'go': go_df, 'kegg': kegg_df}


def plot_enrichment_network(go_df, kegg_df, top_n=30):
    """绘制富集网络图"""
    # 合并GO和KEGG结果
    combined = pd.concat([
        go_df.head(top_n).assign(Category='GO'),
        kegg_df.head(top_n).assign(Category='KEGG')
    ])
    
    # 筛选显著条目
    sig_combined = combined[combined['p.adjust'] < 0.05].head(top_n)
    
    if len(sig_combined) == 0:
        print("没有显著的富集结果")
        return None
    
    # 创建网络图
    G = nx.DiGraph()
    
    # 添加节点
    for _, row in sig_combined.iterrows():
        G.add_node(
            row['Description'],
            category=row['Category'],
            pvalue=row['p.adjust'],
            count=row['Count']
        )
    
    # 布局
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # 绘图
    fig, ax = plt.subplots(figsize=(14, 10))
    
    colors = ['#e74c3c' if n == 'GO' else '#3498db' for n in [G.nodes[n].get('category') for n in G.nodes()]]
    
    nx.draw_networkx_nodes(G, pos, node_color=colors, alpha=0.7, node_size=[G.nodes[n]['count']*50 for n in G.nodes()])
    nx.draw_networkx_edges(G, pos, alpha=0.3, arrows=False)
    nx.draw_networkx_labels(G, pos, font_size=6)
    
    plt.title('GO/KEGG Enrichment Network')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/3_enrichment_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/3_enrichment_network.png'


# ============== 4. PPI网络构建 ==============
def build_ppi_network(gene_list, score_threshold=0.7):
    """
    构建PPI网络
    使用STRING数据库
    """
    try:
        import requests
    except ImportError:
        print("请安装requests库")
        return None
    
    # STRING API批量查询
    genes_str = '%0d'.join(gene_list[:min(len(gene_list), 200)])  # STRING限制200个
    
    url = f"https://string-db.org/api/tsv/network?identifiers={genes_str}&species={CONFIG['species_taxid']}&score_threshold={score_threshold}"
    
    try:
        response = requests.get(url, timeout=30)
        lines = response.text.strip().split('\n')
        
        ppi_edges = []
        for line in lines[1:]:  # 跳过header
            parts = line.split('\t')
            if len(parts) >= 2:
                ppi_edges.append({
                    'protein1': parts[0],
                    'protein2': parts[1],
                    'score': float(parts[2]) if len(parts) > 2 and parts[2] else 0
                })
        
        ppi_df = pd.DataFrame(ppi_edges)
        ppi_df.to_csv('results/4_PPI_network.csv', index=False)
        
        return ppi_df
    except Exception as e:
        print(f"PPI查询失败: {e}")
        return None


def visualize_ppi_network(ppi_df, highlight_genes=None):
    """可视化PPI网络"""
    if ppi_df is None or len(ppi_df) == 0:
        print("没有PPI数据")
        return None
    
    G = nx.Graph()
    
    # 添加边
    for _, row in ppi_df.iterrows():
        G.add_edge(row['protein1'], row['protein2'], weight=row['score'])
    
    # 计算网络统计
    degree_dict = dict(G.degree())
    
    # 筛选高度连接的hub蛋白
    hub_genes = [k for k, v in degree_dict.items() if v >= 5]
    
    print(f"=== PPI网络统计 ===")
    print(f"节点数: {G.number_of_nodes()}")
    print(f"边数: {G.number_of_edges()}")
    print(f"Hub蛋白 (>5连接): {len(hub_genes)}")
    
    # 绘图
    fig, ax = plt.subplots(figsize=(16, 12))
    
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    # 节点大小按degree
    node_sizes = [degree_dict[n] * 30 for n in G.nodes()]
    
    # 节点颜色
    highlight_genes = highlight_genes or []
    node_colors = ['#e74c3c' if n in highlight_genes else '#3498db' for n in G.nodes()]
    
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.7)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    nx.draw_networkx_labels(G, pos, font_size=5)
    
    plt.title('PPI Network (Red=Highlighted Genes)')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/4_PPI_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return hub_genes


# ============== 5. 翻译效率分析 ==============
def translation_efficiency_analysis(tpm_data, intensity_data, deg_data=None, dep_data=None):
    """
    翻译效率(TE)分析
    TE = (Protein Intensity) / (mRNA TPM)
    """
    # 找共同基因
    common_genes = list(set(tpm_data.index) & set(intensity_data.index))
    
    # 计算翻译效率
    te_results = []
    
    for gene in common_genes:
        tpm_vals = tpm_data.loc[gene].values
        intensity_vals = intensity_data.loc[gene].values
        
        # 避免除零
        tpm_mean = np.mean(tpm_vals) + 0.1
        intensity_mean = np.mean(intensity_vals) + 1
        
        te = intensity_mean / tpm_mean
        
        # 获取差异信息
        deg_fc = deg_data.loc[gene, 'log2FC'] if deg_data is not None and gene in deg_data.index else np.nan
        dep_fc = dep_data.loc[gene, 'log2FC'] if dep_data is not None and gene in dep_data.index else np.nan
        
        te_results.append({
            'GeneSymbol': gene,
            'TPM_Mean': tpm_mean,
            'Intensity_Mean': intensity_mean,
            'Translation_Efficiency': te,
            'log2TE': np.log2(te),
            'DEG_log2FC': deg_fc,
            'DEP_log2FC': dep_fc
        })
    
    te_df = pd.DataFrame(te_results)
    
    # 分类翻译效率变化
    # 高TE变化: 蛋白变化超出转录变化的部分
    te_df['TE_Category'] = 'Normal'
    te_df.loc[(te_df['DEP_log2FC'] - te_df['DEG_log2FC']) > 1, 'TE_Category'] = 'Enhanced_TE'  # 翻译增强
    te_df.loc[(te_df['DEP_log2FC'] - te_df['DEG_log2FC']) < -1, 'TE_Category'] = 'Reduced_TE'    # 翻译抑制
    
    te_df.to_csv('results/5_translation_efficiency.csv', index=False)
    
    print(f"=== 翻译效率分析 ===")
    print(te_df['TE_Category'].value_counts())
    
    return te_df


def plot_te_analysis(te_df):
    """可视化翻译效率分析"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. TE分布直方图
    axes[0, 0].hist(te_df['log2TE'].dropna(), bins=50, color='#3498db', alpha=0.7)
    axes[0, 0].set_xlabel('log2(Translation Efficiency)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Translation Efficiency Distribution')
    
    # 2. DEG vs DEP scatter
    valid_te = te_df.dropna(subset=['DEG_log2FC', 'DEP_log2FC'])
    
    colors = valid_te['TE_Category'].map({
        'Normal': '#95a5a6',
        'Enhanced_TE': '#e74c3c', 
        'Reduced_TE': '#3498db'
    })
    
    axes[0, 1].scatter(valid_te['DEG_log2FC'], valid_te['DEP_log2FC'], 
                       c=colors, alpha=0.5, s=20)
    axes[0, 1].plot([-6, 6], [-6, 6], 'k--', alpha=0.3, label='TE=1')
    axes[0, 1].set_xlabel('DEG log2FC (mRNA)')
    axes[0, 1].set_ylabel('DEP log2FC (Protein)')
    axes[0, 1].set_title('mRNA vs Protein Change')
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#95a5a6', label='Normal'),
        Patch(facecolor='#e74c3c', label='Enhanced TE'),
        Patch(facecolor='#3498db', label='Reduced TE')
    ]
    axes[0, 1].legend(handles=legend_elements)
    
    # 3. TE category统计
    category_counts = te_df['TE_Category'].value_counts()
    axes[1, 0].bar(category_counts.index, category_counts.values, 
                   color=['#95a5a6', '#e74c3c', '#3498db'])
    axes[1, 0].set_ylabel('Gene Count')
    axes[1, 0].set_title('Translation Efficiency Categories')
    
    # 4. Boxplot of TE by category
    te_valid = te_df[te_df['TE_Category'] != 'Normal']
    if len(te_valid) > 0:
        categories = te_valid['TE_Category'].unique()
        data_by_cat = [te_valid[te_valid['TE_Category'] == c]['log2TE'].values for c in categories]
        axes[1, 1].boxplot(data_by_cat, labels=categories)
        axes[1, 1].set_ylabel('log2(Translation Efficiency)')
        axes[1, 1].set_title('TE Distribution by Category')
    
    plt.tight_layout()
    plt.savefig('results/5_translation_efficiency.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/5_translation_efficiency.png'


# ============== 主流程 ==============
def main():
    """主分析流程"""
    print("=" * 60)
    print("转录组 + 蛋白质组 联合分析流程")
    print("=" * 60)
    
    # 1. 加载数据
    print("\n[1/5] 加载数据...")
    deg_data, dep_data, tpm_data, intensity_data = load_data()
    
    # 2. 差异重叠分析
    print("\n[2/5] 差异重叠分析...")
    overlap_result = diff_overlap_analysis(deg_data, dep_data)
    plot_overlap_venn(deg_data, dep_data, overlap_result)
    
    # 3. mRNA-蛋白相关性
    print("\n[3/5] mRNA-蛋白表达相关性...")
    corr_df = mrna_protein_correlation(tpm_data, intensity_data)
    plot_correlation_heatmap(corr_df)
    
    # 4. GO/KEGG富集分析 (对重叠基因进行)
    print("\n[4/5] GO/KEGG联合富集...")
    overlap_genes = list(overlap_result['overlap_genes'])
    if len(overlap_genes) > 0:
        enrich_result = enrichment_analysis(overlap_genes, output_prefix='overlap')
        if enrich_result:
            plot_enrichment_network(enrich_result['go'], enrich_result['kegg'])
    
    # 5. PPI网络构建
    print("\n[5/5] PPI网络构建...")
    ppi_df = build_ppi_network(overlap_genes)
    if ppi_df is not None:
        visualize_ppi_network(ppi_df, highlight_genes=overlap_genes[:20])
    
    # 6. 翻译效率分析
    print("\n[Extra] 翻译效率分析...")
    te_df = translation_efficiency_analysis(tpm_data, intensity_data, deg_data, dep_data)
    plot_te_analysis(te_df)
    
    print("\n" + "=" * 60)
    print("分析完成! 结果保存在 results/ 目录")
    print("=" * 60)


if __name__ == '__main__':
    main()
