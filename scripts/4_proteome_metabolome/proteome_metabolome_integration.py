#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework 4: Proteome + Metabolome Integration Analysis
蛋白质组 + 代谢组 联合分析流程

分析模块:
1. 相关性热图 + 典型相关分析 (CCA)
2. 通路联合富集
3. 蛋白-代谢物因果调控关系
4. 多组学biomarker发现

Author: Bio-Design Team
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# ============== 配置参数 ==============
CONFIG = {
    'correlation_threshold': 0.6,     # 相关性阈值
    'pvalue_threshold': 0.05,        # P值阈值
    'cca_components': 2,            # CCA组分数量
    'pathway_pvalue': 0.05,          # 通路显著性
    'biomarker_fdr': 0.05,          # Biomarker筛选FDR
    'min_overlap': 3,               # 通路最小重叠基因/代谢物
}

# ============== 数据加载 ==============
def load_metabolomics_data():
    """
    加载代谢组学数据
    """
    proteome_data = pd.read_csv('input/proteome_quantification.csv', index_col=0)
    metabolome_data = pd.read_csv('input/metabolome_quantification.csv', index_col=0)
    sample_info = pd.read_csv('input/sample_info.csv')
    
    return proteome_data, metabolome_data, sample_info


# ============== 1. 相关性热图 + CCA分析 ==============
def correlation_analysis(proteome_data, metabolome_data, method='spearman'):
    """
    蛋白-代谢物相关性分析
    """
    # 获取共同样本
    common_samples = list(set(proteome_data.columns) & set(metabolome_data.columns))
    print(f"共同样本数: {len(common_samples)}")
    
    # 子集化数据
    prot_common = proteome_data[common_samples]
    met_common = metabolome_data[common_samples]
    
    # 计算相关性矩阵
    correlation_results = []
    
    for prot in prot_common.index:
        for met in met_common.index:
            prot_vals = prot_common.loc[prot].values.astype(float)
            met_vals = met_common.loc[met].values.astype(float)
            
            # 过滤低变异
            if np.std(prot_vals) == 0 or np.std(met_vals) == 0:
                continue
            
            if method == 'spearman':
                r, p = spearmanr(prot_vals, met_vals)
            else:
                r, p = pearsonr(prot_vals, met_vals)
            
            if not np.isnan(r):
                correlation_results.append({
                    'Protein': prot,
                    'Metabolite': met,
                    'Correlation': r,
                    'P_Value': p,
                    'Abs_Correlation': abs(r)
                })
    
    corr_df = pd.DataFrame(correlation_results)
    corr_df['FDR'] = stats.false_discovery_control(corr_df['P_Value'])
    corr_df['Significance'] = corr_df['FDR'] < 0.05
    
    # 保存完整相关性矩阵
    corr_df.to_csv('results/1_protein_metabolite_correlation.csv', index=False)
    
    print(f"=== 相关性分析 ===")
    print(f"总配对数: {len(corr_df)}")
    print(f"显著相关: {corr_df['Significance'].sum()}")
    
    return corr_df


def plot_correlation_heatmap(corr_df, top_n=50):
    """绘制相关性热图"""
    # 筛选强相关
    sig_corr = corr_df[corr_df['Significance']].copy()
    
    if len(sig_corr) == 0:
        print("没有显著相关对")
        return None
    
    # 获取top N
    top_corr = sig_corr.nlargest(top_n, 'Abs_Correlation')
    
    # 创建相关性矩阵
    proteins = top_corr['Protein'].unique()
    metabolites = top_corr['Metabolite'].unique()
    
    corr_matrix = pd.DataFrame(index=proteins, columns=metabolites, dtype=float)
    
    for _, row in top_corr.iterrows():
        corr_matrix.loc[row['Protein'], row['Metabolite']] = row['Correlation']
    
    corr_matrix = corr_matrix.fillna(0)
    
    # 绘制热图
    fig, ax = plt.subplots(figsize=(max(12, len(metabolites)*0.8), max(8, len(proteins)*0.3)))
    
    sns.heatmap(corr_matrix, annot=False, cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, ax=ax,
                cbar_kws={'label': 'Correlation'})
    
    ax.set_title(f'Protein-Metabolite Correlation Heatmap (Top {top_n})')
    plt.xticks(rotation=45, ha='right', fontsize=6)
    plt.yticks(fontsize=6)
    
    plt.tight_layout()
    plt.savefig('results/1_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/1_correlation_heatmap.png'


def cca_analysis(proteome_data, metabolome_data, n_components=2):
    """
    典型相关分析 (CCA)
    找出蛋白组和代谢组之间相关性最大的线性组合
    """
    # 获取共同样本
    common_samples = list(set(proteome_data.columns) & set(metabolome_data.columns))
    
    X = proteome_data[common_samples].T.values  # 蛋白组 (样本 x 蛋白)
    Y = metabolome_data[common_samples].T.values  # 代谢组 (样本 x 代谢物)
    
    # 标准化
    X_scaled = StandardScaler().fit_transform(X)
    Y_scaled = StandardScaler().fit_transform(Y)
    
    # CCA
    cca = CCA(n_components=n_components)
    X_c, Y_c = cca.fit_transform(X_scaled, Y_scaled)
    
    # 结果
    sample_names = common_samples
    cca_scores = pd.DataFrame({
        'Sample': sample_names,
        f'CC1_X': X_c[:, 0],
        f'CC1_Y': Y_c[:, 0],
        f'CC2_X': X_c[:, 1] if n_components > 1 else None,
        f'CC2_Y': Y_c[:, 1] if n_components > 1 else None,
    })
    
    # 相关性载荷
    loadings_prot = pd.DataFrame(
        cca.x_loadings_,
        index=proteome_data.index,
        columns=[f'CC{i+1}' for i in range(n_components)]
    )
    
    loadings_met = pd.DataFrame(
        cca.y_loadings_,
        index=metabolome_data.index,
        columns=[f'CC{i+1}' for i in range(n_components)]
    )
    
    # 保存结果
    cca_scores.to_csv('results/1_CCA_scores.csv', index=False)
    loadings_prot.to_csv('results/1_CCA_protein_loadings.csv')
    loadings_met.to_csv('results/1_CCA_metabolite_loadings.csv')
    
    # CCA相关系数
    corr_cc = []
    for i in range(n_components):
        r, p = pearsonr(X_c[:, i], Y_c[:, i])
        corr_cc.append({'Component': f'CC{i+1}', 'Correlation': r, 'P_Value': p})
    
    corr_cc_df = pd.DataFrame(corr_cc)
    corr_cc_df.to_csv('results/1_CCA_correlations.csv', index=False)
    
    print(f"\n=== CCA分析 ===")
    print(f"CC1 vs CC1相关性: {corr_cc_df.iloc[0]['Correlation']:.3f}")
    print(f"CC2 vs CC2相关性: {corr_cc_df.iloc[1]['Correlation']:.3f}")
    
    return {
        'scores': cca_scores,
        'loadings_prot': loadings_prot,
        'loadings_met': loadings_met,
        'corr': corr_cc_df
    }


def plot_cca_results(cca_result, sample_info=None):
    """可视化CCA结果"""
    scores = cca_result['scores']
    loadings_prot = cca_result['loadings_prot']
    loadings_met = cca_result['loadings_met']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. CCA score plot
    ax = axes[0]
    
    if sample_info is not None:
        # 按分组着色
        scores_merged = scores.merge(sample_info, on='Sample')
        groups = scores_merged['Group'].unique()
        
        for g in groups:
            group_data = scores_merged[scores_merged['Group'] == g]
            ax.scatter(group_data['CC1_X'], group_data['CC1_Y'], 
                      label=g, s=80, alpha=0.7)
    else:
        ax.scatter(scores['CC1_X'], scores['CC1_Y'], s=80, alpha=0.7)
    
    ax.set_xlabel('CC1 (Protein)')
    ax.set_ylabel('CC1 (Metabolite)')
    ax.set_title('CCA Score Plot')
    ax.legend()
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # 2. Loading plot
    ax = axes[1]
    
    # 获取top loadings
    top_prot = loadings_prot.nlargest(20, 'CC1')
    top_met = loadings_met.nlargest(20, 'CC1')
    
    ax.scatter(top_prot['CC1'], [0]*len(top_prot), c='blue', s=100, marker='s', label='Proteins')
    ax.scatter([0]*len(top_met), top_met['CC1'], c='red', s=100, marker='o', label='Metabolites')
    
    # 添加标签
    for prot in top_prot.index[:10]:
        ax.annotate(prot[:10], (top_prot.loc[prot, 'CC1'], 0), fontsize=6)
    for met in top_met.index[:10]:
        ax.annotate(met[:10], (0, top_met.loc[met, 'CC1']), fontsize=6)
    
    ax.set_xlabel('CC1 Loading (Protein)')
    ax.set_ylabel('CC1 Loading (Metabolite)')
    ax.set_title('CCA Loading Plot')
    ax.legend()
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('results/1_CCA_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/1_CCA_visualization.png'


# ============== 2. 通路联合富集 ==============
def pathway_integration_enrichment(proteome_data, metabolome_data, deg_prot=None, deg_met=None):
    """
    蛋白组和代谢组通路联合富集分析
    
    使用不同的富集方法:
    - 蛋白: GO/KEGG
    - 代谢物: MetaboAnalyst, MSEA
    - 联合: 寻找共同的通路
    """
    try:
        from clusterProfiler import enrichKEGG
    except ImportError:
        print("请安装clusterProfiler")
        return None
    
    # 获取差异蛋白和差异代谢物
    if deg_prot is None:
        deg_prot = proteome_data[proteome_data.get('PValue', 1) < 0.05].index.tolist()
    if deg_met is None:
        deg_met = metabolome_data[metabolome_data.get('PValue', 1) < 0.05].index.tolist()
    
    print(f"差异蛋白: {len(deg_prot)}")
    print(f"差异代谢物: {len(deg_met)}")
    
    # KEGG富集 - 蛋白
    kegg_prot = enrichKEGG(
        gene=deg_prot,
        organism='hsa',
        pvalueCutoff=CONFIG['pathway_pvalue']
    )
    kegg_prot_df = pd.DataFrame(kegg_prot) if kegg_prot else pd.DataFrame()
    
    # 保存结果
    kegg_prot_df.to_csv('results/2_protein_KEGG_enrichment.csv', index=False)
    
    # 代谢组富集需要代谢物特异性分析
    # 简化版: 使用差异代谢物进行通路分析
    met_pathways = metabolomics_pathway_enrichment(deg_met)
    
    # 联合通路分析
    # 找出蛋白和代谢物共同富集的通路
    integrated_pathways = integrate_pathways(kegg_prot_df, met_pathways)
    
    return {
        'protein_pathways': kegg_prot_df,
        'metabolite_pathways': met_pathways,
        'integrated_pathways': integrated_pathways
    }


def metabolomics_pathway_enrichment(metabolites):
    """
    代谢物通路富集分析
    使用简化版映射到KEGG代谢通路
    """
    # KEGG代谢通路ID列表 (常见通路)
    metabolic_pathways = {
        'hsa00010': 'Glycolysis / Gluconeogenesis',
        'hsa00020': 'Citrate cycle (TCA cycle)',
        'hsa00030': 'Pentose phosphate pathway',
        'hsa00040': 'Pentose and glucuronate interconversions',
        'hsa00051': 'Fructose and mannose metabolism',
        'hsa00052': 'Galactose metabolism',
        'hsa00053': 'Ascorbate and aldarate metabolism',
        'hsa00061': 'Fatty acid biosynthesis',
        'hsa00071': 'Fatty acid degradation',
        'hsa00072': 'Synthesis and degradation of ketone bodies',
        'hsa00100': 'Steroid biosynthesis',
        'hsa00120': 'Primary bile acid biosynthesis',
        'hsa00220': 'Arginine biosynthesis',
        'hsa00230': 'Purine metabolism',
        'hsa00240': 'Pyrimidine metabolism',
        'hsa00250': 'Alanine, aspartate and glutamate metabolism',
        'hsa00270': 'Cysteine and methionine metabolism',
        'hsa00330': 'Arginine and proline metabolism',
        'hsa00340': 'Histidine metabolism',
        'hsa00350': 'Tyrosine metabolism',
        'hsa00360': 'Phenylalanine metabolism',
        'hsa00380': 'Tryptophan metabolism',
        'hsa00400': 'Phenylalanine, tyrosine and tryptophan biosynthesis',
        'hsa00410': 'beta-Alanine metabolism',
        'hsa00430': 'Taurine and hypotaurine metabolism',
        'hsa00450': 'Selenocompound metabolism',
        'hsa00480': 'Glutathione metabolism',
        'hsa00500': 'Starch and sucrose metabolism',
        'hsa00510': 'Amino sugar and nucleotide sugar metabolism',
        'hsa00520': 'Amino sugar and nucleotide sugar metabolism',
        'hsa00564': 'Glycerophospholipid metabolism',
        'hsa00565': 'Ether lipid metabolism',
        'hsa00590': 'Arachidonic acid metabolism',
        'hsa00591': 'Linoleic acid metabolism',
        'hsa00592': 'alpha-Linolenic acid metabolism',
        'hsa00600': 'Sphingolipid metabolism',
        'hsa00620': 'Pyruvate metabolism',
        'hsa00630': 'Glyoxylate and dicarboxylate metabolism',
        'hsa00640': 'Propanoate metabolism',
        'hsa00650': 'Butanoate metabolism',
        'hsa00670': 'One carbon pool by folate',
        'hsa00700': 'Riboflavin metabolism',
        'hsa00710': 'Carbon fixation in photosynthetic organisms',
        'hsa00730': 'Thiamine metabolism',
        'hsa00740': 'Riboflavin metabolism',
        'hsa00750': 'Vitamin B6 metabolism',
        'hsa00760': 'Nicotinate and nicotinamide metabolism',
        'hsa00770': 'Pantothenate and CoA biosynthesis',
        'hsa00780': 'Biotin metabolism',
        'hsa00790': 'Folate biosynthesis',
        'hsa00830': 'Retinol metabolism',
        'hsa00860': 'Porphyrin and chlorophyll metabolism',
        'hsa00900': 'Terpenoid backbone biosynthesis',
        'hsa00910': 'Nitrogen metabolism',
        'hsa00920': 'Sulfur metabolism',
        'hsa00930': 'Caprolactam degradation',
        'hsa00940': 'Monoterpenoid biosynthesis',
        'hsa00980': 'Metabolism of xenobiotics by cytochrome P450',
        'hsa00982': 'Drug metabolism - cytochrome P450',
        'hsa00983': 'Drug metabolism - other enzymes',
        'hsa01040': 'Biosynthesis of unsaturated fatty acids',
        'hsa01100': 'Metabolic pathways',
        'hsa01200': 'Carbon metabolism',
        'hsa01230': 'Biosynthesis of amino acids',
    }
    
    # 简化版: 直接返回通路列表供后续使用
    results = []
    for pathway_id, pathway_name in metabolic_pathways.items():
        results.append({
            'Pathway': pathway_name,
            'KEGG_ID': pathway_id,
            'N_Hits': 0,  # 需要代谢物-通路映射数据库
            'P_Value': 1
        })
    
    return pd.DataFrame(results)


def integrate_pathways(kegg_prot_df, met_pathways):
    """
    整合蛋白和代谢物的通路富集结果
    找出共同调控的代谢通路
    """
    # 简化版整合
    integrated = []
    
    if len(kegg_prot_df) > 0:
        for _, prot_path in kegg_prot_df.iterrows():
            pathway_name = prot_path.get('Description', '')
            
            # 查找代谢物通路中是否有同名
            met_match = met_pathways[
                met_pathways['Pathway'].str.contains(
                    pathway_name.split(' - ')[0] if ' - ' in pathway_name else pathway_name[:20],
                    case=False, na=False
                )
            ]
            
            if len(met_match) > 0:
                integrated.append({
                    'Pathway': pathway_name,
                    'Protein_Genes': prot_path.get('geneID', ''),
                    'Protein_Pvalue': prot_path.get('pvalue', 1),
                    'Metabolite_Hits': met_match.iloc[0]['N_Hits'],
                    'Integrated': True
                })
            else:
                integrated.append({
                    'Pathway': pathway_name,
                    'Protein_Genes': prot_path.get('geneID', ''),
                    'Protein_Pvalue': prot_path.get('pvalue', 1),
                    'Metabolite_Hits': 0,
                    'Integrated': False
                })
    
    int_df = pd.DataFrame(integrated)
    int_df.to_csv('results/2_integrated_pathways.csv', index=False)
    
    return int_df


def plot_pathway_integration(integrated_pathways):
    """可视化通路整合结果"""
    if integrated_pathways is None or len(integrated_pathways) == 0:
        return None
    
    fig, ax = plt.subplots(figsize=(12, max(8, len(integrated_pathways)*0.4)))
    
    # 按整合状态和P值排序
    int_sorted = integrated_pathways.sort_values('Protein_Pvalue')
    
    colors = ['#e74c3c' if i else '#95a5a6' for i in int_sorted['Integrated']]
    
    ax.barh(range(len(int_sorted)), -np.log10(int_sorted['Protein_Pvalue']+1e-10),
           color=colors)
    ax.set_yticks(range(len(int_sorted)))
    ax.set_yticklabels([p[:50] for p in int_sorted['Pathway']], fontsize=7)
    ax.set_xlabel('-log10(P-value)')
    ax.set_title('Integrated Pathway Analysis\n(Red=Protein+Metabolite, Gray=Protein Only)')
    ax.axvline(x=-np.log10(0.05), color='blue', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('results/2_pathway_integration.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/2_pathway_integration.png'


# ============== 3. 蛋白-代谢物因果调控关系 ==============
def causal_network_analysis(corr_df, proteome_data, metabolome_data):
    """
    构建蛋白-代谢物因果调控网络
    
    调控关系推断:
    - 蛋白上调 + 代谢物上调 -> 蛋白可能调控代谢物
    - 蛋白下调 + 代谢物下调 -> 蛋白可能调控代谢物
    - 蛋白上调 + 代谢物下调 -> 代谢物可能反馈调控蛋白
    """
    # 获取差异蛋白和代谢物
    deg_prots = set(proteome_data[proteome_data.get('PValue', 1) < 0.05].index)
    deg_mets = set(metabolome_data[metabolome_data.get('PValue', 1) < 0.05].index)
    
    # 筛选显著相关的差异蛋白-代谢物对
    sig_pairs = corr_df[corr_df['Significance']].copy()
    
    # 构建因果网络
    G = nx.DiGraph()
    
    regulation_results = []
    
    for _, row in sig_pairs.iterrows():
        prot = row['Protein']
        met = row['Metabolite']
        corr = row['Correlation']
        
        # 判断调控方向
        prot_up = prot in deg_prots and proteome_data.loc[prot, 'log2FC'] > 0 if 'log2FC' in proteome_data.columns else False
        prot_down = prot in deg_prots and proteome_data.loc[prot, 'log2FC'] < 0 if 'log2FC' in proteome_data.columns else False
        met_up = met in deg_mets and metabolome_data.loc[met, 'log2FC'] > 0 if 'log2FC' in metabolome_data.columns else False
        met_down = met in deg_mets and metabolome_data.loc[met, 'log2FC'] < 0 if 'log2FC' in metabolome_data.columns else False
        
        # 推断因果关系
        if prot_up and met_up and corr > 0:
            direction = 'Protein→Metabolite'
            G.add_edge(prot, met, relation='activates', weight=corr)
        elif prot_down and met_down and corr > 0:
            direction = 'Protein→Metabolite'
            G.add_edge(prot, met, relation='inhibits', weight=abs(corr))
        elif prot_up and met_down and corr < 0:
            direction = 'Metabolite↔Protein (Feedback)'
            G.add_edge(prot, met, relation='feedback', weight=abs(corr))
        elif prot_down and met_up and corr < 0:
            direction = 'Metabolite↔Protein (Feedback)'
            G.add_edge(prot, met, relation='feedback', weight=abs(corr))
        else:
            direction = 'Unknown'
        
        regulation_results.append({
            'Protein': prot,
            'Metabolite': met,
            'Correlation': corr,
            'Direction': direction
        })
    
    reg_df = pd.DataFrame(regulation_results)
    reg_df.to_csv('results/3_causal_regulation.csv', index=False)
    
    # 保存网络
    nx.write_graphml(G, 'results/3_protein_metabolite_network.graphml')
    
    print(f"\n=== 因果网络分析 ===")
    print(f"总调控关系: {len(reg_df)}")
    print(reg_df['Direction'].value_counts())
    
    return G, reg_df


def plot_causal_network(G, reg_df, top_n=30):
    """可视化因果网络"""
    if G.number_of_nodes() == 0:
        print("网络为空")
        return None
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # 取top N高度相关对
    top_pairs = reg_df.nlargest(top_n, 'Abs_Correlation' if 'Abs_Correlation' in reg_df.columns else 'Correlation')
    
    # 重建子网络
    H = nx.DiGraph()
    for _, row in top_pairs.iterrows():
        H.add_edge(row['Protein'], row['Metabolite'], 
                  relation=row.get('Direction', 'Unknown'))
    
    # 布局
    pos = nx.spring_layout(H, k=1, iterations=50, seed=42)
    
    # 区分蛋白和代谢物节点
    proteins = [n for n in H.nodes() if n in reg_df['Protein'].values]
    metabolites = [n for n in H.nodes() if n in reg_df['Metabolite'].values]
    
    # 绘制
    nx.draw_networkx_nodes(H, pos, nodelist=proteins,
                          node_color='#3498db', node_size=500, alpha=0.8)
    nx.draw_networkx_nodes(H, pos, nodelist=metabolites,
                          node_color='#e74c3c', node_size=300, alpha=0.8)
    
    nx.draw_networkx_edges(H, pos, alpha=0.4, arrows=True, arrowsize=15)
    nx.draw_networkx_labels(H, pos, font_size=6)
    
    plt.title('Protein-Metabolite Regulatory Network')
    
    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3498db', label='Protein'),
        Patch(facecolor='#e74c3c', label='Metabolite')
    ]
    plt.legend(handles=legend_elements)
    
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/3_causal_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/3_causal_network.png'


# ============== 4. 多组学biomarker发现 ==============
def biomarker_discovery(proteome_data, metabolome_data, corr_df, sample_info):
    """
    多组学biomarker发现
    整合蛋白和代谢物标志物
    """
    # 合并样本信息
    sample_info = sample_info.set_index('Sample')
    
    # 获取分组信息
    groups = sample_info['Group'].unique()
    if len(groups) != 2:
        print("需要两组样本进行biomarker筛选")
        return None
    
    group1, group2 = groups
    
    # 筛选差异蛋白
    prot_case = proteome_data[sample_info[sample_info['Group'] == group1].index].mean(axis=1)
    prot_control = proteome_data[sample_info[sample_info['Group'] == group2].index].mean(axis=1)
    prot_fc = np.log2(prot_case / prot_control)
    
    # 筛选差异代谢物
    met_case = metabolome_data[sample_info[sample_info['Group'] == group1].index].mean(axis=1)
    met_control = metabolome_data[sample_info[sample_info['Group'] == group2].index].mean(axis=1)
    met_fc = np.log2(met_case / met_control)
    
    # 计算联合biomarker score
    biomarker_results = []
    
    # 蛋白biomarkers
    for prot in proteome_data.index:
        if prot in corr_df['Protein'].values:
            related_mets = corr_df[corr_df['Protein'] == prot]
            
            # 找与该蛋白相关的代谢物
            for _, row in related_mets.iterrows():
                met = row['Metabolite']
                if met in metabolome_data.index:
                    combined_score = abs(row['Correlation']) * abs(prot_fc.get(prot, 0))
                    
                    biomarker_results.append({
                        'Feature': prot,
                        'Type': 'Protein',
                        'Related_Metabolite': met,
                        'log2FC': prot_fc.get(prot, 0),
                        'Correlation': row['Correlation'],
                        'Combined_Score': combined_score
                    })
    
    # 代谢物biomarkers
    for met in metabolome_data.index:
        if met in corr_df['Metabolite'].values:
            related_prots = corr_df[corr_df['Metabolite'] == met]
            
            for _, row in related_prots.iterrows():
                prot = row['Protein']
                if prot in proteome_data.index:
                    combined_score = abs(row['Correlation']) * abs(met_fc.get(met, 0))
                    
                    biomarker_results.append({
                        'Feature': met,
                        'Type': 'Metabolite',
                        'Related_Protein': prot,
                        'log2FC': met_fc.get(met, 0),
                        'Correlation': row['Correlation'],
                        'Combined_Score': combined_score
                    })
    
    bio_df = pd.DataFrame(biomarker_results)
    bio_df = bio_df.sort_values('Combined_Score', ascending=False)
    bio_df.to_csv('results/4_multiomics_biomarkers.csv', index=False)
    
    # Top biomarkers
    top_biomarkers = bio_df.head(50)
    
    print(f"\n=== Biomarker发现 ===")
    print(f"总候选数: {len(bio_df)}")
    print(f"蛋白biomarkers: {len(bio_df[bio_df['Type'] == 'Protein'])}")
    print(f"代谢物biomarkers: {len(bio_df[bio_df['Type'] == 'Metabolite'])}")
    
    return bio_df, top_biomarkers


def plot_biomarker_analysis(bio_df, top_biomarkers):
    """可视化biomarker分析"""
    if bio_df is None or len(bio_df) == 0:
        return None
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Biomarker score排名
    ax = axes[0, 0]
    top20 = top_biomarkers.head(20)
    colors = ['#3498db' if t == 'Protein' else '#e74c3c' for t in top20['Type']]
    ax.barh(range(len(top20)), top20['Combined_Score'], color=colors)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels([f"{r['Feature'][:15]} ({r['Type'][:3]})" 
                       for _, r in top20.iterrows()], fontsize=6)
    ax.set_xlabel('Combined Biomarker Score')
    ax.set_title('Top 20 Multi-omics Biomarkers')
    ax.invert_yaxis()
    
    # 2. 相关性 vs 表达变化
    ax = axes[0, 1]
    ax.scatter(bio_df['log2FC'], bio_df['Correlation'], 
              c=bio_df['Combined_Score'], cmap='viridis', alpha=0.5, s=20)
    ax.set_xlabel('log2(Fold Change)')
    ax.set_ylabel('Correlation')
    ax.set_title('Expression Change vs Correlation')
    plt.colorbar(ax.collections[0], ax=ax, label='Combined Score')
    
    # 3. 类型分布
    ax = axes[1, 0]
    type_counts = bio_df['Type'].value_counts()
    ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%',
          colors=['#3498db', '#e74c3c'])
    ax.set_title('Biomarker Type Distribution')
    
    # 4. Score分布
    ax = axes[1, 1]
    ax.hist(bio_df['Combined_Score'], bins=30, color='#9b59b6', alpha=0.7)
    ax.set_xlabel('Combined Biomarker Score')
    ax.set_ylabel('Frequency')
    ax.set_title('Biomarker Score Distribution')
    
    plt.tight_layout()
    plt.savefig('results/4_biomarker_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/4_biomarker_analysis.png'


# ============== 主流程 ==============
def main():
    """主分析流程"""
    print("=" * 60)
    print("蛋白质组 + 代谢组 联合分析流程")
    print("=" * 60)
    
    # 1. 加载数据
    print("\n[1/4] 加载数据...")
    proteome_data, metabolome_data, sample_info = load_metabolomics_data()
    
    # 2. 相关性分析 + CCA
    print("\n[2/4] 相关性热图 + CCA分析...")
    corr_df = correlation_analysis(proteome_data, metabolome_data)
    plot_correlation_heatmap(corr_df)
    cca_result = cca_analysis(proteome_data, metabolome_data, n_components=2)
    plot_cca_results(cca_result, sample_info)
    
    # 3. 通路联合富集
    print("\n[3/4] 通路联合富集分析...")
    pathway_result = pathway_integration_enrichment(proteome_data, metabolome_data)
    if pathway_result:
        plot_pathway_integration(pathway_result['integrated_pathways'])
    
    # 4. 因果网络分析
    print("\n[4/4] 蛋白-代谢物因果调控网络...")
    G, reg_df = causal_network_analysis(corr_df, proteome_data, metabolome_data)
    plot_causal_network(G, reg_df)
    
    # 5. Biomarker发现
    print("\n[Extra] 多组学biomarker发现...")
    bio_df, top_biomarkers = biomarker_discovery(proteome_data, metabolome_data, corr_df, sample_info)
    plot_biomarker_analysis(bio_df, top_biomarkers)
    
    print("\n" + "=" * 60)
    print("分析完成! 结果保存在 results/ 目录")
    print("=" * 60)


if __name__ == '__main__':
    main()
