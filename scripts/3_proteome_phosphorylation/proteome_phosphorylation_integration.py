#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework 3: Proteome + Phosphorylation Integration Analysis
蛋白质组 + 蛋白磷酸化 联合分析流程

分析模块:
1. 磷酸化位点定量 + 全局蛋白关联
2. 激酶活性推断 (KSEA)
3. 信号通路网络构建
4. 磷酸化-only蛋白分析 (信号转导核心)

Author: Bio-Design Team
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact, spearmanr
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# ============== 配置参数 ==============
CONFIG = {
    'fc_threshold': 1.5,              # 差异倍数阈值
    'pvalue_threshold': 0.05,        # P值阈值
    'phosphorylation_fc': 1.5,        # 磷酸化位点差异阈值
    'ksea_min_sites': 3,             # KSEA分析最少位点数
    'kinase_window': 7,              # 激酶底物分析窗口
    'ppi_threshold': 0.7,            # PPI阈值
    'species': 'human',
    'species_taxid': 9606,
}

# ============== 数据加载 ==============
def load_phosphorylation_data():
    """
    加载磷酸化组学数据
    """
    phospho_sites = pd.read_csv('input/phosphorylation_sites.csv', index_col=0)
    proteome_data = pd.read_csv('input/proteome_quantification.csv', index_col=0)
    sample_info = pd.read_csv('input/sample_info.csv')
    
    return phospho_sites, proteome_data, sample_info


# ============== 1. 磷酸化位点定量 + 全局蛋白关联 ==============
def phospho_proteome_association(phospho_sites, proteome_data):
    """
    磷酸化位点与全局蛋白定量关联
    分析位点水平磷酸化变化与蛋白水平变化的关系
    """
    results = []
    
    # 获取所有蛋白
    all_proteins = set(phospho_sites.index) | set(proteome_data.index)
    
    for protein in all_proteins:
        # 蛋白表达信息
        if protein in proteome_data.index:
            prot_data = proteome_data.loc[protein]
            prot_expr = prot_data.mean() if hasattr(prot_data, 'mean') else prot_data.iloc[0]
            prot_change = prot_data.get('log2FC', 0) if hasattr(prot_data, 'get') else 0
            prot_sig = prot_data.get('PValue', 1) < 0.05 if hasattr(prot_data, 'get') else False
        else:
            prot_expr = 0
            prot_change = 0
            prot_sig = False
        
        # 磷酸化位点信息
        if protein in phospho_sites.index:
            sites = phospho_sites.loc[protein]
            
            if isinstance(sites, pd.DataFrame):
                n_sites = len(sites)
                site_changes = sites['log2FC'].values if 'log2FC' in sites.columns else [0]
                site_ps = sites['PValue'].values if 'PValue' in sites.columns else [1]
                mean_change = np.mean(site_changes)
                n_sig = sum(1 for p in site_ps if p < 0.05)
            else:
                n_sites = 1
                mean_change = sites.get('log2FC', 0)
                n_sig = 1 if sites.get('PValue', 1) < 0.05 else 0
        else:
            n_sites = 0
            mean_change = 0
            n_sig = 0
        
        results.append({
            'ProteinID': protein,
            'ProteinExpression': prot_expr,
            'Protein_log2FC': prot_change,
            'Protein_Significant': prot_sig,
            'N_Phospho_Sites': n_sites,
            'N_Significant_Sites': n_sig,
            'Phospho_log2FC': mean_change
        })
    
    assoc_df = pd.DataFrame(results)
    assoc_df.to_csv('results/1_phospho_proteome_association.csv', index=False)
    
    print(f"=== 磷酸化-蛋白关联 ===")
    print(f"有磷酸化位点的蛋白: {len(assoc_df[assoc_df['N_Phospho_Sites'] > 0])}")
    print(f"平均位点数: {assoc_df['N_Phospho_Sites'].mean():.2f}")
    
    return assoc_df


def classify_phospho_proteins(assoc_df):
    """
    磷酸化蛋白分类:
    - Phospho-up + Protein-up
    - Phospho-up + Protein-down
    - Phospho-up + Protein-unchanged (磷酸化-only)
    - Phospho-down + Protein-down
    - Phospho-down + Protein-up
    - Phospho-down + Protein-unchanged (磷酸化-only)
    - No change
    """
    def classify(row):
        p_sig = row['Protein_Significant']
        p_fc = row['Protein_log2FC']
        ph_fc = row['Phospho_log2FC']
        n_sig = row['N_Significant_Sites']
        
        if n_sig == 0:
            return 'No_Phospho_Change'
        
        if p_sig and abs(p_fc) > 0.5:
            if ph_fc > 0.5 and p_fc > 0:
                return 'PhosphoUp_ProtUp'
            elif ph_fc > 0.5 and p_fc < 0:
                return 'PhosphoUp_ProtDown'
            elif ph_fc < -0.5 and p_fc < 0:
                return 'PhosphoDown_ProtDown'
            elif ph_fc < -0.5 and p_fc > 0:
                return 'PhosphoDown_ProtUp'
            elif ph_fc > 0.5:
                return 'PhosphoUp_ProtNS'
            elif ph_fc < -0.5:
                return 'PhosphoDown_ProtNS'
        
        # 磷酸化-only蛋白 (蛋白无显著变化但磷酸化显著)
        if not p_sig and n_sig > 0:
            if ph_fc > 0.5:
                return 'PhosphoOnly_Up'
            elif ph_fc < -0.5:
                return 'PhosphoOnly_Down'
        
        return 'Other'
    
    assoc_df['PhosphoClass'] = assoc_df.apply(classify, axis=1)
    assoc_df.to_csv('results/1_phospho_classification.csv', index=False)
    
    print(f"\n=== 磷酸化蛋白分类 ===")
    print(assoc_df['PhosphoClass'].value_counts())
    
    return assoc_df


# ============== 2. 激酶活性推断 (KSEA) ==============
def ksea_analysis(phospho_sites, kinase_substrate_db=None):
    """
    Kinase Substrate Enrichment Analysis (KSEA)
    
    使用已知激酶-底物关系数据库推断激酶活性变化
    
    激酶活性评分:
    Z-score = (mean(PhosFC) - mean(all_PhosFC)) / std(all_PhosFC)
    
    需要激酶-底物关系数据库 (如 PhosphoSitePlus, KEA2)
    """
    # 简化版KSEA - 实际应使用完整数据库
    # 这里提供一个框架，实际使用需要加载数据库
    
    if kinase_substrate_db is None:
        # 使用内置简化版
        kinase_substrate_db = get_kinase_substrate_builtin()
    
    # 获取所有磷酸化位点的变化
    site_changes = []
    for protein in phospho_sites.index:
        sites = phospho_sites.loc[protein]
        if isinstance(sites, pd.DataFrame):
            for _, site in sites.iterrows():
                if 'log2FC' in site:
                    site_changes.append({
                        'Protein': protein,
                        'Site': site.get('Site', ''),
                        'log2FC': site['log2FC']
                    })
        else:
            if 'log2FC' in sites:
                site_changes.append({
                    'Protein': protein,
                    'Site': sites.get('Site', ''),
                    'log2FC': sites['log2FC']
                })
    
    site_df = pd.DataFrame(site_changes)
    
    # 计算所有位点的均值和标准差
    global_mean = site_df['log2FC'].mean()
    global_std = site_df['log2FC'].std()
    
    # 对每个激酶计算富集评分
    ksea_results = []
    
    for kinase, substrates in kinase_substrate_db.items():
        kinase_sites = site_df[site_df['Protein'].isin(substrates)]
        
        if len(kinase_sites) >= CONFIG['ksea_min_sites']:
            kinase_mean = kinase_sites['log2FC'].mean()
            kinase_std = kinase_sites['log2FC'].std()
            n_sites = len(kinase_sites)
            
            # Z-score计算
            if kinase_std > 0:
                z_score = (kinase_mean - global_mean) / (global_std / np.sqrt(n_sites))
            else:
                z_score = 0
            
            # P-value from Z-score
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
            
            ksea_results.append({
                'Kinase': kinase,
                'N_Substrates': n_sites,
                'Mean_Substrate_FC': kinase_mean,
                'Z_Score': z_score,
                'P_Value': p_value,
                'Activity': 'Increased' if z_score > 2 else ('Decreased' if z_score < -2 else 'Unchanged')
            })
    
    ksea_df = pd.DataFrame(ksea_results)
    ksea_df['FDR'] = stats.false_discovery_control(ksea_df['P_Value'])
    ksea_df = ksea_df.sort_values('Z_Score', ascending=False)
    ksea_df.to_csv('results/2_KSEA_results.csv', index=False)
    
    print(f"\n=== KSEA激酶活性推断 ===")
    print(f"分析激酶数: {len(ksea_df)}")
    print(f"显著激活 (Z>2): {len(ksea_df[ksea_df['Z_Score'] > 2])}")
    print(f"显著抑制 (Z<-2): {len(ksea_df[ksea_df['Z_Score'] < -2])}")
    
    return ksea_df


def get_kinase_substrate_builtin():
    """
    内置简化版激酶-底物关系
    实际分析应使用PhosphoSitePlus等完整数据库
    """
    # 简化版，仅用于演示
    kinase_db = {
        'AKT1': ['GSK3B', 'MTOR', 'BAD', 'PRAS40'],
        'MAPK1': ['ELK1', 'MYC', 'JUN', 'FOS'],
        'MAPK3': ['ELK1', 'MYC', 'JUN', 'FOS'],
        'CDK1': ['PLK1', 'AURKB', 'HIST1H1B'],
        'PRKA': ['ACC1', 'AMPK'],
        'MTOR': ['RPS6KB1', '4EBP1', 'AKT1'],
        'SRC': ['STAT3', 'FAK1', 'BCR'],
        'EGFR': ['STAT3', 'AKT1', 'MAPK1'],
    }
    return kinase_db


def plot_ksea_results(ksea_df, top_n=20):
    """可视化KSEA结果"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Z-score柱状图
    top_kinases = pd.concat([ksea_df.head(top_n), ksea_df.tail(top_n)])
    top_kinases = top_kinases.sort_values('Z_Score')
    
    colors = ['#e74c3c' if z > 0 else '#3498db' for z in top_kinases['Z_Score']]
    axes[0].barh(range(len(top_kinases)), top_kinases['Z_Score'], color=colors)
    axes[0].set_yticks(range(len(top_kinases)))
    axes[0].set_yticklabels(top_kinases['Kinase'], fontsize=8)
    axes[0].axvline(x=2, color='red', linestyle='--', alpha=0.5)
    axes[0].axvline(x=-2, color='blue', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('Z-Score')
    axes[0].set_title('Kinase Activity (KSEA Z-Score)')
    
    # 2. 激酶活性分类饼图
    activity_counts = ksea_df['Activity'].value_counts()
    colors_pie = {
        'Increased': '#e74c3c',
        'Decreased': '#3498db',
        'Unchanged': '#95a5a6'
    }
    axes[1].pie(activity_counts.values, 
               labels=activity_counts.index,
               colors=[colors_pie.get(a, '#95a5a6') for a in activity_counts.index],
               autopct='%1.1f%%')
    axes[1].set_title('Kinase Activity Distribution')
    
    plt.tight_layout()
    plt.savefig('results/2_KSEA_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/2_KSEA_visualization.png'


# ============== 3. 信号通路网络构建 ==============
def build_signaling_network(phospho_sites, ksea_results, assoc_df):
    """
    构建磷酸化信号通路网络
    整合:
    - 显著变化的磷酸化位点
    - 激酶活性变化
    - 蛋白相互作用
    """
    # 获取磷酸化-only蛋白 (信号转导核心)
    phospho_only = assoc_df[assoc_df['PhosphoClass'].isin(['PhosphoOnly_Up', 'PhosphoOnly_Down'])]
    
    # 构建网络
    G = nx.DiGraph()
    
    # 添加节点
    for _, row in ksea_results.iterrows():
        if abs(row['Z_Score']) > 1.5:
            G.add_node(
                row['Kinase'],
                type='kinase',
                activity=row['Activity'],
                z_score=row['Z_Score']
            )
    
    for _, row in phospho_only.iterrows():
        G.add_node(
            row['ProteinID'],
            type='phospho_protein',
            phospho_class=row['PhosphoClass'],
            change=row['Phospho_log2FC']
        )
    
    # 添加边 (激酶-底物关系)
    # 简化版 - 实际需要完整数据库
    kinase_substrate = {
        'AKT1': ['GSK3B', 'MTOR'],
        'MAPK1': ['ELK1', 'FOS'],
        'MTOR': ['RPS6KB1', '4EBP1'],
    }
    
    for kinase, substrates in kinase_substrate.items():
        if kinase in G.nodes():
            for sub in substrates:
                if sub in G.nodes():
                    G.add_edge(kinase, sub, relation='phosphorylates')
    
    # 保存网络
    nx.write_graphml(G, 'results/3_signaling_network.graphml')
    
    print(f"\n=== 信号网络构建 ===")
    print(f"节点数: {G.number_of_nodes()}")
    print(f"边数: {G.number_of_edges()}")
    
    return G, phospho_only


def visualize_signaling_network(G, phospho_only):
    """可视化信号网络"""
    if G.number_of_nodes() == 0:
        print("网络为空")
        return None
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # 布局
    pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    
    # 节点类型分类
    kinases = [n for n, d in G.nodes(data=True) if d.get('type') == 'kinase']
    phos_prots = [n for n, d in G.nodes(data=True) if d.get('type') == 'phospho_protein']
    
    # 绘制节点
    nx.draw_networkx_nodes(G, pos, nodelist=kinases, 
                          node_color='#e74c3c', node_size=800, alpha=0.8,
                          node_shape='s')
    nx.draw_networkx_nodes(G, pos, nodelist=phos_prots,
                          node_color='#3498db', node_size=400, alpha=0.8)
    
    # 绘制边
    nx.draw_networkx_edges(G, pos, alpha=0.4, arrows=True,
                         edge_color='gray', arrowsize=15)
    
    # 绘制标签
    nx.draw_networkx_labels(G, pos, font_size=7)
    
    plt.title('Phosphorylation Signaling Network')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/3_signaling_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/3_signaling_network.png'


# ============== 4. 磷酸化-only蛋白分析 ==============
def phospho_only_analysis(assoc_df, proteome_data):
    """
    分析磷酸化-only蛋白 (信号转导核心)
    这些蛋白的表达没有显著变化，但磷酸化状态发生了显著改变
    可能是信号转导的核心调控节点
    """
    # 筛选磷酸化-only蛋白
    phospho_only = assoc_df[assoc_df['PhosphoClass'].isin(['PhosphoOnly_Up', 'PhosphoOnly_Down'])]
    
    print(f"\n=== 磷酸化-only蛋白分析 ===")
    print(f"磷酸化激活蛋白: {len(phospho_only[phospho_only['PhosphoClass'] == 'PhosphoOnly_Up'])}")
    print(f"磷酸化抑制蛋白: {len(phospho_only[phospho_only['PhosphoClass'] == 'PhosphoOnly_Down'])}")
    
    # 获取磷酸化-only蛋白的详细信息
    phospho_only_genes = phospho_only['ProteinID'].tolist()
    
    # 导出列表
    phospho_only.to_csv('results/4_phospho_only_proteins.csv', index=False)
    
    return phospho_only, phospho_only_genes


def phospho_only_functional_enrichment(phospho_only_genes):
    """
    对磷酸化-only蛋白进行功能富集
    重点关注信号通路
    """
    try:
        from clusterProfiler import enrichGO, enrichKEGG
        import org.Hs.eg.db as orgdb
    except ImportError:
        print("请安装clusterProfiler")
        return None
    
    # GO富集
    go_results = enrichGO(
        gene=phospho_only_genes,
        OrgDb=orgdb,
        keyType='SYMBOL',
        ont='ALL',
        pvalueCutoff=0.05
    )
    
    go_df = pd.DataFrame(go_results)
    go_df.to_csv('results/4_phospho_only_GO_enrichment.csv', index=False)
    
    # KEGG富集
    kegg_results = enrichKEGG(
        gene=phospho_only_genes,
        organism='hsa',
        pvalueCutoff=0.05
    )
    
    kegg_df = pd.DataFrame(kegg_results)
    kegg_df.to_csv('results/4_phospho_only_KEGG_enrichment.csv', index=False)
    
    # 筛选信号通路相关
    signal_pathways = kegg_df[kegg_df['Description'].str.contains(
        'signal|signaling|PI3K|AKT|MAPK|EGFR|insulin|rapamycin', 
        case=False, na=False
    )]
    
    print(f"\n信号通路相关: {len(signal_pathways)}")
    
    return {'go': go_df, 'kegg': kegg_df, 'signal_pathways': signal_pathways}


def plot_phospho_only_analysis(phospho_only, enrich_result):
    """可视化磷酸化-only蛋白分析"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 磷酸化only蛋白分类
    phos_up = phospho_only[phospho_only['PhosphoClass'] == 'PhosphoOnly_Up']
    phos_down = phospho_only[phospho_only['PhosphoClass'] == 'PhosphoOnly_Down']
    
    labels = ['Phospho-Only Up', 'Phospho-Only Down']
    sizes = [len(phos_up), len(phos_down)]
    colors = ['#e74c3c', '#3498db']
    
    axes[0, 0].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%')
    axes[0, 0].set_title('Phospho-Only Protein Categories')
    
    # 2. 磷酸化变化分布
    axes[0, 1].hist(phospho_only['Phospho_log2FC'], bins=30, 
                    color='#3498db', alpha=0.7, edgecolor='black')
    axes[0, 1].axvline(x=0, color='black', linestyle='--')
    axes[0, 1].set_xlabel('Phosphorylation log2FC')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Phospho-Only Protein Change Distribution')
    
    # 3. 位点数分布
    site_counts = phospho_only['N_Phospho_Sites']
    axes[1, 0].hist(site_counts, bins=range(1, site_counts.max()+2),
                    color='#9b59b6', alpha=0.7, edgecolor='black')
    axes[1, 0].set_xlabel('Number of Phosphorylation Sites')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Site Count Distribution (Phospho-Only)')
    
    # 4. 信号通路富集
    if enrich_result is not None and len(enrich_result['kegg']) > 0:
        kegg_sig = enrich_result['kegg'].head(15)
        axes[1, 1].barh(range(len(kegg_sig)), -np.log10(kegg_sig['pvalue']))
        axes[1, 1].set_yticks(range(len(kegg_sig)))
        axes[1, 1].set_yticklabels([t[:35] for t in kegg_sig['Description']], fontsize=7)
        axes[1, 1].set_xlabel('-log10(P-value)')
        axes[1, 1].set_title('KEGG Pathway Enrichment (Phospho-Only)')
        axes[1, 1].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig('results/4_phospho_only_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/4_phospho_only_analysis.png'


# ============== 主流程 ==============
def main():
    """主分析流程"""
    print("=" * 60)
    print("蛋白质组 + 磷酸化组 联合分析流程")
    print("=" * 60)
    
    # 1. 加载数据
    print("\n[1/4] 加载数据...")
    phospho_sites, proteome_data, sample_info = load_phosphorylation_data()
    
    # 2. 磷酸化-蛋白关联
    print("\n[2/4] 磷酸化位点与蛋白表达关联...")
    assoc_df = phospho_proteome_association(phospho_sites, proteome_data)
    assoc_df = classify_phospho_proteins(assoc_df)
    
    # 3. KSEA激酶活性推断
    print("\n[3/4] KSEA激酶活性推断...")
    ksea_df = ksea_analysis(phospho_sites)
    plot_ksea_results(ksea_df)
    
    # 4. 信号通路网络
    print("\n[4/4] 信号通路网络构建...")
    G, phospho_only = build_signaling_network(phospho_sites, ksea_df, assoc_df)
    visualize_signaling_network(G, phospho_only)
    
    # 5. 磷酸化-only蛋白分析
    print("\n[Extra] 磷酸化-only蛋白分析...")
    phospho_only, phospho_only_genes = phospho_only_analysis(assoc_df, proteome_data)
    enrich_result = phospho_only_functional_enrichment(phospho_only_genes)
    plot_phospho_only_analysis(phospho_only, enrich_result)
    
    print("\n" + "=" * 60)
    print("分析完成! 结果保存在 results/ 目录")
    print("=" * 60)


if __name__ == '__main__':
    main()
