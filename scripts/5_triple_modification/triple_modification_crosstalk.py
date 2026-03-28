#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework 5: Triple Modification Crosstalk Analysis (Advanced)
蛋白质组 + 磷酸化组 + 乳酸化组 三组学修饰串扰分析

分析模块:
1. 修饰位点-蛋白表达关联分类
2. 磷酸化-乳酸化位点级串扰分析
3. 酶-底物网络构建
4. 激酶/乳酸化酶活性推断

Author: Bio-Design Team
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact, chi2_contingency, spearmanr
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# ============== 配置参数 ==============
CONFIG = {
    'fc_threshold': 1.5,              # 差异倍数阈值
    'pvalue_threshold': 0.05,        # P值阈值
    'crosstalk_window': 10,           # 串扰分析窗口 (氨基酸)
    'min_co_occurrence': 2,          # 最小共现次数
    'ksea_min_sites': 3,             # KSEA最少位点
    'species': 'human',
    'species_taxid': 9606,
}

# ============== 数据加载 ==============
def load_triple_modification_data():
    """
    加载三组学数据
    """
    proteome_data = pd.read_csv('input/proteome_quantification.csv', index_col=0)
    phospho_sites = pd.read_csv('input/phosphorylation_sites.csv', index_col=0)
    lactylation_sites = pd.read_csv('input/lactylation_sites.csv', index_col=0)
    sample_info = pd.read_csv('input/sample_info.csv')
    
    return proteome_data, phospho_sites, lactylation_sites, sample_info


# ============== 1. 修饰位点-蛋白表达关联分类 ==============
def modification_protein_association(proteome_data, phospho_sites, lactylation_sites):
    """
    三组学修饰位点与蛋白表达关联分析
    
    对每个蛋白，关联其:
    - 蛋白表达变化
    - 磷酸化位点变化
    - 乳酸化位点变化
    """
    results = []
    
    # 获取所有蛋白
    all_proteins = set(proteome_data.index) | set(phospho_sites.index) | set(lactylation_sites.index)
    
    for protein in all_proteins:
        # 蛋白表达
        if protein in proteome_data.index:
            prot_data = proteome_data.loc[protein]
            prot_expr = prot_data.mean() if hasattr(prot_data, 'mean') else prot_data.iloc[0]
            prot_fc = prot_data.get('log2FC', 0) if hasattr(prot_data, 'get') else 0
            prot_p = prot_data.get('PValue', 1) if hasattr(prot_data, 'get') else 1
            prot_sig = prot_p < 0.05
        else:
            prot_expr, prot_fc, prot_p, prot_sig = 0, 0, 1, False
        
        # 磷酸化位点
        if protein in phospho_sites.index:
            phos = phospho_sites.loc[protein]
            if isinstance(phos, pd.DataFrame):
                n_phos = len(phos)
                phos_fc = phos['log2FC'].mean() if 'log2FC' in phos.columns else 0
                phos_sig = (phos['PValue'] < 0.05).sum() if 'PValue' in phos.columns else 0
                phos_sites_list = list(phos.index)
            else:
                n_phos = 1
                phos_fc = phos.get('log2FC', 0)
                phos_sig = 1 if phos.get('PValue', 1) < 0.05 else 0
                phos_sites_list = [phos.index[0]] if hasattr(phos, 'index') else []
        else:
            n_phos, phos_fc, phos_sig = 0, 0, 0
            phos_sites_list = []
        
        # 乳酸化位点
        if protein in lactylation_sites.index:
            lact = lactylation_sites.loc[protein]
            if isinstance(lact, pd.DataFrame):
                n_lact = len(lact)
                lact_fc = lact['log2FC'].mean() if 'log2FC' in lact.columns else 0
                lact_sig = (lact['PValue'] < 0.05).sum() if 'PValue' in lact.columns else 0
                lact_sites_list = list(lact.index)
            else:
                n_lact = 1
                lact_fc = lact.get('log2FC', 0)
                lact_sig = 1 if lact.get('PValue', 1) < 0.05 else 0
                lact_sites_list = [lact.index[0]] if hasattr(lact, 'index') else []
        else:
            n_lact, lact_fc, lact_sig = 0, 0, 0
            lact_sites_list = []
        
        results.append({
            'ProteinID': protein,
            'Protein_log2FC': prot_fc,
            'Protein_Significant': prot_sig,
            'N_Phospho_Sites': n_phos,
            'N_Phospho_Significant': phos_sig,
            'Phospho_log2FC': phos_fc,
            'N_Lactylation_Sites': n_lact,
            'N_Lactylation_Significant': lact_sig,
            'Lactylation_log2FC': lact_fc,
            'Total_Modifications': n_phos + n_lact
        })
    
    assoc_df = pd.DataFrame(results)
    assoc_df.to_csv('results/1_triple_modification_association.csv', index=False)
    
    print(f"\n=== 三组学修饰关联 ===")
    print(f"总蛋白数: {len(assoc_df)}")
    print(f"有磷酸化位点: {len(assoc_df[assoc_df['N_Phospho_Sites'] > 0])}")
    print(f"有乳酸化位点: {len(assoc_df[assoc_df['N_Lactylation_Sites'] > 0])}")
    print(f"同时有磷酸化和乳酸化: {len(assoc_df[(assoc_df['N_Phospho_Sites'] > 0) & (assoc_df['N_Lactylation_Sites'] > 0)])}")
    
    return assoc_df


def classify_modification_patterns(assoc_df):
    """
    修饰模式分类
    
    基于蛋白表达变化和修饰位点变化，对蛋白进行分类:
    
    1. Triple-Active: 蛋白↑ + 磷酸化↑ + 乳酸化↑
    2. Triple-Suppressed: 蛋白↓ + 磷酸化↓ + 乳酸化↓
    3. Phospho-Specific: 仅磷酸化变化
    4. Lact-Specific: 仅乳酸化变化  
    5. Modification-Coordinated: 磷酸化和乳酸化同向，但与蛋白反向
    6. Modification-Opposite: 磷酸化和乳酸化反向
    7. No-Modification: 无修饰变化
    """
    
    def classify(row):
        p_sig = row['Protein_Significant']
        p_fc = row['Protein_log2FC']
        phos_sig = row['N_Phospho_Significant'] > 0
        phos_fc = row['Phospho_log2FC']
        lact_sig = row['N_Lactylation_Significant'] > 0
        lact_fc = row['Lactylation_log2FC']
        
        # 判断修饰变化方向
        phos_up = phos_sig and phos_fc > 0.5
        phos_down = phos_sig and phos_fc < -0.5
        lact_up = lact_sig and lact_fc > 0.5
        lact_down = lact_sig and lact_fc < -0.5
        
        if phos_up and lact_up:
            if p_sig and p_fc > 0:
                return 'Triple-Active'
            elif p_sig and p_fc < 0:
                return 'Coordinated-Up'
        elif phos_down and lact_down:
            if p_sig and p_fc < 0:
                return 'Triple-Suppressed'
            elif p_sig and p_fc > 0:
                return 'Coordinated-Down'
        elif phos_up and lact_down:
            return 'Modification-Opposite'
        elif phos_down and lact_up:
            return 'Modification-Opposite'
        elif phos_up or phos_down:
            return 'Phospho-Specific'
        elif lact_up or lact_down:
            return 'Lact-Specific'
        else:
            return 'No-Modification'
    
    assoc_df['Modification_Pattern'] = assoc_df.apply(classify, axis=1)
    assoc_df.to_csv('results/1_modification_pattern_classification.csv', index=False)
    
    print(f"\n=== 修饰模式分类 ===")
    print(assoc_df['Modification_Pattern'].value_counts())
    
    return assoc_df


def plot_modification_patterns(assoc_df):
    """可视化修饰模式分类"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 修饰模式分布
    pattern_counts = assoc_df['Modification_Pattern'].value_counts()
    colors = plt.cm.Set3(np.linspace(0, 1, len(pattern_counts)))
    axes[0, 0].pie(pattern_counts.values, labels=pattern_counts.index,
                  autopct='%1.1f%%', colors=colors)
    axes[0, 0].set_title('Modification Pattern Distribution')
    
    # 2. 磷酸化 vs 乳酸化散点图
    ax = axes[0, 1]
    has_mod = assoc_df[(assoc_df['N_Phospho_Sites'] > 0) | (assoc_df['N_Lactylation_Sites'] > 0)]
    scatter = ax.scatter(has_mod['Phospho_log2FC'], has_mod['Lactylation_log2FC'],
                        c=has_mod['Total_Modifications'], cmap='viridis',
                        alpha=0.6, s=30)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Phosphorylation log2FC')
    ax.set_ylabel('Lactylation log2FC')
    ax.set_title('Phosphorylation vs Lactylation Changes')
    plt.colorbar(scatter, ax=ax, label='Total Modifications')
    
    # 3. 各模式蛋白的表达和修饰变化
    ax = axes[1, 0]
    patterns = assoc_df['Modification_Pattern'].unique()[:6]  # 取前6个
    data_for_box = []
    labels = []
    for p in patterns:
        subset = assoc_df[assoc_df['Modification_Pattern'] == p]
        data_for_box.append(subset['Protein_log2FC'].values)
        labels.append(p[:10])
    
    bp = ax.boxplot(data_for_box, labels=labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], plt.cm.Set3(np.linspace(0, 1, len(patterns)))):
        patch.set_facecolor(color)
    ax.set_ylabel('Protein log2FC')
    ax.set_title('Protein Expression by Modification Pattern')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 4. 蛋白-修饰位点数量关系
    ax = axes[1, 1]
    ax.scatter(assoc_df['N_Phospho_Sites'], assoc_df['N_Lactylation_Sites'],
              alpha=0.3, s=20, c='#3498db')
    ax.set_xlabel('Number of Phosphorylation Sites')
    ax.set_ylabel('Number of Lactylation Sites')
    ax.set_title('Site Count Correlation')
    
    plt.tight_layout()
    plt.savefig('results/1_modification_patterns.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/1_modification_patterns.png'


# ============== 2. 磷酸化-乳酸化位点级串扰分析 ==============
def site_level_crosstalk(phospho_sites, lactylation_sites, protein_sequence=None):
    """
    位点级磷酸化-乳酸化串扰分析
    
    分析同蛋白上不同修饰位点之间的:
    1. 位置关系 (是否在邻近区域)
    2. 表达相关性
    3. 潜在竞争关系
    """
    crosstalk_results = []
    
    # 找到同时有磷酸化和乳酸化位点的蛋白
    common_proteins = set(phospho_sites.index) & set(lactylation_sites.index)
    print(f"同时有磷酸化和乳酸化位点的蛋白: {len(common_proteins)}")
    
    for protein in common_proteins:
        phos = phospho_sites.loc[protein]
        lact = lactylation_sites.loc[protein]
        
        # 处理多行情况
        if isinstance(phos, pd.DataFrame):
            phos_positions = phos.index.tolist() if hasattr(phos.index, 'tolist') else list(phos.index)
            phos_fcs = phos['log2FC'].values if 'log2FC' in phos.columns else [0]*len(phos)
            phos_ps = phos['PValue'].values if 'PValue' in phos.columns else [1]*len(phos)
        else:
            phos_positions = [phos.index[0] if hasattr(phos.index, '__iter__') else 0]
            phos_fcs = [phos.get('log2FC', 0)]
            phos_ps = [phos.get('PValue', 1)]
        
        if isinstance(lact, pd.DataFrame):
            lact_positions = lact.index.tolist() if hasattr(lact.index, 'tolist') else list(lact.index)
            lact_fcs = lact['log2FC'].values if 'log2FC' in lact.columns else [0]*len(lact)
            lact_ps = lact['PValue'].values if 'PValue' in lact.columns else [1]*len(lact)
        else:
            lact_positions = [lact.index[0] if hasattr(lact.index, '__iter__') else 0]
            lact_fcs = [lact.get('log2FC', 0)]
            lact_ps = [lact.get('PValue', 1)]
        
        # 分析每个磷酸化-乳酸化位点对
        for i, (pp, pf, ppv) in enumerate(zip(phos_positions, phos_fcs, phos_ps)):
            for j, (lp, lf, lpv) in enumerate(zip(lact_positions, lact_fcs, lact_ps)):
                # 位置接近度
                try:
                    pos_dist = abs(int(pp) - int(lp)) if str(pp).isdigit() and str(lp).isdigit() else 999
                except:
                    pos_dist = 999
                
                # 相关性 (如果有重复样本)
                # 这里简化处理，实际应计算同一样本的相关性
                
                crosstalk_results.append({
                    'Protein': protein,
                    'Phospho_Site': pp,
                    'Lact_Site': lp,
                    'Position_Distance': pos_dist,
                    'Phospho_log2FC': pf,
                    'Lact_log2FC': lf,
                    'Phospho_PValue': ppv,
                    'Lact_PValue': lpv,
                    'Co_Localized': pos_dist < CONFIG['crosstalk_window'],
                    'Both_Significant': ppv < 0.05 and lpv < 0.05,
                    'Same_Direction': (pf > 0) == (lf > 0)
                })
    
    crosstalk_df = pd.DataFrame(crosstalk_results)
    crosstalk_df.to_csv('results/2_site_level_crosstalk.csv', index=False)
    
    print(f"\n=== 位点级串扰分析 ===")
    print(f"总位点对数: {len(crosstalk_df)}")
    print(f"共定位 (<{CONFIG['crosstalk_window']}AA): {crosstalk_df['Co_Localized'].sum()}")
    print(f"同时显著: {crosstalk_df['Both_Significant'].sum()}")
    print(f"同向变化: {crosstalk_df['Same_Direction'].sum()}")
    
    return crosstalk_df


def plot_crosstalk_analysis(crosstalk_df):
    """可视化串扰分析"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 位置距离分布
    ax = axes[0, 0]
    ax.hist(crosstalk_df['Position_Distance'], bins=50, color='#3498db', alpha=0.7)
    ax.axvline(x=CONFIG['crosstalk_window'], color='red', linestyle='--', 
              label=f'Window={CONFIG["crosstalk_window"]}')
    ax.set_xlabel('Position Distance (AA)')
    ax.set_ylabel('Count')
    ax.set_title('Phospho-Lact Site Distance Distribution')
    ax.legend()
    
    # 2. 共定位位点对的FC相关性
    ax = axes[0, 1]
    co_local = crosstalk_df[crosstalk_df['Co_Localized']]
    ax.scatter(co_local['Phospho_log2FC'], co_local['Lact_log2FC'],
              alpha=0.6, s=40, c='#e74c3c')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Phosphorylation log2FC')
    ax.set_ylabel('Lactylation log2FC')
    ax.set_title(f'Co-Localized Sites (n={len(co_local)})')
    
    # 3. 同向 vs 反向变化
    ax = axes[1, 0]
    same_dir = crosstalk_df['Same_Direction'].value_counts()
    colors = ['#2ecc71', '#e74c3c']
    ax.pie(same_dir.values, labels=['Same Direction', 'Opposite Direction'],
          autopct='%1.1f%%', colors=colors)
    ax.set_title('Direction Consistency')
    
    # 4. 串扰热图 (按蛋白)
    ax = axes[1, 1]
    # 取位点对最多的蛋白
    top_proteins = crosstalk_df.groupby('Protein').size().nlargest(20).index
    top_crosstalk = crosstalk_df[crosstalk_df['Protein'].isin(top_proteins)]
    
    # 创建简化热图矩阵
    protein_pairs = top_crosstalk.groupby('Protein').agg({
        'Phospho_log2FC': 'mean',
        'Lact_log2FC': 'mean',
        'Co_Localized': 'sum'
    }).head(20)
    
    x = range(len(protein_pairs))
    width = 0.35
    ax.bar([i - width/2 for i in x], protein_pairs['Phospho_log2FC'], 
          width, label='Phosphorylation', color='#3498db')
    ax.bar([i + width/2 for i in x], protein_pairs['Lact_log2FC'],
          width, label='Lactylation', color='#e74c3c')
    ax.set_xticks(x)
    ax.set_xticklabels([p[:10] for p in protein_pairs.index], rotation=45, ha='right', fontsize=6)
    ax.set_ylabel('log2FC')
    ax.set_title('Top 20 Proteins with Most Crosstalk Sites')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('results/2_crosstalk_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/2_crosstalk_analysis.png'


# ============== 3. 酶-底物网络构建 ==============
def enzyme_substrate_network(assoc_df, phospho_sites, lactylation_sites, ksea_results=None):
    """
    构建酶-底物调控网络
    
    整合:
    - 激酶-磷酸化底物
    - Writer/Eraser-乳酸化底物 (简化处理，使用假设的writer/eraser)
    - 蛋白表达变化
    """
    G = nx.DiGraph()
    
    # 简化版酶-底物关系
    # 实际应使用磷酸化位点数据库 (PhosphoSitePlus, dbPTM等)
    
    known_kinases = {
        'AKT1': ['GSK3B', 'MTOR', 'BAD', 'PRAS40', 'ELK1', 'FOS'],
        'MAPK1': ['ELK1', 'MYC', 'JUN', 'FOS', 'ATF2'],
        'MAPK3': ['ELK1', 'MYC', 'JUN', 'FOS'],
        'CDK1': ['PLK1', 'AURKB', 'CCNB1'],
        'PRKA': ['ACC1', 'AMPK'],
        'MTOR': ['RPS6KB1', '4EBP1', 'AKT1'],
        'SRC': ['STAT3', 'FAK1'],
        'EGFR': ['STAT3', 'AKT1', 'MAPK1'],
        'PKCA': ['RAF1', 'MEK1'],
    }
    
    # 假设的乳酸化酶 (简化版)
    # 实际乳酸化酶尚未完全确定，这里使用假设
    known_lactylases = {
        'HDAC1': ['P53', 'STAT3'],  # 去乳酸化酶候选
        'SIRT1': ['P53', 'HIF1A'],  # 去乳酸化酶候选
    }
    
    known_lactylwriters = {
        'EP300': ['HIF1A', 'STAT3', 'P53'],  # 写入酶候选
        'CREBBP': ['HIF1A', 'STAT3'],
    }
    
    # 添加激酶节点
    active_kinases = []
    if ksea_results is not None:
        active_kinases = ksea_results[ksea_results['Z_Score'] > 1.5]['Kinase'].tolist()
    
    for kinase in known_kinases.keys():
        substrates = known_kinases[kinase]
        
        # 判断激酶活性
        is_active = kinase in active_kinases
        
        G.add_node(
            kinase,
            type='kinase',
            is_active=is_active,
            n_substrates=len(substrates)
        )
        
        # 添加底物边
        for sub in substrates:
            if sub in assoc_df['ProteinID'].values:
                G.add_edge(
                    kinase, sub,
                    relation='phosphorylates',
                    weight=1.0
                )
    
    # 添加乳酸化酶节点
    for enzyme in set(list(known_lactylases.keys()) + list(known_lactylwriters.keys())):
        G.add_node(
            enzyme,
            type='lactylase',
            is_active=True,
            n_substrates=len(known_lactylwriters.get(enzyme, [])) + len(known_lactylases.get(enzyme, []))
        )
    
    # 保存网络
    nx.write_graphml(G, 'results/3_enzyme_substrate_network.graphml')
    
    print(f"\n=== 酶-底物网络 ===")
    print(f"节点数: {G.number_of_nodes()}")
    print(f"边数: {G.number_of_edges()}")
    
    return G


def visualize_enzyme_network(G, assoc_df):
    """可视化酶-底物网络"""
    if G.number_of_nodes() == 0:
        return None
    
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # 布局
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # 节点分类
    kinases = [n for n, d in G.nodes(data=True) if d.get('type') == 'kinase']
    lactylases = [n for n, d in G.nodes(data=True) if d.get('type') == 'lactylase']
    substrates = [n for n in G.nodes() if n not in kinases and n not in lactylases]
    
    # 绘制节点
    nx.draw_networkx_nodes(G, pos, nodelist=kinases,
                          node_color='#e74c3c', node_size=1000, alpha=0.8,
                          node_shape='s')
    nx.draw_networkx_nodes(G, pos, nodelist=lactylases,
                          node_color='#9b59b6', node_size=800, alpha=0.8,
                          node_shape='s')
    nx.draw_networkx_nodes(G, pos, nodelist=substrates,
                          node_color='#3498db', node_size=300, alpha=0.6)
    
    # 绘制边
    nx.draw_networkx_edges(G, pos, alpha=0.3, arrows=True, arrowsize=10)
    
    # 标签
    nx.draw_networkx_labels(G, pos, font_size=6)
    
    plt.title('Enzyme-Substrate Regulatory Network')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label='Kinase'),
        Patch(facecolor='#9b59b6', label='Lactylase'),
        Patch(facecolor='#3498db', label='Substrate')
    ]
    plt.legend(handles=legend_elements, loc='upper left')
    
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/3_enzyme_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/3_enzyme_network.png'


# ============== 4. 激酶/乳酸化酶活性推断 ==============
def enzyme_activity_inference(assoc_df, ksea_results=None):
    """
    综合激酶和乳酸化酶活性推断
    
    基于底物磷酸化/乳酸化变化推断酶活性
    """
    activity_results = []
    
    # 激酶活性 (基于KSEA结果)
    if ksea_results is not None:
        for _, row in ksea_results.iterrows():
            activity_results.append({
                'Enzyme': row['Kinase'],
                'Enzyme_Type': 'Kinase',
                'N_Substrates': row['N_Substrates'],
                'Mean_Substrate_FC': row['Mean_Substrate_FC'],
                'Z_Score': row['Z_Score'],
                'P_Value': row['P_Value'],
                'Activity_Status': row['Activity']
            })
    
    # 乳酸化酶活性推断 (简化版)
    # 由于乳酸化酶尚未完全确定，这里基于已知底物变化进行推断
    lactylase_substrates = {
        'EP300': ['HIF1A', 'STAT3', 'P53'],
        'CREBBP': ['HIF1A', 'STAT3'],
        'HDAC1': ['P53', 'STAT3'],
        'SIRT1': ['P53', 'HIF1A'],
    }
    
    for enzyme, substrates in lactylase_substrates.items():
        enzyme_substrates = assoc_df[assoc_df['ProteinID'].isin(substrates)]
        
        if len(enzyme_substrates) > 0:
            mean_fc = enzyme_substrates['Lactylation_log2FC'].mean()
            n_sig = (enzyme_substrates['N_Lactylation_Significant'] > 0).sum()
            
            activity_results.append({
                'Enzyme': enzyme,
                'Enzyme_Type': 'Lactylase',
                'N_Substrates': len(substrates),
                'Mean_Substrate_FC': mean_fc,
                'Z_Score': mean_fc * np.sqrt(n_sig) if n_sig > 0 else 0,
                'P_Value': 0.05 if n_sig > 0 else 1,
                'Activity_Status': 'Increased' if mean_fc > 0.5 else ('Decreased' if mean_fc < -0.5 else 'Unchanged')
            })
    
    activity_df = pd.DataFrame(activity_results)
    activity_df = activity_df.sort_values('Z_Score', ascending=False)
    activity_df.to_csv('results/4_enzyme_activity_inference.csv', index=False)
    
    print(f"\n=== 酶活性推断 ===")
    if len(activity_df) > 0:
        print(activity_df[activity_df['Enzyme_Type'] == 'Kinase']['Activity_Status'].value_counts())
        print(activity_df[activity_df['Enzyme_Type'] == 'Lactylase']['Activity_Status'].value_counts())
    
    return activity_df


def plot_enzyme_activity(activity_df):
    """可视化酶活性推断结果"""
    if activity_df is None or len(activity_df) == 0:
        return None
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 激酶活性
    ax = axes[0]
    kinases = activity_df[activity_df['Enzyme_Type'] == 'Kinase'].sort_values('Z_Score')
    
    if len(kinases) > 0:
        colors = ['#e74c3c' if z > 0 else '#3498db' for z in kinases['Z_Score']]
        ax.barh(range(len(kinases)), kinases['Z_Score'], color=colors)
        ax.set_yticks(range(len(kinases)))
        ax.set_yticklabels(kinases['Enzyme'], fontsize=8)
        ax.axvline(x=2, color='red', linestyle='--', alpha=0.5)
        ax.axvline(x=-2, color='blue', linestyle='--', alpha=0.5)
        ax.set_xlabel('Z-Score')
        ax.set_title('Kinase Activity Inference')
    
    # 乳酸化酶活性
    ax = axes[1]
    lactylases = activity_df[activity_df['Enzyme_Type'] == 'Lactylase'].sort_values('Z_Score')
    
    if len(lactylases) > 0:
        colors = ['#e74c3c' if z > 0 else '#3498db' for z in lactylases['Z_Score']]
        ax.barh(range(len(lactylases)), lactylases['Z_Score'], color=colors)
        ax.set_yticks(range(len(lactylases)))
        ax.set_yticklabels(lactylases['Enzyme'], fontsize=8)
        ax.axvline(x=1, color='red', linestyle='--', alpha=0.5)
        ax.axvline(x=-1, color='blue', linestyle='--', alpha=0.5)
        ax.set_xlabel('Z-Score')
        ax.set_title('Lactylase Activity Inference')
    
    plt.tight_layout()
    plt.savefig('results/4_enzyme_activity.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/4_enzyme_activity.png'


# ============== 5. 综合串扰网络可视化 ==============
def comprehensive_crosstalk_network(assoc_df, crosstalk_df, G):
    """
    创建综合串扰网络图
    
    展示:
    - 蛋白节点 (按修饰模式着色)
    - 修饰位点 (磷酸化/乳酸化)
    - 酶-底物边
    - 串扰边
    """
    fig, ax = plt.subplots(figsize=(18, 14))
    
    # 创建子网络
    # 取关键蛋白
    key_proteins = assoc_df[
        (assoc_df['Total_Modifications'] >= 2) & 
        ((assoc_df['N_Phospho_Significant'] > 0) | (assoc_df['N_Lactylation_Significant'] > 0))
    ]['ProteinID'].tolist()[:50]
    
    H = G.subgraph(key_proteins).copy()
    
    if len(H.nodes()) == 0:
        print("没有足够的节点构建网络")
        return None
    
    # 布局
    pos = nx.spring_layout(H, k=1.5, iterations=50, seed=42)
    
    # 节点颜色按修饰模式
    patterns = assoc_df.set_index('ProteinID')['Modification_Pattern'].to_dict()
    
    color_map = {
        'Triple-Active': '#e74c3c',
        'Triple-Suppressed': '#3498db',
        'Phospho-Specific': '#f39c12',
        'Lact-Specific': '#9b59b6',
        'Modification-Opposite': '#e67e22',
        'No-Modification': '#95a5a6'
    }
    
    node_colors = [color_map.get(patterns.get(n, 'No-Modification'), '#95a5a6') for n in H.nodes()]
    
    # 节点大小
    degrees = dict(H.degree())
    node_sizes = [max(200, degrees[n] * 100) for n in H.nodes()]
    
    nx.draw_networkx_nodes(H, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.8)
    nx.draw_networkx_edges(H, pos, alpha=0.3, arrows=True)
    nx.draw_networkx_labels(H, pos, font_size=6)
    
    plt.title('Comprehensive Modification Crosstalk Network')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=c, label=l) for l, c in color_map.items()
    ]
    plt.legend(handles=legend_elements, loc='upper left', fontsize=8)
    
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('results/5_comprehensive_crosstalk_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return 'results/5_comprehensive_crosstalk_network.png'


# ============== 主流程 ==============
def main():
    """主分析流程"""
    print("=" * 60)
    print("蛋白质组 + 磷酸化组 + 乳酸化组 三组学修饰串扰分析")
    print("=" * 60)
    
    # 1. 加载数据
    print("\n[1/5] 加载数据...")
    proteome_data, phospho_sites, lactylation_sites, sample_info = load_triple_modification_data()
    
    # 2. 修饰位点-蛋白关联分类
    print("\n[2/5] 修饰位点-蛋白表达关联分类...")
    assoc_df = modification_protein_association(proteome_data, phospho_sites, lactylation_sites)
    assoc_df = classify_modification_patterns(assoc_df)
    plot_modification_patterns(assoc_df)
    
    # 3. 位点级串扰分析
    print("\n[3/5] 磷酸化-乳酸化位点级串扰分析...")
    crosstalk_df = site_level_crosstalk(phospho_sites, lactylation_sites)
    plot_crosstalk_analysis(crosstalk_df)
    
    # 4. 酶-底物网络构建
    print("\n[4/5] 酶-底物网络构建...")
    ksea_results = None  # 如果有KSEA结果可以传入
    G = enzyme_substrate_network(assoc_df, phospho_sites, lactylation_sites, ksea_results)
    visualize_enzyme_network(G, assoc_df)
    
    # 5. 激酶/乳酸化酶活性推断
    print("\n[5/5] 激酶/乳酸化酶活性推断...")
    activity_df = enzyme_activity_inference(assoc_df, ksea_results)
    plot_enzyme_activity(activity_df)
    
    # 6. 综合串扰网络
    print("\n[Extra] 综合串扰网络可视化...")
    comprehensive_crosstalk_network(assoc_df, crosstalk_df, G)
    
    print("\n" + "=" * 60)
    print("分析完成! 结果保存在 results/ 目录")
    print("=" * 60)


if __name__ == '__main__':
    main()
