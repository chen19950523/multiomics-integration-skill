# 框架3：蛋白质组 + 磷酸化 联合分析

## 背景

磷酸化（Phosphorylation）是最重要的PTM之一，调控几乎所有信号通路。磷酸化蛋白组学（Phosphoproteomics）可鉴定数千个磷酸化位点的动态变化。

**核心应用：**
- 信号通路激活状态
- 激酶-底物调控网络
- 药物作用机制（激酶抑制剂）
- 细胞响应外界刺激

## 分析流程

### Step 1: 磷酸化蛋白组分析

**富集方法：** TiO2、IMAC、Fe-NTA
**定量软件：** MaxQuant (PhosphoRS)、Pirat、PhosphoPlex

```python
phospho_df = pd.DataFrame({
    'ProteinID': ['P00533', 'P04637'],
    'Gene': ['EGFR', 'TP53'],
    'Site': ['S1068', 'S15'],
    'Sequence': ['SASLDNPD', 'SRAHSSEL'],
    'log2FC': [2.1, -1.5],
    'P-value': [0.001, 0.003],
    'Localization_prob': [0.95, 0.88],  # 定位概率
    'Regulation': ['Up', 'Down']
})
```

### Step 2: 全局蛋白组关联

```python
def phospho_protein_association(phospho_df, proteome_df):
    """
    磷酸化位点数据与全局蛋白表达关联
    """
    merged = phospho_df.merge(
        proteome_df[['ProteinID', 'Gene', 'log2FC_protein', 'P-value_protein']],
        on=['ProteinID', 'Gene'],
        how='left',
        suffixes=('_phospho', '_protein')
    )
    
    # 标记磷酸化-only蛋白（蛋白不变但磷酸化变化）
    merged['Phospho_only'] = (
        abs(merged['log2FC_protein'].fillna(0)) < 0.5
    ) & (abs(merged['log2FC_phospho']) > 1)
    
    return merged
```

### Step 3: 磷酸化-only分析（关键！）

**磷酸化-only蛋白 = 信号转导的核心节点（蛋白不变但磷酸化激活）**

```python
def phospho_only_analysis(merged_df):
    """
    筛选磷酸化-only蛋白
    这些是纯粹的信号转导分子
    """
    phospho_only = merged_df[merged_df['Phospho_only'] == True]
    
    print(f"Total phospho sites: {len(merged_df)}")
    print(f"Phospho-only sites (protein unchanged): {len(phospho_only)}")
    print(f"Phospho-only genes: {phospho_only['Gene'].nunique()}")
    
    return phospho_only
```

### Step 4: 激酶活性推断（KSEA）

```python
def ksea_analysis(phospho_df, kinase_substrate_db, fc_thresh=1, p_thresh=0.05):
    """
    Kinase Substrate Enrichment Analysis (KSEA)
    基于磷酸化底物丰度推断激酶活性
    
    数据库：PhosphoSitePlus, PSP
    """
    # 激酶-底物对应关系
    # kinase_substrate_db = { kinase: [substrate_genes] }
    
    results = []
    for kinase, substrates in kinase_substrate_db.items():
        # 获取该激酶的底物在数据中的磷酸化变化
        kinase_sites = phospho_df[
            (phospho_df['Gene'].isin(substrates)) & 
            (abs(phospho_df['log2FC']) > fc_thresh)
        ]
        
        if len(kinase_sites) >= 3:  # 至少3个底物
            # 计算平均log2FC（激酶活性 proxy）
            mean_fc = kinase_sites['log2FC'].mean()
            # 超几何检验评估富集显著性
            from scipy.stats import hypergeometric
            background = len(phospho_df)
            substrates_in_data = len(kinase_sites)
            total_substrates = len(substrates)
            regulated = len(phospho_df[abs(phospho_df['log2FC']) > fc_thresh])
            
            pval = hypergeometric.sf(
                substrates_in_data - 1, 
                background, 
                total_substrates, 
                regulated
            )
            
            results.append({
                'Kinase': kinase,
                'Substrates_measured': substrates_in_data,
                'Mean_log2FC': mean_fc,
                'P-value': pval,
                'Activity': 'Active' if mean_fc > 0 else 'Inhibited'
            })
    
    return pd.DataFrame(results).sort_values('P-value')
```

### Step 5: 信号通路网络构建

```python
import networkx as nx

def build_signaling_network(phospho_df, kinase_activity_df):
    """
    构建"激酶 → 底物蛋白 → 下游效应"网络
    """
    G = nx.DiGraph()
    
    # 添加激活的激酶
    active_kinases = kinase_activity_df[
        (kinase_activity_df['P-value'] < 0.05) & 
        (abs(kinase_activity_df['Mean_log2FC']) > 1)
    ]
    
    for _, kin in active_kinases.iterrows():
        G.add_node(kin['Kinase'], type='Kinase', activity=kin['Activity'])
        
        # 添加底物（简化版）
        substrates = phospho_df[phospho_df['Gene'].isin(get_kinase_substrates(kin['Kinase']))]
        for _, sub in substrates.iterrows():
            G.add_edge(kin['Kinase'], sub['Gene'], type='phosphorylation')
    
    return G

def visualize_signaling_network(G, output='signaling_network.png'):
    plt.figure(figsize=(16, 12))
    pos = nx.spring_layout(G, k=2)
    
    # 按类型着色
    kinases = [n for n, d in G.nodes(data=True) if d.get('type') == 'Kinase']
    substrates = list(set(G.nodes()) - set(kinases))
    
    nx.draw_networkx_nodes(G, pos, nodelist=kinases, node_color='red', node_size=500)
    nx.draw_networkx_nodes(G, pos, nodelist=substrates, node_color='lightblue', node_size=200)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    nx.draw_networkx_labels(G, pos, font_size=6)
    
    plt.savefig(output, dpi=300)
```

### Step 6: 关键信号通路分析

```python
SIGNALING_PATHWAYS = {
    'PI3K/AKT': ['AKT1', 'AKT2', 'AKT3', 'MTOR', 'GSK3B', 'PRKCA'],
    'MAPK/ERK': ['MAPK1', 'MAPK3', 'RAF1', 'MEK1', 'MEK2', 'ERK1', 'ERK2'],
    'JAK/STAT': ['JAK1', 'JAK2', 'STAT1', 'STAT3', 'STAT5A', 'STAT5B'],
    'p53': ['TP53', 'MDM2', 'CHEK1', 'CHEK2', 'ATM', 'ATR'],
    'Cell_Cycle': ['CDK1', 'CDK2', 'CDK4', 'CDK6', 'CCNB1', 'CCND1'],
    'Apoptosis': ['BCL2', 'BAX', 'CASP3', 'CASP9', 'PARP1'],
    'TGF_beta': ['TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3', 'SMAD4'],
    'WNT': ['CTNNB1', 'APC', 'AXIN1', 'GSK3B', 'LEF1'],
}

def pathway_activity_inference(phospho_df, pathways):
    """
    基于磷酸化位点丰度推断信号通路激活状态
    """
    results = []
    for pathway_name, genes in pathways.items():
        pathway_sites = phospho_df[phospho_df['Gene'].isin(genes)]
        
        if len(pathway_sites) > 0:
            results.append({
                'Pathway': pathway_name,
                'Phospho_sites': len(pathway_sites),
                'Mean_log2FC': pathway_sites['log2FC'].mean(),
                'Up_sites': len(pathway_sites[pathway_sites['log2FC'] > 0]),
                'Down_sites': len(pathway_sites[pathway_sites['log2FC'] < 0]),
                'Activity': 'Activated' if pathway_sites['log2FC'].mean() > 0 else 'Inhibited'
            })
    
    return pd.DataFrame(results)
```

## 激酶-底物数据库

| 数据库 | 描述 | 链接 |
|--------|------|------|
| PhosphoSitePlus | 最大的PTM数据库 | phosphosite.org |
| KinomeXpress | 激酶-底物关系 | kinomexpress.com |
| dbPTM | PTM数据库 | dbPTM.bmcnctu.org |
| UniProt | 蛋白注释 | uniprot.org |

## 典型应用场景

| 研究场景 | 关键发现 | 典型通路 |
|---------|---------|---------|
| 药物处理 | 激酶抑制剂效果 | p-ERK, p-AKT变化 |
| 时间序列 | 信号通路动态 | MAPK级联 |
| 耐药机制 | 旁路激活 | bypass signaling |
| 刺激响应 | 外界信号转导 | calcium signaling |

## 可视化

```python
def visualize_phosphorylation_results(kinase_activity_df, pathway_activity_df):
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # 1. 激酶活性条形图
    ax1 = axes[0, 0]
    kinase_activity_df = kinase_activity_df.sort_values('Mean_log2FC')
    colors = ['red' if x > 0 else 'blue' for x in kinase_activity_df['Mean_log2FC']]
    ax1.barh(kinase_activity_df['Kinase'], kinase_activity_df['Mean_log2FC'], color=colors)
    ax1.set_xlabel('Mean log2FC of substrates')
    ax1.set_title('Kinase Activity Inference (KSEA)')
    ax1.axvline(0, color='black', linestyle='-', alpha=0.3)
    
    # 2. 通路激活热图
    ax2 = axes[0, 1]
    pathway_matrix = pathway_activity_df.set_index('Pathway')[['Up_sites', 'Down_sites']]
    sns.heatmap(pathway_matrix, annot=True, fmt='d', cmap='RdBu_r', ax=ax2)
    ax2.set_title('Signaling Pathway Activity')
    
    # 3. 磷酸化位点分布
    ax3 = axes[1, 0]
    ax3.hist(phospho_df['log2FC'], bins=50, edgecolor='black', alpha=0.7)
    ax3.axvline(0, color='red', linestyle='--')
    ax3.set_xlabel('log2FC (Phosphorylation)')
    ax3.set_ylabel('Count')
    ax3.set_title('Distribution of Phosphorylation Changes')
    
    # 4. 蛋白表达 vs 磷酸化
    ax4 = axes[1, 1]
    merged = phospho_df.merge(proteome_df, on='Gene', suffixes=('_phospho', '_protein'))
    sns.scatterplot(data=merged, x='log2FC_protein', y='log2FC_phospho', alpha=0.3, ax=ax4)
    ax4.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax4.axvline(0, color='red', linestyle='--', alpha=0.3)
    ax4.set_xlabel('Protein log2FC')
    ax4.set_ylabel('Phosphorylation log2FC')
    ax4.set_title('Protein Expression vs Phosphorylation')
    
    plt.tight_layout()
    plt.savefig('phosphorylation_analysis.png', dpi=300)
```

## 参考文献

- Bodenmiller et al., Nat Methods 2007 (Phosphoproteomics)
- Craig et al., MCP 2014 (KSEA)
- Needham et al., Nat Rev Cancer 2019 (Kinase signaling networks)
