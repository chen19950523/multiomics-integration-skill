# 框架5：三组学修饰串扰分析（蛋白组 + 磷酸化 + 乳酸化）

## 背景

磷酸化和乳酸化是两种重要的PTM，它们可能：
1. **协同调控**——同一蛋白被双重激活
2. **竞争关系**——共享相同或邻近位点
3. **单向调控**——磷酸化影响乳酸化或反过来
4. **独立调控**——各自独立作用于不同位点

三组学联合分析是当前PTM研究的前沿，适合发表高水平文章（STTT、NC级别）。

## 分析流程

### Step 1: 三组学数据整合

```python
import pandas as pd
import numpy as np

def load_triple_omics(proteome_file, phospho_file, lactyl_file):
    """
    加载三组学数据
    """
    proteome = pd.read_csv(proteome_file, sep='\t')
    phospho = pd.read_csv(phospho_file, sep='\t')
    lactyl = pd.read_csv(lactyl_file, sep='\t')
    
    return proteome, phospho, lactyl

# 数据格式
proteome_df = pd.DataFrame({
    'ProteinID': ['P12345', 'P67890'],
    'Gene': ['GAPDH', 'AKT1'],
    'log2FC': [1.2, -0.5],
    'P-value': [0.001, 0.02]
})

phospho_df = pd.DataFrame({
    'ProteinID': ['P12345', 'P12345', 'P67890'],
    'Gene': ['GAPDH', 'GAPDH', 'AKT1'],
    'Site': ['S200', 'T201', 'S473'],
    'log2FC': [2.1, 0.8, -1.5],
    'P-value': [0.001, 0.05, 0.002]
})

lactyl_df = pd.DataFrame({
    'ProteinID': ['P12345', 'P67890'],
    'Gene': ['GAPDH', 'AKT1'],
    'Site': ['K184', 'K310'],
    'log2FC': [1.8, -0.3],
    'P-value': [0.003, 0.1]
})
```

### Step 2: 蛋白-修饰表达关联分类

```python
def classify_modification_patterns(proteome_df, phospho_df, lactyl_df, fc_thresh=0.58):
    """
    对蛋白-磷酸化-乳酸化进行9类分类
    
    基于：蛋白表达变化、磷酸化变化、乳酸化变化
    """
    # 合并数据
    merged = proteome_df[['ProteinID', 'Gene', 'log2FC']].rename(
        columns={'log2FC': 'Protein_FC'}
    )
    
    # 磷酸化汇总（按蛋白，取最大变化位点）
    phospho_agg = phospho_df.groupby(['ProteinID', 'Gene']).agg({
        'log2FC': ['max', 'mean'],
        'Site': 'count'
    }).reset_index()
    phospho_agg.columns = ['ProteinID', 'Gene', 'Phospho_FC_max', 'Phospho_FC_mean', 'Phospho_site_count']
    
    # 乳酸化汇总
    lactyl_agg = lactyl_df.groupby(['ProteinID', 'Gene']).agg({
        'log2FC': ['max', 'mean'],
        'Site': 'count'
    }).reset_index()
    lactyl_agg.columns = ['ProteinID', 'Gene', 'Lactyl_FC_max', 'Lactyl_FC_mean', 'Lactyl_site_count']
    
    # 三者合并
    triple = merged.merge(phospho_agg, on=['ProteinID', 'Gene'], how='outer')
    triple = triple.merge(lactyl_agg, on=['ProteinID', 'Gene'], how='outer')
    triple = triple.fillna(0)
    
    # 分类
    def get_pattern(row):
        p = row['Protein_FC']
        ph = row['Phospho_FC_max']
        la = row['Lactyl_FC_max']
        
        # 简化为协同/竞争/独立分类
        if abs(ph) > fc_thresh and abs(la) > fc_thresh:
            if (ph > 0 and la > 0) or (ph < 0 and la < 0):
                return 'Synergistic'
            else:
                return 'Competitive'
        elif abs(ph) > fc_thresh and abs(la) <= fc_thresh:
            return 'Phospho_only'
        elif abs(la) > fc_thresh and abs(ph) <= fc_thresh:
            return 'Lactyl_only'
        elif abs(ph) > fc_thresh and abs(la) > fc_thresh:
            return 'Dual_modification'
        else:
            return 'Unchanged'
    
    triple['Pattern'] = triple.apply(get_pattern, axis=1)
    
    return triple

def visualize_pattern_distribution(triple_df):
    """
    可视化修饰模式分布
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    pattern_counts = triple_df['Pattern'].value_counts()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 饼图
    axes[0].pie(pattern_counts.values, labels=pattern_counts.index, autopct='%1.1f%%')
    axes[0].set_title('Distribution of Modification Patterns')
    
    # 条形图
    axes[1].bar(pattern_counts.index, pattern_counts.values)
    axes[1].set_ylabel('Count')
    axes[1].set_title('Modification Pattern Counts')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig('modification_patterns.png', dpi=300)
```

### Step 3: 位点级串扰分析

```python
def site_level_crosstalk(phospho_df, lactyl_df, distance_threshold=10):
    """
    位点级串扰分析
    检测磷酸化和乳酸化位点是否在空间上接近
    
    同一蛋白上：磷酸化位点 vs 乳酸化位点
    如果距离 < threshold AA，认为可能存在串扰
    """
    # 找出同时有磷酸化和乳酸化数据的蛋白
    common_proteins = set(phospho_df['ProteinID']) & set(lactyl_df['ProteinID'])
    
    crosstalk_pairs = []
    for prot in common_proteins:
        ph_sites = phospho_df[phospho_df['ProteinID'] == prot]
        la_sites = lactyl_df[lactyl_df['ProteinID'] == prot]
        
        for _, ph in ph_sites.iterrows():
            for _, la in la_sites.iterrows():
                # 计算位点距离（简化版，需要实际蛋白序列）
                # 这里用氨基酸位置差估算
                pos_diff = abs(ph['Position'] - la['Position'])
                
                if pos_diff <= distance_threshold:
                    crosstalk_pairs.append({
                        'Protein': prot,
                        'Phospho_site': ph['Site'],
                        'Lactyl_site': la['Site'],
                        'Position_diff': pos_diff,
                        'Phospho_FC': ph['log2FC'],
                        'Lactyl_FC': la['log2FC'],
                        'Potential_crosstalk': True
                    })
    
    return pd.DataFrame(crosstalk_pairs)
```

### Step 4: 酶-底物调控网络

```python
import networkx as nx

def build_modification_network(triple_df, kinase_activity, lactylase_activity):
    """
    构建"激酶-磷酸化-乳酸化酶"联合调控网络
    """
    G = nx.DiGraph()
    
    # 添加蛋白节点
    for _, row in triple_df.iterrows():
        if row['Pattern'] in ['Synergistic', 'Competitive', 'Dual_modification']:
            G.add_node(
                row['Gene'], 
                type='Hub_protein',
                pattern=row['Pattern'],
                protein_fc=row['Protein_FC'],
                phospho_fc=row['Phospho_FC_max'],
                lactyl_fc=row['Lactyl_FC_max']
            )
    
    # 添加激酶边
    for kin, data in kinase_activity.items():
        if data['Activity'] == 'Active':
            G.add_node(kin, type='Kinase')
            # 连接到磷酸化底物
            for substrate in data['substrates']:
                if substrate in G.nodes():
                    G.add_edge(kin, substrate, type='Kinase-Phospho')
    
    # 添加乳酸化酶边
    for enzyme, data in lactylase_activity.items():
        if data['Activity'] == 'Active':
            G.add_node(enzyme, type='Lactylase')
            for substrate in data['substrates']:
                if substrate in G.nodes():
                    G.add_edge(enzyme, substrate, type='Lactylase-Lactyl')
    
    return G

def visualize_modification_network(G, output='modification_network.png'):
    """
    可视化修饰调控网络
    """
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(16, 12))
    pos = nx.spring_layout(G, k=2)
    
    # 按类型着色
    node_colors = []
    for node in G.nodes():
        if G.nodes[node].get('type') == 'Kinase':
            node_colors.append('red')
        elif G.nodes[node].get('type') == 'Lactylase':
            node_colors.append('green')
        else:
            node_colors.append('lightblue')
    
    nx.draw_networkx(G, pos, node_color=node_colors, node_size=500, font_size=6, alpha=0.8)
    
    plt.savefig(output, dpi=300)
```

### Step 5: 通路共调控分析

```python
PATHWAY_MODIFICATION_PATTERNS = {
    'Glycolysis': {
        'genes': ['HK1', 'HK2', 'PGK1', 'PKM2', 'LDHA', 'ENO1'],
        'expected': 'Phospho+Lactyl synergistic (metabolic adaptation)'
    },
    'mTORC1': {
        'genes': ['MTOR', 'AKT1S1', 'ULK1', 'RPS6KB1'],
        'expected': 'Phospho dominant (signaling activation)'
    },
    'Mitochondrial_function': {
        'genes': ['MT-ATP6', 'NDUFA9', 'SDHA', 'COX5B'],
        'expected': 'Lactyl dominant (metabolic sensing)'
    },
    'Epigenetic_regulation': {
        'genes': ['H1', 'H2A', 'H2B', 'H3', 'H4', 'KAT2A'],
        'expected': 'Lactyl dominant (lactate-histone crosstalk)'
    },
    'Inflammatory_signaling': {
        'genes': ['NFKB1', 'RELA', 'STAT3', 'JUN', 'FOS'],
        'expected': 'Phospho+Lactyl dual (signal+metabolism)'
    }
}

def pathway_modification_analysis(triple_df, pathway_definitions):
    """
    分析每个通路中的修饰调控模式
    """
    results = []
    for pathway_name, pathway_info in pathway_definitions.items():
        genes = pathway_info['genes']
        pathway_data = triple_df[triple_df['Gene'].isin(genes)]
        
        if len(pathway_data) > 0:
            results.append({
                'Pathway': pathway_name,
                'Genes_detected': len(pathway_data),
                'Pattern_distribution': pathway_data['Pattern'].value_counts().to_dict(),
                'Mean_Phospho_FC': pathway_data['Phospho_FC_max'].mean(),
                'Mean_Lactyl_FC': pathway_data['Lactyl_FC_max'].mean(),
                'Interpretation': pathway_info['expected']
            })
    
    return pd.DataFrame(results)
```

### Step 6: 高阶串扰矩阵

```python
def crosstalk_matrix(triple_df, top_proteins):
    """
    创建磷酸化-乳酸化串扰矩阵
    用于热图展示
    """
    matrix_data = []
    for prot in top_proteins:
        prot_data = triple_df[triple_df['Gene'] == prot].iloc[0]
        matrix_data.append({
            'Protein': prot,
            'Phospho': prot_data['Phospho_FC_max'],
            'Lactyl': prot_data['Lactyl_FC_max'],
            'Protein': prot_data['Protein_FC']
        })
    
    matrix_df = pd.DataFrame(matrix_data).set_index('Protein')
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 10))
    sns.heatmap(
        matrix_df[['Phospho', 'Lactyl']], 
        annot=True, 
        fmt='.2f', 
        cmap='RdBu_r', 
        center=0,
        ax=ax
    )
    ax.set_title('Phospho-Lactyl Crosstalk Matrix')
    plt.savefig('crosstalk_matrix.png', dpi=300)
    
    return matrix_df
```

## 串扰类型总结

| 串扰类型 | 磷酸化 | 乳酸化 | 机制解释 |
|---------|--------|--------|---------|
| **协同激活** | ↑ | ↑ | 双重激活，信号+代谢双重放大 |
| **协同抑制** | ↓ | ↓ | 双重抑制 |
| **竞争互斥** | ↑ | ↓ | 磷酸化激活，乳酸化抑制（反式激活）|
| **磷酸化主导** | ↑ | - | 纯信号转导 |
| **乳酸化主导** | - | ↑ | 代谢感知调控 |
| **独立调控** | ↑ | ↑ (不同位点) | 同一蛋白，不同功能域 |

## 高分文章逻辑模板

```
标题：Phosphorylation-lactylation crosstalk reveals novel mechanism of metabolic reprogramming in XX cancer

1. 临床问题
   - XX癌症中代谢重编程的机制不清

2. 三组学筛选
   - 蛋白组+磷酸化组+乳酸化组
   - 鉴定"双修饰蛋白"（协同/竞争模式）

3. 机制探索
   - 激酶-底物、乳酸化酶-底物网络
   - 关键位点突变体验证
   - 细胞功能实验

4. 临床关联
   - 双修饰signature的预后价值

5. 总结
   - 首次揭示磷酸化-乳酸化串扰在XX中的机制
```

## 参考文献

- Zhang et al., Nature 2019 (Lactylation discovery)
- Liu et al., Nat Cell Biol 2022 (Lactylation in stemness)
- Wang et al., Cell Metab 2021 (Lactylation in tumor immunity)
- Pan et al., STTT 2023 (Phospho-lactyl crosstalk in cancer)
