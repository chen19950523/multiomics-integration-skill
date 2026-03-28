# 框架1：转录组 + 蛋白质组 联合分析

## 核心逻辑

mRNA与蛋白质的表达相关性，探究"转录后调控"机制。
同一基因的mRNA与蛋白表达不一致 = 存在翻译后调控或降解。

## 分析流程

### Step 1: 数据预处理

**转录组**
```python
# 使用Salmon/kallisto定量，或者从他公司获取FPKM/TPM矩阵
# 差异分析：DESeq2/edgeR
import pandas as pd
from scipy import stats

def differential_expression_rna(count_matrix, control, treatment):
    log2fc = count_matrix.mean(axis=1).apply(lambda x: np.log2(x+1))
    # 简化示例
    return log2fc, pvalue
```

**蛋白组**
```python
# 蛋白组通常由公司提供差异分析结果
# 格式：ProteinID, log2FC, P-value, Regulation
proteome_df = pd.read_csv('proteome_diff.tsv', sep='\t')
```

### Step 2: 差异重叠分析（4象限分类）

| 象限 | mRNA | Protein | 含义 |
|------|------|---------|------|
| Q1 | ↑ | ↑ | 转录激活（一致） |
| Q2 | ↑ | ↓ | 转录后抑制/降解 |
| Q3 | ↓ | ↑ | 翻译激活/稳定 |
| Q4 | ↓ | ↓ | 转录抑制（一致） |

```python
import matplotlib.pyplot as plt
import numpy as np

def quadrant_analysis(deg_df, dep_df, fc_thresh=1, pval_thresh=0.05):
    merged = deg_df.merge(dep_df, on='gene', suffixes=('_mRNA', '_protein'))
    
    conditions = [
        (merged['log2FC_mRNA'] > fc_thresh) & (merged['log2FC_protein'] > fc_thresh),  # Q1
        (merged['log2FC_mRNA'] > fc_thresh) & (merged['log2FC_protein'] < -fc_thresh),  # Q2
        (merged['log2FC_mRNA'] < -fc_thresh) & (merged['log2FC_protein'] > fc_thresh),  # Q3
        (merged['log2FC_mRNA'] < -fc_thresh) & (merged['log2FC_protein'] < -fc_thresh),  # Q4
    ]
    labels = ['Q1: Co-up', 'Q2: mRNA-up/Protein-down', 'Q3: mRNA-down/Protein-up', 'Q4: Co-down']
    
    for condition, label in zip(conditions, labels):
        genes = merged[condition]['gene'].tolist()
        print(f"{label}: {len(genes)} genes")
    
    return merged
```

### Step 3: mRNA-蛋白相关性分析

```python
import seaborn as sns
from scipy.stats import spearmanr, pearsonr

def correlation_analysis(rna_df, protein_df):
    merged = rna_df.merge(protein_df, on='gene')
    
    # Spearman相关性（更稳健，不假设线性关系）
    rho, pval = spearmanr(merged['rna_expression'], merged['protein_expression'])
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='rna_expression', y='protein_expression', data=merged, alpha=0.5)
    plt.xlabel('mRNA Expression (TPM)')
    plt.ylabel('Protein Expression (iBAQ)')
    plt.title(f'Spearman Correlation: ρ={rho:.3f}, p={pval:.2e}')
    plt.savefig('mrna_protein_correlation.png', dpi=300)
    
    return rho, pval
```

### Step 4: 翻译效率分析

```python
def translation_efficiency(rna_df, protein_df):
    """
    翻译效率 = 蛋白表达量 / mRNA表达量
    高TE：蛋白翻译增强或蛋白更稳定
    低TE：翻译抑制或蛋白降解加速
    """
    merged = rna_df.merge(protein_df, on='gene')
    merged['translation_efficiency'] = merged['protein_expression'] / (merged['rna_expression'] + 1)
    
    # 差异翻译效率基因
    merged['TE_zscore'] = (merged['translation_efficiency'] - merged['translation_efficiency'].mean()) / merged['translation_efficiency'].std()
    
    high_te = merged[merged['TE_zscore'] > 2]['gene'].tolist()
    low_te = merged[merged['TE_zscore'] < -2]['gene'].tolist()
    
    print(f"High TE genes (z>2): {len(high_te)}")
    print(f"Low TE genes (z<-2): {len(low_te)}")
    
    return merged, high_te, low_te
```

### Step 5: GO/KEGG联合富集

```python
from scipy.stats import hypergeometric

def enrichment_analysis(gene_list, background, go_annotations, kegg_pathways):
    """
    对基因列表进行GO/KEGG富集分析
    使用超几何检验
    """
    results = []
    for term, genes in go_annotations.items():
        overlap = len(set(gene_list) & set(genes))
        if overlap > 0:
            pval = hypergeometric.sf(overlap-1, len(background), len(genes), len(gene_list))
            results.append({'Term': term, 'Overlap': overlap, 'P-value': pval})
    
    return pd.DataFrame(results).sort_values('P-value')
```

### Step 6: PPI网络构建

```python
import networkx as nx

def build_ppi_network(genes, string_db_path='string_interactions.tsv'):
    """
    使用STRING数据库构建PPI网络
    """
    string_df = pd.read_csv(string_db_path, sep='\t')
    interactions = string_df[
        (string_df['gene1'].isin(genes)) & 
        (string_df['gene2'].isin(genes)) &
        (string_df['combined_score'] > 700)
    ]
    
    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row['gene1'], row['gene2'], weight=row['combined_score'])
    
    return G

def visualize_ppi(G, output='ppi_network.png'):
    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G, k=2)
    nx.draw_networkx(G, pos, node_size=100, font_size=6, alpha=0.7)
    plt.savefig(output, dpi=300)
```

## 典型结果解读

### 案例：某肿瘤vs正常组织

| 类型 | 基因数 | 生物学解释 |
|------|--------|-----------|
| Q1 (Co-up) | 324 | 核心激活通路（转录+翻译双重激活）|
| Q2 (mRNA↑/Protein↓) | 89 | 转录后抑制/泛素化降解 |
| Q3 (mRNA↓/Protein↑) | 56 | 翻译激活/蛋白稳定化 |
| Q4 (Co-down) | 287 | 核心抑制通路 |

**重点关注：**
- Q2中下调的蛋白可能是被降解的肿瘤抑制因子
- Q3中上调的蛋白可能是翻译激活的癌基因

## 参考文献

- GSE156632 (mRNA-protein correlation in breast cancer)
- Vogel et al., Mol Cell Proteomics 2010
