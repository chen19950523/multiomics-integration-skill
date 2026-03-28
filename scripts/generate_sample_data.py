#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Sample Input Data for Testing
生成示例输入数据用于测试

运行此脚本生成测试数据：
python generate_sample_data.py
"""

import pandas as pd
import numpy as np
import os

np.random.seed(42)


def generate_sample_data():
    """生成示例测试数据"""
    
    # 创建输出目录
    os.makedirs('input', exist_ok=True)
    
    # ============ Framework 1: 转录组 + 蛋白质组 ============
    
    # 基因列表
    genes = [f'GENE_{i:03d}' for i in range(1, 201)]
    
    # 样本名
    samples = ['Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'Treat_1', 'Treat_2', 'Treat_3']
    
    # 生成DEGs数据
    deg_data = pd.DataFrame({
        'GeneSymbol': genes[:50],
        'log2FC': np.random.uniform(-3, 3, 50),
        'PValue': np.random.uniform(0, 0.05, 50)
    })
    deg_data.set_index('GeneSymbol', inplace=True)
    deg_data.to_csv('input/DEGs.csv')
    
    # 生成DEPs数据
    dep_data = pd.DataFrame({
        'ProteinID': genes[20:70],
        'log2FC': np.random.uniform(-2.5, 2.5, 50),
        'PValue': np.random.uniform(0, 0.05, 50)
    })
    dep_data.set_index('ProteinID', inplace=True)
    dep_data.to_csv('input/DEPs.csv')
    
    # 生成TPM数据
    tpm_data = pd.DataFrame(
        np.random.lognormal(mean=5, sigma=1, size=(200, 6)),
        index=genes,
        columns=samples
    )
    tpm_data.to_csv('input/TPM.csv')
    
    # 生成蛋白intensity数据
    intensity_data = pd.DataFrame(
        np.random.lognormal(mean=10, sigma=1.5, size=(200, 6)),
        index=genes,
        columns=samples
    )
    intensity_data.to_csv('input/protein_intensity.csv')
    
    print("[完成] Framework 1 测试数据已生成")
    
    # ============ Framework 2: 蛋白质组 + 乳酸化组 ============
    
    # 蛋白列表
    proteins = [f'PROT_{i:03d}' for i in range(1, 151)]
    
    # 乳酸化位点数据
    lactylation_sites = []
    for prot in proteins[:80]:
        n_sites = np.random.randint(1, 5)
        for site in range(n_sites):
            lactylation_sites.append({
                'ProteinID': prot,
                'Site': np.random.randint(10, 500),
                'Sequence': 'MKLAVLGAAG[AST]'[:15],  # 简化的序列
                'Intensity': np.random.lognormal(8, 2),
                'log2FC': np.random.uniform(-2, 2),
                'PValue': np.random.uniform(0, 0.1)
            })
    
    lact_df = pd.DataFrame(lactylation_sites)
    lact_df.to_csv('input/lactylation_sites.csv', index=False)
    
    # 蛋白定量数据
    proteome_data = pd.DataFrame(
        np.random.lognormal(mean=10, sigma=1.5, size=(150, 6)),
        index=proteins,
        columns=samples
    )
    proteome_data.to_csv('input/proteome_quantification.csv')
    
    print("[完成] Framework 2 测试数据已生成")
    
    # ============ Framework 3: 蛋白质组 + 磷酸化组 ============
    
    # 磷酸化位点数据
    phospho_sites = []
    for prot in proteins[:100]:
        n_sites = np.random.randint(1, 8)
        for site in range(n_sites):
            phospho_sites.append({
                'ProteinID': prot,
                'Site': np.random.randint(5, 800),
                'Position': f"S{np.random.randint(50, 700)}",  # 磷酸化位点格式
                'log2FC': np.random.uniform(-2.5, 2.5),
                'PValue': np.random.uniform(0, 0.1),
                'LocalizationProb': np.random.uniform(0.5, 1)
            })
    
    phos_df = pd.DataFrame(phospho_sites)
    phos_df.to_csv('input/phosphorylation_sites.csv', index=False)
    
    print("[完成] Framework 3 测试数据已生成")
    
    # ============ Framework 4: 蛋白质组 + 代谢组 ============
    
    # 代谢物列表
    metabolites = [f'MET_{i:03d}' for i in range(1, 101)]
    
    # 代谢物定量数据
    metabolome_data = pd.DataFrame(
        np.random.lognormal(mean=6, sigma=2, size=(100, 6)),
        index=metabolites,
        columns=samples
    )
    metabolome_data.to_csv('input/metabolome_quantification.csv')
    
    print("[完成] Framework 4 测试数据已生成")
    
    # ============ 样本信息 ============
    
    sample_info = pd.DataFrame({
        'Sample': samples,
        'Group': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment']
    })
    sample_info.to_csv('input/sample_info.csv', index=False)
    
    print("[完成] 样本信息已生成")
    
    print("\n" + "=" * 50)
    print("所有测试数据已生成到 input/ 目录")
    print("=" * 50)


if __name__ == '__main__':
    generate_sample_data()
