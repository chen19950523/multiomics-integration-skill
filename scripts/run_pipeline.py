#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-omics Integration Pipeline - Main Orchestrator
多组学联合分析流程 - 主协调器

支持四种联合分析框架的自动化运行

Usage:
    python run_pipeline.py --framework 1 --input data/
    python run_pipeline.py --framework all --config config.yaml
"""

import argparse
import os
import sys
import subprocess
from datetime import datetime

# 工作目录
WORK_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(WORK_DIR)


def setup_directories():
    """创建必要的目录结构"""
    dirs = ['input', 'results', 'logs', 'config']
    for d in dirs:
        path = os.path.join(PROJECT_DIR, d)
        os.makedirs(path, exist_ok=True)
    print(f"目录结构已创建")


def check_dependencies():
    """检查必要的Python包"""
    required_packages = [
        'pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn',
        'scikit-learn', 'networkx', 'clusterProfiler'
    ]
    
    missing = []
    for pkg in required_packages:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    
    if missing:
        print(f"缺少以下包: {', '.join(missing)}")
        print("请运行: pip install -r requirements.txt")
        return False
    return True


def run_framework(framework_num, input_dir='input', output_dir='results'):
    """
    运行指定框架的分析
    
    Parameters:
    -----------
    framework_num : int or str
        框架编号 1-5 或 'all'
    input_dir : str
        输入数据目录
    output_dir : str
        输出结果目录
    """
    framework_scripts = {
        '1': 'scripts/1_transcriptome_proteome/transcriptome_proteome_integration.py',
        '2': 'scripts/2_proteome_lactylation/proteome_lactylation_integration.py',
        '3': 'scripts/3_proteome_phosphorylation/proteome_phosphorylation_integration.py',
        '4': 'scripts/4_proteome_metabolome/proteome_metabolome_integration.py',
        '5': 'scripts/5_triple_modification/triple_modification_crosstalk.py',
    }
    
    framework_names = {
        '1': 'Transcriptome-Proteome',
        '2': 'Proteome-Lactylation',
        '3': 'Proteome-Phosphorylation',
        '4': 'Proteome-Metabolome',
        '5': 'Triple-Modification Crosstalk'
    }
    
    if framework_num == 'all':
        frameworks = list(framework_scripts.keys())
    else:
        frameworks = [str(framework_num)]
    
    results_summary = {}
    
    for fw in frameworks:
        if fw not in framework_scripts:
            print(f"错误: 无效的框架编号 {fw}")
            continue
        
        script_path = os.path.join(PROJECT_DIR, framework_scripts[fw])
        
        if not os.path.exists(script_path):
            print(f"错误: 脚本不存在 {script_path}")
            continue
        
        print("\n" + "=" * 60)
        print(f"开始运行: {framework_names[fw]} 框架")
        print(f"脚本: {script_path}")
        print("=" * 60)
        
        start_time = datetime.now()
        
        try:
            # 运行分析脚本
            # 实际部署时替换为实际的subprocess调用
            print(f"[INFO] 分析脚本路径: {script_path}")
            print(f"[INFO] 输入目录: {input_dir}")
            print(f"[INFO] 输出目录: {output_dir}")
            
            # 记录运行时间
            elapsed = (datetime.now() - start_time).total_seconds()
            
            results_summary[fw] = {
                'status': 'Ready',
                'framework': framework_names[fw],
                'elapsed_seconds': elapsed,
                'script': script_path
            }
            
            print(f"[完成] {framework_names[fw]} 框架准备就绪")
            
        except Exception as e:
            print(f"[错误] {framework_names[fw]} 运行失败: {str(e)}")
            results_summary[fw] = {
                'status': 'Error',
                'framework': framework_names[fw],
                'error': str(e)
            }
    
    return results_summary


def generate_report(results_summary):
    """生成分析报告"""
    report_path = os.path.join(PROJECT_DIR, 'results', 'pipeline_report.txt')
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("Multi-omics Integration Pipeline Report\n")
        f.write("多组学联合分析流程报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("-" * 40 + "\n")
        f.write("框架运行状态:\n")
        f.write("-" * 40 + "\n")
        
        for fw, result in results_summary.items():
            status = result.get('status', 'Unknown')
            framework = result.get('framework', 'Unknown')
            f.write(f"  {fw}. {framework}: {status}\n")
        
        f.write("\n" + "-" * 40 + "\n")
        f.write("输出目录结构:\n")
        f.write("-" * 40 + "\n")
        f.write("  results/\n")
        f.write("    1_transcriptome_proteome/\n")
        f.write("    2_proteome_lactylation/\n")
        f.write("    3_proteome_phosphorylation/\n")
        f.write("    4_proteome_metabolome/\n")
        f.write("    5_triple_modification/\n")
    
    print(f"\n报告已生成: {report_path}")
    return report_path


def main():
    parser = argparse.ArgumentParser(
        description='Multi-omics Integration Pipeline / 多组学联合分析流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py --framework 1
  python run_pipeline.py --framework all
  python run_pipeline.py --framework 2 --input ./data --output ./results

Supported Frameworks:
  1 - Transcriptome + Proteome
  2 - Proteome + Lactylation  
  3 - Proteome + Phosphorylation
  4 - Proteome + Metabolome
  5 - Triple Modification Crosstalk (Phospho + Lacto + Proteome)
  all - Run all frameworks
        """
    )
    
    parser.add_argument('--framework', '-f', default='all',
                       help='分析框架编号 (1-5) 或 all (默认: all)')
    parser.add_argument('--input', '-i', default='input',
                       help='输入数据目录 (默认: input)')
    parser.add_argument('--output', '-o', default='results',
                       help='输出结果目录 (默认: results)')
    parser.add_argument('--config', '-c', default='config/pipeline_config.yaml',
                       help='配置文件路径 (默认: config/pipeline_config.yaml)')
    parser.add_argument('--check-deps', action='store_true',
                       help='检查Python依赖包')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Multi-omics Integration Pipeline")
    print("多组学联合分析流程")
    print("=" * 60)
    
    # 检查依赖
    if args.check_deps or not check_dependencies():
        check_dependencies()
        sys.exit(0 if check_dependencies() else 1)
    
    # 创建目录
    setup_directories()
    
    # 运行分析
    results = run_framework(args.framework, args.input, args.output)
    
    # 生成报告
    generate_report(results)
    
    print("\n" + "=" * 60)
    print("流程执行完成!")
    print("=" * 60)


if __name__ == '__main__':
    main()
