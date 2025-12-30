import pandas as pd
import os

# 读取进化树数据
data = {}
fa = open(r'D:\\projects\\HZ\\data\\1059sam-Results\\mag.tax.xls')
falines = fa.readlines()
fa.close()
headline = falines[0].strip('\n\r').split('\t')
for i in range(1, len(falines)):
    falines[i] = falines[i].strip('\n\r').split('\t')
    taxonomy = falines[i][-1]
    taxonomy = taxonomy.split(';')
    if len(taxonomy) < 6:
        continue
    for o in range(len(taxonomy)):
        taxonomy[o] = taxonomy[o].split('__')[1]
    if taxonomy[5] not in data:
        data[taxonomy[5]] = {'D': taxonomy[0], 'P': taxonomy[1], 'C': taxonomy[2], 'O': taxonomy[3], 'F': taxonomy[4], 'G': taxonomy[5]}

print(f"已加载 {len(data)} 个genus的分类信息")

# 读取Genus相关性结果
genus_corr_file = 'corr_tls_bac_Genus.txt'
genus_df = pd.read_csv(genus_corr_file, sep='\t')

# 筛选pearson显著的结果
# genus_sig = genus_df[genus_df['pearson_sig'] == True].copy()
genus_sig = genus_df[genus_df['pearson_p'] < 0.05].copy()
print(f"Genus显著结果: {len(genus_sig)} 条")

# 读取Species相关性结果
species_corr_file = 'corr_tls_bac_Species.txt'
species_df = pd.read_csv(species_corr_file, sep='\t')

# 筛选pearson显著的结果
species_sig = species_df[species_df['pearson_sig'] == True].copy()
# species_sig = species_df[species_df['pearson_p'] < 0.05].copy()
print(f"Species显著结果: {len(species_sig)} 条")

# 处理Genus结果，补全进化树信息
genus_results = []
for idx, row in genus_sig.iterrows():
    genus_name = row['compared_list']
    result = {
        'taxonomy_level': 'Genus',
        'name': genus_name,
        'pearson_corr': row['pearson_corr'],
        'pearson_p': row['pearson_p'],
        'pearson_fdr': row['pearson_fdr'],
        'pearson_sig': row['pearson_sig']
    }
    
    # 查找分类信息
    if genus_name in data:
        result.update(data[genus_name])
    else:
        # 如果找不到，使用未知值
        result.update({'D': 'Unknown', 'P': 'Unknown', 'C': 'Unknown', 
                      'O': 'Unknown', 'F': 'Unknown', 'G': genus_name})
    genus_results.append(result)

# 处理Species结果，补全进化树信息
species_results = []
for idx, row in species_sig.iterrows():
    species_name = row['compared_list']
    # 从species名称中提取genus（通常是第一个单词）
    genus_name = species_name.split()[0] if ' ' in species_name else species_name

    result = {
        'taxonomy_level': 'Species',
        'name': species_name,
        'pearson_corr': row['pearson_corr'],
        'pearson_p': row['pearson_p'],
        'pearson_fdr': row['pearson_fdr'],
        'pearson_sig': row['pearson_sig']
    }
    
    # 查找分类信息（基于genus）
    if genus_name in data:
        result.update(data[genus_name])
        # 更新genus为species名称（完整分类）
        result['G'] = genus_name  # 保持genus名称
        result['S'] = species_name  # 添加species名称
    else:
        # 如果找不到，使用未知值
        result.update({'D': 'Unknown', 'P': 'Unknown', 'C': 'Unknown', 
                      'O': 'Unknown', 'F': 'Unknown', 'G': genus_name, 'S': species_name})
    species_results.append(result)

# 合并结果并输出
all_results = genus_results + species_results

# 转换为DataFrame以便输出
output_df = pd.DataFrame(all_results)

# 重新排列列的顺序，使其更易读
column_order = [
    'taxonomy_level',
    'name',
    'D',
    'P',
    'C',
    'O',
    'F',
    'G',
    'S',
    'pearson_corr',
    'pearson_p',
    'pearson_fdr',
    'pearson_sig'
]

# 确保所有列都存在
for col in column_order:
    if col not in output_df.columns:
        output_df[col] = ''

output_df = output_df[column_order]

# 输出到文件
output_file = 'significant_bacteria_phylogeny_correlation.txt'
output_df.to_csv(output_file, sep='\t', index=False)
print(f"\n结果已保存到: {output_file}")
print(f"总共 {len(all_results)} 个显著细菌（{len(genus_results)} 个Genus, {len(species_results)} 个Species）")

# 同时输出分类统计
print("\n=== 统计信息 ===")
print(f"Genus显著结果数: {len(genus_results)}")
print(f"Species显著结果数: {len(species_results)}")
print(f"找到分类信息的Genus数: {sum(1 for r in genus_results if r['D'] != 'Unknown')}")
print(f"找到分类信息的Species数: {sum(1 for r in species_results if r['D'] != 'Unknown')}")

# 显示前几个结果作为预览
print("\n=== 前10个结果预览 ===")
print(output_df.head(10).to_string())