import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# 使用方法: python plot_identity.py A_identity.tsv
file_path = sys.argv[1]

print(f"正在读取文件: {file_path} ...")

# 读取 Diamond/BLAST 输出结果 (outfmt 6)
# 列名: qseqid, sseqid, pident, ...
df = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2], names=['qseqid', 'sseqid', 'pident'])

# --- 关键步骤：数据清洗 ---
# 1. 去除自身比对 (Query == Subject)
df = df[df['qseqid'] != df['sseqid']]

# 2. 只要高相似度区域的数据
df_high = df[df['pident'] >= 70]

print(f"数据读取完毕。用于绘图的数据点数量: {len(df_high)}")

# --- 绘图 ---
plt.figure(figsize=(10, 6))

# 绘制直方图 (Histogram) 和 密度曲线 (KDE)
sns.histplot(df_high['pident'], bins=100, kde=True, color='skyblue', edgecolor='black')

plt.title('Sequence Identity Distribution (70% - 100%)', fontsize=15)
plt.xlabel('Identity (%)', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 标记常见的 CD-HIT 阈值位置
for threshold in [90, 93, 94, 95, 96, 98]:
    plt.axvline(x=threshold, color='red', linestyle='--', alpha=0.5)
    plt.text(threshold, plt.gca().get_ylim()[1]*0.9, f'{threshold}%', color='red', ha='center')

output_img = file_path.replace('.tsv', '_dist.png')
plt.savefig(output_img, dpi=300)
print(f"图表已保存为: {output_img}")
plt.show()
