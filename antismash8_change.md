# antismash8更新变化

## 1 BGC检测

检测BGC规则：81增加到101条，检测更多类型BGC

新增命令行选项：
- `--hmmdetection-limit-to-rule-names`: 限制只检测指定规则名称的原始簇protoclusters   
- `--hmmdetection-limit-to-rule-categories`: 限制只检测指定规则类别的原始簇protoclusters   

## 2 萜烯分析

萜烯合酶/环化酶（terpene synthases/cyclases）的分析

## 3 基因功能分析

新增定制化注释`Tailoring functions`，肽骨架后修饰相关的功能注释，尤其是糖肽类？
可以从“卤代吡咯（halogenated pyrrole）”到具有特定位点选择性的色氨酸 卤化酶（如哪一个碳原子被卤化）的具体预测。
来源于MITE数据库

## 4 NRPS/PKS改进

新增了针对铁载体相关 β-羟化酶（siderophore-associated β-hydroxylases）和界面结构域（interface domains）的识别模型，以及一个更通用的 α/β 水解酶（α/β-hydrolase）识别模型

antiSMASH 现在还会检查 NRPS 的缩合（condensation, C）结构域和外消旋化（epimerization, E）结构域的活性位点是否具备催化残基；如果这些关键残基缺失，则会被标记为“失活”。

已有的 NRPS 腺苷酸化（adenylation, A）结构域的底物特异性预测功能：新增指向外部 PARAS 底物特异性预测软件的链接

## 5 其他

overview中similarity评分 百分数 变 高中低。高75+，中50+，低15+，15以下不相似

移除真菌基因预测功能，GlimmerHMM

移除 ClusterBlast 结果在 GenBank 文件注释中的展示