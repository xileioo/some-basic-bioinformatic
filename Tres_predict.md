# **A T cell resilience model associated with response to immunotherapy in multiple tumor types**


https://baijiahao.baidu.com/s?id=1731779825270948210&wfr=spider&for=pc
https://www.nature.com/articles/s41591-022-01799-y

2022年5 月 2号，美国NIH的姜鹏团队在Nature Medicine上发表了文章A T cell resilience model associated with response to immunotherapy in multiple tumor types，发布了美国国家癌症中心的 Tres (Tumor resilience T cell, https://resilience.ccr.cancer.gov) 平台，用来寻找在实体肿瘤中受到各种抑制因子压力下的依然坚韧的T细胞的特征，和预测T细胞在免疫疗法中的效率。Tres还揭示了 FIBP 为新的T细胞免疫代谢检查点和可能的免疫疗法新靶点。


近几年来，单细胞测序技术在肿瘤研究中的应用产生了大量的单细胞基因表达谱，描绘了来自肿瘤的T细胞群体的各种状态。T细胞的生长快慢，感知周围环境中的信号，都会在基因表达谱上得到体现。首先基于 CytoSig 数据库（https://cytosig.ccr.cancer.gov），Tres 平台评估肿瘤环境中每个 T细胞所感知到的细胞因子。比如，TGF-beta 和PGE2是常见的免疫抑制因子。TRAIL 是T细胞自杀机制触发因子。如果这些细胞因子的下游通路启动了，这说明这个T细胞处在一个非常不利的环境中。同时T细胞的健康状态可以用细胞周期和DNA复制通路基因的表达量来衡量。被压制的或者要死亡的细胞这些通路的活性往往会很低。基于以上计算评估的变量，Tres 来寻找哪个T细胞处于各种抑制因子的压力下，但依然保持了健康状态。这些T细胞被定义为肿瘤韧性T细胞（tumor resilient T cells简称Tres）。这些Tres 特征展示了重要的临床应用。

姜鹏团队系统收集这几年发表的具有免疫疗法临床效果的T细胞基因表达谱研究数据。这些数据包括anti PD1/CTLA4 在治疗前肿瘤中分离的T细胞基因表达，CAR T 给病人注射的产品细胞，和进行CAR T 或者肿瘤浸润T细胞（TIL）转移疗法中用来制造T细胞产品的病人血液或者肿瘤样本（pre-manufacture）。基于最简单的相关系数计算，如果这些T细胞样本和韧性T细胞特征谱正向相关，那么相对应的免疫治疗就会取得很好的疗效。如果负向相关，相对应的免疫疗法效果就会非常不好。

特别指出的是 Tres 平台在只用制造前样本（pre-manufacture）对于细胞疗法中疗效不好的病人预测精度几乎达到了完全正确。目前，CAR T 和 TIL 治疗的费用非常昂贵，而且可能造成致命的毒性。此结果预示着Tres 平台也许可以帮助明显对于细胞治疗没有希望的病人来避免不必要的巨额费用和副作用风险。

姜鹏团队还分析了Tres 的基因特征。在来自于19种癌症168个肿瘤单细胞表达数据中，FIBP 基因的T细胞高表达，几乎都预示了T细胞的低韧性。张渝博士利用基因编辑技术，在人体和小鼠的T细胞中敲除了这个基因，并发现在FIBP敲除后，T 细胞对肿瘤细胞的杀伤效果大大增加了。在跟进了机制研究后，张渝博士发现了FIBP 敲除显著降低了CD8 T 细胞的胆固醇代谢和吸收。很多实体肿瘤中胆固醇的浓度都很高。尽管适量的胆固醇会保障T细胞活性，但过度的胆固醇浓度会大大降低T细胞的肿瘤杀伤能力和导致T细胞耗竭。FIBP 敲除可以有效的下调T细胞对环境中胆固醇的吸收率和降低自身的胆固醇合成。

综上所述，Tres 平台利用了大数据建模为癌症免疫治疗，特别是T细胞疗法提供了一个重要的研发工具，并且对指导细胞疗法临床决策提供了一个重要的依据。

姜鹏博士为论文通讯作者。张渝医生博士和 吴创博士为共同第一作者。关新元教授，戴池平，吴船和陈祚珈博士也做出了关键贡献。

方法：<br/>
Variable interaction test in multivariate linear regression <br/>
多元线性回归中的变量交互测试 <br/>

1. Calculate Suppression and Proliferation score for each cell <br/>
    Supperesion score: CytoSig <br/>
    Proliferation score: Proliferation indicates the T cell proliferation score, computed through a linear regression approach. The output variable is the single-cell RNA-seq transcriptome. The explanatory variable is a binary vector with value 1 for all genes in the cell cycle and DNA replication pathways from the KEGG database; and value 0 for all other genes in KEGG. The proliferation score is computed as the t value (Coefficient / Stderr) of the explanatory variable, representing whether a T cell is proliferative <br/>
  
  我的推测： (1,1,1,0,0,0...0) * coefficient = (3.1, 7.5, 12.6, 1.5, -3, -6...2)  <br/>
              m\*1 X 1\*1 = m\*1 <br/>

2.  Calulate Tres Score for each genes in each tumor sample <br/>
    G represents the expression level of gene G across CD8+ T cells in a tumor. <br/>
    For each gene G, we performed the following regression for each cytokine: <br/>

d + a × suppression + b × G + c × suppression × G = proliferation <br/>

我的推测：input: G1, suppression --> output: proliferation 对单一肿瘤样本中的每个细胞的gene1 和 score进行线性拟合 <br/>
            最后的Tres score 是相对gene1 在 某个tumor中的一个打分 <br/>

  The Tres score is defined as the t value: c / StdErr(c) <br/>
  To further understand the variable interaction test, we can rewrite the model as follows: <br/>
  d + b × G + (a + c × G) × suppression = proliferation <br/>

