[AUCell](https://www.bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#calculate-enrichment-for-the-gene-signatures-auc) <br/>
[AUCell package](https://bioconductor.org/packages/devel/bioc/manuals/AUCell/man/AUCell.pdf)

step1：Build gene-expression rankings for each cell。根据表达对所有基因在单个细胞中的排序，进而推导所有基因在所有细胞中进行统一排序。得到基因排序后的基因X细胞表达矩阵。top 5%的基因用于后续AUC（Area Under the Curve）计算。

step2:  Calculate enrichment for the gene signatures (AUC)。对于每个细胞，按照AUC滑动计算top 5%基因中属于我们特征geneset的基因。曲线下的面积就是AUC。得到对于我们特征geneset在每个细胞中的AUC score.

Step3: Determine the cells with the given gene signatures or active gene sets. 用AUC score - total cells，画密度分布图。选择曲线的双峰拐点，或者1% outliers的值确定AUC threshold。

优势：genesets中的基因越多，这个方法越准确。且可以对不同的细胞，不同的测序水平，有针对性的计算不同的阈值，减少组织差异和测序带来的误差。实验室的决策树方法，只是用了10个基因作为signature，不适用该方法。

不足：该方法的threshold是自动计算的，生物学意义可能需要针对具体问题，手动调试，确定最终阈值。
