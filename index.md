## 芯片数据分析
### 1、数据预处理
- **探针——基因名替换**：将探针的名字替换为基因的名称
- **去重**
  - 重复行：探针/基因名和行表达数据完全相同，直接去重，保留一行即可。
  ```R
  #两个函数都行
  data = unique(data)
  data = data[!duplicated(data),]
  ```
  - 多个探针对应同一个基因名称：1）对相同基因名称的行，按列（样本）**取平均（mean）,取最大值（max）**，作为该基因的在样本中的表达值；2）1）对相同基因名称的行，按列（样本）**求和（sum）**，作为该基因的在样本中的表达值；3）挑选相同基因中，行平均值最大的一行。
  ```R
  # expr_mat第一列是"Symbol":基因名称
  # 1）
  expr_mat = aggregate(formula = .~Symbol,data=expr,FUN = mean)
  expr_mat = aggregate(formula = .~Symbol,data=expr,FUN = max) # 不合理，一般不用
  # dpylr 应用，同aggregate
  library(dplyr)
  expr_mat %>% group_by(ID_REF) %>% summarise_all(mean)
  # 2）
  expr_max = aggregate(formula = .~Symbol,data=expr, FUN = sum)
  # 取平均和求和在后续计算两组间 fold change时，结果一样（可理论推导验证一下）
  # 3）
  #计算行平均值，按降序排列
  index = order(rowMeans(expr[,-1]),decreasing = T)
  #调整表达谱的基因顺序
  expr_ordered = expr[index,]
  #对于有重复的基因，保留第一次出现的那个，即行平均值大的那个;
  expr_max = expr[!duplicated(expr_ordered$Symbol),]
  expr_max
  ```
    > 基因名去重复的一点思考：
    > 这两种思路的差别在于,1）只取表达量最高的基因，认为只有这个基因有意义，其余表达量靠后的相同基因不重要。2）合并所有具有相同基因名的基因，考虑到了所有基因的表达情况，考虑更全面，因此个人更推荐。 实际分析中，由于我们一般差异分析是对不同样品中的同一个基因进行的比较，因此这两种方法实际差别并不大，按自己需求选择即可。（其实分析时有些人嫌麻烦，还会直接合并symbol和表达矩阵后随缘保留一个重复基因，这种方法就见仁见智了。
    > For probesets that map to identical Entrez gene names, select the one with highest IQR (for Affy, select mean for Agilent)，也就是四分位间距IQR，这个概念主要是在boxplot图表里面显示出来。当然了，不同芯片平台也是有一些细微的差别。一代Array探针可以这么做，RNA seq会出现一个gene symbol对应多个isform的数据（有点类似array的这种情况吧），1）一个md Anderson 的老师说他们用**最长的CCDS的那个transcript**作为这个基因的代表，2）另一个ucla的老师说他们是**将所有的isform表达量加起来**作为这个基因的表达量。

    >是否过滤基因：要过滤
    [ref link](https://www.jieandze1314.com/post/cnposts/249/)
- **标准化**
  [reference link](https://blog.csdn.net/tommyhechina/article/details/80356110)

### 2. 差异分析
- 样本数较少时
- [reference link](https://www.jianshu.com/p/a5196698ba98)
### 3. 假设检验
- p value计算 
$\log_2 \frac{mean(\vec{x_A})}{mean(\vec{x_B})}$

#### 3.1 单边检验还是双边检验？
[reference link](https://towardsdatascience.com/permutation-test-in-r-77d551a9f891)

双侧的原假设用的是*a ≠ a~0~*，讨论的是相等性质的问题（有无差异）；单侧检验的假设是 *a > a~0~* 或者 *a < a~0~*，讨论的是不等性质的问题（有差异，且优于/劣于）。

- 单变量置换检验
  Permutation test 可以称作是置换检验，Fisher于20世纪30年代提出的一种基于大量计算（computationally intensive），利用样本数据的全（或随机）排列，进行统计推断的方法。其属于一种**非参数检验**，**对样本总体分布情况无要求**，特别适用于**总体分布未知的<u>小样本数据</u>**，即使样本数据小到无法使用比如说t检验（这一点是说Permutation test可以用于样本量非常小的数据）。当然，Permutation test也可以用于**分布未知的大样本量的数据**，以及某些难以用常规方法分析资料的假设检验问题。因此，其应用非常广泛。在具体使用上它和Bootstrap Methods类似，通过对样本进行顺序上的置换，重新计算统计检验量，构造经验分布，然后在此基础上求出P-value进行推断。

  **Strategy**: 对两组样本进行顺序上的随机置换，并重新计算统计检验量（一般是两组均值差），把上述过程重复多遍（比如说1000遍），就可以构造出统计检验量的经验分布，然后对比两组样本的统计检验量和构造出的统计检验量经验分布，就可以计算求出P值。最终也是利用p-value的值来判断假设是否成立的，看p值的大小，p值小于0.05时，是说明拒绝H0，大于0.05.则是说明服从0假设。
```R
# 已有的轮子
tmp_data = data.frame(exp=expr, group = factor(rep(c("case","control"),c(25,20)),levels=c("case","control")))
tmp_res = coin::oneway_test(exp~group, data = tmp_data, alternative = "two.sided")
p_value = coin::pvalue(tmp_res)
```

```R
# 自己造轮子，仅示例
d <- as.data.frame(cbind(rnorm(1:20, 500, 50), c(rep(0, 10), rep(1, 10))))
treatment <- d$V2
outcome <- d$V1

# Difference in means
# diff(x, lags = 1)函数计算 x中间隔lags-1的两项差值，所以lags=1表示计算相邻两项差值（后项-前项）
original <- diff(tapply(outcome, treatment, mean))
mean(outcome[treatment==1])-mean(outcome[treatment==0])

#Permutation test
permutation.test <- function(treatment, outcome, n){
    distribution=c()
    result=0
    distribution=unlist(lapply(1:n,function(x) 
        diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))    ))
    result=sum(abs(distribution) >= abs(original))/(n)
    return(list(result, distribution))
}
test1 <- permutation.test(treatment, outcome, 10000)
hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")
```
- 最终结果
```R
final_res = data.frame(log2FC = log2FC_res, p.value = p_value)
final_res = dpylr::arrange(final_res, p.value)
# Benjamini Hochberg方法计算FDR，假阳性率，方法大概：
final_res[["p.adj"]] = p.adjust(final_res$p.value,method = "BH")
```
 - - - 
[Principle and implementation of enrichment analysis](https://www.fatalerrors.org/a/principle-and-implementation-of-enrichment-analysis.html)

[富集分析实现原理](https://zhuanlan.zhihu.com/p/426988267)

[disgenet](https://www.disgenet.org/home/)
