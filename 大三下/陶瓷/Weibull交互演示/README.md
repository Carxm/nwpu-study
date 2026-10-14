# Weibull 交互演示

[打开或下载 weibull_demo.html](weibull_demo.html)。这是从 2026 年陶瓷课程展示整理的教学页面，用于观察 selected stress 与 Weibull strength distribution 的关系。

下载后用浏览器打开 HTML 即可，无需安装 Python，也不需要联网。GitHub 文件页一般显示 HTML 源码，需下载后打开。

## 怎么操作

- `microsep. scenario`：在 0.00、0.10、0.20、0.30、0.50 mm 五个情景中选择。页面读取固定表中的应力值，**不会运行或重新求解 finite element analysis（FEA）**。
- `m`：改变 shape parameter，观察强度分散程度及 Weibull plot 的斜率。
- `σ₀`：改变 characteristic strength（MPa），观察整条分布的位置。
- `n samples`：改变红色 synthetic points 的数量，帮助理解有限样本的散点表现；不会改变理论分布或用样本重新拟合 m、σ₀。

建议先固定 σ₀，只改变 m；再固定 m 改变 σ₀；最后改变所选应力，比较同一分布在不同应力下的 CDF 数值。

## 如何理解显示值

蓝线采用 two-parameter Weibull CDF：

\[
P_f(\sigma)=1-\exp[-(\sigma/\sigma_0)^m].
\]

红点由页面生成，用于展示 low-strength tail 与样本散布，**不是实验数据或文献逐点数据**。`Pf index` 是给定分布和所选应力下的示意指标。它没有纳入完整应力场、有效体积、复杂载荷历史、疲劳与环境退化等因素，不能解释成临床植入物失效率、二十年寿命预测或 ISO 合格判定。

固定表保留原课程页面的五组案例值，来源为 CW2 目录中保存的 Weibull 交互展示资料。它们没有在本页面中重新验证；尤其不要把改变下拉选项理解为新完成了一次 FEA。学习使用时可对照[课程项目](../课程项目/README.md)中的相关讨论。
