# Weibull 分析工具

这个脚本读取 specimen-level failure data，使用 two-parameter Weibull model 进行 maximum likelihood estimation（MLE），再用 bootstrap 给出参数、分位数及曲线的不确定性。原脚本来自 CW1 性能分析资料；归档版补充了输入检查与明确标为 synthetic 的使用示例。

它适用于每个试样都有完整失效值的 **uncensored data（未删失数据）**。不能把未失效的 run-out 直接写成失效值，也不能直接处理 right-censored data、频数汇总、混合失效机制或重复测量相关性。

## 先看示例

[synthetic_strength.csv](examples/synthetic_strength.csv) 有两组各 60 个模拟强度值，单位 MPa。生成时设定 scale = 1200 MPa，shape 分别为 4 和 12，随机种子为 20261014。它是教学模拟，**不是本人实验结果，也不是从论文提取的数据**；有限样本拟合值不会精确等于生成参数。

可直接查看[示例统计表](examples/example_output/weibull_summary.csv)与[示例概率图](examples/example_output/weibull_probability_plot.png)。如果需要重新生成 CSV，运行 `examples/generate_synthetic_data.py`，它会覆盖本目录内的这份 synthetic 示例文件。

## 运行方法

电脑需要 Python 3。在本 README 所在目录打开终端。如果已经有可用环境，直接使用它；没有时可建立单独环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

运行随包示例，将新结果保存到 `weibull_output/`：

```bash
python3 weibull_analysis.py examples/synthetic_strength.csv \
  --label "Synthetic flexural strength (MPa)" \
  --bootstrap 1000 \
  --confidence 0.95 \
  --seed 42 \
  --outdir weibull_output
```

会生成 `weibull_summary.csv` 和 `weibull_probability_plot.png`。完整参数可通过 `python3 weibull_analysis.py --help` 查看。

## 输入格式

CSV 默认数值列叫 `value`，分组列叫 `group`。没有 `group` 列时，所有数据作为一组。可以通过 `--value-column`、`--group-column` 指定自己的列名。`unit` 和 `data_kind` 是示例的说明列，脚本不会自动识别或换算单位；同一次比较应统一单位和试验条件。

每一行应代表一个独立试样。失效值必须有限且严格大于零；缺失、非数值、无穷值、非正值和缺少的分组标签会报错，不会被静默删除。每组至少需要 3 条且不能全部相同；**3 条只是软件运行下限，不代表足以支持可靠的统计结论**。

`--facet-column scenario` 可按 `scenario` 列分图。同一个 group 必须只属于一个 facet；不同情景下请使用唯一的组名。`--x-reference` 与 `--x-reference-label` 只添加参考线，不会自动进行标准合格判定。

## 输出怎么看

| 字段 | 含义 |
| --- | --- |
| `n`、`mean`、`std` | 试样数、均值与样本标准差 |
| `shape_m` | Weibull shape parameter；在该模型内描述分散程度 |
| `scale_eta` | Scale parameter；CDF 为 1−exp(−1)，约 63.2% 时的失效值，与输入同单位 |
| `b10`、`b50` | 模型的 10% 和 50% failure quantiles，与输入同单位 |
| `*_ci_low`、`*_ci_high` | Bootstrap percentile interval 的上下限 |

概率图用 median ranks 放置数据点，参数来自 MLE，不是对这些点作简单直线回归。阴影是 **pointwise bootstrap band**，不能当作整条曲线同时成立的置信带。默认 confidence 为 0.95、bootstrap 为 1000、seed 为 42；改变 confidence 会同时改变表格区间和曲线阴影。

拟合假定同组试样可用同一种 two-parameter Weibull distribution 描述，location 固定为零。使用结果前仍需检查样本代表性、测试方式、尺寸/体积差异、失效机制及模型适合程度；低概率尾部尤其依赖样本量与分布假设。脚本没有进行完整的 brittle-fracture effective-volume 修正，也没有把试样强度转换为组件或人体内的寿命预测。
