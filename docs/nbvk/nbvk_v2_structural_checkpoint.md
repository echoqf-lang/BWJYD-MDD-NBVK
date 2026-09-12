# NBVK 2.0 结构层检查点报告

日期：2026-07-26

## 结论

本次计算**不支持** ARG1+LCN2+LTF 在控制组合结构后仍具有异常的全局网络效率损伤。

主要原因有两层：

1. 严格联合匹配只得到 599 个合法三元组，低于预注册的 1,000 个最低支持量，因此主要条件检验按规则停止。
2. 即使透明地查看预设放宽分析，观察损伤也未进入右尾显著区。

## 结果

观察三节点全局效率损伤：

- `D_E = 0.001627784`，即约 0.1628%。

| 匹配层级 | 合法组合数 | 零分布中位数 | 97.5%分位数 | 观察百分位 | 右尾P（+1） | 解释 |
|---|---:|---:|---:|---:|---:|---|
| 严格：诱导结构、距离、社区、共同邻域、邻域并集 | 599 | 0.000443439 | 0.001976556 | 93.82% | 0.06333 | 支持量不足，仅描述 |
| 放宽1：删除社区占用约束 | 999 | 0.000445031 | 0.002110119 | 92.99% | 0.07100 | 仍低于最低支持量，仅探索 |
| 放宽2：仅保留诱导结构与距离类别 | 6,546 | 0.000488500 | 0.002334230 | 89.18% | 0.10829 | 预设探索性敏感性，不显著 |

观察值在三层中均低于零分布 97.5% 分位数。

## 与 NBVK 1.0 的关系

观察 `D_E` 与 NBVK 1.0 保存值之差为 `3.11e-15`，达到数值精度一致，说明全局损伤引擎回归核对通过。

NBVK 1.0 的 Bliss 型剩余信号在未控制组合内部结构时出现；NBVK 2.0 控制三节点诱导结构与距离后，未获得右尾支持。这与此前红队发现“LCN2–LTF 直接边及共同邻域主导旧信号”的风险一致。

## 预注册执行状态

- 严格主要检验：因 `N=599<1000` 停止。
- 放宽1：`N=999<1000`，仍不足。
- 放宽2：完成全部 6,546 个唯一组合精确枚举。
- 未进行可选停止或按结果增加抽样。
- 未把任何放宽结果标记为确认性。

## 后续门控

预注册要求结构检验先通过，模块选择性才可进入“候选药理网络选择性扰动”的主要判定。结构门未通过，因此本轮不继续用 GO 模块结果挽救主要机制主张。

模块注册表仍可作为独立的描述性/图件准备工作另行生成，但不能改变当前结构层结论。

## 可用于论文的限定表述

> After conditioning on the induced three-node structure and pairwise distance categories, the observed global-efficiency loss caused by joint removal of ARG1, LCN2, and LTF was not greater than that of comparable triplets (exploratory relaxed null: D_E=0.001628; empirical P=0.108). The fully specified joint-matching set contained only 599 triplets and therefore failed the prespecified minimum-support criterion. These results do not support a robust combination-specific global network perturbation.

## 复现与校验

- Python 单元测试：4/4 通过。
- NBVK 1.0 原有核心测试：通过。
- 输出范围、组合数、集合嵌套和 P 值已独立复算。
- 分析计划 SHA-256：`3f3f34828f5771eba3cc3089c2b1988ae4ba2bd476a2bf1e6969833f358058cc`
- 预注册 SHA-256：`4250e6560f8c9383a9adc98525f05467a91e04174376f8bd3b75218bea397ec0`
- 结构推断表 SHA-256：`6ed73ddb28b479fb1760f6d538b70868f91653a6684d668ecfa54a773904e131`

## 证据边界

该分析仍为同一网络上的事后方法开发。非显著结果不是“证明三靶点没有作用”，而是说明现有 PPI 网络数据不支持“超出组合结构的异常全局扰动”这一特定命题。
