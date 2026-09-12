# 八味解郁汤–MDD候选网络NBVK分析报告

## 结论摘要

预注册的三节点**总网络损伤假设不获支持**。ARG1、LCN2和LTF联合删除造成的保留节点对通信效率损失为0.001628，即约0.1628%。该值高于拓扑匹配零分布均值0.000752，但低于其97.5%分位数0.002556；右尾经验P=0.1153，7项损伤检验BH FDR=0.2018。

预设Bliss型组合偏离量在score≥0.700主网络中为正：

- 三节点`Syn_Bliss = 3.083×10^-5`；
- 拓扑匹配右尾经验P=0.00530；
- 4项协同检验BH FDR=0.0106。

但该结果不能解释为稳健的“拓扑协同”。它约占三节点总损伤的1.89%，约95.9%来自LCN2+LTF组合，并对网络阈值及组合内部结构敏感。最准确的表述是：

> 在STRING score≥0.700网络和预设的节点边际degree+betweenness匹配零模型下，ARG1/LCN2/LTF呈现一个较小的、阈值特异性的Bliss型组合剩余项；总损伤并未超过拓扑匹配基线，且该剩余项尚不能排除LCN2–LTF直接关联和共享邻域所造成的组合结构效应。

## 冻结设计

- 主网络：STRING human functional network，score≥0.700，最大连通分量1,226节点、13,416条边。
- 观察集合：3个单节点、3个双节点和1个三节点。
- 主要指标：相同保留节点对在删除前后的全局效率损失。
- 主要基线：每个候选节点分别匹配`log1p(degree)`和betweenness百分位，50节点匹配池。
- 每组10,000次Monte Carlo draws，种子20260726。
- 损伤与协同分别构成7项和4项BH检验族。

## 预注册主要结果

### 总损伤

| 扰动集合 | D_E | 拓扑匹配P | BH FDR | 判定 |
|---|---:|---:|---:|---|
| ARG1 | 0.001473 | 0.0796 | 0.2018 | 不支持 |
| LCN2 | 0.000068 | 0.4236 | 0.4942 | 不支持 |
| LTF | 0.000057 | 0.5991 | 0.5991 | 不支持 |
| ARG1+LCN2 | 0.001542 | 0.1029 | 0.2018 | 不支持 |
| ARG1+LTF | 0.001530 | 0.1134 | 0.2018 | 不支持 |
| LCN2+LTF | 0.000154 | 0.3686 | 0.4942 | 不支持 |
| ARG1+LCN2+LTF | 0.001628 | 0.1153 | 0.2018 | 不支持 |

三节点删除后新增的断连主要对应ARG1删除后SLC7A2孤立；联合删除没有产生与LCN2/LTF相对应的额外网络碎裂。

### Bliss型组合偏离量

| 组合 | Syn_Bliss | 拓扑匹配P | BH FDR |
|---|---:|---:|---:|
| ARG1+LCN2 | 1.023×10^-6 | 0.2889 | 0.3332 |
| ARG1+LTF | 9.003×10^-7 | 0.3332 | 0.3332 |
| LCN2+LTF | 2.956×10^-5 | 0.00170 | 0.00680 |
| ARG1+LCN2+LTF | 3.083×10^-5 | 0.00530 | 0.0106 |

这里的“Bliss型”是预定义的网络残差统计量，而不是药理学Bliss独立性证明。不同删除集合的效率分母对应不同保留节点对集合，也限制其概率独立解释。

## 预注册敏感性结果

### 网络定义

| 网络版本 | 三节点D_E | 损伤FDR | 三节点Syn_Bliss | 协同FDR |
|---|---:|---:|---:|---:|
| score≥0.700，无权主网络 | 0.001628 | 0.2018 | 3.083×10^-5 | 0.0106 |
| score≥0.400，无权网络 | 0.000214 | 0.2410 | -9.771×10^-7 | 0.9925 |
| score≥0.700，加权网络 | 0.001672 | 0.2780 | 2.361×10^-5 | 0.0242 |

三节点总损伤在三个网络版本中均未达到预注册标准。Bliss型剩余项在主网络及加权网络为正，但在score≥0.400网络变为负值，说明其依赖网络阈值。

值得注意的是，score≥0.400网络中的LCN2+LTF双节点剩余项仍为正且FDR=0.0236，但三节点项为负。这进一步说明信号主要由特定节点对而非稳定的三节点整体结构驱动。

### 匹配方式和随机种子

- degree-only匹配：三节点总损伤FDR=0.0548；三节点协同FDR=0.0028。
- 第二随机种子20260727：拓扑匹配三节点总损伤FDR=0.2033；协同FDR=0.00980。
- 加权网络和第二种子结果说明数值不是单次Monte Carlo偶然，但不能解决组合内部拓扑未匹配的问题。

## 事后探索性压力测试

LCN2与LTF：

- 在score≥0.700网络中存在直接STRING功能关联边；
- 共享7个一阶邻居；
- 其匹配池重叠39/50；
- LCN2+LTF的协同剩余占三节点剩余项约95.9%。

条件于随机组合具有与观察组合相同的内部边数后：

| 组合 | 条件化P | BH FDR |
|---|---:|---:|
| LCN2+LTF | 0.0540 | 0.1079 |
| ARG1+LCN2+LTF | 0.0420 | 0.1079 |

该分析为看到结果后的压力测试，不能作为新的确认性检验；它说明原协同显著性可能由未匹配的组合内部边结构解释。

LCN2–LTF边的STRING combined score为0.996，但主要来自text-mining分量0.992；实验分量仅0.072，数据库分量为0。因此它不能直接解释为物理蛋白互作。

## 方法学价值的合理边界

当前NBVK仍展示了区别于中心性排序的分析框架：

1. 直接量化删除后保留节点对的通信效率损失；
2. 区分总损伤与组合偏离量；
3. 用随机和节点边际拓扑匹配基线检验观察结果；
4. 能识别“组合显著但总损伤不显著”以及“结果依赖网络定义”的情况。

但当前结果不足以证明NBVK已经解决组合层级拓扑混杂。若要主张方法学创新，下一版本需在独立预注册中加入组合层级匹配，包括内部边数、节点间距离、共同邻居/Jaccard、社区关系及可枚举组合的精确零分布。

## 可用于论文的结论

可以报告：

> The prespecified triple-node knockout did not produce greater overall communication-efficiency loss than topology-matched node sets. A small positive Bliss-form residual was detected in the high-confidence STRING network, but it was threshold-sensitive and largely attributable to the directly associated LCN2–LTF pair. Accordingly, the result was interpreted as a conditional network-interaction signal rather than evidence of pharmacological synergy or a validated treatment mechanism.

不应报告：

- “三靶点产生强网络协同”；
- “NBVK验证了八味解郁汤治疗抑郁症的机制”；
- “协同超过所有拓扑匹配组合”；
- “LCN2与LTF具有已验证的物理互作”。

## 可复现性与偏离

- 预注册、冻结输入SHA、全部观察值、完整零分布、代码、测试和运行环境均保存在本公开仓库的`docs/nbvk/`、`data/network_inputs/`、`results/nbvk_v1/`、`code/nbvk_v1/`和`environment/nbvk/`目录中。
- 所有真实结论均从已保存的完整零分布重新复算。
- 实现修复和未完全输出的预注册项目见`docs/nbvk/nbvk_v1_deviations.md`和`docs/nbvk/nbvk_v1_implementation_audit.md`。
- 组合内部边条件化结果明确存放于`results/exploratory_posthoc`，不得与预注册结果混合。
