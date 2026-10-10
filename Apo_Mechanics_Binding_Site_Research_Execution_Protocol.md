# Apo 状态蛋白质全局力学响应与配体结合位点研究
## 冻结版研究执行规划、项目审计与验收标准

**版本：v1.0 Frozen Protocol**  
**研究定位：计算生物物理 / 结构生物信息学 / 配体结合位点预测**  
**核心原则：先证伪、后扩库；先物理问题、后预测应用；开发集与最终测试集严格隔离。**

---

# 0. 项目一句话定义

本项目只回答一个核心科学问题：

> **蛋白质在无配体 apo 状态下的全局弹性响应，是否包含超越常规局部结构描述符的配体结合位点信息，并能否推广到 cryptic binding sites？**

英文工作表述：

> **Does intrinsic mechanical response in ligand-free apo structures provide ligand-binding-site information beyond conventional local structural descriptors?**

这里的“beyond”严格限定为：

> **beyond conventional local structural descriptors**，而不是 beyond structure itself。

因为 ANM/PRS 的力学量本身仍由三维坐标决定。本文要证明的是：将静态结构映射为弹性网络并求取全局扰动响应后得到的物理表征，是否能提供普通局部几何、packing 和邻域上下文未充分表达的信息。

---

# 1. 冻结后的论文逻辑

全文只保留四个连续结论：

1. **Ligand-binding regions exhibit characteristic mechanical-response profiles.**
2. **Mechanical response adds binding-site information beyond conventional local structural descriptors and neighborhood context.**
3. **The signal is already detectable in ligand-free apo structures and includes contributions from spatially distant perturbation sources.**
4. **Apo-state mechanical information generalizes to cryptic binding sites.**

其中：

- Result 1 为现象建立；
- Result 2 为“增量信息”证明；
- Result 3 为全文最重要的机制与 apo→holo 证据；
- Result 4 为最终外部泛化验证。

任何新增分析只有在直接服务以上四点时才允许进入主文。

---

# 2. 项目整体架构

```text
                    Protocol lock
                         │
                         ▼
            Seal CryptoBench test set
                         │
                         ▼
              Build AHoJ dev pool
                         │
                         ▼
              Pilot 200–300 pairs
                         │
           ┌─────────────┴─────────────┐
           │                           │
        No-Go                         Go
   停止扩库/收缩主张                    │
                                       ▼
                          Full AHoJ development study
                                       │
                 ┌─────────────────────┼─────────────────────┐
                 ▼                     ▼                     ▼
              Holo test           Apo→Holo             Source-distance
          描述 mechanical       核心前瞻验证           机制贡献分析
                 │                     │                     │
                 └─────────────────────┴─────────────────────┘
                                       │
                                       ▼
                         Freeze final model/parameters
                                       │
                                       ▼
                           CryptoBench sealed test
                                       │
                                       ▼
                                  Final paper
```

---

# 3. 不再改变的核心设计

## 3.1 研究单位

研究基本单位定义为：

> **protein–pocket pair**

而不是 UniProt，也不是 PDB structure。

一个 UniProt 可以有多个物理上不同的口袋：

```text
UniProt X
├── Pocket A：orthosteric ligand site
├── Pocket B：allosteric site
└── Pocket C：cofactor site
```

不同 pocket 必须保持独立身份。

### 同一 pocket 的多个 holo

若 AHoJ 映射表明多个 holo structures 对应同一个物理 pocket，可在该 pocket 内整合 binding-residue 标签：

\[
B_{pocket}=\bigcup_{h\in same\ pocket} B_h
\]

但不同 pocket 绝不能跨 pocket 取并集。

---

## 3.2 正标签、负标签和屏蔽标签

针对目标 pocket \(P\)：

### Positive

由对应 holo structure 定义并成功映射到 apo 的目标 pocket lining/binding residues：

\[
y_i=1
\]

### Masked

同一个 UniProt 中属于**其他已知物理口袋**的残基：

\[
y_i=\mathrm{MASK}
\]

这些残基：

- 不作为当前 pocket 的 positive；
- 也绝不能作为 negative；
- 不参与当前 pocket 的 AP/AUROC 评价。

### Negative

当前 target pocket 之外、且不属于任何其他已知 pocket、并满足结构可评价条件的残基：

\[
y_i=0
\]

这条规则防止“模型正确识别了另一个真实口袋，却被评价系统当成错误预测”。

---

# 4. 数据体系与泄漏控制

## 4.1 Dataset A：AHoJ Development Pool

AHoJ-DB 是本项目的主要数据来源，用于：

- Pilot；
- holo 描述性实验；
- apo→holo 主实验；
- GAM / LightGBM 开发；
- source-distance analysis；
- 参数和模型冻结。

AHoJ-DB 原论文报告其针对 468,293 个 biologically relevant protein–ligand interactions、27,983 个 proteins 建立 apo/holo pocket mapping，并报告约 24% binding sites 位于多链界面、不到一半位点具有可用 apo form。实际项目使用时必须以下载版本和原始元数据重新统计，不允许直接把论文总体数字当成最终样本量。

## 4.2 Dataset B：CryptoBench Sealed Test

CryptoBench 只用于最终 cryptic-site validation。

已核实的正式论文信息：

- 总计 1,107 个 apo structures；
- 1,361 个 cryptic pockets；
- benchmark train/test = 885 / 222 apo structures；
- benchmark 构建阶段另外采用了 10% sequence identity clustering 以降低 split leakage。

因此 CryptoBench **不是与 AHoJ 天然独立的数据源**，因为其由 AHoJ apo–holo 数据筛选而来。

### 强制隔离规则

在任何 Pilot 或模型训练开始之前：

1. 下载 CryptoBench 官方 split；
2. 建立 `cryptobench_test_sealed.tsv`；
3. 从 AHoJ development pool 删除 CryptoBench test 中出现的蛋白；
4. 若官方提供 cluster membership，删除与 test 属于同一官方 leakage-control cluster 的条目；
5. 额外使用 MMseqs2 对 AHoJ development sequences 与 CryptoBench test sequences 做同源检查；默认辅助排除规则为 **sequence identity ≥30% 且 alignment coverage ≥80%**；
6. 所有删除均记录 exclusion reason。

> CryptoBench test 在模型、超参数、特征定义全部冻结前不得用于任何调参、特征选择或阈值选择。

---

# 5. Protocol Lock：项目开始前必须完成

正式计算前生成：

```text
protocol/
├── analysis_protocol_v1.0.md
├── feature_dictionary.tsv
├── inclusion_exclusion_rules.md
├── cryptobench_test_sealed.tsv
├── sequence_cluster_rules.md
├── go_nogo_criteria.md
└── environment_lock.yml
```

并保存：

- Git commit hash；
- 日期；
- Conda/pip environment；
- ANM cutoff；
- 固定随机种子列表；
- 数据版本或下载日期；
- CryptoBench exclusion list checksum。

### Protocol Lock 验收标准

必须同时满足：

- [ ] CryptoBench test 已封存；
- [ ] AHoJ 开发集排除规则已写死；
- [ ] 正/负/MASK 标签规则已写死；
- [ ] ANM primary cutoff 固定为 10 Å；
- [ ] Static 和 Mechanics 特征定义已写死；
- [ ] CV 策略已写死；
- [ ] Go/No-Go 判据已写死；
- [ ] Pilot 开始前已提交 Git hash 或 OSF 时间戳版本。

若任一项未满足，不允许进入 Pilot。

---

# 6. AHoJ 数据清洗

## 6.1 Primary ligand 范围

主分析优先研究 biologically relevant small organic ligands。

建议固定排除：

- water；
- 单原子离子；
- PEG、glycerol 等常见 crystallization additives；
- nucleic acid；
- peptide/protein ligand；
- 明显非生理性 buffer molecule。

建议 small-molecule 最低 heavy-atom 数设为 ≥6；如 AHoJ 已提供 biological relevance 标注，以 AHoJ 标注为第一层过滤，再进行上述二次过滤。

## 6.2 结构质量

每个 protein–pocket pair 至少要求：

- holo pocket 可明确映射到 UniProt residues；
- apo 与 holo residue mapping 唯一且可追踪；
- pocket 关键 residues 无严重缺失；
- target pocket 映射后至少存在 5 个 positive residues；
- apo/holo biological assembly state 相同或可明确构建为等价 assembly；
- primary analysis 排除膜蛋白；
- interface pocket 必须保留形成 pocket 所需的全部 chains。

## 6.3 每个 pocket 只选择一个 primary apo

为了避免伪重复：

> 一个 target pocket 在 primary analysis 中只保留一个 primary apo structure。

选择规则必须预先固定，例如依次优先：

1. biological assembly 与 holo 一致；
2. pocket mapping 完整；
3. 分辨率/结构质量更高；
4. 缺失残基更少；
5. 若仍并列，按 PDB ID 字典序或固定随机种子选择。

其他 apo structures 留给 reproducibility analysis，不进入 primary per-pocket evaluation。

## 6.4 Strict apo sensitivity subset

Primary apo：target pocket 无对应配体。

Sensitivity subset：更加严格地要求：

- target pocket 无 ligand；
- 结构中无可能显著改变整体构象的大分子/小分子结合物，或按照预先定义的 ligand whitelist/blacklist 排除。

Primary 和 strict-apo 结果方向应一致，但 strict subset 不承担主要功效。

---

# 7. 数据审计文件

每一个 protein–pocket pair 必须在 manifest 中保留：

```text
pair_id
uniprot_id
pocket_id
apo_pdb
holo_pdb
apo_chains
holo_chains
assembly_id
ligand_id
ligand_heavy_atoms
positive_residue_count
masked_other_pocket_residue_count
core_residue_count
membrane_flag
missing_pocket_residue_flag
strict_apo_flag
cryptobench_overlap_flag
sequence_cluster_id
inclusion_status
exclusion_reason
```

### 数据清洗阶段审计问题

1. 是否存在同一个 PDB/chain 被错误重复计数？
2. 是否存在同一个物理 pocket 被拆成多个伪独立样本？
3. 是否存在不同 pocket 被错误合并？
4. 其他已知 pocket 是否正确 MASK？
5. interface pocket 是否保留完整 assembly？
6. CryptoBench final-test homolog 是否全部排除？
7. exclusion reason 是否 100% 可追踪？

### 数据阶段验收标准

- [ ] 100% 样本具有唯一 `pair_id`；
- [ ] 100% positive residues 可从 holo 映射回 UniProt/apo；
- [ ] 每个 target pocket ≥5 positive residues；
- [ ] 其他已知 pocket residue 不作为 negative；
- [ ] CryptoBench sealed-test contamination = 0；
- [ ] 主分析膜蛋白 contamination = 0；
- [ ] assembly mismatch = 0；
- [ ] 每一条排除记录均有 reason code。

---

# 8. Residue-level Static 特征

主模型 Static baseline 必须足够强，避免 mechanics 只是补充最基本的局部结构上下文。

## 8.1 单残基特征

固定为：

1. **AA identity**：20 类氨基酸；
2. **RSA**：relative solvent accessibility；
3. **Depth**：residue depth；
4. **WCN**：weighted contact number；
5. **SS**：secondary structure；
6. **Radial position**：残基到结构质心距离，按蛋白尺度归一化；
7. **Terminal distance**：到最近 sequence terminus 的归一化序列距离。

其中 terminal distance 专门用于排除末端天然高柔性产生的伪相关。

## 8.2 邻域上下文特征

固定以 Cα 距离 **10 Å** 为邻域，计算：

- neighborhood mean RSA；
- neighborhood mean Depth；
- neighborhood mean WCN。

这三项进入 Static baseline。

它们的作用是回答：

> mechanics 的增量是否只是因为 mechanical features 天然进行了空间平滑和邻域整合？

## 8.3 Static spatial-smoothing control

GAM/LightGBM 得到 Static prediction \(p_i^{static}\) 后，再计算 10 Å 邻域平滑版本：

\[
\tilde p_i^{static}=\frac{1}{|N_i|}\sum_{j\in N_i} p_j^{static}
\]

主结果同时比较：

- Static；
- Static-smoothed；
- Static+Mechanics。

如果 Static+Mechanics 只优于 Static，却不优于 Static-smoothed，则不能声称 mechanical response 本身提供了独立物理增量，只能表述为“non-local neighborhood context improves prediction”。

---

# 9. B-factor 低成本替代对照

仅在适合比较的实验结构子集进行，例如 X-ray structures。

对每个蛋白内部将 B-factor 转换为 percentile 或 z-score，避免不同晶体实验的绝对量纲直接比较。

比较：

```text
Static
Static + normalized B-factor
Static + Mechanics
Static + normalized B-factor + Mechanics
```

该分析不属于 Pilot 硬 Gate，但在正式论文中用于回答：

> ANM-derived mechanics 是否只是实验 B-factor 的替代？

验收标准：若加入 B-factor 后 Mechanics 的 ΔAP 仍保持正向，则可写为 mechanics provides information beyond normalized experimental flexibility proxy；否则相应收缩表述。

---

# 10. ANM/PRS 主计算

## 10.1 固定参数

Primary ANM：

- coarse graining：Cα；
- cutoff：**10 Å**；
- spring constant：统一量纲下固定值；
- Hessian 去除 rigid-body zero modes；
- primary calculation 使用完整非刚体伪逆 \(H^+\)，不得在 Pilot 结果出来后调整 cutoff 追求最佳性能。

Cutoff sensitivity 仅在主分析完成后测试：

\[
8,10,12,15\ \text{Å}
\]

其中 10 Å 始终为预注册 primary result。

## 10.2 核心残基与末端处理

必须先定义 `core_residue_mask`，再计算 protein-internal percentile。

最低要求：

- 未解析的 terminal residues 不进入分析；
- 明确的长无序/低接触 terminal tails 按固定算法裁剪；
- 同一算法对全部数据统一执行；
- 原始 residue index 与裁剪后 index 必须可双向映射。

不得针对单个蛋白人工调整裁剪范围。

---

# 11. Mechanical Feature 1：%DFI

对扰动源 \(j\) 对目标残基 \(i\) 的响应定义 \(R_{ij}\)。对于各向同性单位力，可利用 \(H^+\) 的 3×3 block 范数计算响应，而不需要大量 Monte Carlo 力方向采样。

Residue-level raw DFI 按预先固定实现计算，之后只在 `core_residue_mask` 中转为蛋白内 percentile：

\[
\%DFI_i=\mathrm{PercentileRank}(DFI_i)
\]

### %DFI 项目审计

- residue 数是否与 mapping 一致；
- rigid modes 是否正确剔除；
- Hessian 是否对称；
- \(H H^+ H\approx H\) 数值误差是否在容许范围；
- 是否存在 NaN/Inf；
- percentile 是否只在核心残基计算；
- apo/holo 使用完全相同实现。

### 验收标准

- [ ] Hessian symmetry error < 预设数值容差；
- [ ] 非刚体有效 eigenvalues >0；
- [ ] 每个样本 mechanics vector 与 residue manifest 一一对应；
- [ ] NaN/Inf rate = 0；
- [ ] 结果可由固定 seed/environment 复现。

---

# 12. Mechanical Feature 2：Response Dispersion

排除 self-response：

\[
p_{ij}=\frac{R_{ij}}{\sum_{k\neq i}R_{ik}},\quad j\neq i
\]

定义归一化 response dispersion：

\[
D_i=-\frac{\sum_{j\neq i}p_{ij}\log p_{ij}}{\log(N-1)}
\]

之后在核心残基上转换为 protein-internal percentile。

物理解释严格限定为：

> **mechanical response concentration/dispersion**

不得解释为 directional signaling、sensor/effector asymmetry 或确定的 allosteric pathway。

## 12.1 共线性预检

在 Pilot 中必须报告：

- Spearman(%DFI, WCN)；
- Spearman(Dispersion, WCN)；
- Spearman(%DFI, Dispersion)；
- GAM/线性辅助模型的 VIF 或等价冗余指标。

若 Dispersion 与 WCN 出现极高冗余，例如 \(|\rho|>0.8\)，则：

- Dispersion 不再作为“独立机制”叙述；
- 可保留为 secondary feature；
- 主论文允许最终简化为 Static + %DFI。

这一变化不视为 protocol failure，因为该降级规则已预先写入。

---

# 13. Primary 模型一：GAM

## 13.1 Static GAM

形式示意：

\[
\mathrm{logit}[P(B_i=1)] =
\beta_{AA}+s(RSA)+s(Depth)+s(WCN)+s(Radial)+s(Terminal)
+\beta_{SS}
+s(NbRSA)+s(NbDepth)+s(NbWCN)
+te(RSA,Depth)
\]

其中 `te(RSA, Depth)` 作为最低成本的非线性交互防线。

## 13.2 Static + Mechanics GAM

增加：

\[
s(\%DFI)+s(Dispersion)
\]

若 Dispersion 在预定义共线性审计中被降级，则只增加：

\[
s(\%DFI)
\]

## 13.3 训练权重

每个 protein–pocket pair 总权重相等。

若 pair \(p\) 有 \(N_p\) 个参与训练的 residues：

\[
w_{pi}=1/N_p
\]

避免大蛋白主导 loss。

---

# 14. Confirmatory 模型二：LightGBM

LightGBM 不承担物理解释，但必须承担模型形式确认。

使用：

- 与 GAM 完全相同的 CV folds；
- 完全相同的 Static feature set；
- 完全相同的 Mechanics feature set；
- 仅在 development fold 内进行有限超参数选择；
- 不查看 CryptoBench test。

目的：

> 证明机械增量不是 GAM 函数形式或 smoothing specification 的产物。

---

# 15. Sequence-cluster Grouped CV

严禁 residue-random split。

## 15.1 Development clustering

对 AHoJ development proteins 使用 MMseqs2 聚类。

Primary development grouping 推荐：

- sequence identity ≤30% cluster separation；
- coverage rule 预先固定，例如 80%。

所有属于同一 cluster 的 protein–pocket pairs 必须进入同一 fold。

## 15.2 Fold 约束

同一：

- UniProt；
- sequence cluster；
- 高同源蛋白；

不得横跨 train/test fold。

## 15.3 审计

每次运行生成：

```text
splits/run_x/
├── train_clusters.tsv
├── validation_clusters.tsv
├── test_clusters.tsv
└── leakage_check.tsv
```

验收标准：

- [ ] UniProt overlap = 0；
- [ ] cluster overlap = 0；
- [ ] CryptoBench sealed-test overlap = 0；
- [ ] 所有模型使用完全相同 fold assignment。

---

# 16. Primary 评价指标

## 16.1 每个 protein–pocket pair 内单独评价

对每个 pair \(p\)：

\[
AP_p^{Static}
\]

和：

\[
AP_p^{Static+Mechanics}
\]

定义：

\[
\Delta AP_p=AP_p^{Static+Mechanics}-AP_p^{Static}
\]

主汇总：

- median \(\Delta AP\)；
- cluster-bootstrap 95% CI；
- win rate。

AUROC 为 secondary metric。

## 16.2 Win / Tie / Loss

设置预先固定 tolerance：

\[
\epsilon=0.005
\]

- Win：\(\Delta AP>0.005\)
- Tie：\(|\Delta AP|\le0.005\)
- Loss：\(\Delta AP<-0.005\)

win rate 在排除 ties 后计算，并附 cluster-bootstrap CI；可补充 sign test。

## 16.3 Cluster bootstrap

bootstrap 单位不是 residue，也不是单个 pocket，而是：

> **sequence cluster**

每次整体抽取 cluster，再汇总其中所有 protein–pocket pairs。

---

# 17. 全残基评价与困难评价

## 17.1 Primary：可评价残基范围

在 target pocket 的 apo protein 中：

- positive = target pocket residues；
- masked = 其他已知 pocket residues；
- negative = 其余可评价 residues。

## 17.2 Secondary：surface-only evaluation

为了避免“排除大量埋藏残基”本身使任务过于容易，增加：

> **surface-accessible residues only**

例如按固定 RSA threshold 定义表面残基，threshold 必须在 protocol 中提前固定。

## 17.3 Secondary：pocket-like candidate evaluation

如果可稳定运行 fpocket/P2Rank 等几何工具，则只在其候选 cavity/pocket lining residues 内重新评价。

目的不是把 P2Rank 当特征，而是问：

> 在都已经具有 pocket-like geometry 的候选区域中，mechanics 是否还能进一步重排真实口袋？

这部分属于 secondary analysis，不影响 Pilot 主 Gate。

---

# 18. 管线阴性对照 1：Mechanical Profile Swap

目的：验证 pipeline 不会因为 mechanics 的分布形式、蛋白长度或 percentile 本身产生伪增益。

## 方法

1. 按 core residue length 将 proteins 分 bin；
2. 在 bin 内随机交换 protein 的 %DFI/Dispersion profile；
3. donor profile 若长度略有差异，按归一化 sequence coordinate 插值到 target 长度；
4. 保留 target 的 Static features 和 labels；
5. 重复多个固定 seed，例如 100 次；
6. 完整重跑 Static + swapped Mechanics。

## 预期

\[
\Delta AP_{swap}\approx 0
\]

## 验收标准

- observed real-mechanics median \(\Delta AP\) 必须高于 profile-swap null distribution 的 95th percentile；
- swap-null 的中心应接近 0；
- 如果 swapped mechanics 也产生系统性正增益，则 pipeline 存在伪信号，立即 No-Go，先排查流程。

---

# 19. 管线阴性对照 2：Static Spatial Smoothing

比较：

```text
Static
Static-smoothed
Static+Mechanics
```

### 审计解释

- 若 Static+Mechanics > Static 且 > Static-smoothed：支持 mechanics 不只是邻域平滑；
- 若 Static+Mechanics ≈ Static-smoothed：机械增益可能主要来自非局部邻域上下文；
- 若 Static-smoothed > Static+Mechanics：不能主张 mechanics 有独立优势。

该结果决定论文措辞强度。

---

# 20. Pilot 0：200–300 个 protein–pocket pairs

Pilot 是整个项目最重要的 fail-fast 阶段。

## 20.1 Pilot 输入

从已剔除 CryptoBench final-test contamination 的 AHoJ development pool 中选择：

- 200–300 个 protein–pocket pairs；
- soluble proteins；
- small-molecule ligands；
- 每 pocket ≥5 positive residues；
- assembly matched；
- pocket mapping 完整；
- sequence-cluster grouped。

Pilot 不要求覆盖全部 AHoJ 分布，目标是验证主效应是否存在。

## 20.2 Pilot 只做以下项目

### Static

- AA；
- RSA；
- depth；
- WCN；
- SS；
- radial position；
- terminal distance；
- 10 Å neighborhood mean RSA/depth/WCN。

### Mechanics

- %DFI；
- response dispersion。

### Models

- GAM；
- LightGBM。

### Controls

- Static smoothing；
- mechanical profile swap；
- mechanics/WCN 共线性审计。

暂时不做：

- GNN；
- ProteinMPNN/ESM embeddings；
- DCI；
- MD；
- membrane imANM；
- web server；
- archetype clustering；
- 大规模 BioLiP2 扩展。

---

# 21. Pilot Go / No-Go 判据

## 21.1 Primary GAM Gate

必须同时满足：

### G1：统计方向

\[
95\%\ CI[median(\Delta AP_{apo})]\ \text{lower bound}>0
\]

### G2：最小实际效应

预注册 SESOI：

\[
median(\Delta AP_{apo})\ge0.01
\]

0.01 是本项目内部 minimum practically relevant effect，不宣称为领域通用标准。

### G3：蛋白级一致性

排除 ties 后：

\[
WinRate>50\%
\]

其中 tie tolerance = 0.005。

## 21.2 Confirmatory LightGBM Gate

要求：

\[
median(\Delta AP_{apo,LightGBM})>0
\]

且 cluster-bootstrap 95% CI 的下界 >0。

LightGBM 不强制达到 GAM 的 0.01 SESOI，但方向必须稳定为正。

## 21.3 Pipeline Integrity Gate

- real mechanics 的增益必须超过 profile-swap 95% null；
- Static+Mechanics 不得只表现为明显等价于简单 Static spatial smoothing；
- CryptoBench leakage audit 必须为 0。

## 21.4 Go

全部主要 Gate 通过：进入完整 AHoJ 扩展。

## 21.5 Conditional Go

如果 GAM 达到 SESOI，但 LightGBM CI 横跨 0，或 mechanics 与 Static-smoothed 差异很弱：

- 暂不扩大论文主张；
- 先检查样本量、fold variance、特征冗余；
- 允许增加 AHoJ development sample 做 confirmatory analysis，但不得改变特征定义和 CryptoBench test。

## 21.6 No-Go

以下任一情况出现：

- GAM median ΔAP <0.01；
- GAM 95% CI 下界 ≤0 且扩大合理样本量后仍然如此；
- LightGBM 方向与 GAM 相反；
- swapped mechanics 同样产生类似增益；
- 泄漏无法排除；

则停止“mechanics improves apo→holo prediction”的主故事。

可转为较弱问题：

- 描述 binding-site mechanics；
- mechanics 与 packing 的冗余；
- 某些特定 pocket classes 的力学异质性。

不得继续堆模型强行追求正结果。

---

# 22. Pilot 后功效分析

Pilot 结束后，使用实际获得的：

\[
SD(\Delta AP)
\]

和 cluster-level variance 估计完整 development study 所需样本量。

原则：

- 不使用简单 residue 数当独立样本量；
- 以 sequence cluster / protein–pocket 为统计单位；
- 保留 pilot 的所有数据和结果；
- Pilot 可继续作为 development 数据，但不能影响 CryptoBench sealed test。

所有功效估计写入：

```text
analysis/power_analysis_after_pilot.md
```

并记录此时的 Git commit。

---

# 23. Experiment 1：Holo 描述性证据

## 问题

> 已知 ligand-binding pocket 在去除 ligand 后，是否表现出 characteristic mechanical-response profile？

## 输入

AHoJ development pool 的 holo structures。

## 操作

1. 删除 ligand；
2. 使用蛋白 biological assembly；
3. 计算 Static；
4. 计算 %DFI/Dispersion；
5. 使用同一个 grouped-CV framework；
6. 比较 Static vs Static+Mechanics。

## 注意

这是描述性证据。

不得由 holo 结果直接推导：

> mechanics pre-exists before ligand binding。

因为 holo structure 已经历 ligand-associated structural selection/imprint。

## 验收标准

- median ΔAP >0；
- cluster-bootstrap CI 报告；
- 结果方向与 apo 主实验一致；
- 若 holo 强而 apo 不成立，只能解释为 binding-associated mechanical signature，不能写 pre-existing apo information。

---

# 24. Experiment 2：Apo → Holo 核心实验

这是全文最重要实验。

## 24.1 输入

- apo structure；
- 所有输入特征完全从 apo 计算；
- positive labels 完全由对应 holo pocket 映射得到。

即：

\[
X_{apo}\rightarrow y_{holo}
\]

## 24.2 主比较

\[
Static_{apo}
\]

vs

\[
Static_{apo}+Mechanics_{apo}
\]

主结果：

\[
\Delta AP_{apo}
\]

## 24.3 允许的核心结论

若通过：

> Ligand-free apo structures contain mechanical-response information associated with holo-defined binding regions beyond conventional local structural descriptors.

避免直接写：

- mechanics causes ligand binding；
- mechanics determines pocket formation；
- apo state proves causal preorganization。

---

# 25. Apo→Holo 构象变化分层

只保留一个简单分层，不做 archetype clustering。

对 target pocket 计算 apo→holo pocket structural change，例如 pocket Cα RMSD：

\[
RMSD_{pocket}^{apo\rightarrow holo}
\]

按 development dataset median 或预注册 quantile 分为：

- Low-change；
- High-change。

分别计算：

\[
median(\Delta AP_{Low})
\]

\[
median(\Delta AP_{High})
\]

如果 High-change 中 mechanics 增益仍然存在，尤其是高于 Low-change，则会形成重要机制证据：

> mechanics retains predictive information when apo local geometry differs substantially from the holo binding configuration.

该结果自然连接 cryptic-pocket validation。

---

# 26. 同一 pocket 多个 apo 的重复性审计

从 development pool 中抽取具有 ≥2 个独立 apo structures 的 protein–pocket pairs。

计算：

1. residue-wise Spearman(%DFI_apo1, %DFI_apo2)；
2. residue-wise Spearman(Dispersion_apo1, Dispersion_apo2)；
3. 两个 apo 对同一 holo pocket 的 \(\Delta AP\) 是否方向一致。

不设置人为 ICC 门槛。

目的只回答：

> 主效应是否严重依赖某一个特定 PDB snapshot？

若重复性很低，应报告其与 apo conformational heterogeneity 的关系，而不是简单删除不稳定样本。

---

# 27. Source-distance-resolved mechanical response

这一分析不称为 local-vs-global ANM。

因为 \(H^+\) 始终来自完整蛋白网络，任何 \(R_{ij}\) 都已经包含全局 elastic propagation。

它回答：

> **不同空间距离的 perturbation sources 对 binding-site mechanical information 贡献多少？**

## 27.1 壳层定义

建议固定：

```text
0–10 Å
10–20 Å
20–30 Å
30–40 Å
>40 Å
```

排除 self term。

## 27.2 壳层 Mechanical contribution

对 target residue \(i\)：

\[
M_i^{(a,b)}=\sum_{j:a<d_{ij}\le b}R_{ij}
\]

这里使用**求和**，不使用 shell mean，因为全 DFI 本身就是各 perturbation-source contribution 的累加，shell sum 才能形成可解释的分解。

## 27.3 同壳层 packing control

定义：

\[
WCN_i^{(a,b)}=\sum_{j:a<d_{ij}\le b}\frac{1}{d_{ij}^2}
\]

每个 shell 单独测试：

```text
Static + shell-WCN
vs
Static + shell-WCN + shell-Mechanics
```

不同 shell 不同时进入同一模型，避免壳层共线性。

## 27.4 壳层占据

若某 residue 在某 shell 的 source residue 数低于预设最小值，例如 <5：

- 该 residue 对该 shell 记 NA；
- 不进行人工补零；
- 每个 shell 报告有效 residue/pocket 数。

## 27.5 允许的结论

若中远距离 shell 在控制同距离 WCN 后仍有正 ΔAP：

> spatially distant perturbation sources contribute binding-site information beyond distance-dependent packing density.

不要写：

> proves a global allosteric network。

---

# 28. ANM cutoff sensitivity

Primary cutoff 始终：

\[
r_c=10\ \text{Å}
\]

正式主结果完成后再运行：

\[
8,10,12,15\ \text{Å}
\]

检查：

- %DFI residue rank correlation；
- median ΔAP；
- GAM partial-effect shape；
- apo→holo 主结论方向。

验收：

> 主要科学结论不应依赖单一 10 Å cutoff 才成立。

无需要求每个 cutoff 数值完全一致。

---

# 29. Experiment 3：CryptoBench 最终外部验证

CryptoBench 在此阶段之前不得用于调参。

## 29.1 Test A：Zero-shot transfer

训练：

\[
AHoJ_{development}
\]

直接测试：

\[
CryptoBench_{test}
\]

不进行再训练。

这是最重要的泛化测试之一。

问题：

> 在普通 AHoJ apo→holo relationships 中学习到的 mechanics increment，能否直接迁移到 cryptic sites？

## 29.2 Test B：Official benchmark setting

在 CryptoBench 官方 train split 上重新训练完全相同的：

- Static GAM；
- Static+Mechanics GAM；
- Static LightGBM；
- Static+Mechanics LightGBM。

只在官方 test split 上评价。

该设置用于与 published CryptoBench baselines 处于更接近的评价环境。

## 29.3 外部 baseline

P2Rank、PocketMiner 等作为：

> published / off-the-shelf reference baselines

而不是本项目模型输入。

在正式写作中必须明确这些工具的预训练数据与本研究 development set 并不完全一致，因此只作为参考，不把不同训练数据下的绝对分数差异解释成严格的模型优劣。

## 29.4 评估管线验证

在打开本项目最终结果之前，先使用公开输出或官方工具尽量复现至少一个 CryptoBench published baseline 的评价指标，以确认：

- residue mapping；
- labels；
- test subset；
- metric implementation；

没有错误。

---

# 30. CryptoBench 最终验收

核心目标不是必须打败 SOTA，而是确认主规律外推。

最低成功标准：

\[
median(\Delta AP_{CryptoBench})>0
\]

且：

- GAM 与 LightGBM 方向一致；
- AP 为主指标；
- AUROC 为 secondary；
- 同样报告 win/tie/loss；
- 所有结果明确区分 zero-shot 与 official-train setting。

若 Static+Mechanics 未击败某个复杂 published model，但稳定优于同架构 Static baseline，论文核心结论仍成立。

---

# 31. 完整阶段门与验收表

| Gate | 阶段 | 核心问题 | 必须输出 | 验收标准 | 未通过处理 |
|---|---|---|---|---|---|
| G0 | Protocol lock | 是否已防止分析后改规则 | 协议、Git hash、sealed list | 所有规则冻结 | 不得开始 Pilot |
| G1 | Data QC | 标签/assembly/泄漏是否正确 | pair manifest、exclusion log | leakage=0；mapping 100% 可追踪 | 修复数据，不训练 |
| G2 | Mechanics QC | ANM/DFI 是否数值可靠 | Hessian/eigen/QC report | 无 NaN；mapping 一致；矩阵检查通过 | 修复计算 |
| G3 | Pilot GAM | mechanics 是否有实际增量 | ΔAP、cluster-bootstrap | CI lower>0；median ΔAP≥0.01 | No-Go/收缩假设 |
| G4 | Pilot LightGBM | 是否依赖 GAM 形式 | LightGBM ΔAP | median>0 且 CI lower>0 | Conditional Go/扩样本确认 |
| G5 | Pipeline null | 是否存在伪信号 | profile-swap null | real > null 95th percentile | No-Go，查管线 |
| G6 | Spatial context | mechanics 是否只等于邻域平滑 | Static-smoothed 对照 | mechanics 仍有增量则支持强表述 | 收缩为 neighborhood-context 结论 |
| G7 | Full apo→holo | 主结论是否在扩展数据成立 | held-out grouped CV | ΔAP 正且稳定 | 不进入强结论 |
| G8 | Source distance | 中远距离扰动是否有贡献 | shell ΔAP + shell WCN | 至少部分远距离 shell 独立增量 | 不宣称 long-range contribution |
| G9 | CryptoBench | 是否外推 cryptic sites | zero-shot + official setting | Static+Mech > Static | 若失败，cryptic 结论删除 |

---

# 32. 项目审计体系

整个项目必须可以从最终 figure 反向追踪到原始 PDB/label。

建议目录：

```text
project/
├── protocol/
├── raw/
│   ├── ahoj/
│   └── cryptobench_sealed/
├── manifests/
│   ├── source_manifest.tsv
│   ├── pair_manifest.tsv
│   ├── residue_mapping.tsv
│   ├── exclusion_log.tsv
│   └── split_manifest.tsv
├── features/
│   ├── static/
│   └── mechanics/
├── models/
│   ├── gam/
│   └── lightgbm/
├── results/
│   ├── pilot/
│   ├── holo/
│   ├── apo_holo/
│   ├── source_distance/
│   └── cryptobench/
├── audit/
│   ├── run_manifest.tsv
│   ├── leakage_report.md
│   ├── numerical_qc.md
│   ├── deviations.md
│   └── reproducibility_report.md
└── figures/
```

## 32.1 每次 run 必须记录

```text
run_id
date
git_commit
environment_hash
data_manifest_hash
split_manifest_hash
feature_version
model_version
random_seed
command/config
output_path
status
notes
```

## 32.2 Protocol deviation

任何 protocol lock 后发生的变化都写入：

`audit/deviations.md`

包含：

- 改了什么；
- 为什么改；
- 是否在看见对应结果之前提出；
- 是否影响 final test；
- 哪些结果需重新计算。

不得静默修改。

---

# 33. 论文主图规划

## Figure 1 — Study design and data integrity

展示：

- protein–pocket pair 定义；
- apo/holo mapping；
- other-pocket MASK；
- AHoJ development vs CryptoBench sealed test；
- Static / Mechanics pipeline。

## Figure 2 — Mechanical signatures of ligand-binding regions

展示：

- %DFI；
- response dispersion；
- binding vs nonbinding；
- WCN/packing relationship；
- GAM partial-effect curves。

## Figure 3 — Mechanics beyond static structural context

展示：

- Static vs Static+Mechanics per-pair AP；
- ΔAP distribution；
- win/tie/loss；
- GAM + LightGBM；
- Static-smoothed control；
- profile-swap negative control。

## Figure 4 — Apo mechanics identifies holo-defined binding regions

展示：

- apo structure 上的 mechanics；
- holo ligand site overlay；
- apo→holo ΔAP；
- Low-change vs High-change pocket；
- repeated apo reproducibility。

## Figure 5 — Source-distance and cryptic-site generalization

A：source-distance shell ΔAP 与 same-shell WCN control；

B：CryptoBench zero-shot；

C：CryptoBench official train/test；

D：2–3 个 representative cryptic pocket cases。

---

# 34. Supplementary 最低内容

建议放入 Supplementary：

1. ANM cutoff = 8/10/12/15 Å sensitivity；
2. B-factor comparator；
3. strict-apo subset；
4. multiple-apo Spearman；
5. surface-only evaluation；
6. pocket-like candidate evaluation；
7. Dispersion/WCN 共线性；
8. LightGBM 详细参数；
9. inclusion/exclusion flowchart；
10. ligand class distribution；
11. shell occupancy statistics。

---

# 35. 不允许在主项目中途加入的内容

除非 Pilot 明确提出必须解决的问题，否则第一篇论文不加入：

- GNN；
- ProteinMPNN embedding；
- ESM/ESM-IF embedding；
- DCI；
- all-atom MD；
- imANM membrane model；
- web server；
- archetype clustering；
- 大规模 agent workflow；
- de novo pocket predictor。

这些均属于后续研究，不得用于“救”没有通过 Go/No-Go 的主假设。

---

# 36. 最短实际执行清单

## Step 1 — 封存 CryptoBench

输出：

- `cryptobench_test_sealed.tsv`
- test sequence FASTA；
- test cluster mapping；
- checksum。

**验收：** Development pool 与 sealed test overlap = 0。

---

## Step 2 — 构建 AHoJ protein–pocket pairs

输出：

- `pair_manifest.tsv`
- apo/holo mapping；
- ligand metadata。

**验收：** 每个 pair 有唯一 pocket identity。

---

## Step 3 — 生成 residue labels

输出：

- positive；
- negative；
- MASK；
- residue mapping table。

**验收：** 其他真实 pocket 不进入 negative。

---

## Step 4 — 生成 Static features

输出：

- AA/RSA/Depth/WCN/SS/Radial/Terminal；
- 10 Å neighborhood averages。

**验收：** missing rate 和异常值报告完成。

---

## Step 5 — ANM/PRS

固定：

\[
r_c=10\ \text{Å}
\]

输出：

- %DFI；
- dispersion；
- ANM QC。

**验收：** 数值 QC 全部通过。

---

## Step 6 — Pilot grouped CV

固定 sequence clusters。

运行：

- GAM Static；
- GAM Static+Mechanics；
- LightGBM Static；
- LightGBM Static+Mechanics。

**验收：** fold leakage = 0。

---

## Step 7 — Pipeline controls

运行：

- Static smoothing；
- mechanical profile swap；
- WCN/mechanics 共线性。

**验收：** real mechanics 不能被 null/smoothing 完全解释。

---

## Step 8 — Go/No-Go

核心数字：

\[
\boxed{\Delta AP_{apo}}
\]

Go 要求：

- GAM 95% CI lower >0；
- GAM median ΔAP ≥0.01；
- win rate >50%；
- LightGBM median ΔAP >0 且 CI lower >0；
- profile-swap null 不产生同等增益。

---

## Step 9 — Full AHoJ

只有 Go 后才扩展。

运行：

- holo descriptive；
- apo→holo primary；
- high/low conformational change；
- cutoff sensitivity；
- B-factor sensitivity；
- repeated apo reproducibility。

---

## Step 10 — Source-distance analysis

运行 shell：

```text
0–10
10–20
20–30
30–40
>40 Å
```

每个 shell：

```text
Static + shell WCN
vs
Static + shell WCN + shell mechanics
```

---

## Step 11 — Freeze final model

生成：

```text
final_model_config.yaml
final_feature_list.tsv
final_git_hash.txt
```

之后不得再使用 CryptoBench 结果调整模型。

---

## Step 12 — CryptoBench final test

先复核官方评价管线，再执行：

1. zero-shot transfer；
2. official train/test setting。

保存所有结果后才允许查看 final summary figure。

---

# 37. 项目最终验收标准

项目达到“完整 SCI 一区论文故事”至少需要：

## Level A：核心发现

- apo Static+Mechanics 在 strict grouped evaluation 下稳定优于 Static；
- median ΔAP 达到预注册实际效应阈值；
- GAM 与 LightGBM 方向一致；
- profile-swap null 不产生相同效果；
- Static spatial smoothing 无法完全解释 mechanics 增益。

## Level B：物理解释

- %DFI 的作用不能完全由 WCN/packing 解释；
- 至少部分中/远距离 perturbation-source shell 在 same-shell WCN 控制后仍提供信息；
- cutoff sensitivity 不改变主要结论。

## Level C：外部泛化

- CryptoBench sealed test 中 Static+Mechanics > Static；
- zero-shot 和/or official-train setting 中至少一种获得稳定正增量；
- 结果不依赖单一 pocket 或单一蛋白家族。

同时达到 A+B+C，可支撑完整标题级主张。

若只有 A：可发表为 apo mechanical signature / binding-site prediction study，但不强调 long-range/cryptic。

若 A+B：可形成较强 computational biophysics 论文，但删除 cryptic generalization。

若 A+C 但 B 不成立：可强调 predictive mechanical representation，但避免强物理机制表述。

---

# 38. 写作时允许和禁止的表述

## 可以写

- mechanical-response information beyond conventional local structural descriptors；
- apo-state mechanical signatures associated with holo-defined binding regions；
- spatially distant perturbation sources contribute predictive information；
- mechanics improves within-protein ranking of ligand-binding residues；
- mechanics provides complementary information to static pocket descriptors。

## 不应直接写

- mechanics contains information not present in structure；
- mechanics causes ligand binding；
- ANM proves the allosteric pathway；
- DFI uniquely determines ligandability；
- apo mechanics proves future pocket opening；
- the method is state of the art，除非严格同条件 benchmark 确实证明。

---

# 39. 项目冻结后的核心决策树

```text
Pilot apo ΔAP
    │
    ├── <0.01 / CI crosses 0
    │       │
    │       └── STOP expansion
    │           检查 packing redundancy / label / subgroup
    │
    └── ≥0.01 and CI>0
            │
            ▼
      LightGBM confirms?
            │
       ┌────┴────┐
       │         │
      No        Yes
       │         │
 Conditional    ▼
   validation   Null controls pass?
                 │
            ┌────┴────┐
            │         │
           No        Yes
            │         │
       Audit pipeline ▼
                  Full AHoJ
                      │
                      ▼
               Source-distance
                      │
                      ▼
               Freeze final model
                      │
                      ▼
               CryptoBench test
```

---

# 40. 最终研究思想

这项工作的价值不在于“DFI 能不能预测 binding site”这一单独现象。

真正需要建立的是一个更加严格的证据链：

\[
\text{Local static structure}
\]

并不能完全解释：

\[
\text{where ligand-binding residues occur}
\]

而完整蛋白弹性网络中的扰动响应：

\[
\mathbf X
\rightarrow
H
\rightarrow
H^+
\rightarrow
R
\rightarrow
DFI
\]

为这些位点提供了额外、可重复、可外推的表征。

最关键的证据不是 holo 中看到 DFI 差异，而是：

\[
\boxed{
Structure_{apo}
\rightarrow
Mechanics_{apo}
\rightarrow
BindingSite_{holo}
}
\]

并且这一增量在：

- sequence-cluster-disjoint evaluation；
- 强 Static baseline；
- neighborhood context control；
- LightGBM confirmatory model；
- profile-swap null；
- sealed CryptoBench cryptic-site test；

下仍然存在。

如果这条证据链成立，论文可以非常清楚地回答：

> **蛋白质 apo 构象中的弹性扰动响应并不是一个简单的局部柔性指标，而是一种能够补充传统静态几何描述、并与未来 holo binding region 相关联的物理表征。**

这就是本项目应该追求的最终 SCI 一区级结论。

---

# 41. 已核实的关键数据来源

1. **AHoJ-DB**  
   Feidakis CP, Krivak R, Hoksza D, Novotny M. *AHoJ-DB: A PDB-wide Assignment of apo & holo Relationships Based on Individual Protein–Ligand Interactions.* Journal of Molecular Biology, 2024, 436(17):168545.  
   DOI: 10.1016/j.jmb.2024.168545  
   PubMed: https://pubmed.ncbi.nlm.nih.gov/38508305/  
   论文报告：468,293 biologically relevant protein–ligand interactions、27,983 proteins；约 24% binding sites 位于多个 chains 界面；不到 50% processed binding sites 有 apo form。

2. **CryptoBench**  
   *CryptoBench: cryptic protein–ligand binding sites dataset and benchmark.* Bioinformatics, 2025, 41(1):btae745.  
   正式论文: https://pmc.ncbi.nlm.nih.gov/articles/PMC11725321/  
   正式论文报告：1,107 apo structures、1,361 cryptic pockets；正式 benchmark train/test = 885/222；构建 benchmark split 时采用额外 10% sequence-identity clustering 以降低同源泄漏。

3. **CryptoBench data/code**  
   数据和模型： https://osf.io/pz4a9/  
   代码和 tutorial： https://github.com/skrhakv/CryptoBench/  
   执行时以正式发布版本的 split 文件为唯一准绳，不采用二手论文中因过滤条件不同出现的 179、193 等子集数字替代官方 222 test structures。

---

# 42. 最终冻结声明

从本版本开始，不再增加新的主要方法。

项目下一步唯一任务是跑通：

```text
CryptoBench seal
    ↓
AHoJ pair construction
    ↓
Pocket-specific label + other-pocket MASK
    ↓
Static features
    ↓
ANM/PRS (%DFI + dispersion)
    ↓
Sequence-cluster grouped CV
    ↓
GAM + LightGBM
    ↓
Negative controls
    ↓
Go / No-Go
```

**真正决定项目是否继续的核心结果只有一个：**

\[
\boxed{
\Delta AP_{apo}
=AP_{Static+Mechanics}-AP_{Static}
}
\]

在严格无泄漏、蛋白内评价、cluster bootstrap 和阴性对照下，该量若具有稳定且达到预设最低实际效应的正增量，则进入完整研究；否则优先收缩科学主张，而不是增加更复杂模型追求正结果。

