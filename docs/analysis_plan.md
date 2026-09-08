# 统计分析计划与机器学习建模方案

## 1. 研究目标

### 中文

在 ICU 急性脑卒中患者中，评估代谢、炎症相关血液学指标、生命体征及衰弱相关共病负担联合预测院内 30 天死亡的能力，并比较模型在 MIMIC-IV v3.1 与 eICU-CRD v2.0 中的可迁移性。

### English

To evaluate the ability of combined metabolic, inflammation-related hematologic, physiologic, and frailty-related comorbidity features to predict 30-day in-hospital mortality among critically ill patients with acute stroke, and to assess model transportability between MIMIC-IV v3.1 and eICU-CRD v2.0.

## 2. 研究队列与结局

| 分析 | 数据集 | 样本量 | 事件数 | 预测窗口 | 结局 |
|---|---|---:|---:|---|---|
| 主模型开发 | MIMIC-IV | 7,245 | 1,194 | ICU 入科后 0–24 h | 24 h landmark 后至 30 天院内死亡 |
| 外部验证 | eICU-CRD | 4,264 | 509 | ICU 入科后 0–24 h | 24 h landmark 后至 30 天院内死亡 |
| 早期敏感性模型 | MIMIC-IV/eICU-CRD | 7,592/4,502 | 1,474/591 | ICU 入科后 0–6 h | ICU 入科后 30 天院内死亡 |

0–24 小时 landmark 合格队列：MIMIC 7,247 例、eICU 4,264 例；排除结局缺失记录后，实际进入主模型的记录分别为 7,245 例和 4,264 例，对应结局事件分别为 1,194 和 509 例。MIMIC 中的 1,204 例是全队列 `death_24h_to_30d_inhospital` 事件数，其中 10 例不属于满足 24 小时 landmark 条件的记录；另有 2 例结局缺失，因此均不能进入主 landmark 模型。早期敏感性模型的可分析记录数为 MIMIC 7,592 例和 eICU 4,502 例。

主要建模方向为“ MIMIC-IV 开发 → eICU-CRD 外部验证”。反向“eICU 开发 → MIMIC 验证”作为预先规定的可迁移性敏感性分析，不把两个数据库随机混合后称作外部验证。

## 3. 数据划分与验证策略

### 3.1 MIMIC-IV 开发集

- 按患者而不是按住院记录随机划分，避免同一患者的多次住院跨越训练集和测试集。
- 80% 作为模型开发集，20% 作为一次性内部测试集；按结局和卒中表型分层。
- 在 80% 开发集中使用 5 折分组交叉验证进行超参数选择、插补参数估计和概率校准。
- 内部测试集只在模型和参数完全固定后使用一次。
- 固定随机种子并在代码、补充材料和分析日志中记录。

### 3.2 eICU-CRD 外部验证集

- 使用全部符合条件的 eICU landmark 队列进行外部验证，不重新训练模型。
- 直接应用 MIMIC 开发阶段固定的变量转换、插补参数、标准化参数和分类阈值。
- 允许进行“仅截距再校准”作为敏感性分析，但必须同时报告未经再校准的原始外部验证结果。

### 3.3 24 小时模型的严格时间逻辑

24 小时模型仅纳入 `landmark_24h_eligible = 1` 的患者，并将 `death_24h_to_30d_inhospital` 作为结局。不能将 0–24 小时特征与从 ICU 入科起计算的完整 30 天死亡直接配对作为主分析，否则早期死亡患者会因无法完成 24 小时观察而被系统排除，且预测时点与结局起点不一致。

## 4. 预先规定的特征集合

### 4.1 跨数据库主模型

- 人口学：年龄、性别。
- 卒中表型：AIS、ICH、SAH；混合/模糊表型合并为 `mixed_or_ambiguous`，避免极少数 CVST 混合记录造成不稳定估计。
- 衰弱相关代理：9 项慢性共病构成的 `frailty_comorbidity_count`。
- 代谢/肾功能及血液学：glucose、sodium、potassium、chloride、bicarbonate、BUN、creatinine、total calcium、WBC、platelets、hemoglobin、hematocrit。
- 生命体征：heart rate、respiratory rate、SpO2、SBP、DBP、MAP。

### 4.2 时间窗口对应的特征

- 0–6 小时模型使用各变量的 `first_6h` 及对应的 `missing_6h`。
- 0–24 小时 landmark 模型使用各变量的 `mean_24h` 及对应的 `missing_24h`。
- 0–24 小时变量只在 landmark 队列中使用；对未达 landmark 的记录，相关值在分析数据中保持为空。
- `subject_id`、`hadm_id`、`stay_id`、`patientunitstayid`、数据库来源及任何结局字段只用于分组、审计或验证，不作为预测特征。

### 4.3 不进入主模型的变量

CRP、乳酸、白蛋白、胆红素、镁、磷、温度、中性粒细胞比例/绝对计数、淋巴细胞比例/绝对计数不进入跨数据库主模型。它们仅用于预先规定的敏感性分析，原因是覆盖率低、单位或构造定义不完全一致，或存在明显的数据库记录差异。

## 5. 数据预处理

1. 在正式建模前完成单位审计；所有跨数据库共同变量必须转换为同一临床单位。
2. 根据变量字典预先规定生理合理范围；超出范围的值设为缺失，不根据结局进行删值或截尾。
3. 连续变量使用 MIMIC 开发集内训练折估计的中位数进行插补，并保留缺失指示变量；验证集和 eICU 只能使用训练阶段得到的插补参数。
4. 对弹性网逻辑回归，在训练折内完成中心化和尺度标准化；XGBoost 使用相同的已固定输入版本以便比较。
5. 性别设置为统一的二元类别并保留未知类别；年龄异常或缺失采用训练集规则处理。
6. 不在主分析中使用 SMOTE 或随机过采样。两库事件比例并非极端失衡，过采样还可能破坏概率校准；如需要，类别权重仅作为敏感性分析并重新校准概率。
7. 插补、标准化、特征筛选和超参数调优全部限制在训练折内，防止数据泄漏。

## 6. 模型设置

### 6.1 主机器学习模型

采用 gradient-boosted decision trees（XGBoost）作为主机器学习模型，原因是其适合结构化临床数据、可处理非线性和变量交互，并能输出 SHAP 解释。

### 6.2 比较模型

- 弹性网逻辑回归：作为可解释的临床统计基线。
- 随机森林：作为非线性树模型敏感性比较。

只比较上述预先规定的三类模型，不根据测试集性能临时增加或删除模型。

### 6.3 调参与校准

- XGBoost 和随机森林的超参数在 MIMIC 开发集的 5 折交叉验证中选择。
- 评价函数以 AUROC 与 AUPRC 的综合表现为主，同时监测 Brier score，避免只优化区分度而损害校准。
- 当前实现不进行额外的 Platt scaling 或 isotonic calibration；内部测试集和 eICU 外部验证报告模型原始概率的 Brier score、校准截距、校准斜率及校准图。若后续加入再校准，必须在训练折内拟合，并在内部测试集和 eICU 上固定应用，不能根据外部验证结果重新选择。
- 不使用测试集或 eICU 外部验证集选择特征、超参数或阈值。

## 7. 性能评价

### 7.1 主要评价指标

- 区分度：AUROC 和 AUPRC。
- 校准：Brier score、校准截距、校准斜率和校准图。
- 临床效用：decision-curve analysis，在预先报告的风险阈值范围内计算 net benefit。

### 7.2 补充评价指标

在固定于开发阶段的阈值下报告 sensitivity、specificity、PPV、NPV 和 F1 score。若需要二分类阈值，使用训练/验证集确定的 Youden 阈值，不能在内部测试集或外部验证集上重新寻找最佳阈值。

### 7.3 置信区间

使用患者级 bootstrap 计算 95% 置信区间，建议 1,000 次重抽样。eICU 外部验证时以 `patienthealthsystemstayid` 为重抽样单位，MIMIC 时以 `subject_id` 为分组单位。

## 8. 解释性与亚组分析

- 对 XGBoost 报告全局 SHAP 重要性、SHAP 分布图和主要变量的非线性关系。
- 对弹性网逻辑回归报告标准化系数及方向，但不将预测模型系数解释为因果效应。
- 预先规定 AIS、ICH、SAH 三个严格亚型；混合/模糊表型纳入总体模型但不纳入互斥亚型模型。
- 补充按数据库、年龄（例如 <65/≥65）和性别分层的性能评价；亚组结果以估计值和置信区间为主，不进行无计划的结果导向筛选。

## 9. 敏感性分析

1. 0–6 小时早期模型与 0–24 小时 landmark 模型比较。
2. 严格亚型队列与包含混合/模糊表型的总体队列比较。
3. MIMIC 开发/eICU 验证与反向 eICU 开发/MIMIC 验证比较。
4. 训练集插补、完整案例分析和不同缺失指示策略比较。
5. 仅 MIMIC 中加入 GCS、temperature 和 HFRS/FI-Lab 的增量价值分析；不把这些结果描述为跨数据库外部验证。
6. 用单个共病指示变量替代 `frailty_comorbidity_count`，评估衰弱相关代理构造方式对性能和校准的影响。

## 10. 关键偏倚控制

- 排除或单独报告 `outcome_missing_flag = 1` 的记录。
- 不使用出院状态、住院死亡标志、死亡时间、住院时长或任何结局字段作为特征。
- 不使用 24 小时之后的数据预测 24 小时 landmark 结局。
- 不将 MIMIC 的出院诊断共病编码解释为完全确定的入 ICU 前基线状态；该指标在论文中称为衰弱相关共病负担代理，并在局限性中说明诊断时序不确定性。
- 以开发集固定规则处理异常值、单位和缺失，不根据外部验证集表现调整规则。

## 11. 可直接用于论文 Methods 的英文段落

The primary model will be developed in MIMIC-IV v3.1 and externally validated in eICU-CRD v2.0. For the 24-hour landmark analysis, only patients who remained under observation for at least 24 hours after ICU admission will be included, and the outcome will be death occurring from 24 hours after ICU admission through 30 days. The prespecified predictors will include age, sex, stroke phenotype, a frailty-related comorbidity burden proxy, 12 harmonized laboratory variables, and six vital-sign variables. For the early prediction analysis, the first valid measurement within 6 hours will be used; for the landmark analysis, the mean value within the first 24 hours will be used. Missingness indicators will be retained, and continuous variables will be imputed using training-fold medians only. An XGBoost model will serve as the primary machine-learning model, with elastic-net logistic regression and random forest as comparator models. Hyperparameters and calibration procedures will be selected using grouped five-fold cross-validation within the MIMIC development set. Model performance will be assessed using the AUROC, AUPRC, Brier score, calibration intercept and slope, calibration plots, and decision-curve analysis. All preprocessing and model-selection steps will be performed without access to the internal test set or the external validation data.

## 12. 下一步

下一步为编写可复现的 R/Python 建模脚本，读取两库导出的分析数据集，完成患者级划分、训练折内插补与调参、内部测试、eICU 外部验证、校准和 SHAP 分析。
