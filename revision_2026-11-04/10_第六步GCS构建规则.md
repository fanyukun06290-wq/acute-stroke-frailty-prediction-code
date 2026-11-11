# 第六步：严格时间窗的GCS记录构建与质量门槛

本规则在新评分提取、构建及任何相关模型拟合之前保存。当前步骤只构建评分及检查可用性，不根据外部结局选择评分算法，不训练或替换主模型。

## 时间与来源

- 保持既有6小时和24小时合格队列、患者分组及结局定义不变。
- MIMIC使用原始chartevents的220739（眼）、223900（言语）、223901（运动）。charttime及storetime均须在landmark之前，charttime不早于入ICU；拒绝storetime早于charttime的异常记录。
- eICU使用physicalExam的GCS路径，physicalExamOffset须在[0,landmark)。该offset表示录入时间，不能证明实际检查时间完全同步。只从严格数字分量路径解析分数；不从自由文字猜分。
- 仅使用实际ICU停留内的记录。记录均来自先前已通过ICU和医院landmark资格检查的首个stay，仍核对有效ICU结束时间。
- 不使用MIMIC首日汇总表或eICU的APACHE日汇总分量代替严格时窗。

## 分量、总分与选择

- 眼1–4、言语1–5、运动1–6，须为整数。相同记录时间和分量的重复相同值去重；多个不同值或同时有不可评估标志时，拒绝相应分量。
- 完整总分只相加相同MIMIC charttime或相同eICU录入offset的三种完整、有效分量。不跨时间前向填充，不使用各分量在整个窗口的最小值拼总分。
- 同一eICU时间点若有直接记录总分且与分量和冲突，拒绝该时间点的完整总分；直接总分不单独补足缺失分量。相同时间点可能合并多个检查，是剩余局限，不能称已证明同一次床旁检查。
- 保留窗口内最低有效完整GCS及最低有效运动分量；同分取更晚记录。另保留首次和最后一次有效记录，便于审计及后续敏感性分析。窗口内最低值是修回中明确规定的构建规则，不宣称为原研究预设。

## 不可评估、气道与镇静

- MIMIC的No Response-ETT即使数值为1，也不能作为有效言语1分；ETT、Trach、unable、unassessable、non-testable等显式文字标志使该分量不可用，不替换为1或15。
- eICU同一录入offset标有estimated due to meds、unable to score due to meds或sedated时，拒绝该时点的总分与运动分量。
- 额外提取eICU respiratoryCare中在landmark之前录入的非空airwayType。只允许用检查之前或同一分钟、最近120分钟以内的气道记录。相同最新offset若有矛盾气道类型，视为未知；不能使用未来记录。
- eICU同期记录人工气道时，不使用完整GCS总分。言语1分而气道未知时，也不使用完整总分；明确No Artificial Airway时才允许。该保守规则不能消除所有未记录的插管影响。
- 运动分量不因言语不可评估而自动失效，但显式药物/镇静无法评分标志仍排除。它不能称为完整GCS或无镇静影响的神经功能测量。
- 两库对镇静、插管的记录能力不同：不能把没有显式标志解释为未镇静，不能声称完整GCS已实现完全等价的无镇静评估。因此优先考察运动分量的神经功能对照；完整总分仅作为严格受限的补充候选，是否进入后续模型要依据质量审计，不以外部性能挑选。

## 输出与范围

患者级事件、候选检查和stay级评分仅写入private_data。公开可读结果只包含汇总可用性、拒绝原因、规则、核验和文件哈希。评分缺失保持缺失，不按正常或无反应处理，不改变原队列。

该步骤不代表已经完成NIHSS、SOFA、APACHE等验证预测模型比较。下一步才能按冻结训练划分建立神经功能参考并评价增益，仍须区分参考变量构建、概率模型拟合和标准评分验证。

来源：[GCS官方评估说明](https://www.glasgowcomascale.org/downloads/GCS-Assessment-Aid-English.pdf)、[eICU physicalExam](https://eicu.mit.edu/eicutables/physicalexam/)、[eICU respiratoryCare](https://eicu.mit.edu/eicutables/respiratorycare/)。
