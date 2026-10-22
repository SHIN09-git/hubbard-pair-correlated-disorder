**PHY-PAIR-02 结项总结**

结项日期：2026-10-22。作者：JIN HANZHE。

项目按后续明确的“核心贡献审查＋最小研究稿＋可移植公开交付”范围结项。研究、限定范围的文献对照、论文与复现材料均已交付；成果现为公开预印本。J. Phys. A 投稿已被拒，期刊录用未实现，不作为本次结项的完成项。

**主要成果**

- 在一维 Bose–Hubbard 双玻色子、共同弱静态相关势模型内，给出二阶解离核的闭式表达及正性证明，整理态归一化、交换对称、动量折叠、连续谱端点和特殊极限的核验。
- 得到解离与弹性反散射的关闭条件、参数分类，以及“反散射关闭但解离仍开放”的窗口；推导两种谱归一化下的端点幂律和固定方差平台。所有关闭结论均限于无序强度的二阶系数，不等于全阶稳定。
- 完成含证明与检查附录的 18 页英文稿。有限环数值支持固定方差与固定通带高度可使解离响应朝相反方向变化；与既有公式的共有因子和本模型的具体结果已区分，不将黄金律、谱筛选或平方根边缘本身归为首创。

公开发布前，保存数据复现、24 项单元检查、25 点散射态／预解式对照及论文重编译通过；25 点最大相对差约为 6.90×10⁻⁵。公开 PDF 和源码 ZIP 与本地定稿的文件哈希一致。上述检查的范围见[公开核验记录](https://github.com/SHIN09-git/hubbard-pair-correlated-disorder/blob/main/VALIDATION.md)。

**公开链接**

- 论文预印本（版本 1.0）：[10.5281/zenodo.23552352](https://doi.org/10.5281/zenodo.23552352)。
- GitHub 仓库：[hubbard-pair-correlated-disorder](https://github.com/SHIN09-git/hubbard-pair-correlated-disorder)。论文 DOI 已写入 README、CITATION.cff 和版本下载页。
- 数据与代码存档：[10.5281/zenodo.23390681](https://doi.org/10.5281/zenodo.23390681)。该 DOI 与论文 DOI 分开引用。

**复现入口**

从 [v1.0.0 下载页](https://github.com/SHIN09-git/hubbard-pair-correlated-disorder/releases/tag/v1.0.0)下载完整仓库源码，先阅读[复现说明](https://github.com/SHIN09-git/hubbard-pair-correlated-disorder/blob/main/reproducibility/README.md)和[论文—文件—命令对应表](https://github.com/SHIN09-git/hubbard-pair-correlated-disorder/blob/main/reproducibility/PAPER_TO_FILES.md)。使用 Python 3.12 的新虚拟环境，从解压后的仓库根目录执行：

```text
cd reproducibility
python -m pip install -r requirements.txt
python reproduce.py --target verify --output reproductions/closeout_integrity
python reproduce.py --target all --output reproductions/closeout_saved_data
```

每次使用新的输出目录。默认入口核验文件并从保存数据重建表格、统计、检查和图，不重跑全部原始动力学，也不产生新的独立样本。完整动力学重放另见复现说明中的 `REPRODUCE_RESEARCH.py`。本次结项整理未新增科学计算。

**仍未解决的问题**

- 未求得非零高阶解离率、有限无序强度下的长期寿命或端点联合极限；高阶动量支撑允许过程发生，不能证明其振幅非零。
- 接近关闭边界时，有限环再现妨碍渐近率的稳定估计。有限尺寸、离散谱与场实现的影响尚未完全分离；事后指数替代模型尚无受控重求和推导。无限链主例的二阶时间窗偏差虽小于 0.01%，不能据此排除有限环的谱与时间窗耦合。
- Li–He（2018）保留摘要级全文缺口；Dufour–Orso 投影方程到本文核的完整推导对照未完成。先行文献审查有明确边界，不能据此断言穷尽性首创。

以上事项作为后续研究问题归档，不自动扩大本次已收窄的项目范围。当前没有待完成的公开交付；后续仅在发现具体错误或作者决定启动新研究、再投稿时重开。

公开核验记录：[发布核验说明](VALIDATION.md)；论文及源码的文件校验值可在 [Zenodo 预印本记录](https://zenodo.org/records/23552352)查看。
