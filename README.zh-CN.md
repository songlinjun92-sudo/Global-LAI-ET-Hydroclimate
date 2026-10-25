# 全球 LAI–ET 水热背景研究｜GitHub 开源代码

这是面向公开发布的模块化 Python 研究代码库。**GLEAM 是主蒸散产品，PML-V2.2b 用于 ET 产品独立验证**。核心研究期为 **1982–2020 年**，空间分辨率为 **1°×1°**。

核心模块包含：

1. **数据预处理接口**：已年度化、已对齐的 1° GeoTIFF 制作年度面板；也可以直接导入现有年度 CSV。
2. **PML 数据空间汇总**：读取官方年度 NetCDF 的 ET、Ec、Ei、Es，将 0.1° 有效像元聚合为 1° 年度数据。单位是 mm/year，不额外乘天数。
3. **RF 建模与验证**：LAI、Ta、P、Rs、VPD 五变量；气候模型对照、10° 地理块交叉验证、时间留出。
4. **局地 ICE**：历史 LAI 的 q05–q95、20 个点、同像元多年平均气候固定；只估计模型响应斜率，不作因果解释。
5. **AI/PFT 分析**：三套 AI，连续 AI 梯度、10° 地理块 Spearman/Bootstrap、PFT 内部分层、半干旱与湿润区的对称 Kitagawa 分解。
6. **PML 与 GLEAM 同样本比较**：ET↔ET、Ec↔T、Ei↔Ei、Es↔Eb，后者仅为最相近组分而非完全相同定义；使用匹配样本重新训练两套产品模型。
7. **结果制图与记录**：配对梯度图、PFT 热图、统计 CSV 和运行参数清单。

## 安装与测试

首先在 GitHub 下载或克隆本仓库，在仓库目录打开 Anaconda Prompt：

```bat
conda env create -f environment.yml
conda activate lai-et-reproducible
python -m pytest -q
python -m lai_et demo --root examples/demo_run
```

上面的 `demo` 自动生成纯合成数据，不使用你的真实数据，**其任何数值都不能写进论文**。测试成功后可删除 `examples/demo_run` 文件夹。

## 使用自己的论文数据

将 `configs/config.example.yaml` 复制为 `configs/config.yaml`。配置内**全部是相对仓库根目录的路径**，没有个人电脑专属盘符。

你至少需要准备：

- `data/processed/annual_panel.csv.gz`：一年一像元，含 `year,lon,lat,LAI_gs,Ta_ann,P_ann,Rs_ann,VPD_ann,ET_ann,T_ann,Ei_ann,Eb_ann`。
- `data/processed/pixel_attributes.csv`：一像元一行，含 `lon,lat,PFT,AI_A,AI_B,AI_C`。
- `data/processed/pml_panel.csv.gz`：独立 ET 验证使用，含 `year,lon,lat,PML_ET,PML_Ec,PML_Ei,PML_Es`。

如果原来的 GLEAM、ERA5、LAI 等数据已经有年度 1° TIFF，可以在配置中填写 `raster_inputs.annual`，然后运行：

```bat
python -m lai_et make-panel --config configs/config.yaml
```

如果下载了 PML 的 39 年年度 NetCDF，在数据文件夹下分别建立 `ET/Ec/Ei/Es` 子目录，并运行：

```bat
python -m lai_et prepare-pml --config configs/config.yaml
```

最后开展全部分析：

```bat
python -m lai_et run --config configs/config.yaml
```

程序会写出模型、ICE、AI 相关系数与 CI、PFT 对比、Kitagawa 分解、跨产品配对差异、图件和运行参数。完整字段约定见 `docs/DATA_CONTRACT.md`。

## 特别注意：本代码不等于已经无差异复现全部论文数值

当前已经通过自动化测试、合成数据全流程运行、原研究格式的小规模真实数据结构试运行。**没有在本次制作过程中用你完整的全部第三方原始栅格、39 年 PML NetCDF 及论文最终全部模型配置进行重跑**，因此不能承诺结果数值与已归档论文每一项统计完全相同。

对于正式投稿的代码可用性，建议先根据真实 V9/V16 研究版本确定软件/产品版本、LAI 生长季口径、模型训练规模、样本掩膜与参数，再用完整数据运行，并逐项核对论文的表格与图形。你公开上传的是可复现的算法代码，数据必须由使用者按数据源规定自行取得。

## 上传 GitHub 前检查

- 更换 `CITATION.cff` 中暂用的贡献者名称，并经实际合作者同意再公布代码许可。
- 检查配置文件有无账号信息、邮箱、个人路径，**不上传第三方受限原始数据**。
- 不上传 `outputs/`、`data/raw/` 和 `data/processed/`，仓库已提供 `.gitignore`。
- 阅读英文 `README.md`（面向国际读者）及 `docs/RELEASE_CHECKLIST.md`。
- 先上传后再在论文 `Code availability` 中写 GitHub URL；必要时通过 Zenodo 创建 DOI。
