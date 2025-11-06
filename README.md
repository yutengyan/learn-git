# PtSnO/Al₂O₃ Trajectory Analysis Tool

## 📖 简介

这是一个专门用于分析负载型 PtSnO 团簇在 Al₂O₃ 载体上的分子动力学轨迹的 Python 工具。

### ✨ 主要功能

1. **自动识别原子组** - 基于 z 坐标自动区分载体（Support）和团簇（Cluster）
2. **自动 Unwrap** - 处理周期性边界条件导致的坐标跳跃
3. **团簇重组** - 修复被周期性边界"撕裂"的团簇结构
4. **双模式居中** - 生成两种不同居中方式的轨迹用于不同类型的分析
5. **灵活的文件命名** - 支持完全自定义所有输出文件名
6. **原始轨迹保存** - 可选保存 wrapped 原始坐标作为参考

---

## 🚀 快速开始

### 安装依赖

```bash
conda activate vscode-1  # 或你的环境名称
# 需要的包：numpy, scipy
```

### 基本用法

```bash
# 最简单的用法（自动检测团簇原子数）
python analyze_ptsn_trajectory.py trajectory.xyz

# 指定载体原子数（默认 240）
python analyze_ptsn_trajectory.py trajectory.xyz --support 240

# 手动指定团簇原子数
python analyze_ptsn_trajectory.py trajectory.xyz --support 240 --cluster 14

# 自定义输出文件名前缀
python analyze_ptsn_trajectory.py trajectory.xyz --output my_run_

# 完全自定义文件名 + 保存原始轨迹
python analyze_ptsn_trajectory.py sampling-simply.xyz \
    --surface-out sampling-simply_ase_unwrapped.xyz \
    --cluster-out test_cluster_centered.xyz \
    --save-wrapped \
    --wrapped-out sampling-simply_wrapped_original.xyz
```

---

## 📝 命令行参数

### 必需参数

| 参数 | 说明 |
|------|------|
| `trajectory` | 输入的 XYZ 轨迹文件 |

### 可选参数

| 参数 | 缩写 | 类型 | 默认值 | 说明 |
|------|------|------|--------|------|
| `--support` | `-s` | 整数 | 240 | 载体原子数量 |
| `--cluster` | `-c` | 整数 | 自动计算 | 团簇原子数量（默认 = 总原子数 - 载体数）|
| `--output` | `-o` | 字符串 | 无 | 输出文件名前缀（当不使用下面的自定义名称时）|
| `--surface-out` | - | 字符串 | `surface_centered.xyz` | 载体居中轨迹的文件名 |
| `--cluster-out` | - | 字符串 | `cluster_centered.xyz` | 团簇居中轨迹的文件名 |
| `--index-out` | - | 字符串 | `index_zsplit.ndx` | 索引文件的文件名 |
| `--save-wrapped` | - | 开关 | False | **是否保存原始 wrapped 轨迹** |
| `--wrapped-out` | - | 字符串 | `wrapped_original.xyz` | 原始 wrapped 轨迹的文件名 |
| `--help` | `-h` | - | - | 显示帮助信息 |

---

## 📂 输出文件

### 1️⃣ **索引文件** (`index_zsplit.ndx`)

GROMACS 风格的原子组索引文件，包含三个组：

```
[ Support ]         - 载体原子（Al₂O₃，240 个原子）
[ PtSnOCluster ]    - 团簇原子（PtₓSnᵧOᵧ，自动检测数量）
[ System ]          - 全部原子
```

#### 🔍 **索引文件如何生成？**

1. **读取第一帧**的所有原子坐标
2. **按 z 坐标从低到高排序**所有原子
3. **前 `--support` 个原子**（z 坐标最低）→ Support 组
4. **剩余的原子**（z 坐标最高）→ PtSnOCluster 组
5. **原子 ID 从 1 开始**（符合 GROMACS 约定）

**为什么按 z 坐标分组？**
- 载体（Al₂O₃）在下方，z 坐标较小
- 团簇（PtSnO）负载在表面，z 坐标较大
- 这是基于物理结构的自然分组方式

**用途：**
- ✅ 配合 GROMACS/MDAnalysis 等工具进行分析
- ✅ 定义原子选择
- ✅ 计算分组的物理量（如质心、RMSD 等）

---

### 2️⃣ **载体居中轨迹** (`surface_centered.xyz`)

将载体（Support）的质心固定在盒子中心的轨迹。

**特点：**
- 载体位置固定
- 可观察团簇相对于载体的运动

**用途：**
- ✅ 分析团簇在表面的**迁移行为**
- ✅ 计算团簇质心的扩散系数
- ✅ 观察团簇在载体表面的吸附位点变化
- ✅ 可视化团簇的宏观运动

---

### 3️⃣ **团簇居中轨迹** (`cluster_centered.xyz`)

将团簇（Cluster）的质心固定在盒子中心的轨迹。

**特点：**
- 团簇质心位置固定
- 可观察团簇内部原子的相对运动

**用途：**
- ✅ 分析团簇内部的**原子重排**
- ✅ 计算团簇的 RMSD 和形状变化
- ✅ 观察 Pt-Sn-O 键长、键角变化
- ✅ 研究团簇的结构稳定性

---

### 4️⃣ **原始 Wrapped 轨迹** (`wrapped_original.xyz`，可选）

保存的原始坐标（周期性边界条件限制在盒子内）。

#### ❓ **`--save-wrapped` 是原始的 xyz 文件吗？**

**✅ 是的！但有细节：**

- **内容：** 与你输入的 `trajectory.xyz` **完全相同**（除非你的输入文件已经部分 unwrapped）
- **坐标状态：** **Wrapped** - 所有原子坐标都在 `[0, L]` 范围内（L 是盒子尺寸）
- **何时生成：** 仅当使用 `--save-wrapped` 参数时生成

**流程：**
```
输入文件 → 读取 → [如果 --save-wrapped] → 立即保存原始坐标
              ↓
           Unwrap + 团簇重组 + 居中
              ↓
           输出处理后的轨迹
```

**用途：**
- ✅ 作为**对照参考** - 对比 unwrap 前后的差异
- ✅ **验证处理正确性** - 确认 unwrap 没有引入错误
- ✅ **备份原始数据** - 保留未处理的原始轨迹

**举例：**
```bash
python analyze_ptsn_trajectory.py sampling-simply.xyz \
    --save-wrapped \
    --wrapped-out original_backup.xyz
```
生成的 `original_backup.xyz` 与 `sampling-simply.xyz` 内容相同。

---

### 📋 **文件命名方式**

#### 方式 1：使用默认文件名
```bash
python analyze_ptsn_trajectory.py trajectory.xyz
```
生成：
- `index_zsplit.ndx`
- `surface_centered.xyz`
- `cluster_centered.xyz`

#### 方式 2：使用前缀
```bash
python analyze_ptsn_trajectory.py trajectory.xyz --output run1_
```
生成：
- `run1_index_zsplit.ndx`
- `run1_surface_centered.xyz`
- `run1_cluster_centered.xyz`

#### 方式 3：完全自定义文件名
```bash
python analyze_ptsn_trajectory.py sampling-simply.xyz \
    --surface-out unwrapped_ase.xyz \
    --cluster-out centered_cluster.xyz \
    --index-out atoms_groups.ndx \
    --save-wrapped \
    --wrapped-out original.xyz
```
生成：
- `atoms_groups.ndx`
- `unwrapped_ase.xyz`
- `centered_cluster.xyz`
- `original.xyz`

---

## 🔧 工作原理

### 1. 读取轨迹
- 支持扩展 XYZ 格式（包含 Lattice 信息）
- 自动提取盒子尺寸和原子坐标

### 2. 原子分组
```
按 z 坐标排序：
┌─────────────────┐
│  z 最大的原子    │ → Cluster (PtSnO 团簇)
├─────────────────┤
│  z 最小的原子    │ → Support (Al₂O₃ 载体)
└─────────────────┘
```

### 3. 团簇重组（Make Whole）
使用迭代质心法修复第一帧被周期性边界撕裂的团簇：
```python
迭代过程：
1. 计算团簇质心
2. 对每个原子应用最小镜像约定
3. 将远离质心的原子"拉回"
4. 重复直到收敛（通常 2-4 次迭代）
```

### 4. 轨迹 Unwrap
```python
对每一帧：
  如果原子位移 > 盒子尺寸/2：
    → 这是跨盒跳跃，平移回来
```

### 5. 居中处理
- **载体居中**：`shift = box_center - support_COM`
- **团簇居中**：`shift = box_center - cluster_COM`

---

## 📊 使用示例

### 示例 1：基本分析
```bash
python analyze_ptsn_trajectory.py sampling-simply.xyz
```

**输出：**
```
======================================================================
  Sn7Pt4O3/Al2O3 Trajectory Processing - Dual-Mode Centering
======================================================================
Reading trajectory: sampling-simply.xyz
Auto-detected: 254 total atoms = 240 support + 14 cluster
Successfully read 1751 frames, 254 atoms per frame

Identifying atom groups based on z-coordinates...
   Support (240 atoms):
      Al: 96
      O: 144
   Cluster (14 atoms):
      Pt: 4
      Sn: 7
      O: 3

Reassembling cluster in first frame (making whole)...
   Converged after 4 iterations
   First frame cluster max internal distance: 8.71 A
   OK: Cluster is now compact

...处理完成...
```

### 示例 2：批量处理多个轨迹
```bash
# 处理多个轨迹，使用不同的输出前缀
python analyze_ptsn_trajectory.py run1.xyz --output run1_
python analyze_ptsn_trajectory.py run2.xyz --output run2_
python analyze_ptsn_trajectory.py run3.xyz --output run3_
```

### 示例 3：不同尺寸的系统
```bash
# 200 个载体原子的系统
python analyze_ptsn_trajectory.py small_system.xyz --support 200

# 300 个载体原子 + 20 个团簇原子
python analyze_ptsn_trajectory.py large_system.xyz --support 300 --cluster 20
```

---

## ⚠️ 注意事项

### 输入文件要求
1. **XYZ 格式**：必须是扩展 XYZ 格式，包含 `Lattice` 信息
2. **原子顺序**：假设 z 坐标低的是载体，高的是团簇
3. **完整性**：所有帧的原子数必须相同

### 常见问题

#### 1. 轨迹未正确 unwrap？
**现象：**
```
ERROR: Detected 555 frames with possible box-crossing jumps!
WARNING: Trajectory may not be properly unwrapped
```

**解决：**
- 脚本会**自动 unwrap**，无需手动处理
- 如果仍有问题，检查 LAMMPS 输出是否正确

#### 2. 团簇仍然分散？
**现象：**
```
WARNING: Cluster still appears fragmented (max dist > 10 A)
```

**可能原因：**
- 团簇本身就是扩展结构（非球形）
- 原子分组错误（检查 z 坐标范围）
- 载体原子数设置错误

#### 3. Support 和 Cluster z 坐标重叠？
**现象：**
```
WARNING: Support and Cluster z-coordinates overlap!
```

**说明：**
- 这可能表示团簇部分嵌入载体
- 检查是否需要调整 `--support` 参数

---

## 🔬 后续分析建议

### 使用载体居中轨迹 (`surface_centered.xyz`)
```python
# 示例：计算团簇质心相对载体的位移
import numpy as np

# 读取轨迹...
# cluster_com = ...  # 团簇质心坐标

# 计算 XY 平面的扩散
xy_displacement = cluster_com[:, :2] - cluster_com[0, :2]
msd = (xy_displacement**2).sum(axis=1).mean()
```

### 使用团簇居中轨迹 (`cluster_centered.xyz`)
```python
# 示例：计算团簇 RMSD
from scipy.spatial.distance import cdist

# 读取团簇坐标...
# cluster_coords = ...  # [n_frames, 14, 3]

# 计算相对第一帧的 RMSD
rmsd = []
for frame in cluster_coords:
    rmsd.append(np.sqrt(((frame - cluster_coords[0])**2).sum() / 14))
```

---

## 📚 技术细节

### 最小镜像约定（Minimum Image Convention）
```
对于坐标差 Δ：
  if Δ > L/2:   Δ = Δ - L  (原子在右边界，拉回左边)
  if Δ < -L/2:  Δ = Δ + L  (原子在左边界，拉回右边)
```

### 迭代质心法收敛条件
- 当一次迭代中没有原子需要移动时，算法收敛
- 通常 2-4 次迭代即可收敛
- 最多迭代 10 次

---

## 🤝 贡献与反馈

如有问题或建议，请通过 GitHub Issues 反馈。

---

## 📄 许可证

MIT License

---

## 🔖 版本历史

### v1.0.0 (2025-11-06)
- ✅ 自动检测团簇原子数
- ✅ 自动 unwrap 轨迹
- ✅ 迭代质心法重组团簇
- ✅ 双模式居中输出
- ✅ 支持自定义输出前缀
- ✅ 完整的命令行参数支持

---

## 📧 联系方式

作者：[Your Name]  
邮箱：[Your Email]  
GitHub: [Your GitHub]
