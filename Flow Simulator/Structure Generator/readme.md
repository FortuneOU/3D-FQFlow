# V-System 2.0 血管网络生成器

## 项目概述

V-System 2.0 是一个用于生成合成血管网络的工具。它使用 Lindenmayer 系统（L-系统）生成血管树结构，并将其转换为 3D 体素表示。这个项目已经被重构，所有功能现在都整合到一个单独的 Python 文件中，便于使用和维护。

## 特点

-   使用 L-系统生成逼真的血管树结构
-   支持 3D 体素渲染
-   可定制的血管网络参数
-   支持生成单个或多个血管网络
-   提供最大强度投影预览功能

## 安装

只需确保您安装了所需的依赖项：

```bash
pip install numpy scipy tifffile
```

## 使用方法

### 基本用法

```bash
python generate_vessel_network.py
```

这将使用默认参数生成一个血管网络，并将其保存在`./output/`目录中。

### 命令行参数

-   `--output`: 输出目录路径（默认：`./output/`）
-   `--count`: 要生成的网络数量（默认：1）
-   `--iterations`: L-系统迭代次数（默认：7）
-   `--diameter`: 基础血管直径（默认：20.0）
-   `--volume`: 组织体积尺寸，格式为"x,y,z"（默认："512,512,140"）

### 示例

生成 3 个不同的血管网络，并显示预览：

```bash
python generate_vessel_network.py --count 3
```

指定输出目录和组织体积大小：

```bash
python generate_vessel_network.py --output ./MyNetworks/ --volume 600,600,200
```

## 编程接口

您也可以在自己的 Python 脚本中使用这个生成器：

```python
from generate_vessel_network import VesselNetworkGenerator

# 创建生成器
generator = VesselNetworkGenerator(output_path="./output/", tissue_volume=(512, 512, 140))

# 生成单个网络
image = generator.generate_network(d0=20.0, iterations=7)

# 保存网络
generator.save_network(image, "my_vessel_network.tiff")

# 或者生成多个网络
file_paths = generator.generate_multiple_networks(count=5)
```

## 工作原理

1. **L-系统生成**: 使用分形规则生成血管树结构
2. **海龟解析**: 将 L-系统命令转换为 3D 坐标
3. **插值**: 使用贝塞尔曲线对坐标进行插值
4. **体素化**: 将 3D 坐标转换为体素表示
5. **渲染**: 生成最终的 3D 血管网络图像

## 注意事项

-   生成的图像默认以负值保存，这对某些医学图像处理软件更为友好
-   较高的迭代次数会生成更复杂的血管网络，但也会增加计算时间
-   所有功能已整合到单个文件中，无需额外的模块

## 许可证

此项目遵循 GPL 许可证。

## 参考链接

[V-System: Vascular Lindenmayer Systems](https://github.com/psweens/V-System)
