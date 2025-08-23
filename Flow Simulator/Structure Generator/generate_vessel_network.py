#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
血管网络生成器 V-System 2.0

这个脚本整合了V-System项目的功能，用于生成合成血管网络。
它使用Lindenmayer系统（L-系统）生成血管树结构，并将其转换为3D体素表示。
所有功能已经整合到单个文件中，便于使用和维护。
"""

import argparse
import os
import random
import time

import numpy as np
import tifffile
from scipy.spatial.transform import Rotation as R

# ----------------- 工具函数 -----------------


def normalise(x, scale=1.0):
    """向量归一化"""
    return scale * (x / np.linalg.norm(x)) if np.linalg.norm(x) > 0 else x


def rotate(pitch_angle=0.0, roll_angle=0.0, yaw_angle=0.0, vector=None):
    """旋转向量 - 使用scipy的旋转矩阵"""
    # 创建旋转对象（按照y-x-z顺序应用旋转）
    rotation = R.from_euler("yxz", [pitch_angle, roll_angle, yaw_angle], degrees=False)
    # 应用旋转到向量
    return rotation.apply(vector)


def check_boundary(image):
    """设置边界像素为零"""
    # 更简洁的边界设置方法
    image[0, :, :] = image[-1, :, :] = 0
    image[:, 0, :] = image[:, -1, :] = 0
    image[:, :, 0] = image[:, :, -1] = 0
    return image


def generate_points(newNodes=None, points=None, usenan=True):
    """生成新点"""
    if newNodes is None or newNodes.size == 0:
        newNodes = points
    else:
        newNodes = np.hstack((newNodes, points))

    if usenan:
        nanVec = np.full((4, 1), np.nan)  # 使用np.full替代空数组+赋值
        newNodes = np.hstack((newNodes, nanVec))

    return newNodes


def bezier_interpolation(coords=None):
    """贝塞尔曲线插值"""
    if coords is None:
        return np.array([])

    X, Y, Z, diam = [], [], [], []

    try:
        for c in coords:
            X.append(float(c[0]))
            Y.append(float(c[1]))
            Z.append(float(c[2]))
            diam.append(float(c[5]))
    except:
        pass

    nodes = np.vstack((X, Y, Z, diam))

    vessel = np.array([])
    newNodes = np.array([])

    for i in range(0, nodes.shape[1] - 1):
        if np.isnan(nodes[0, i]):
            _, idx = np.unique(vessel, axis=1, return_index=True)
            last = vessel[:, -1]
            vessel = vessel[:, np.sort(idx)]

            if vessel.shape[1] > 2:
                daughters = np.arange(0, vessel.shape[1] - 5, 5)
                for j in daughters:
                    theChosenOne = vessel[:, j : (j + 6)]
                    newNodes = generate_points(newNodes, theChosenOne)

            vessel = np.array([])
        else:
            if vessel.size == 0:
                vessel = nodes[:, i].reshape(-1, 1)
            else:
                vessel = np.hstack((vessel, nodes[:, i].reshape(-1, 1)))

    return newNodes


# ----------------- L系统属性 -----------------

# 默认L系统属性
properties = {"k": 3, "epsilon": 10, "randmarg": 3, "sigma": 5, "stochparams": True}


def set_properties(props=None):
    """设置L系统属性"""
    global properties
    if props is not None:
        properties.update(props)
    return properties


def get_length(d0):
    """根据父分支直径计算分支长度"""
    c0 = d0 * properties["epsilon"]
    return np.random.uniform(c0 - properties["randmarg"], c0 + properties["randmarg"])


def cal_param(text, params):
    """计算括号内的值"""
    txt = text[:]
    for i in params:
        txt = txt.replace(i, str(params[i]))
    return str(params["co"] / eval(txt))


def cal_bifurcation(d0):
    """计算给定直径的分叉直径和角度"""
    resp = {}
    k = properties["k"]
    sigma = properties["sigma"]
    stochparams = properties["stochparams"]

    dOpti = d0 / 2 ** (1.0 / k)
    if stochparams:
        d1 = abs(np.random.normal(dOpti, dOpti / sigma))
    else:
        d1 = dOpti

    if d1 >= d0:
        d1 = dOpti

    d2 = (d0**k - d1**k) ** (1.0 / k)
    alpha = d2 / d1

    # 模拟人体分叉角度的方程
    xtmp = (1 + alpha**3) ** (4.0 / 3) + 1 - alpha**4
    xtmpb = 2 * ((1 + alpha**3) ** (2.0 / 3))
    a1 = np.arccos(np.clip(xtmp / xtmpb, -1.0, 1.0))  # 使用np.arccos和np.clip避免数值问题

    xtmp = (1 + alpha**3) ** (4.0 / 3) + (alpha**4) - 1
    xtmpb = 2 * alpha**2 * ((1 + alpha**3) ** (2.0 / 3))
    a2 = np.arccos(np.clip(xtmp / xtmpb, -1.0, 1.0))  # 使用np.arccos和np.clip避免数值问题

    resp["d1"] = d1
    resp["d2"] = d2
    resp["d0"] = d0
    resp["th1"] = np.degrees(a1)  # 使用np.degrees替代手动计算
    resp["th2"] = np.degrees(a2)  # 使用np.degrees替代手动计算
    resp["co"] = get_length(d0)

    return resp


# ----------------- L系统函数 -----------------


def F(n, d0):
    """F分形模式生成函数"""
    if n > 0:
        params = cal_bifurcation(d0)
        theta1 = params["th1"]
        theta2 = params["th2"]
        return (
            S(n - 1, d0)
            + "["
            + "+("
            + str(theta1)
            + ")"
            + "/("
            + str(np.random.uniform(22.5, 27.5) * random.randint(-1, 1))
            + ")"
            + F(n - 1, params["d1"])
            + "]"
            + "["
            + "-("
            + str(theta2)
            + ")"
            + "/("
            + str(np.random.uniform(22.5, 27.5) * random.randint(-1, 1))
            + ")"
            + F(n - 1, params["d2"])
            + "]"
        )
    else:
        return "F"


def S(n, d0, val=5, margin=0.5):
    """S函数：生成树枝形状"""
    r = random.random()
    if r >= 0.0 and r < margin:
        return "{" + S1(n, d0, val) + "}"
    if r >= margin and r < 1.0:
        return "{" + S2(n, d0, val) + "}"


def S1(n, d0, val=5):
    """S1递归函数：生成树结构的字符串表示"""
    if n > 0:
        params = cal_bifurcation(d0)
        randInt = random.randint(-1, 1)
        theta1 = params["th1"] * randInt
        theta2 = params["th2"] * abs(randInt)
        if random.random() < 0.1:
            growth = np.sort(np.random.uniform(1.0, 1.5, 3))
            decay = -np.sort(-np.random.uniform(1.0, 1.5, 2))
        else:
            growth = np.ones(3)
            decay = np.ones(2)
        descrip = (
            D(n - 1, growth[0] * params["d0"], val)
            + "+("
            + str(theta1)
            + ")"
            + D(n - 1, growth[1] * params["d0"], val)
            + "-("
            + str(theta2)
            + ")"
            + D(n - 1, growth[2] * params["d0"], val)
            + "-("
            + str(theta1)
            + ")"
            + D(n - 1, decay[0] * params["d0"], val)
            + "+("
            + str(theta2)
            + ")"
            + D(n - 1, decay[1] * params["d0"], val)
        )
        return descrip
    else:
        return "S"


def S2(n, d0, val=5):
    """S2递归函数：生成树结构的字符串表示"""
    if n > 0:
        params = cal_bifurcation(d0)
        randInt = random.randint(-1, 1)
        theta1 = params["th1"] * randInt
        theta2 = params["th2"] * abs(randInt)
        if random.random() < 0.1:
            growth = np.sort(np.random.uniform(1.0, 1.5, 3))
            decay = -np.sort(-np.random.uniform(1.0, 1.5, 2))
        else:
            growth = np.ones(3)
            decay = np.ones(2)
        descrip = (
            D(n - 1, growth[0] * params["d0"], val)
            + "-("
            + str(theta1)
            + ")"
            + D(n - 1, growth[1] * params["d0"], val)
            + "+("
            + str(theta2)
            + ")"
            + D(n - 1, growth[2] * params["d0"], val)
            + "+("
            + str(theta1)
            + ")"
            + D(n - 1, decay[0] * params["d0"], val)
            + "-("
            + str(theta2)
            + ")"
            + D(n - 1, decay[1] * params["d0"], val)
        )
        return descrip
    else:
        return "S"


def D(n, d0, val=5):
    """D函数：生成L系统符号字符串"""
    if n > 0:
        params = cal_bifurcation(d0)
        p1 = cal_param("co/" + str(int(val)), params)
        return "f(" + p1 + "," + str(params["d0"]) + ")"
    else:
        return "D"


# ----------------- 海龟解析 -----------------


def eval_brackets(index, turtle):
    """提取海龟程序中括号内的值并计算"""
    a = ""
    b = ""
    double = 0
    neg1 = 0
    neg2 = 0

    for i in range(index + 2, len(turtle)):
        if turtle[i] == ")":
            break
        elif turtle[i] == ",":
            double = 1
        elif turtle[i] == "-":
            if double == 0:
                neg1 = 1
            else:
                neg2 = 1
        elif not (turtle[i] == "(" or turtle[i] == "+"):
            if double == 0:
                a = a + str(turtle[i])
            else:
                b = b + str(turtle[i])

    if double == 1:
        aa = eval(a)
        bb = eval(b)
        if neg1 == 1:
            aa *= -1
        if neg2 == 1:
            bb *= -1
        return aa, bb
    else:
        aa = eval(a)
        if neg1 == 1:
            aa *= -1
        return aa, 0.0


def branching_turtle_to_coords(turtle_program, d0, theta=20.0, phi=20.0):
    """将海龟命令转换为坐标"""
    DEGREES_TO_RADIANS = np.pi / 180
    saved_states = list()
    stateSize = 10
    dx = 0
    dy = 0
    dz = 0
    lseg = 1.0

    startidx = 3  # random.randint(1,3)
    if startidx == 1:
        state = (1.0, 0.1, 0.1, 0, 0, d0, lseg, dx, dy, dz)
    elif startidx == 2:
        state = (0.1, 1.0, 0, 0, 0, d0, lseg, dx, dy, dz)
    else:
        state = (0.1, 0.1, 1.0, 0, 0, d0, lseg, dx, dy, dz)

    yield state

    index = 0

    for command in turtle_program:
        x, y, z, alpha, beta, diam, lseg, dx, dy, dz = state

        if command.lower() in "abcdefghijs":
            if command.islower():
                lseg, tdiam = eval_brackets(index, turtle_program)
                dx, dy, dz = rotate(
                    pitch_angle=beta * DEGREES_TO_RADIANS,
                    roll_angle=alpha * DEGREES_TO_RADIANS,
                    vector=normalise(np.array([x, y, z]), lseg),
                )

                if tdiam > 0.0:
                    diam = tdiam

                x += dx
                y += dy
                z += dz

            state = (x, y, z, alpha, beta, diam, lseg, dx, dy, dz)
            yield state

        elif command == "+":
            phi, _ = eval_brackets(index, turtle_program)
            state = (x, y, z, alpha + phi, beta, diam, lseg, dx, dy, dz)

        elif command == "-":
            phi, _ = eval_brackets(index, turtle_program)
            state = (x, y, z, alpha - phi, beta, diam, lseg, dx, dy, dz)

        elif command == "/":
            theta, _ = eval_brackets(index, turtle_program)
            state = (x, y, z, alpha, beta + theta, diam, lseg, dx, dy, dz)

        elif command == "[":
            saved_states.append(state)

        elif command == "]":
            state = saved_states.pop()
            # 生成nan值序列
            yield tuple(np.full(stateSize, np.nan))
            x, y, z, alpha, beta, diam, lseg, dx, dy, dz = state
            yield state

        index += 1


# ----------------- 体素计算 -----------------


def diam_voxels(x, y, z, idx, r, img_stack, tVol, res=(1, 1, 1)):
    """创建以(x, y, z)为中心、半径为r的球体中体素的二值图像"""
    t0, t1, t2 = tVol[0], tVol[1], tVol[2]
    r0, r1, r2 = res[0], res[1], res[2]

    # 坐标范围检查提前
    x_int, y_int, z_int = int(x - 1), int(y - 1), int(z - 1)

    if idx == 0:
        for i in range(-r, r + 1):
            y_i = int(y + i - 1)
            if 0 <= y_i < t1:
                for j in range(-r, r + 1):
                    z_j = int(z + j - 1)
                    if 0 <= z_j < t2:
                        if np.sqrt((i * r1) ** 2 + (j * r2) ** 2) <= r:
                            img_stack[x_int, y_i, z_j] = 1
    elif idx == 1:
        for i in range(-r, r + 1):
            x_i = int(x + i - 1)
            if 0 <= x_i < t0:
                for j in range(-r, r + 1):
                    z_j = int(z + j - 1)
                    if 0 <= z_j < t2:
                        if np.sqrt((i * r0) ** 2 + (j * r2) ** 2) <= r:
                            img_stack[x_i, y_int, z_j] = 1
    else:
        for i in range(-r, r + 1):
            x_i = int(x + i - 1)
            if 0 <= x_i < t0:
                for j in range(-r, r + 1):
                    y_j = int(y + j - 1)
                    if 0 <= y_j < t1:
                        if np.sqrt((i * r0) ** 2 + (j * r1) ** 2) <= r:
                            img_stack[x_i, y_j, z_int] = 1

    return img_stack


def transversal(x, y, z, dx, dy, dz, diam, img_stack=None, tVol=(600, 600, 700), resolution=(1, 1, 1)):
    """计算两点之间线段横切的体素 - Bresenham的3D线算法"""
    # 计算坐标变化
    dx = dx - x
    dy = dy - y
    dz = dz - z

    # 计算变化方向和绝对值
    sx = np.sign(dx)
    sy = np.sign(dy)
    sz = np.sign(dz)

    ax = abs(dx)
    ay = abs(dy)
    az = abs(dz)

    bx = 2 * ax
    by = 2 * ay
    bz = 2 * az

    exy = ay - ax
    exz = az - ax
    ezy = ay - az

    n = ax + ay + az

    # 计算半径
    r = int(np.floor(0.5 * diam))

    # 寻找管/血管内的体素
    while n > 0:
        if exy < 0:
            if exz < 0:
                x += sx
                exy += by
                exz += bz
                if 0 <= x < tVol[0]:
                    img_stack = diam_voxels(x, y, z, 0, r, img_stack, tVol, res=resolution)
            else:
                z += sz
                exz -= bx
                ezy += by
                if 0 <= z < tVol[2]:
                    img_stack = diam_voxels(x, y, z, 2, r, img_stack, tVol, res=resolution)
        else:
            if ezy < 0:
                z += sz
                exz -= bx
                ezy += by
                if 0 <= z < tVol[2]:
                    img_stack = diam_voxels(x, y, z, 2, r, img_stack, tVol, res=resolution)
            else:
                y += sy
                exy -= bx
                ezy -= bz
                if 0 <= y < tVol[1]:
                    img_stack = diam_voxels(x, y, z, 1, r, img_stack, tVol, res=resolution)

        n -= 1

    return img_stack


def discretisation(p1, dmin, dmax, tVol):
    """将3D点从实际坐标转换为离散体素坐标"""
    # 向量化操作
    coords = np.zeros(3)
    for i in range(3):
        coords[i] = np.floor(tVol[i] * (p1[i] - dmin[i]) / (dmax[i] - dmin[i]))
    return np.append(coords, p1[3]).reshape(-1, 1)


def process_network(data, tVol):
    """处理3D医学图像数据，生成点数据"""
    # 找到数据中非NaN的最小/最大值，留出10%的边界
    dmin = np.nanmin(data, axis=1) * 1.1
    dmax = np.nanmax(data, axis=1) * 1.1

    # 初始化新数组
    newarray = np.array([])
    img_stack = np.zeros(tVol, dtype=int)

    # 处理数据
    for i in range(data.shape[1]):
        if not np.isnan(data[0, i]):
            tempvec = discretisation(data[:, i], dmin, dmax, tVol)
        else:
            tempvec = data[:, i].reshape(-1, 1)
        newarray = generate_points(newarray, tempvec, usenan=False)

    # 处理网络
    j = 0
    for i in range(newarray.shape[1] - 1):
        j += 1
        if not np.isnan(newarray[0, i]) and j < newarray.shape[1] and not np.isnan(newarray[0, j]):
            diam = 0.5 * (newarray[3, i] + newarray[3, j])
            transversal(
                newarray[0, i],
                newarray[1, i],
                newarray[2, i],
                newarray[0, j],
                newarray[1, j],
                newarray[2, j],
                diam,
                img_stack,
                tVol,
            )

    return check_boundary(img_stack)


# ----------------- 血管网络生成器类 -----------------


class VesselNetworkGenerator:
    """血管网络生成器类"""

    def __init__(self, output_path="./output/", tissue_volume=(512, 512, 140)):
        """初始化血管网络生成器"""
        self.output_path = output_path
        self.tissue_volume = tissue_volume

        # 确保输出目录存在
        os.makedirs(output_path, exist_ok=True)
        print(f"确保输出目录存在: {output_path}")

    def set_properties(self, props=None):
        """设置L-系统的属性"""
        if props is None:
            props = {
                "k": 3,
                "epsilon": random.uniform(4, 10),  # 长度与直径的比例
                "randmarg": random.uniform(2, 4),  # 长度与直径之间的随机边界
                "sigma": 5,  # 高斯分布的类型偏差
                "stochparams": True,  # 生成的参数是否也是随机的
            }

        return set_properties(props)

    def generate_network(self, d0=20.0, iterations=7, props=None):
        """生成血管网络"""
        # 设置V-System属性
        props = self.set_properties(props)

        print(f"生成血管网络，基础直径: {d0:.2f}，迭代次数: {iterations}")

        # 运行V-System语法进行n次迭代
        turtle_program = F(iterations, d0)

        # 将语法转换为坐标
        coords = list(branching_turtle_to_coords(turtle_program, d0))

        # 分析/排序坐标数据
        update = bezier_interpolation(coords)

        # 运行3D体素遍历以生成V-System网络的二进制掩码
        image = process_network(update, tVol=self.tissue_volume)

        return image

    def save_network(self, image, filename):
        """保存生成的血管网络"""
        # 转换为8位格式
        image_8bit = (255 * image).astype("uint8")

        # 保存图像体积
        full_path = os.path.join(self.output_path, filename)
        tifffile.imwrite(
            full_path,
            np.transpose(image_8bit, (2, 0, 1)),
            bigtiff=False,
        )
        print(f"保存血管网络到: {full_path}")

        return full_path

    def generate_multiple_networks(
        self,
        count=5,
        d0_mean=20.0,
        d0_std=5.0,
        min_iterations=6,
        max_iterations=14,
    ):
        """生成多个血管网络"""
        file_paths = []

        for i in range(count):
            # 随机分配基础直径
            d0 = np.random.normal(d0_mean, d0_std)

            # 随机分配V-System递归循环的次数
            iterations = random.randint(min_iterations, max_iterations)

            # 设置V-System属性
            props = self.set_properties()

            # 生成网络
            image = self.generate_network(d0, iterations, props)

            # 保存网络
            filename = f"VesselNetwork_i{iterations}_{i}.tiff"
            file_path = self.save_network(image, filename)
            file_paths.append(file_path)

        return file_paths


def main():
    """主函数"""
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="生成合成血管网络")
    parser.add_argument("--output", type=str, default="./output/", help="输出目录")
    parser.add_argument("--count", type=int, default=1, help="要生成的网络数量")
    parser.add_argument("--iterations", type=int, default=7, help="L-系统迭代次数")
    parser.add_argument("--diameter", type=float, default=20.0, help="基础血管直径")
    parser.add_argument("--volume", type=str, default="512,512,140", help="组织体积尺寸 (x,y,z)")
    args = parser.parse_args()

    # 开始计时
    start_time_total = time.time()

    # 解析组织体积
    tissue_volume = tuple(map(int, args.volume.split(",")))

    # 创建生成器
    generator = VesselNetworkGenerator(output_path=args.output, tissue_volume=tissue_volume)

    if args.count == 1:
        # 生成单个网络
        print("开始生成血管网络...")
        start_time = time.time()

        image = generator.generate_network(d0=args.diameter, iterations=args.iterations)
        generator.save_network(image, f"VesselNetwork_i{args.iterations}.tiff")

        elapsed_time = time.time() - start_time
        print(f"血管网络生成完成！耗时: {elapsed_time:.2f} 秒")
    else:
        # 生成多个网络
        print(f"开始生成 {args.count} 个血管网络...")

        for i in range(args.count):
            start_time = time.time()
            print(f"生成第 {i+1}/{args.count} 个网络...")

            # 为每个网络生成随机参数
            d0 = np.random.normal(args.diameter, 5.0)
            iterations = random.randint(6, 14)
            props = generator.set_properties()

            # 生成并保存网络
            image = generator.generate_network(d0=d0, iterations=iterations, props=props)
            generator.save_network(image, f"VesselNetwork_i{iterations}_{i}.tiff")

            elapsed_time = time.time() - start_time
            print(f"第 {i+1} 个网络生成完成！耗时: {elapsed_time:.2f} 秒")

    # 计算总耗时
    total_elapsed_time = time.time() - start_time_total
    print(f"所有血管网络生成完成！总耗时: {total_elapsed_time:.2f} 秒")


if __name__ == "__main__":
    main()
