import os
import glob
import numpy as np
import pyvista as pv
from scipy.io import savemat, loadmat
import time

# ====================== 工具函数 ======================

def preload_meshes(files):
    """预加载所有网格到内存"""
    print("🔄 预加载网格数据...")
    meshes = []
    for fname in files:
        meshes.append(pv.read(fname))
    print(f"✅ 共加载 {len(meshes)} 个网格文件")
    return meshes


def get_velocities_batch(positions, mesh):
    """批量速度插值优化版本"""
    if len(positions) == 0:
        return np.zeros((0, 3))
    pts = pv.PolyData(positions)
    sampled = pts.sample(mesh)
    if "Velocity" not in sampled.point_data:
        return np.zeros_like(positions)
    return sampled.point_data["Velocity"] * 10


def rk2_step_batch(positions, dt, mesh):
    """二阶 Runge-Kutta 积分"""
    k1 = get_velocities_batch(positions, mesh)
    predicted_positions = positions + dt * k1
    k2 = get_velocities_batch(predicted_positions, mesh)
    return positions + (dt / 2) * (k1 + k2)

def load_latest_positions(filename):
    """从MAT文件中加载最近一次保存的粒子坐标（支持 cell + struct 层级）"""
    import numpy as np
    from scipy.io import loadmat

    data = loadmat(filename, struct_as_record=False, squeeze_me=False)
    scatters = data["scatters"]

    # scatters: (N,1) cell array，每个 cell 是 (1,1) struct
    try:
        last_cell = scatters[-1, 0]      # cell 内容
        struct_obj = last_cell[0, 0]     # 取 (1,1) 的 struct
        pos_field = struct_obj.pos if hasattr(struct_obj, "pos") else struct_obj["pos"]
    except Exception as e:
        raise ValueError(f"❌ 无法解析 {filename}，建议打印 scatters[-1,0] 看实际结构") from e

    # 提取 pos 数据
    if isinstance(pos_field, np.ndarray) and pos_field.ndim == 2:
        last_snapshot = pos_field
    elif isinstance(pos_field, np.ndarray) and pos_field.size == 1:
        last_snapshot = pos_field.item()
    elif isinstance(pos_field, np.void):
        # 万一是 structure
        last_snapshot = pos_field[0, 0]
    else:
        raise ValueError(f"❌ 未识别 pos 数据结构类型: {type(pos_field)}")

    last_snapshot = np.array(last_snapshot)
    print(f"✅ 已从 {filename} 读取 {last_snapshot.shape[0]} 个粒子位置")
    return last_snapshot

# ====================== 主程序 ======================

def main_optimized(resume_file=None):
    folder = "VelocityField"
    file_pattern = os.path.join(folder, "result_*.vtu")
    files = sorted(glob.glob(file_pattern))
    meshes = preload_meshes(files)

    # 基本参数
    n_substeps = 5
    dt = 0.0025
    stagnant_tol = 1e-6
    n_cycles = 3

    # ============ 初始化粒子或加载断点 ============
    if resume_file is not None:
        current_positions = load_latest_positions(resume_file)
        import re
        m = re.search(r"cycle(\d+)_part(\d+)", resume_file)
        if not m:
            raise ValueError("❌ 无法从文件名解析出循环/分段编号，请确认文件名格式")
        start_cycle = int(m.group(1)) - 1  # e.g. cycle02 -> index 1
        start_part = int(m.group(2))       # e.g. part10 -> 下一次应从 part11
        print(f"🔁 从循环 cycle={start_cycle+1}, part={start_part} 继续计算…")
    else:
        # 初始粒子
        n_particles1, n_particles2 = 25000, 25000
        particles1 = np.column_stack([
            np.random.uniform(-5.0, -4.5, n_particles1),
            np.random.uniform(-5.2, -4.7, n_particles1),
            np.full(n_particles1, -4.8)
        ])
        particles2 = np.column_stack([
            np.random.uniform(0.5, 1.0, n_particles2),
            np.random.uniform(-2.1, -1.6, n_particles2),
            np.full(n_particles2, -10.5)
        ])
        particles = np.vstack([particles1])
        current_positions = particles.copy()
        start_cycle = 0
        start_part = 0

    # ====================== 主循环 ======================
    for cycle in range(start_cycle, n_cycles):
        print(f"\n===== 第 {cycle + 1} 次循环 =====")
        start_time = time.time()
        trajectories = []

        for step, mesh in enumerate(meshes):
            # 如果是恢复模式，跳过已完成的段
            if resume_file and cycle == start_cycle and step < start_part * 5:
                continue

            if step % 1 == 0:
                print(f"处理步骤 {step + 1}/{len(meshes)}")

            for sub in range(n_substeps):
                trajectories.append(current_positions.copy())
                current_positions = rk2_step_batch(current_positions, dt / n_substeps, mesh)

                # 停滞粒子判断与重置
                vel = get_velocities_batch(current_positions, mesh)
                speed = np.linalg.norm(vel, axis=1)
                stagnant_mask = speed < stagnant_tol
                n_reseed = stagnant_mask.sum()
                if n_reseed > 0:
                    reseed_indices = np.where(stagnant_mask)[0]
                    rand_vals = np.random.rand(len(reseed_indices))
                    mask1 = rand_vals < 1
                    if mask1.any():
                        idx1 = reseed_indices[mask1]
                        current_positions[idx1, 0] = np.random.uniform(-5.0, -4.5, len(idx1))
                        current_positions[idx1, 1] = np.random.uniform(-5.2, -4.7, len(idx1))
                        current_positions[idx1, 2] = -4.8
                    if (~mask1).any():
                        idx2 = reseed_indices[~mask1]
                        current_positions[idx2, 0] = np.random.uniform(0.5, 1.0, len(idx2))
                        current_positions[idx2, 1] = np.random.uniform(-2.1, -1.6, len(idx2))
                        current_positions[idx2, 2] = -10.5
                    print(f"  >> 重置 {n_reseed} 个粒子")

            # === 每5个文件保存一次 ===
            if (step + 1) % 5 == 0 or (step + 1) == len(meshes):
                print(f"💾 保存中：已处理 {step + 1} 个文件...")
                trajectories = np.array(trajectories)
                scatters = [{'pos': snapshot} for snapshot in trajectories]
                scatters = np.array(scatters, dtype=object).reshape(-1, 1)
                save_filename = (
                    f"particle_trajectories_cycle{cycle + 1:02d}_part{(step // 5) + 1:02d}.mat"
                )
                savemat(save_filename, {"scatters": scatters})
                print(f"✅ 已保存轨迹到 {save_filename}")
                trajectories = []

        cycle_time = time.time() - start_time
        print(f"✅ 循环 {cycle + 1} 完成，耗时 {cycle_time:.2f} 秒")

    print("🎉 全部计算完成！")


# ====================== 程序入口 ======================

if __name__ == "__main__":
    import sys
    resume_file = None
    if len(sys.argv) > 2 and sys.argv[1].lower() == "resume":
        resume_file = sys.argv[2]
    trajectories = main_optimized(resume_file)
