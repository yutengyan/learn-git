#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sn7Pt4O3/Al2O3 Dual-Mode Trajectory Processing Script

Features:
1. Auto-identify Support (240 atoms) and Cluster (14 atoms) based on z-coordinate
2. Verify unwrap correctness (detect box-crossing jumps)
3. Generate two centered trajectories:
   - surface_centered.xyz (support-centered, for surface migration)
   - cluster_centered.xyz (cluster-centered, for internal rearrangement)
"""

import numpy as np
import sys
from pathlib import Path

class TrajectoryAnalyzer:
    def __init__(self, xyz_file, support_n=240, cluster_n=14):
        """
        Parameters:
            xyz_file: Input unwrapped XYZ trajectory file
            support_n: Number of Support atoms (default 240)
            cluster_n: Number of Cluster atoms (default 14)
        """
        self.xyz_file = Path(xyz_file)
        self.support_n = support_n
        self.cluster_n = cluster_n
        self.total_atoms = support_n + cluster_n
        
        # Store data
        self.frames = []
        self.box_vectors = []
        self.atom_types = []
        self.support_ids = []
        self.cluster_ids = []
        
    def read_xyz(self):
        """Read extended XYZ format trajectory"""
        print(f"Reading trajectory: {self.xyz_file}")
        
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
        
        i = 0
        frame_count = 0
        
        while i < len(lines):
            # 读取原子数
            try:
                n_atoms = int(lines[i].strip())
            except:
                break
            
            if n_atoms != self.total_atoms:
                print(f"Warning: Frame {frame_count} has {n_atoms} atoms, expected {self.total_atoms}")
                break
            
            # Read Lattice info
            header = lines[i + 1].strip()
            box = self._parse_lattice(header)
            
            # Read atom coordinates
            coords = []
            types = []
            for j in range(n_atoms):
                parts = lines[i + 2 + j].split()
                types.append(parts[0])  # Element symbol
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            
            self.frames.append(np.array(coords))
            self.box_vectors.append(box)
            
            # Record atom types only for first frame
            if frame_count == 0:
                self.atom_types = types
            
            i += 2 + n_atoms
            frame_count += 1
        
        self.frames = np.array(self.frames)
        self.box_vectors = np.array(self.box_vectors)
        
        print(f"Successfully read {frame_count} frames, {self.total_atoms} atoms per frame")
        print(f"   Box size: {self.box_vectors[0]}")
        
    def _parse_lattice(self, header):
        """Parse Lattice information"""
        # Format: Lattice="a 0 0 0 b 0 0 0 c" ...
        parts = header.split('"')
        if len(parts) >= 2:
            values = [float(x) for x in parts[1].split()]
            return np.array([values[0], values[4], values[8]])  # a, b, c
        else:
            return np.array([30.0, 30.0, 30.0])  # Default value
    
    def identify_groups(self):
        """Identify Support and Cluster based on initial frame z-coordinates"""
        print("\nIdentifying atom groups based on z-coordinates...")
        
        # Sort by z-coordinate of first frame
        z_coords = self.frames[0, :, 2]
        sorted_ids = np.argsort(z_coords)
        
        # First 240 are Support, last 14 are Cluster
        self.support_ids = sorted_ids[:self.support_n]
        self.cluster_ids = sorted_ids[self.support_n:self.support_n + self.cluster_n]
        
        # Count element composition
        support_elements = [self.atom_types[i] for i in self.support_ids]
        cluster_elements = [self.atom_types[i] for i in self.cluster_ids]
        
        print(f"   Support ({self.support_n} atoms):")
        for elem in set(support_elements):
            count = support_elements.count(elem)
            print(f"      {elem}: {count}")
        
        print(f"   Cluster ({self.cluster_n} atoms):")
        for elem in set(cluster_elements):
            count = cluster_elements.count(elem)
            print(f"      {elem}: {count}")
        
        # Verify z range
        support_z_range = [z_coords[self.support_ids].min(), z_coords[self.support_ids].max()]
        cluster_z_range = [z_coords[self.cluster_ids].min(), z_coords[self.cluster_ids].max()]
        
        print(f"\n   Support z range: [{support_z_range[0]:.2f}, {support_z_range[1]:.2f}] A")
        print(f"   Cluster z range: [{cluster_z_range[0]:.2f}, {cluster_z_range[1]:.2f}] A")
        
        if support_z_range[1] > cluster_z_range[0]:
            print(f"   WARNING: Support and Cluster z-coordinates overlap!")
    
    def verify_unwrap(self):
        """Verify trajectory is properly unwrapped (detect abnormal cross-frame displacement)"""
        print("\nVerifying unwrap correctness...")
        
        issues = []
        
        for i in range(1, len(self.frames)):
            # 计算帧间位移
            displacement = self.frames[i] - self.frames[i-1]
            max_disp = np.abs(displacement).max(axis=1)
            
            # 检测异常位移（阈值设为盒子长度的 1/3）
            threshold = self.box_vectors[0].min() / 3.0
            bad_atoms = np.where(max_disp > threshold)[0]
            
            if len(bad_atoms) > 0:
                issues.append({
                    'frame': i,
                    'atoms': bad_atoms,
                    'max_disp': max_disp[bad_atoms].max()
                })
        
        if issues:
            print(f"   ❌ 检测到 {len(issues)} 帧存在疑似跨盒跳跃！")
            print(f"   前 5 个问题帧:")
            for issue in issues[:5]:
                print(f"      帧 {issue['frame']}: {len(issue['atoms'])} 原子，最大位移 {issue['max_disp']:.2f} Å")
            print(f"\n   ⚠️  轨迹可能未正确 unwrap，建议检查 LAMMPS 输出！")
            return False
        else:
            print(f"   ✅ 未检测到异常跨盒跳跃，轨迹 unwrap 正确")
            
            # 额外验证：检查团簇整体漂移
            cluster_com = self.frames[:, self.cluster_ids, :].mean(axis=1)
            total_drift = np.linalg.norm(cluster_com[-1] - cluster_com[0])
            print(f"   ✅ 团簇质心总漂移: {total_drift:.2f} Å")
            
            return True
    
    def center_trajectory(self, ref_group, output_file):
        """
        将轨迹以指定原子组的质心居中
        
        参数:
            ref_group: 'support' 或 'cluster'
            output_file: 输出文件名
        """
        if ref_group == 'support':
            ref_ids = self.support_ids
            desc = "基底 (Support)"
        elif ref_group == 'cluster':
            ref_ids = self.cluster_ids
            desc = "团簇 (Cluster)"
        else:
            raise ValueError("ref_group 必须是 'support' 或 'cluster'")
        
        print(f"\n📐 生成 {desc} 居中轨迹...")
        
        centered_frames = []
        
        for frame_idx, frame in enumerate(self.frames):
            # 计算参考组质心
            com = frame[ref_ids].mean(axis=0)
            
            # 将盒子中心设为原点
            box_center = self.box_vectors[frame_idx] / 2.0
            
            # 平移整个体系
            shift = box_center - com
            centered_frame = frame + shift
            
            centered_frames.append(centered_frame)
        
        centered_frames = np.array(centered_frames)
        
        # 写入 XYZ 文件
        self._write_xyz(centered_frames, output_file, f"{desc} 居中")
        
        print(f"   ✅ 已保存至: {output_file}")
        
        return centered_frames
    
    def _write_xyz(self, frames, filename, description=""):
        """写入 XYZ 格式轨迹"""
        with open(filename, 'w') as f:
            for frame_idx, frame in enumerate(frames):
                # 写入原子数
                f.write(f"{self.total_atoms}\n")
                
                # 写入 Lattice 信息
                box = self.box_vectors[min(frame_idx, len(self.box_vectors)-1)]
                lattice_str = f'Lattice="{box[0]:.8f} 0.00000000 0.00000000 0.00000000 {box[1]:.8f} 0.00000000 0.00000000 0.00000000 {box[2]:.8f}"'
                f.write(f'{lattice_str} pbc="1 1 1" Properties=species:S:1:pos:R:3 # {description} Frame {frame_idx}\n')
                
                # 写入原子坐标
                for atom_idx, (atom_type, coords) in enumerate(zip(self.atom_types, frame)):
                    f.write(f"{atom_type} {coords[0]:.8f} {coords[1]:.8f} {coords[2]:.8f}\n")
    
    def generate_index_file(self, output_file="index_zsplit.ndx"):
        """生成 GROMACS 风格的索引文件"""
        print(f"\n📋 生成索引文件: {output_file}")
        
        with open(output_file, 'w') as f:
            # Support 组（原子 ID 从 1 开始）
            f.write("[ Support ]\n")
            for i, atom_id in enumerate(self.support_ids + 1):
                f.write(f"{atom_id:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n\n")
            
            # Cluster 组
            f.write("[ PtSnCluster ]\n")
            for i, atom_id in enumerate(self.cluster_ids + 1):
                f.write(f"{atom_id:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n\n")
            
            # System 组
            f.write("[ System ]\n")
            for i in range(self.total_atoms):
                f.write(f"{i+1:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n")
        
        print(f"   ✅ 索引文件已保存")
    
    def run(self):
        """执行完整分析流程"""
        print("=" * 70)
        print("  Sn₇Pt₄O₃/Al₂O₃ 轨迹处理 - 双模式居中分析")
        print("=" * 70)
        
        # 1. 读取轨迹
        self.read_xyz()
        
        # 2. 识别原子组
        self.identify_groups()
        
        # 3. 验证 unwrap
        is_valid = self.verify_unwrap()
        
        if not is_valid:
            response = input("\n❓ 检测到可能的问题，是否继续处理？(y/n): ")
            if response.lower() != 'y':
                print("❌ 用户取消操作")
                return
        
        # 4. 生成索引文件
        self.generate_index_file()
        
        # 5. 生成基底居中轨迹
        self.center_trajectory('support', 'surface_centered.xyz')
        
        # 6. 生成团簇居中轨迹
        self.center_trajectory('cluster', 'cluster_centered.xyz')
        
        print("\n" + "=" * 70)
        print("✅ 处理完成！生成文件:")
        print("   - index_zsplit.ndx       (原子组索引)")
        print("   - surface_centered.xyz   (基底居中，用于表面迁移分析)")
        print("   - cluster_centered.xyz   (团簇居中，用于内部重排分析)")
        print("=" * 70)


def main():
    if len(sys.argv) < 2:
        print("用法: python analyze_ptsn_trajectory.py <trajectory.xyz> [support_n] [cluster_n]")
        print("\n示例:")
        print("  python analyze_ptsn_trajectory.py sampling-simply.xyz")
        print("  python analyze_ptsn_trajectory.py traj.xyz 240 14")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    support_n = int(sys.argv[2]) if len(sys.argv) > 2 else 240
    cluster_n = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    
    analyzer = TrajectoryAnalyzer(xyz_file, support_n, cluster_n)
    analyzer.run()


if __name__ == "__main__":
    main()