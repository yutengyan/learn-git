#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sn7Pt4O3/Al2O3 Dual-Mode Trajectory Processing Script

Features:
1. Auto-identify Support (240 atoms) and Cluster (14 atoms) based on z-coordinate
2. Auto-unwrap wrapped trajectories using minimum image convention
3. Verify unwrap correctness (detect box-crossing jumps)
4. Generate two centered trajectories:
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
    
    def make_cluster_whole(self):
        """
        Reassemble cluster atoms that are split across periodic boundaries
        Uses the first cluster atom as reference and brings others close via minimum image
        """
        print("\nReassembling cluster (making whole)...")
        
        for frame_idx in range(len(self.frames)):
            frame = self.frames[frame_idx]
            box = self.box_vectors[frame_idx]
            
            # Use first cluster atom as reference
            ref_pos = frame[self.cluster_ids[0]].copy()
            
            # Bring all other cluster atoms close to reference using minimum image
            for cluster_idx in self.cluster_ids[1:]:
                atom_pos = frame[cluster_idx]
                delta = atom_pos - ref_pos
                
                # Apply minimum image convention
                for dim in range(3):
                    if delta[dim] > box[dim] / 2.0:
                        frame[cluster_idx, dim] -= box[dim]
                    elif delta[dim] < -box[dim] / 2.0:
                        frame[cluster_idx, dim] += box[dim]
        
        # Verify cluster is now compact
        cluster_coords = self.frames[0, self.cluster_ids, :]
        from scipy.spatial.distance import pdist
        max_dist = pdist(cluster_coords).max()
        
        print(f"   First frame cluster max internal distance: {max_dist:.2f} A")
        
        if max_dist > 10.0:
            print(f"   WARNING: Cluster still appears fragmented (max dist > 10 A)")
        else:
            print(f"   OK: Cluster is now compact")
    
    def unwrap_trajectory(self):
        """
        Unwrap trajectory using minimum image convention
        Particularly important for cluster atoms that may cross periodic boundaries
        """
        print("\nUnwrapping trajectory using minimum image convention...")
        
        unwrapped_frames = [self.frames[0].copy()]  # First frame as reference
        
        for i in range(1, len(self.frames)):
            prev_frame = unwrapped_frames[-1]
            curr_frame = self.frames[i].copy()
            box = self.box_vectors[i]
            
            # Calculate displacement from previous frame
            delta = curr_frame - prev_frame
            
            # Apply minimum image convention
            # If displacement > box/2, atom crossed boundary
            for dim in range(3):
                # Atoms that jumped forward (crossed right boundary)
                mask_forward = delta[:, dim] > box[dim] / 2.0
                curr_frame[mask_forward, dim] -= box[dim]
                
                # Atoms that jumped backward (crossed left boundary)
                mask_backward = delta[:, dim] < -box[dim] / 2.0
                curr_frame[mask_backward, dim] += box[dim]
            
            unwrapped_frames.append(curr_frame)
        
        self.frames = np.array(unwrapped_frames)
        print(f"   Successfully unwrapped {len(self.frames)} frames")
        
        # Verify unwrap quality for cluster
        cluster_com = self.frames[:, self.cluster_ids, :].mean(axis=1)
        total_drift = np.linalg.norm(cluster_com[-1] - cluster_com[0])
        print(f"   Cluster center of mass total drift: {total_drift:.2f} A")
    
    def verify_unwrap(self):
        """Verify trajectory is properly unwrapped (detect abnormal cross-frame displacement)"""
        print("\nVerifying unwrap correctness...")
        
        issues = []
        
        for i in range(1, len(self.frames)):
            # Calculate inter-frame displacement
            displacement = self.frames[i] - self.frames[i-1]
            max_disp = np.abs(displacement).max(axis=1)
            
            # Detect abnormal displacement (threshold = 1/3 of box length)
            threshold = self.box_vectors[0].min() / 3.0
            bad_atoms = np.where(max_disp > threshold)[0]
            
            if len(bad_atoms) > 0:
                issues.append({
                    'frame': i,
                    'atoms': bad_atoms,
                    'max_disp': max_disp[bad_atoms].max()
                })
        
        if issues:
            print(f"   ERROR: Detected {len(issues)} frames with possible box-crossing jumps!")
            print(f"   First 5 problem frames:")
            for issue in issues[:5]:
                print(f"      Frame {issue['frame']}: {len(issue['atoms'])} atoms, max displacement {issue['max_disp']:.2f} A")
            print(f"\n   WARNING: Trajectory may not be properly unwrapped, please check LAMMPS output!")
            return False
        else:
            print(f"   OK: No abnormal box-crossing jumps detected, trajectory is properly unwrapped")
            
            # Additional verification: check cluster overall drift
            cluster_com = self.frames[:, self.cluster_ids, :].mean(axis=1)
            total_drift = np.linalg.norm(cluster_com[-1] - cluster_com[0])
            print(f"   OK: Cluster center of mass total drift: {total_drift:.2f} A")
            
            return True
    
    def center_trajectory(self, ref_group, output_file):
        """
        Center trajectory on the center of mass of specified atom group
        
        Parameters:
            ref_group: 'support' or 'cluster'
            output_file: output filename
        """
        if ref_group == 'support':
            ref_ids = self.support_ids
            desc = "Support"
        elif ref_group == 'cluster':
            ref_ids = self.cluster_ids
            desc = "Cluster"
        else:
            raise ValueError("ref_group must be 'support' or 'cluster'")
        
        print(f"\nGenerating {desc}-centered trajectory...")
        
        centered_frames = []
        
        for frame_idx, frame in enumerate(self.frames):
            # Calculate reference group center of mass
            com = frame[ref_ids].mean(axis=0)
            
            # Set box center as origin
            box_center = self.box_vectors[frame_idx] / 2.0
            
            # Translate entire system
            shift = box_center - com
            centered_frame = frame + shift
            
            centered_frames.append(centered_frame)
        
        centered_frames = np.array(centered_frames)
        
        # Write XYZ file
        self._write_xyz(centered_frames, output_file, f"{desc}-centered")
        
        print(f"   Saved to: {output_file}")
        
        return centered_frames
    
    def _write_xyz(self, frames, filename, description=""):
        """Write XYZ format trajectory"""
        with open(filename, 'w') as f:
            for frame_idx, frame in enumerate(frames):
                # Write atom count
                f.write(f"{self.total_atoms}\n")
                
                # Write Lattice information
                box = self.box_vectors[min(frame_idx, len(self.box_vectors)-1)]
                lattice_str = f'Lattice="{box[0]:.8f} 0.00000000 0.00000000 0.00000000 {box[1]:.8f} 0.00000000 0.00000000 0.00000000 {box[2]:.8f}"'
                f.write(f'{lattice_str} pbc="1 1 1" Properties=species:S:1:pos:R:3 # {description} Frame {frame_idx}\n')
                
                # Write atom coordinates
                for atom_idx, (atom_type, coords) in enumerate(zip(self.atom_types, frame)):
                    f.write(f"{atom_type} {coords[0]:.8f} {coords[1]:.8f} {coords[2]:.8f}\n")
    
    def generate_index_file(self, output_file="index_zsplit.ndx"):
        """Generate GROMACS-style index file"""
        print(f"\nGenerating index file: {output_file}")
        
        with open(output_file, 'w') as f:
            # Support group (atom ID starts from 1)
            f.write("[ Support ]\n")
            for i, atom_id in enumerate(self.support_ids + 1):
                f.write(f"{atom_id:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n\n")
            
            # Cluster group
            f.write("[ PtSnCluster ]\n")
            for i, atom_id in enumerate(self.cluster_ids + 1):
                f.write(f"{atom_id:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n\n")
            
            # System group
            f.write("[ System ]\n")
            for i in range(self.total_atoms):
                f.write(f"{i+1:5d} ")
                if (i + 1) % 15 == 0:
                    f.write("\n")
            f.write("\n")
        
        print(f"   Index file saved")
    
    def run(self):
        """Execute complete analysis workflow"""
        print("=" * 70)
        print("  Sn7Pt4O3/Al2O3 Trajectory Processing - Dual-Mode Centering")
        print("=" * 70)
        
        # 1. Read trajectory
        self.read_xyz()
        
        # 2. Identify atom groups
        self.identify_groups()
        
        # 3. Reassemble cluster if it's fragmented in first frame
        self.make_cluster_whole()
        
        # 4. Check if unwrap is needed
        is_valid = self.verify_unwrap()
        
        if not is_valid:
            print("\n" + "=" * 70)
            print("  AUTO-UNWRAPPING TRAJECTORY")
            print("=" * 70)
            self.unwrap_trajectory()
            
            # Verify again after unwrapping
            print("\nRe-verifying after unwrap...")
            is_valid_after = self.verify_unwrap()
            
            if not is_valid_after:
                print("\n   WARNING: Still detected issues after unwrapping")
                response = input("\nContinue processing anyway? (y/n): ")
                if response.lower() != 'y':
                    print("User cancelled operation")
                    return
        
        # 4. Generate index file
        self.generate_index_file()
        
        # 5. Generate support-centered trajectory
        self.center_trajectory('support', 'surface_centered.xyz')
        
        # 6. Generate cluster-centered trajectory
        self.center_trajectory('cluster', 'cluster_centered.xyz')
        
        print("\n" + "=" * 70)
        print("Processing complete! Generated files:")
        print("   - index_zsplit.ndx       (atom group indices)")
        print("   - surface_centered.xyz   (support-centered, for surface migration)")
        print("   - cluster_centered.xyz   (cluster-centered, for internal rearrangement)")
        print("=" * 70)


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_ptsn_trajectory.py <trajectory.xyz> [support_n] [cluster_n]")
        print("\nExamples:")
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