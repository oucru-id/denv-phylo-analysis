import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
import re
from Bio import Phylo
import matplotlib.patches as mpatches
import math

def generate_plots(df, output_prefix):
    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    distances = df.where(mask).stack().values
    
    # Heatmap
    n_samples = df.shape[0]
    fig_dim = max(12, n_samples * 0.6)
    plt.figure(figsize=(fig_dim, fig_dim))
    font_size = 10 if n_samples < 20 else 8
    
    sns.heatmap(
        df, 
        cmap="viridis", 
        annot=True, 
        fmt="d", 
        square=True,
        cbar_kws={"shrink": 0.8},
        annot_kws={"size": font_size}
    )
    
    plt.title("SNP Distance Heatmap")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap.png")
    plt.close()

    # Violin Plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(y=distances, color="lightblue", inner="quartile")
    plt.title("Plot of SNP Distances")
    plt.ylabel("SNP Distance")
    plt.savefig(f"{output_prefix}_violin.png")
    plt.close()

def get_serotype(conclusion):
    if pd.isna(conclusion): return "Unknown"
    conclusion = str(conclusion).upper()
    
    if "DENV-1" in conclusion or "DENV1" in conclusion: return "DENV-1"
    if "DENV-2" in conclusion or "DENV2" in conclusion: return "DENV-2"
    if "DENV-3" in conclusion or "DENV3" in conclusion: return "DENV-3"
    if "DENV-4" in conclusion or "DENV4" in conclusion: return "DENV-4"
    
    if "NC_001477" in conclusion: return "DENV-1"
    if "NC_001474" in conclusion: return "DENV-2"
    if "NC_001475" in conclusion: return "DENV-3"
    if "NC_002640" in conclusion: return "DENV-4"
    
    return "Unknown"

def get_lineage_colors(tree, meta_df):
    lineage_map = {}
    sample_serotype = {}
    for _, row in meta_df.iterrows():
        sample_serotype[row['sample_id']] = get_serotype(row.get('conclusion'))

    for clade in tree.get_terminals():
        lineage_map[clade.name] = sample_serotype.get(clade.name, "Unknown")

    unique_clades = ["DENV-1", "DENV-2", "DENV-3", "DENV-4", "Unknown"]
    palette = ["#4C78A8", "#59A14F", "#EDC948", "#F28E2B", "#BAB0AC"]
    color_map_mpl = {clade: color for clade, color in zip(unique_clades, palette)}
    
    return lineage_map, color_map_mpl, unique_clades

def get_legend_labels():
    return {
        "DENV-1": "DENV-1 (Ref: NC_001477.1)",
        "DENV-2": "DENV-2 (Ref: NC_001474.2)",
        "DENV-3": "DENV-3 (Ref: NC_001475.3)",
        "DENV-4": "DENV-4 (Ref: NC_002640.4)",
        "Unknown": "Unknown"
    }

def plot_rectangular_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    def to_bio_color(mpl_color):
        if isinstance(mpl_color, str) and mpl_color.startswith('#'):
            h = mpl_color.lstrip('#')
            return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
        return tuple(int(x * 255) for x in mpl_color)
    
    color_map_bio = {k: to_bio_color(v) for k, v in color_map_mpl.items()}
    gray_bio = (128, 128, 128)

    def color_clade(clade):
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_bio.get(l, gray_bio)
            clade.color = c
            return c
        else:
            child_colors = [color_clade(c) for c in clade]
            first_color = child_colors[0]
            if all(c == first_color for c in child_colors):
                clade.color = first_color
                return first_color
            else:
                clade.color = gray_bio 
                return gray_bio

    color_clade(tree.root)

    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    def get_branch_label(clade):
        if clade.branch_length and clade.branch_length > 0.001:
            return f"{clade.branch_length:.3f}"
        return None

    Phylo.draw(
        tree, 
        axes=ax, 
        do_show=False, 
        show_confidence=False,
        label_func=lambda x: x.name if x.is_terminal() else "",
        branch_labels=get_branch_label,
    )
    
    terminals = tree.get_terminals()
    for i, clade in enumerate(terminals):
        y_pos = i + 1
        x_pos = tree.distance(tree.root, clade)
        l = lineage_map.get(clade.name, "Unknown")
        c_mpl = color_map_mpl.get(l, "gray")
        ax.scatter(x_pos, y_pos, color=c_mpl, s=80, zorder=10, edgecolors='white', linewidth=0.5)

    legend_labels = get_legend_labels()
    handles = [mpatches.Patch(color=color_map_mpl[c], label=legend_labels.get(c, c)) for c in unique_clades]
    plt.legend(handles=handles, title="Serotype", loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    
    plt.title("Phylogenetic Tree (Rectangular)", fontsize=14)
    plt.xlabel("Genetic Distance", fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def get_coords(tree):
    coords = {}
    leaves = tree.get_terminals()
    total_leaves = len(leaves)
    
    for i, leaf in enumerate(leaves):
        angle = (2 * math.pi * i) / total_leaves
        coords[leaf] = {'theta': angle}

    for clade in tree.get_nonterminals(order='postorder'):
        children_angles = [coords[c]['theta'] for c in clade.clades]
        if children_angles:
            avg_angle = sum(children_angles) / len(children_angles)
            coords[clade] = {'theta': avg_angle}

    coords[tree.root]['r'] = 0
    for clade in tree.get_nonterminals(order='preorder'):
        parent_r = coords[clade]['r']
        for child in clade.clades:
            length = child.branch_length if child.branch_length else 0.01
            coords[child] = coords.get(child, {})
            coords[child]['r'] = parent_r + length
            
    return {k: (v['r'], v['theta']) for k,v in coords.items()}

def plot_circular_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    coords = get_coords(tree)
    
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111, projection='polar')
    
    ax.set_frame_on(False)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    max_r = max(r for r, t in coords.values())
    
    parents = {c: p for p in tree.find_clades() for c in p.clades}

    for clade in tree.find_clades(order='level'):
        if clade == tree.root: continue
        parent = parents.get(clade)
        if not parent: continue

        r1, t1 = coords[parent]
        r2, t2 = coords[clade]
        
        if abs(t1 - t2) > 0:
            theta_range = np.linspace(t1, t2, num=20)
            r_range = [r1] * len(theta_range)
            ax.plot(theta_range, r_range, color='gray', linewidth=0.5)
            
        ax.plot([t2, t2], [r1, r2], color='gray', linewidth=0.5)

    for clade, (r, theta) in coords.items():
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_mpl.get(l, "gray")
            
            ax.scatter(theta, r, color=c, s=40, zorder=10, edgecolors='white', linewidth=0.5)
            
            rot = math.degrees(theta)
            if 90 < rot < 270:
                rot += 180
                ha = 'right'
                label_r = r + (max_r * 0.02)
            else:
                ha = 'left'
                label_r = r + (max_r * 0.01)
                
            ax.text(theta, label_r, clade.name, rotation=rot, ha=ha, va='center', fontsize=8)

    legend_labels = get_legend_labels()
    handles = [mpatches.Patch(color=color_map_mpl[c], label=legend_labels.get(c, c)) for c in unique_clades]
    fig.legend(handles=handles, title="Serotype", loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False)
    
    plt.title("Phylogenetic Tree (Circular)", fontsize=16, y=1.05)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def plot_unrooted_tree(tree, lineage_map, color_map_mpl, unique_clades, output_file):
    coords = get_coords(tree)
    
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.set_aspect('equal')
    ax.axis('off')
    
    cart_coords = {}
    max_r = 0
    for clade, (r, theta) in coords.items():
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        cart_coords[clade] = (x, y)
        if r > max_r: max_r = r

    parents = {c: p for p in tree.find_clades() for c in p.clades}

    for clade in tree.find_clades(order='level'):
        if clade == tree.root: continue
        parent = parents.get(clade)
        if not parent: continue

        x1, y1 = cart_coords[parent]
        x2, y2 = cart_coords[clade]
        
        ax.plot([x1, x2], [y1, y2], color='gray', linewidth=0.5)

    for clade, (x, y) in cart_coords.items():
        if clade.is_terminal():
            l = lineage_map.get(clade.name, "Unknown")
            c = color_map_mpl.get(l, "gray")
            
            ax.scatter(x, y, color=c, s=40, zorder=10, edgecolors='white', linewidth=0.5)
            
            r, theta = coords[clade]
            rot = math.degrees(theta)
            
            if 90 < rot < 270:
                rot += 180
                ha = 'right'
                lx = x + (max_r * 0.02) * math.cos(theta)
                ly = y + (max_r * 0.02) * math.sin(theta)
            else:
                ha = 'left'
                lx = x + (max_r * 0.01) * math.cos(theta)
                ly = y + (max_r * 0.01) * math.sin(theta)
            
            ax.text(lx, ly, clade.name, rotation=rot, ha=ha, va='center', fontsize=8, rotation_mode='anchor')

    legend_labels = get_legend_labels()
    handles = [mpatches.Patch(color=color_map_mpl[c], label=legend_labels.get(c, c)) for c in unique_clades]
    fig.legend(handles=handles, title="Serotype", loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False)
    
    plt.title("Phylogenetic Tree (Radial)", fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_phylo_trees(tree_file, meta_df, output_prefix):
    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"Error reading tree file: {e}")
        return

    target_root = None
    for leaf in tree.get_terminals():
        if "NC_002640" in leaf.name:
            target_root = leaf.name
            break
    
    if target_root:
        print(f"Re-rooting tree on {target_root} (DENV-4)...")
        tree.root_with_outgroup(target_root)
    else:
        print("DENV-4 reference not found, using midpoint rooting...")
        tree.root_at_midpoint()
    
    tree.ladderize()

    lineage_map, color_map, unique_clades = get_lineage_colors(tree, meta_df)
    
    plot_rectangular_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_rectangular.png")
    plot_unrooted_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_unrooted.png")
    plot_circular_tree(tree, lineage_map, color_map, unique_clades, f"{output_prefix}_circular.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help="Path to distance_matrix.tsv")
    parser.add_argument('--metadata', required=True, help="Path to metadata.tsv")
    parser.add_argument('--tree', required=False, help="Path to phylo_tree.nwk")
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', index_col=0)
    df.index.name = "Sample"
    
    meta_df = pd.read_csv(args.metadata, sep='\t')

    generate_plots(df, "stats")
    
    if args.tree:
        generate_phylo_trees(args.tree, meta_df, "phylo_tree")

if __name__ == "__main__":
    main()