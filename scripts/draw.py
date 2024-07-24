import argparse
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import to_hex
from copy import deepcopy
import seaborn as sns
import numpy as np

special_nodes=["A", "E", "J", "Q", "L0", "L1", "L2", "L3", "Z", "H", "F", "D", "M", "B", "C", "Z", "G"]

special_nodes += [
    'Node5515', 'I1a1e', 'Node5517', 'Node5516', 'I1a*2', 'I1a', 'Node5519', 'Node5518', 
 'I1a1d', 'Node5525', 'I1a1a3a', 'Node5520', 'I1a1', 'I1a1c', 'Node5524', 'I1a1a1', 
 'I1a1a', 'Node5523', 'I1a1b', 'Node5521', 'Node5514', 'Node5522', 'I1a1a3', 'I1a1a2', 
 'I1a*1', 'Node5481', 'Node5480', 'Node5482', 'Node5483', 'Node5484', 'Node5545', 'I2e', 
 'Node5534', 'Node5533', 'Node5542', 'Node5541', 'Node5543', 'Node5535', 'Node5536', 
 'Node5531', 'I7', 'Node5526', 'Node5527', 'I', 'Node5479', 'Node5547', 'Node5551', 
 'I2d', 'I2', 'I2f', 'Node5503', 'Node5501', 'Node5544', 'I3b', 'Node5537', 'I2b', 
 'Node5509', 'Node5538', 'I2c', 'Node5546', 'I3', 'Node5540', "I2'3", "I3c", "Node5502", 
 'I1e', 'I3d1', 'Node5539', 'I3d', 'Node5506', 'I1f', 'Node5510'
                 ]

def hardcode_prune(tree, clades_to_prune):
    for clade_name in clades_to_prune:
        clade = next(tree.find_clades(name=clade_name), None)
        if clade:
            try:
                tree.prune(clade)
                print(f"Pruned hardcoded clade: {clade_name}")
            except ValueError as e:
                pass
                #print(f"Error pruning clade {clade_name}: {e}")
    return tree

def prune_internal_node(tree, node_name, special_nodes, data_nodes):
    internal_node = next(tree.find_clades(name=node_name), None)
    if internal_node is None:
        return tree

    if node_name in special_nodes:
        return tree  # Skip pruning for special internal nodes

    terminal_nodes = [terminal.name for terminal in internal_node.get_terminals()]

    for terminal_node in terminal_nodes:
        print(terminal_node)
        terminal_clade = next(tree.find_clades(name=terminal_node), None)
        if terminal_clade and (terminal_clade.name not in data_nodes) and (terminal_clade.name not in special_nodes):
            tree.collapse(terminal_clade)

    return tree

def collapse_clades(tree, special_nodes, data_nodes):
    def is_essential_clade(clade):
        if clade.name in special_nodes or clade.name in data_nodes:
            return True
        for descendant in clade.get_terminals():
            if descendant.name in special_nodes or descendant.name in data_nodes:
                return True
        return False

    essential_clades = set()
    for clade in tree.get_nonterminals(order='postorder'):
        if is_essential_clade(clade):
            essential_clades.add(clade)
            parent = get_parent(tree, clade)
            if parent:
                essential_clades.add(parent)

    for clade in tree.get_nonterminals(order='level'):
        if clade not in essential_clades:
            prune_internal_node(tree, clade.name, special_nodes, data_nodes)

    # Print the essential nonterminal clades and reasons
    for clade in essential_clades:
        if not clade.is_terminal():
            if clade.name in special_nodes:
                reason = "clade name is in special_nodes"
            elif clade.name in data_nodes:
                reason = "clade name is in data_nodes"
            else:
                reason = "contains descendants in special_nodes or data_nodes"
            print(f"Essential clade: {clade.name}, Reason: {reason}")

    return tree

def strip_suffix(node_name):
    if '.' in node_name:
        return node_name.split('.')[0]
    return node_name


def read_mcmc_data(file_path, burnin_fraction=0.1):
    with open(file_path, 'r') as file:
        header = next(file).strip().split()
        source_indices = [i for i, x in enumerate(header) if 'Source' in x]
        source_names = [header[i] for i in source_indices]
        node_data = {source_name: [] for source_name in source_names}
        branch_iterations = {source_name: {} for source_name in source_names}
        internal_nodes = set()

        lines = file.readlines()
        skip_lines = int(burnin_fraction * len(lines))

        for line in lines[skip_lines:]:
            parts = line.strip().split()
            for source_idx, source_name in zip(source_indices, source_names):
                node = strip_suffix(parts[source_idx])
                if node.startswith("Node"):
                    internal_nodes.add(node)
                node_data[source_name].append(node)
                if node not in branch_iterations[source_name]:
                    branch_iterations[source_name][node] = 0
                branch_iterations[source_name][node] += 1

    print("MCMC data read")

    for source_name in source_names:
        print(f"Source: {source_name}, Iterations: {branch_iterations[source_name]}")

    return node_data, source_names, branch_iterations, internal_nodes

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2] if len(node_path) > 1 else None


def prune_tree(tree, source_clades, data_nodes=[], max_nodes=120, special_nodes=special_nodes):
    pruned_tree = deepcopy(tree)
    marked_clades = set(special_nodes) | set(data_nodes)
    missing_nodes = []  # List to store missing nodes

    print("Starting pruning process...")

    # Identify and mark all clades to keep
    total_data_nodes = len(data_nodes)

    for index, node_name in enumerate(data_nodes):
        node = next(pruned_tree.find_clades(name=node_name), None)
        if node:
            if not node.is_terminal():
                marked_clades.update([descendant.name for descendant in node.get_terminals()])
            #print(f"Marked clade: {node_name}")
        else:
            missing_nodes.append(node_name)  # Add missing nodes to the list
            #print(f"Missing clade: {node_name}")
        # Progress tracking
        progress = (index + 1) / total_data_nodes * 100
        print(f"Progress: Marking essential clades {progress:.2f}% complete")

    # Set to store essential clades
    essential_clades = set()

    # Function to mark essential clades recursively
    def mark_essential_clades(clade):
        if clade.name in marked_clades:
            essential_clades.add(clade)
            parent = get_parent(pruned_tree, clade)
            if parent:
                mark_essential_clades(parent)
            #print(f"Essential clade marked: {clade.name}")

    # Mark essential clades starting from data nodes and special nodes
    total_marked_clades = len(marked_clades)
    for index, node_name in enumerate(marked_clades):
        clade = next(pruned_tree.find_clades(name=node_name), None)
        if clade:
            mark_essential_clades(clade)
        # Progress tracking
        progress = (index + 1) / total_marked_clades * 100
        print(f"Progress: Marking essential clades recursively {progress:.2f}% complete")

    print("Initial essential clades marked. Starting pruning of non-essential clades...")

    # Prune non-essential clades
    all_clades = list(pruned_tree.find_clades(order='level'))
    total_clades = len(all_clades)
    for index, clade in enumerate(all_clades):
        if (clade not in essential_clades and clade != pruned_tree.root):
            try:
                pruned_tree.prune(clade)
                print(f"Pruned clade: {clade.name}")
            except ValueError as e:
                #print(f"Error pruning clade {clade.name}: {e}")
                pass
        # Progress tracking
        progress = (index + 1) / total_clades * 100
        print(f"Progress: Pruning non-essential clades {progress:.2f}% complete")

    # If tree size still exceeds max_nodes, collapse further
    total_clades = len(pruned_tree.get_terminals()) + len(pruned_tree.get_nonterminals())
    print(f"Tree size after initial pruning: {total_clades}")

    if total_clades > max_nodes:
        print(f"Tree size exceeds max_nodes ({max_nodes}). Collapsing clades...")
        pruned_tree = collapse_clades(pruned_tree, special_nodes, data_nodes)
    else:
        print(f"No further pruning needed. Tree size is within the limit.")

    # Hardcoded list of clades to prune
    clades_to_prune = []
    pruned_tree = hardcode_prune(pruned_tree, clades_to_prune)

    final_clade_count = len(pruned_tree.get_terminals()) + len(pruned_tree.get_nonterminals())
    print(f"Final tree size: {final_clade_count}")

    if missing_nodes:  # Print warning if there are any missing nodes
        print("Warning: The following node names from the data could not be found in the tree file:")
        print(missing_nodes)

    print("Pruning process completed.")

    return pruned_tree


def visualize_tree_branches(tree, source_clades, source_names, branch_iterations, output_file, title_suffix, max_labels=120):

    def label_filter(clade):
        #if clade.name.startswith("Node"):
        #    return clade.name
        if clade.name in labeled_nodes:
            return clade.name
        return None

    data_nodes = set()
    for source, nodes in branch_iterations.items():
        data_nodes.update(nodes.keys())

    prioritized_nodes = data_nodes.union(set(special_nodes))
    all_terminal_nodes = set([clade.name for clade in tree.get_terminals()])

    # Separate nodes with data and other prioritized nodes
    nodes_with_data = data_nodes.intersection(prioritized_nodes)
    leaf_nodes_with_data = nodes_with_data.intersection(all_terminal_nodes)
    internal_nodes_with_data = nodes_with_data.difference(all_terminal_nodes)
    print(leaf_nodes_with_data)

    # Combine and prioritize nodes with data first
    labeled_nodes = leaf_nodes_with_data.copy()
    remaining_labels = max_labels - len(labeled_nodes)

    if remaining_labels > 0:
        additional_nodes = list(internal_nodes_with_data)[:remaining_labels]
        labeled_nodes.update(additional_nodes)
        remaining_labels -= len(additional_nodes)

    # If there is still room, add other prioritized nodes
    if remaining_labels > 0:
        remaining_prioritized_nodes = prioritized_nodes.difference(labeled_nodes)
        additional_nodes = list(remaining_prioritized_nodes)[:remaining_labels]
        labeled_nodes.update(additional_nodes)

    fig, axes = plt.subplots(figsize=(20, 18))

    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xticklabels([])
    axes.set_yticklabels([])
    axes.set_xlabel('')
    axes.set_ylabel('')

    # Set font properties
    plt.rc('font', size=17, weight='bold')

    source_colors = sns.color_palette("Set1", len(source_names))
    source_colors_hex = [to_hex(color) for color in source_colors]
    source_to_color = {source: color for source, color in zip(source_names, source_colors_hex)}
    # In the visualize_tree_branches function:
    min_iterations = min(min(count for count in branch_iterations[source].values() if count > 0) for source in source_names)
    max_iterations = max(max(branch_iterations[source].values()) for source in source_names)

    def get_color_intensity(iteration, min_iterations, max_iterations):
        if iteration == 0:
            return 0  # No data, no color

        # Directly scale the iteration count to the 0-7 range
        return 7 * (iteration / max_iterations)

    colormap = plt.cm.viridis  # Choose a colormap. Options include viridis, plasma, inferno, magma, etc.


    for clade in tree.find_clades():
        clade.color = 'gray'
        clade.width = 1
        for source, nodes in source_clades.items():
            if clade.name in nodes:
                iterations = branch_iterations[source].get(clade.name, 0)
                clade.width = 1 + 6 * (iterations / max_iterations)
                if len(source_names) == 1:
                    # Gradient coloring by iteration count using a colormap
                    color_intensity = get_color_intensity(iterations, min_iterations, max_iterations)
                    print(clade.name, iterations, color_intensity, clade.width)
                    rgba_color = colormap(color_intensity / 7)
                    clade.color = to_hex(rgba_color)
                else:
                    # Fixed color per source
                    clade.color = source_to_color[source]
                    print(clade.name, iterations, clade.color, clade.width)
                break

    Phylo.draw(tree, axes=axes, label_func=label_filter, do_show=False,
               label_colors={'color': 'black', 'weight': 'bold'})

    modified_source_names = [source.replace('_', ' ') for source in source_names]

    if len(source_names) > 1:  # Only create and display legend if there are multiple sources
        legend_handles = [plt.Line2D([0], [0], color=color, lw=4) for color in source_colors_hex]
        plt.legend(handles=legend_handles,
                   labels=modified_source_names, loc='lower right', fontsize=20, title_fontsize='20')

    # Adding the colorbar if there is only one source
    if len(source_names) == 1:
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=min_iterations, vmax=max_iterations))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Sampling Frequency', fontsize=16)  # Updated label
        tick_count = 5  # Number of ticks on the colorbar
        tick_positions = np.linspace(min_iterations, max_iterations, tick_count)
        tick_labels = [f'{int(val)}' for val in tick_positions]
        cbar.set_ticks(tick_positions)
        cbar.set_ticklabels(tick_labels)

    plt.suptitle(f"MCMC Source Placements, {title_suffix}", fontsize=22, fontweight='bold', y=0.92)
    plt.ylabel("")
    plt.xlabel("")
    plt.savefig(output_file, dpi=400)
    print(f"Tree visualized and saved as {output_file}.")

def main():
    parser = argparse.ArgumentParser(description='Process phylogenetic data.')
    parser.add_argument('mcmc_data_file', type=str, help='MCMC data file path')
    parser.add_argument('newick_file', type=str, help='Newick tree file path')
    parser.add_argument('output_file', type=str, help='Output file path')
    parser.add_argument('title_suffix', type=str, help='Title suffix for the plot')
    parser.add_argument('--max-nodes', type=int, default=120, help='Maximum number of nodes in the pruned tree')
    parser.add_argument('--max-labels', type=int, default=120, help='Maximum number of labels in the plot')
    parser.add_argument('--burnin-fraction', type=float, default=0.0, help='Fraction of iterations to exclude as burn-in')
    args = parser.parse_args()

    source_clades, source_names, branch_iterations, internal_nodes = read_mcmc_data(args.mcmc_data_file, args.burnin_fraction)

    original_tree = Phylo.read(args.newick_file, 'newick')
    small_tree = prune_tree(original_tree, source_clades, internal_nodes, max_nodes=args.max_nodes)

    visualize_tree_branches(small_tree, source_clades, source_names, branch_iterations, args.output_file, args.title_suffix, max_labels=args.max_labels)
    print("Tree visualized and saved.")

if __name__ == '__main__':
    main()
