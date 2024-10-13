import sys
from ete3 import Tree

def reroot_tree(newick_file, branch_id, output_file):
    # Load the tree
    tree = Tree(newick_file, format=1)

    # Find the node or branch to reroot on
    target_node = tree.search_nodes(name=branch_id)
    
    if not target_node:
        print(f"Branch/Node ID '{branch_id}' not found in the tree.")
        sys.exit(1)
    
    # Reroot the tree
    tree.set_outgroup(target_node[0])

    # Write the rerooted tree to the output file
    tree.write(outfile=output_file, format=1)
    print(f"Tree successfully rerooted on '{branch_id}' and saved to '{output_file}'.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python reroot_tree.py <newick_file> <branch_id> <output_file>")
        sys.exit(1)

    newick_file = sys.argv[1]
    branch_id = sys.argv[2]
    output_file = sys.argv[3]

    reroot_tree(newick_file, branch_id, output_file)

