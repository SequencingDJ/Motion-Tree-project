import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram




def trajectory_load(TOPOLOGY, TRAJECTORY, SELECTION):
    print("MOTION TREE ANALYSIS")



    u = mda.Universe(
    TOPOLOGY,
    TRAJECTORY
)

    ca = u.select_atoms(SELECTION)

    n_atoms = len(ca)

    n_total_frames = len(u.trajectory)

    print(
    f"Number of selected C-alpha atoms: {n_atoms}",
)

    print(
    f"Number of trajectory frames: {n_total_frames}",
)
    return u,ca,n_atoms,n_total_frames

def selecting_frames(FRAME_STEP, n_total_frames):
    print(
    "\nSelecting trajectory frames...",
    flush=True
)

    frame_indices = np.arange(
    0,
    n_total_frames,
    FRAME_STEP,
    dtype=int
)

    n_frames = len(frame_indices)

    print(
    f"Frames used: {n_frames}",
    flush=True
)

    print(
    f"Frame step: {FRAME_STEP}",
    flush=True
)
    return frame_indices,n_frames


def residule_info(ca):
    print(
    "\nSelected residue information:",
    flush=True
)

    print(
    f"First residue ID: {ca.resids[0]}",
    flush=True
)

    print(
    f"Last residue ID:  {ca.resids[-1]}",
    flush=True
)

    np.savetxt(
    "residue_numbers.txt",
    ca.resids,
    fmt="%d"
)

    print(
    "Residue numbers saved to residue_numbers.txt",
    flush=True
)

def c_alpha_pairs(n_atoms):
    print(
    "\nGenerating unique C-alpha atom pairs...",
    flush=True
)

# Only calculate the upper triangle:
#
#       N(N-1)/2
#
# unique pairs.

    i_idx, j_idx = np.triu_indices(
    n_atoms,
    k=1
)

    n_pairs = len(i_idx)

    print(
    f"Number of unique C-alpha pairs: {n_pairs:,}",
    flush=True
)
    return i_idx,j_idx,n_pairs

def get_pairs(n_atoms):
    print(
    "\nGenerating unique C-alpha atom pairs...",
    flush=True
)

# Only the upper triangle:


    i_idx, j_idx = np.triu_indices(
    n_atoms,
    k=1
)

    n_pairs = len(i_idx)

    print(
    f"Number of unique C-alpha pairs: {n_pairs:,}",
    flush=True
)
    return i_idx,j_idx,n_pairs



# The following code relates to marking the top 7 nodes on the MT graph, this is for the purpose of creating pretty figures and is specific to ADK.
# This can be put back above 'graph settings' in the actual .py file, if required.

# It was created through the use of AI, however it has little impact on the functionality on the actual MT and can be ignored.
""" # ============================================================
# MARK ONLY THE TOPMOST NODES
# ============================================================

 LABEL_INTERNAL_NODES = True



if LABEL_INTERNAL_NODES:

    print(
        "\nMarking topmost dendrogram nodes...",
        flush=True
    )

    # --------------------------------------------------------
    # Determine the depth of every internal node
    #
    # Root node:
    #     depth 0
    #
    # Its two children:
    #     depth 1
    #
    # Their children:
    #     depth 2
    #
    # This gives a maximum of:
    #
    #     1 + 2 + 4 = 7 nodes
    #
    # matching the nodes shown in the reference figure.
    # --------------------------------------------------------

    node_depth = {}

    # The final/root node is the last internal node
    root_node = 2 * n_atoms - 2

    def calculate_node_depth(node_id):

        # Leaf nodes are not internal nodes
        if node_id < n_atoms:
            return None

        # Already calculated
        if node_id in node_depth:
            return node_depth[node_id]

        # Linkage row corresponding to this internal node
        row_number = node_id - n_atoms

        child_a = int(linkage_matrix[row_number, 0])
        child_b = int(linkage_matrix[row_number, 1])

        # Root has depth 0
        if node_id == root_node:
            depth = 0

        else:
            # Parent depth + 1
            parent_depth = None

            for parent_row_number, parent_row in enumerate(
                linkage_matrix
            ):

                parent_id = n_atoms + parent_row_number

                parent_child_a = int(parent_row[0])
                parent_child_b = int(parent_row[1])

                if (
                    parent_child_a == node_id
                    or parent_child_b == node_id
                ):
                    parent_depth = calculate_node_depth(parent_id)
                    break

            if parent_depth is None:
                raise RuntimeError(
                    f"Could not determine parent of node {node_id}"
                )

            depth = parent_depth + 1

        node_depth[node_id] = depth

        return depth

    # Calculate depth for all internal nodes
    for row_number in range(len(linkage_matrix)):

        node_id = n_atoms + row_number

        calculate_node_depth(node_id)

    # --------------------------------------------------------
    # Keep only the top three levels:
    #
    # depth 0 -> 1 node
    # depth 1 -> 2 nodes
    # depth 2 -> 4 nodes
    #
    # Total = 7 nodes
    # --------------------------------------------------------

    top_nodes = [
        node_id
        for node_id, depth in node_depth.items()
        if depth <= 2
    ]
    # 7 nodes
    # Sort from root downward
    top_nodes.sort(
        key=lambda node_id: (
            node_depth[node_id],
            node_id
        )
    )

    print(
        f"Topmost nodes selected: {len(top_nodes)}",
        flush=True
    )

    print(
        "Node depths:",
        [(node_id, node_depth[node_id]) for node_id in top_nodes],
        flush=True
    )


y_min, y_max = ax.get_ylim()
y_range = y_max - y_min

if y_range <= 0:
        y_range = 1.0

for display_number, node_id in enumerate(
        top_nodes,
        start=1
    ):

        x_position = node_x[node_id]
        y_position = node_y[node_id]

        # Circle at the actual merge point
        ax.scatter(
            x_position,
            y_position,
            s=8,
            zorder=10
        )

        # Number the selected nodes 1–7
        ax.text(
            x_position + 8,
            y_position,
            str(display_number),
            ha="left",
            va="center",
            fontsize=8,
            fontweight="bold",
            
        )

 """
""" # ============================================================
# DETERMINE VISUAL POSITION OF EVERY INTERNAL NODE
# ============================================================

# The dendrogram result gives the actual rendered leaf order.
# We use that order to determine the x-coordinate of every node.

leaves = dendro["leaves"]

leaf_position = {
    int(atom_index): int(position)
    for position, atom_index
    in enumerate(leaves)
}

node_x = {}
node_y = {}

# With scipy's standard dendrogram coordinates:
#
#   first leaf  = x 5
#   second leaf = x 15
#   etc.

for atom_index, position in leaf_position.items():

    node_x[atom_index] = (
        5.0 + 10.0 * position
    )

    node_y[atom_index] = 0.0


# Every linkage row creates one internal node.
for row_number, row in enumerate(linkage_matrix):

    child_a = int(row[0])
    child_b = int(row[1])

    node_id = n_atoms + row_number

    node_x[node_id] = (
        node_x[child_a]
        +
        node_x[child_b]
    ) / 2.0

    node_y[node_id] = float(row[2])

"""