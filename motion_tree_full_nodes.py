import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from modules.calpha import *


# ai was used for:
# the progress bars, 
# outputting cluster information for chimerax & numbering of the topmost nodes under 'topmost nodes' in calpha module.



# ============================================================
# SETTINGS

# ------------------------------------------------------------
# Input files


TOPOLOGY = "protein.gro"
TRAJECTORY = "md_aligned.xtc"

# ------------------------------------------------------------
# Atom selection

SELECTION = "protein and name CA"
# This can be changed to instead computer the AAMT.

# ------------------------------------------------------------
# Motion Tree parameters
# ------------------------------------------------------------
# Number of largest Dmn values used for the cluster (step 3)
N_LARGEST = 20 

FRAME_STEP = 1

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
DISTANCE_FLUCTUATIONS_FILE = "Dmn.npy"

LINKAGE_FILE = "motion_tree_linkage.npy"
NODE_LEAVES_FILE = "motion_tree_node_leaves.txt"
FIGURE_FILE = "motion_tree.png"



#---------------
# trajectory loading

u, ca, n_atoms, n_total_frames = trajectory_load(TOPOLOGY, TRAJECTORY, SELECTION)


frame_indices, n_frames = selecting_frames(FRAME_STEP, n_total_frames)


# ============================================================
# GENERATE UNIQUE C-ALPHA PAIRS
# ============================================================



i_idx, j_idx, n_pairs = get_pairs(n_atoms)

# ============================================================
# Dmn calc
# ============================================================

print(
    "\nCalculating distance fluctuations Dmn...",
    flush=True
)

sum_distances = np.zeros(
    n_pairs,
    dtype=np.float64
)

sum_squared_distances = np.zeros(
    n_pairs,
    dtype=np.float64
)


# ============================================================
# loop trajectory frames
# ============================================================

for frame_number, trajectory_frame in enumerate(
    frame_indices
):

    
    # Load frame
    

    u.trajectory[trajectory_frame]

    
    # C-alpha coordinates EQUATION 1
    

    positions = ca.positions

    
    # Pairwise displacement vectors EQUATION 2
    

    delta = (
        positions[i_idx]
        -
        positions[j_idx]
    )

    
    # Pairwise distances EQUATION 3
     

    current_distances = np.sqrt(
        np.sum(
            delta * delta,
            axis=1
        )
    )

    
    # SUMS EQUATION 4 
    

    sum_distances += current_distances

    sum_squared_distances += (
        current_distances
        *
        current_distances
    )

    
    # Progress
    

    if (
        frame_number % 1000 == 0
        or
        frame_number == n_frames - 1
    ):

        percentage = (
            100.0
            *
            (frame_number + 1)
            /
            n_frames
        )

        print(
            f"Frame "
            f"{frame_number + 1:>8}/{n_frames:<8}"
            f" "
            f"({percentage:6.2f}%)",
            flush=True
        )


# ============================================================
# CALCULATE Dmn
# ============================================================

print(
    "\nAll trajectory frames processed.",
    flush=True
)

print(
    "Calculating mean distances...",
    flush=True
)
# Mean EQUATION 5
mean_distances = (
    sum_distances
    /
    n_frames
)
# Mean squared EQUATION 6
mean_squared_distances = (
    sum_squared_distances
    /
    n_frames
)

# Variance: EQUATION 7


variance = (
    mean_squared_distances
    -
    mean_distances * mean_distances
)

# floating-point precision.

variance = np.maximum(
    variance,
    0.0
)

# Distance fluctuation:


D_values = np.sqrt(
    variance
)


np.save(
    DISTANCE_FLUCTUATIONS_FILE,
    D_values
)

def get_D_value(a, b):

    """
    memory efficient method to extract matrix values
    """

    if a == b:
        return 0.0

    if a > b:
        a, b = b, a

    index = (
        n_atoms * a
        -
        (a * (a + 1)) // 2
        +
        b
        -
        a
        -
        1
    )

    return D_values[index]


# ============================================================
# MOTION TREE CLUSTER DISSIMILARITY
# ============================================================

def cluster_dissimilarity(
    cluster_a,
    cluster_b
):

    """
    Step 2 of the Algorithm 
    Calculate the Motion Tree dissimilarity between two
    clusters.

    The N_LARGEST largest inter-cluster Dmn values are averaged.
    """

    n_cross_pairs = (
        len(cluster_a)
        *
        len(cluster_b)
    )
# for the case where cluster scores are less than the required, the mean of the entire list is taken instead.
# Step 4
    n_required = min(
        N_LARGEST,
        n_cross_pairs
    )

    values = np.empty(
        n_cross_pairs,
        dtype=np.float64
    )

    index = 0

    for a in cluster_a:

        for b in cluster_b:

            values[index] = get_D_value(
                a,
                b
            )

            index += 1



    largest = np.partition(
        values,
        -n_required
    )[-n_required:]

    score = np.mean(
        largest
    )

    return score


# step 1.

clusters = [
    [i]
    for i in range(n_atoms)
]


# ============================================================
#  CLUSTER IDs
# ============================================================

# Original atoms/residues:
#
#     0 ... n_atoms-1
#
# Newly created clusters:
#
#     n_atoms, n_atoms+1, ...

cluster_ids = list(
    range(n_atoms)
)

next_cluster_id = n_atoms


linkage_matrix = []


# ============================================================
# Merging Clusters
# ============================================================

print("\n")
print("=" * 70)
print("CONSTRUCTING MOTION TREE")
print("=" * 70)

merge_number = 0

#Step 7

while len(clusters) > 1:

    merge_number += 1

    best_pair = None

    best_score = np.inf

    n_current_clusters = len(
        clusters
    )

    # ========================================================
    # step 2-6
    # ========================================================

    for a in range(
        n_current_clusters
    ):

        for b in range(
            a + 1,
            n_current_clusters
        ):

            cluster_a = clusters[a]

            cluster_b = clusters[b]

            score = cluster_dissimilarity(
                cluster_a,
                cluster_b
            )

            if score < best_score:

                best_score = score

                best_pair = (
                    a,
                    b
                )


    # 
    # select clusters to merge
    # 

    a, b = best_pair

    cluster_a = clusters[a]

    cluster_b = clusters[b]

    id_a = cluster_ids[a]

    id_b = cluster_ids[b]

    # 
    # merge
    # 

    merged_cluster = (
        cluster_a
        +
        cluster_b
    )

    merged_size = len(
        merged_cluster
    )

    # 
    # residue information is used later for chimera
    # 

    merged_resids = [
        int(ca.resids[i])
        for i in merged_cluster
    ]

    # ========================================================
    # SAVE LINKAGE ROW
    # ========================================================

    linkage_matrix.append(
        [
            id_a,
            id_b,
            best_score,
            merged_size
        ]
    )

    # ========================================================
    # PRINT PROGRESS
    # ========================================================

    print(
        f"Merge "
        f"{merge_number:4d}/{n_atoms - 1}: "
        f"{len(cluster_a):4d} + "
        f"{len(cluster_b):4d}"
        f" -> "
        f"{merged_size:4d}"
        f" | clusters remaining: "
        f"{n_current_clusters - 1:4d}",
        flush=True
    )

    # ========================================================
    # REMOVE OLD CLUSTERS
    # ========================================================

    for index in sorted(
        [a, b],
        reverse=True
    ):

        del clusters[index]

        del cluster_ids[index]

    # ========================================================
    # ADD NEW CLUSTER
    # ========================================================

    clusters.append(
        merged_cluster
    )

    cluster_ids.append(
        next_cluster_id
    )

    next_cluster_id += 1



print("CLUSTERING COMPLETE")



# ============================================================
# linkage matrix
# ============================================================

linkage_matrix = np.asarray(
    linkage_matrix,
    dtype=np.float64
)

print(
    "\nLinkage matrix shape:",
    flush=True
)

print(
    linkage_matrix.shape,
    flush=True
)


print(
    f"\nSaving linkage matrix to "
    f"{LINKAGE_FILE}...",
    flush=True
)

np.save(
    LINKAGE_FILE,
    linkage_matrix
)

print(
    "Linkage matrix saved successfully.",
    flush=True
)


# ============================================================
# CREATE graph
# ============================================================



# Residue labels for the leaves.
labels = [
    str(resid)
    for resid in ca.resids
]


# ============================================================
# data for chimerax
# ============================================================

print(
    "\nCalculating leaves belonging to every internal node...",
    flush=True
)

# SciPy linkage convention:
#
#   leaf nodes:
#       0 ... n_atoms-1
#
#   internal nodes:
#       n_atoms ... 2*n_atoms-2
#
# For each internal node, store the original C-alpha indices
# underneath that node.
node_leaves = {}

for row_number, row in enumerate(linkage_matrix):

    child_a = int(row[0])
    child_b = int(row[1])

    node_id = n_atoms + row_number

    if child_a < n_atoms:
        leaves_a = [child_a]
    else:
        leaves_a = node_leaves[child_a]

    if child_b < n_atoms:
        leaves_b = [child_b]
    else:
        leaves_b = node_leaves[child_b]

    node_leaves[node_id] = leaves_a + leaves_b


# 
# SAVE NODE -> LEAF MAPPING
# 

print(
    f"\nSaving node -> leaf mapping to "
    f"{NODE_LEAVES_FILE}...",
    flush=True
)

with open(NODE_LEAVES_FILE, "w") as f:

    f.write("MOTION TREE INTERNAL NODE -> LEAF MAPPING\n")
    f.write("===========================================\n\n")

    f.write(
        f"Number of C-alpha atoms: {n_atoms}\n"
    )

    f.write(
        f"Number of internal nodes: {len(node_leaves)}\n\n"
    )

    for row_number, row in enumerate(linkage_matrix):

        node_id = n_atoms + row_number

        child_a = int(row[0])
        child_b = int(row[1])

        distance = float(row[2])

        atom_indices = node_leaves[node_id]

        residue_numbers = [
            int(ca.resids[i])
            for i in atom_indices
        ]

        f.write(
            f"Node {node_id}\n"
        )

        f.write(
            f"Distance: {distance:.6f} Å\n"
        )

        f.write(
            f"Children: {child_a}, {child_b}\n"
        )

        f.write(
            f"Number of leaves: {len(residue_numbers)}\n"
        )

        f.write(
            "Leaves: "
            +
            " ".join(
                map(str, residue_numbers)
            )
            +
            "\n\n"
        )

print(
    "Node -> leaf mapping saved successfully.",
    flush=True
)


# Figure

fig, ax = plt.subplots(

    figsize=(8, 6) 
    )



dendro = dendrogram(
    linkage_matrix,
    labels=labels,
    leaf_rotation=90,
    leaf_font_size=7,
    ax=ax
)





# ============================================================
# graph settiungs
# ============================================================



# Code below is specific to ADK, if another molecule is run, alternative values and names will be required to label each domain!
# Values below were taken directly from a previously run MT. These help with labelling each region.
""" ax.text(500, -0.02, "NMP",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top")

ax.text(150, -0.02, "ATP lid",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top")

ax.text(1250, -0.02, "Core domain",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top") """
        
ax.set_xlabel("Residue") 

ax.set_ylabel(
    "MT Score (Å)"
)
plt.ylim(0,6)

ax.set_title(
    "Motion Tree - Free ADK"
)
ax.set_xticks([]) 


plt.tight_layout()


# ============================================================
# SAVE DENDROGRAM
# ============================================================

print(
    f"\nSaving Motion Tree to "
    f"{FIGURE_FILE}...",
    flush=True
)

plt.savefig(
    FIGURE_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close(
    fig
)

print(
    "Motion Tree figure saved successfully.",
    flush=True
)


# ============================================================
# FINAL OUTPUT
# ============================================================

print("\n")
print("=" * 70)
print("FILES GENERATED")
print("=" * 70)

print(
    f"Dmn values:        "
    f"{DISTANCE_FLUCTUATIONS_FILE}",
    flush=True
)

print(
    f"Linkage matrix:    "
    f"{LINKAGE_FILE}",
    flush=True
)

print(
    f"Motion Tree plot:  "
    f"{FIGURE_FILE}",
    flush=True
)


print(
    f"Node -> leaf map:  "
    f"{NODE_LEAVES_FILE}",
    flush=True
)

print(
    "\nMotion Tree analysis finished successfully.",
    flush=True
)

