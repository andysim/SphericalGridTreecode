import numpy as np
import scipy.special
import glob
import time
import itertools

# This is a very simple proof-of-concept code, designed to provide a minimal working example
# for the Spherical Grid Treecode trick for evaluating electrostatics.  This version makes
# a list of remote nodes for each leaf node, avoiding any downward pass.

def first_node(level):
    return (8**level - 1) // 7

class TreeNode(object):
    # A simple placeholder object, representing nodes in the tree
    def __init__(self, depth, nvecs ):    
        self.remote_boxes = { level:[] for level in range(depth+1) }
        self.direct_boxes = []
        self.parent_address = -1
        self.charges = []
        self.coords = []
        self.sin_vecs = { level : np.zeros((nvecs)) for level in range(depth+1) }
        self.cos_vecs = { level : np.zeros((nvecs)) for level in range(depth+1) }


def get_gauss_hermite_weights(tol, rmax, ω):
    # Obtain Gauss-Hermite quadrature nodes and weights, trimming the small values according
    # to the estimates from Limpanuparb, et al. in dx.doi.org/10.1021/ct301110y, eqns. (14) - (18).
    ωR = ω*rmax
    estN = ωR**2 + 2*ωR*(np.sqrt(-np.log10(tol))-1) + 2
    N = int(np.ceil(estN))
    num_gh = 2 * N
    x, w = scipy.special.roots_hermite(num_gh)
    # Pick out only the unfiltered weights
    estN_p = (2/np.pi)*np.sqrt(-(N+1)*np.log(tol)) - 1
    N_p = min(int(np.ceil(estN_p)), N)
    istart = N - N_p
    return -x[istart:N], w[istart:N]


def estimate_t_design_degree(tol, ω, xvals, wvals, rmax):
    # Obtain Gauss-Hermite quadrature nodes and weights, trimming the small values according
    # to the estimates from Limpanuparb, et al. in dx.doi.org/10.1021/ct301110y, eqns. (14) - (18).
    Nprime = len(xvals)
    rhs = tol / Nprime
    degree = 0
    while True:
        max_lhs = 0.0
        for wval,xval in zip(wvals, xvals):
            q_n = 4 * np.sqrt(ω * xval)
            λ_n = 2 * ω * xval
            j_L = scipy.special.spherical_jn(degree, λ_n*rmax)
            lhs = (q_n * j_L)**2 / (4 * np.pi)
            max_lhs = max(max_lhs, lhs)
        if max_lhs < rhs:
            break
        degree += 1
    # round up to odd number
    return degree + 1 - degree % 2

def load_t_design(degree):
    matches = []
    for filename in glob.glob("SS31-Mar-2016/*"):
        if f"ss{degree:03}" in filename:
            matches.append(filename)
    assert len(matches) == 1
    full_Sm = np.loadtxt(matches[0])
    full_dim = full_Sm.shape[0]
    # Keep only the upper hemisphere, and double results later to account for symmetry
    return full_Sm[:full_dim//2,:]

def get_cubature_grid(tol, rmax, ω):
    # Find the radial (Gauss-Hermite) quadrature info first
    x,w = get_gauss_hermite_weights(tol, rmax, ω)
    num_gh = len(x)
    # Find the appropriate spherical t-design degree
    degree = estimate_t_design_degree(tol, ω, x, w, rmax)
    Sm = load_t_design(degree)
    Sm_dim = Sm.shape[0]
    weight_vec = np.repeat(w, Sm_dim) / Sm_dim
    λ = 2 * ω * np.repeat(x, Sm_dim)
    a = λ * np.tile(Sm[:,0], num_gh)
    b = λ * np.tile(Sm[:,1], num_gh)
    c = λ * np.tile(Sm[:,2], num_gh)
    return np.array(weight_vec), a, b, c, num_gh, Sm_dim


def node_address(level, x, y, z):
    dim = 1 << level
    return first_node(level) + dim*dim * x + dim * y + z


def direct_energy_full(chargesA, coordsA, chargesB, coordsB):
    num_atomsA = len(chargesA)
    num_atomsB = len(chargesB)
    assert len(coordsA) == num_atomsA
    assert len(coordsB) == num_atomsB
    energy = 0.0
    for A in range(num_atomsA):
        for B in range(num_atomsB):
            energy += chargesA[A] * chargesB[B] / np.linalg.norm(np.subtract(coordsA[A], coordsB[B]))
    return energy


def direct_energy(charges, coords):
    num_atoms = len(charges)
    assert len(coords) == num_atoms
    energy = 0.0
    for A in range(num_atoms):
        for B in range(A):
            energy += charges[A] * charges[B] / np.linalg.norm(np.subtract(coords[A], coords[B]))
    return energy


def promote_nodes_if_all_adjacent(neighbor_list, depth):
    # This simply looks for the situation where all 8 octants are evaluated on the same (or higher)
    # level, and groups them all together by removing them from the lower, promoting them up instead.
    for level in range(depth, 1, -1):
        leaves_per_dimension = 1 << level
        for this_x in range(0, leaves_per_dimension, 2):
            for this_y in range(0, leaves_per_dimension, 2):
                for this_z in range(0, leaves_per_dimension, 2):
                    this_octant_addresses = [
                        node_address(level, this_x, this_y, this_z),
                        node_address(level, this_x, this_y, this_z+1),
                        node_address(level, this_x, this_y+1, this_z),
                        node_address(level, this_x+1, this_y, this_z),
                        node_address(level, this_x, this_y+1, this_z+1),
                        node_address(level, this_x+1, this_y, this_z+1),
                        node_address(level, this_x+1, this_y+1, this_z),
                        node_address(level, this_x+1, this_y+1, this_z+1)
                    ]
                    if np.count_nonzero(np.isin(neighbor_list, this_octant_addresses)) == 8:
                        # All 8 boxes are in the list; replace them with a single box one level up
                        for element in this_octant_addresses:
                            neighbor_list.remove(element)
                        neighbor_list.append(node_address(level-1, this_x//2, this_y//2, this_z//2))


def run_test(box_width, depth, atoms_per_dimension):
    rng = np.random.default_rng()#123)
    ZERO = 1e-7

    tol = 1e-5
    print_level = 3
    num_atoms = atoms_per_dimension**3
    leaves_per_dimension = 2**depth
    leaf_width = box_width / leaves_per_dimension

    # erfc(ɑ cutoff)
    # -------------- = ZERO
    #     cutoff
    α = scipy.special.erfcinv(leaf_width * ZERO) / leaf_width

    # Add some padding to account for partially intersected remote nodes
    rmax = 2.5 * leaf_width

    weights, pts_x, pts_y, pts_z, num_gh, num_Sm = get_cubature_grid(tol, rmax*(2**depth), α/(2**depth))
    nvecs = len(weights)

    # Start with a regular cube of atoms, but apply a small random perturbation to break symmetry
    coords = np.zeros((num_atoms, 3))
    atom_spacing = box_width / atoms_per_dimension
    coords1d  = np.linspace(atom_spacing/2, box_width-atom_spacing/2, atoms_per_dimension)
    for n,regular_coord in enumerate(itertools.product(coords1d, coords1d, coords1d)):
        # Take the regular grid of positions, and add a small random shift
        x,y,z = regular_coord + 2*rng.random(3)-1
        # Make sure everything is within the bounding box
        x = x if x < box_width else box_width - 1e-5
        y = y if y < box_width else box_width - 1e-5
        z = z if z < box_width else box_width - 1e-5
        x = x if x > 0.0 else 1e-5
        y = y if y > 0.0 else 1e-5
        z = z if z > 0.0 else 1e-5
        coords[n,:] = x, y, z

    # Randomize charges, but make it neutral overall
    charges = 2*rng.random((num_atoms))
    charges -= np.mean(charges)

    # Offset into the node array, indexing the leaves
    leaf_beg = first_node(depth)
    leaf_end = leaf_beg + 8**depth

    tree = [TreeNode(depth, nvecs) for i in range(leaf_end)]
    A = {}

    print(f"num atoms: {num_atoms}")
    print(f"Total quadrature points {nvecs}: {num_gh} radial, {num_Sm} angular")
    print(f"box size = {box_width:.2f}, leaf width = {leaf_width:.2f}")
    print(f"ewald_ɑ = {α:.3f}, test of cutoff = {scipy.special.erfc(leaf_width * α)/leaf_width:.2e}")
    print(f"total number of nodes = {len(tree):d}")

    #
    # These tasks are geometry-independent, and their cost is amortized
    #
    if print_level > 2: print("Figuring out the parent list... ", end="", flush=True)
    for level in range(1, depth+1):
        # Figure out each node's parent.  This is only done once, so we use a terrible
        # algorithm here for simplicity.
        nodes_in_this_dimension = 2**level
        for this_x in range(nodes_in_this_dimension):
            px = this_x//2
            for this_y in range(nodes_in_this_dimension):
                py = this_y//2
                for this_z in range(nodes_in_this_dimension):
                    pz = this_z//2
                    this_node_address = node_address(level, this_x, this_y, this_z)
                    tree[this_node_address].parent_address = node_address(level-1, px, py, pz)
    if print_level > 2: print("Done!")

    if print_level > 2: print("Figuring out the neighbor list...", end="", flush=True)
    # Find the neighbor list for each leaf node
    for this_x in range(leaves_per_dimension):
        for this_y in range(leaves_per_dimension):
            for this_z in range(leaves_per_dimension):
                this_node_address = node_address(depth, this_x, this_y, this_z)
                for neighbor_x in range(leaves_per_dimension):
                    dx = max(abs(this_x - neighbor_x) - 1, 0)
                    for neighbor_y in range(leaves_per_dimension):
                        dy = max(abs(this_y - neighbor_y) - 1, 0)
                        for neighbor_z in range(leaves_per_dimension):
                            dz = max(abs(this_z - neighbor_z) - 1, 0)
                            neighbor_node_address = node_address(depth, neighbor_x, neighbor_y, neighbor_z)
                            if neighbor_node_address == this_node_address:
                                continue
                            distance = np.sqrt(dx**2 + dy**2 + dz**2)  # in units of box width
                            if distance < 1.0:
                                # This is a direct box, which is handled in real space by applying 1/R
                                tree[this_node_address].direct_boxes.append(neighbor_node_address)
                            else:
                                # This is a reciprocal box
                                # Find out which "level" of the tree this pair of boxes interacts on, i.e., which
                                scale = depth - int(np.floor(np.log2(distance)))
                                tree[this_node_address].remote_boxes[scale].append(neighbor_node_address)
    if print_level > 2: print("Done!")

    if print_level > 2: print("Optimizing tree...", end="", flush=True)
    # A slow algorithm to find groupings of adjacent octants that can be consolidated
    # and interacted one level higher on the tree
    for scale in range(depth + 1):
        for leaf in tree[leaf_beg:leaf_end]:
            promote_nodes_if_all_adjacent(leaf.remote_boxes[scale], depth)
    if print_level > 2: print("Done!")

    #
    # The tasks below are needed for each change of geometry
    #

    t_fac_start = time.time()

    # Place the charges in the right leaf
    if print_level > 2: print("Placing atoms into the correct leaf... ", end="", flush=True)
    for (atom,q) in enumerate(charges):
        boxes = np.array(np.floor(np.divide(coords[atom], leaf_width)), dtype=np.int64)
        bin_addr = node_address(depth, *boxes)
        tree[bin_addr].coords.append(coords[atom])
        tree[bin_addr].charges.append(q)
    if print_level > 2: print("Done!")


    if print_level > 2: print("Building the structure factors... ", end="", flush=True)
    old_cos_vec = np.zeros(nvecs)
    old_sin_vec = np.zeros(nvecs)
    new_cos_vec = np.zeros(nvecs)
    new_sin_vec = np.zeros(nvecs)
    for leafnum in range(leaf_beg, leaf_end):
        for scale in range(depth+1):
            tree[leafnum].cos_vecs[scale][:] = 0
            tree[leafnum].sin_vecs[scale][:] = 0
        for atom, q in enumerate(tree[leafnum].charges):
            Rxyz = tree[leafnum].coords[atom]
            cosθ = Rxyz[0] * pts_x + Rxyz[1] * pts_y + Rxyz[2] * pts_z
            new_cos_vec = np.cos(cosθ)
            new_sin_vec = np.sin(cosθ)
            old_cos_vec[:] = new_cos_vec
            old_sin_vec[:] = new_sin_vec
            tree[leafnum].cos_vecs[0][:] += q * new_cos_vec
            tree[leafnum].sin_vecs[0][:] += q * new_sin_vec
            # Now the scale 0 vectors are done, propagate down the tree
            for scale in range(1, depth+1):
                new_cos_vec[:] = 1 - 2 * np.square(old_sin_vec)
                new_sin_vec[:] = 2 * old_sin_vec * old_cos_vec
                old_cos_vec[:] = new_cos_vec
                old_sin_vec[:] = new_sin_vec
                tree[leafnum].cos_vecs[scale][:] += q * new_cos_vec
                tree[leafnum].sin_vecs[scale][:] += q * new_sin_vec
    # Fold the scale factors into the cos/sin vectors to make the energy a simple dot product.
    scalefac = np.sqrt((4*(2**-depth)*α/np.pi) * weights)
    for scale in range(0,depth+1):
        for leafnum in range(leaf_beg, leaf_end):
            tree[leafnum].cos_vecs[scale] *= scalefac
            tree[leafnum].sin_vecs[scale] *= scalefac
        scalefac *= np.sqrt(2)
    if print_level > 2: print("Done!")

    # Now iterate backwards from the end of the leaves, contributing partial structure factors to
    # every node's parent.  At the end, the top level node should have the complete sum
    if print_level > 2: print("Accumulating structure factors... ", end="", flush=True)
    for leafnum in range(leaf_end-1, 0, -1):
        parent = tree[leafnum].parent_address
        for scale in range(0, depth+1):
            tree[parent].cos_vecs[scale][:] += tree[leafnum].cos_vecs[scale][:] 
            tree[parent].sin_vecs[scale][:] += tree[leafnum].sin_vecs[scale][:] 
    if print_level > 2: print("Done!")

    if print_level > 2: print("Computing energy using decomposition...", end="", flush=True)
    Erec = 0.0
    remote_cos_sum = np.zeros((nvecs))
    remote_sin_sum = np.zeros((nvecs))
    for scale in range(depth+1):
        for leafnum in range(leaf_beg, leaf_end):
            remote_cos_sum[...] = 0.0
            remote_sin_sum[...] = 0.0
            for remote_box in tree[leafnum].remote_boxes[scale]:
                remote_cos_sum += tree[remote_box].cos_vecs[scale]
                remote_sin_sum += tree[remote_box].sin_vecs[scale]
            Erec += 0.5 * (tree[leafnum].cos_vecs[scale] @ remote_cos_sum +
                           tree[leafnum].sin_vecs[scale] @ remote_sin_sum)

    Edir = 0.0
    for leaf in tree[leaf_beg:leaf_end]:
        # The direct contribution within each leaf
        Edir += direct_energy(leaf.charges, leaf.coords)
        # The direct contribution between the homebox and its nearest neighbor
        for remote_box in leaf.direct_boxes:
            Edir += 0.5 * direct_energy_full(leaf.charges, leaf.coords, tree[remote_box].charges, tree[remote_box].coords)
    if print_level > 2: print("Done!")

    t_fac_stop = time.time()

    t_reg_start = time.time()
    if print_level > 2: print("Computing direct energy (N^2 algorithm)...", end="", flush=True)
    Ereg = direct_energy(charges, coords)
    if print_level > 2: print("Done!")
    t_reg_stop = time.time()

    Etot = Edir + Erec
    print(f"Energy (Dir) : {Edir:16.12f}")
    print(f"Energy (Rec) : {Erec:16.12f}")
    print(f"Total Ewald Energy :  {Etot:16.12f}")

    print(f"Energy (conventional) {Ereg:16.12f}\n")
    print(f"Relative Error = {(Etot - Ereg)/Ereg:16.12f}\n")

    t_fac = t_fac_stop - t_fac_start
    t_reg = t_reg_stop - t_reg_start
    print(f"Time for factorized algorithm: {t_fac:6.2f}s")
    print(f"Time for regular algorithm:    {t_reg:6.2f}s")

# Run some random boxes
run_test(box_width=20, depth=2, atoms_per_dimension=8)
run_test(box_width=40, depth=3, atoms_per_dimension=16)
