import pandas as pd
import numpy as np
import sys
import glob
import os
import re
import Bio.PDB.PDBParser
import warnings
import math
warnings.filterwarnings("ignore", message="Used element '.' for Atom")

levels = ["class", "arch", "topo", "superfam"]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cath-filename", 
                    help="CATH domain list input file")
parser.add_argument("--pdb-dir", 
                    help="PDB directory")
parser.add_argument("--n-splits", default=10, type=int,
                    help="Number of splits (default: %(default)s)")
parser.add_argument("--atom-selector-regexp", default="CA",
                    help="Atom selector (default: %(default)s)")
parser.add_argument("--max-distance", default=50.0,
                    help="Maximum distance from atom to center of mass (default: %(default)s)")
parser.add_argument("--extract-at-level", default="arch",
                    help="Which CATH-level to use, i.e. class, arch, topo, or superfam (default: %(default)s)")
parser.add_argument("--sub-category-level", default="superfam",
                    help="Which CATH-level to use, i.e. class, arch, topo, or superfam (default: %(default)s)")
parser.add_argument("--min-size", default=500, type=int,
                    help="Minimum number of elements in category in order for it to be included (default: %(default)s)")
parser.add_argument("--min-resolution", default=3.5, type=float,
                    help="The minimum resolution for entries to be included (note that the resolution grows when this number drops) (default: %(default)s)")
parser.add_argument("--print-group-sizes-only", default=False, action="store_true",
                    help="Just print out the unfiltered group sizes - useful for deciding on a min-size value (default: %(default)s)")

args = parser.parse_args()

print("# Arguments")
for key, value in sorted(vars(args).items()):
    print(key, "=", value)

extract_at_level = levels.index(args.extract_at_level)+1

sub_category_col_range = list(range(1, levels.index(args.sub_category_level)+1+1))  # add one offset and one because end is excluded

# Read data into pandas dataframe
data = pd.read_csv(args.cath_filename, sep='\s+', header=None, usecols=[0,1,2,3,4,11], comment="#")

# Iterate over PDB files, and use IDs to filter dataframe
pdb_filenames = glob.glob(os.path.join(args.pdb_dir, "*"))
pdb_ids = [os.path.basename(name) for name in pdb_filenames]
data = data[data[0].isin(pdb_ids)]
data = data.reset_index(drop=True)

print ("Processing PDB files for distance to CM")

# Create Bio.PDB parser object
pdb_parser = Bio.PDB.PDBParser()

# Iterate over PDB files, and calculate distance from center of mass. Add as additional column
data['dist'] = 0
max_distance_to_cm = np.zeros(len(data))
# import pickle
# max_distance_to_cm = pickle.load(open('max_distance_to_cm.pickle', 'rb'))
for index, row in data.iterrows():
    # Extract CATH classification and PDB ID
    cath_id = row[0]
    pdb_filename = os.path.join(args.pdb_dir, cath_id)

    print("%d/%d:" % (index,len(data)), pdb_filename)

    # Parse structure
    structure = pdb_parser.get_structure(pdb_filename, pdb_filename)

    positions_all_atoms = []
    # Retrieve all atom coordinates
    for atom in structure.get_atoms():
        # The all-atom selection is used to calculate distance from center of mass
        positions_all_atoms.append(atom.get_coord())

    # Translate to center of mass
    positions_all_atoms = np.array(positions_all_atoms)
    positions_all_atoms = positions_all_atoms - np.mean(positions_all_atoms, axis=0)

    max_distance_to_cm[index] = np.max(np.linalg.norm(positions_all_atoms, axis=1))
data = data.assign(dist=max_distance_to_cm)
    

# Group dataframe by [class, architecture] levels
minimum_group_len = None
for name, group in data.groupby(list(range(1,extract_at_level+1))):

    # Select elements with at least min_arch_size structures with resolution at least min_resolution (note higher resolution is smaller number)
    # and with specified maximum distance to CM.
    if len(group[np.logical_and(group[11] < args.min_resolution, group['dist'] < args.max_distance)]) > args.min_size:

        # Print entry and size
        if args.print_group_sizes_only:
            print(name, len(group))
            continue

        # Group by all subcategories (i.e., to superfamily level)
        sub_category_groups = group[np.logical_and(group[11] < args.min_resolution, group['dist'] < args.max_distance)].groupby(sub_category_col_range)

        # Calculate group length by summing subgroup lengths
        group_len = np.sum([len(g[1]) for g in sub_category_groups])

        # Update minimum group length
        minimum_group_len = group_len if minimum_group_len is None else min(minimum_group_len, group_len)
                                
groups_reduced = []
for name, group in data.groupby(list(range(1,extract_at_level+1))):

    # Select elements with at least min_arch_size structures with resolution at least min_resolution (note higher resolution is smaller number)
    # and with specified maximum distance to CM.
    if len(group[np.logical_and(group[11] < args.min_resolution, group['dist'] < args.max_distance)]) > args.min_size:

        # Reduce to the minimum number
        n_entries = minimum_group_len
        
        # Reduce to minimum number of elements, but attempting to select evenly from superfamilies
        group_reduced = []

        # Group by all subcategories (i.e., to superfamily level)
        sub_category_groups = group[np.logical_and(group[11] < args.min_resolution, group['dist'] < args.max_distance)].groupby(sub_category_col_range)

        # Skip if there is only one subcategory (this violates the constraint that all splits should have all groups represented)
        if len(sub_category_groups) < args.n_splits:
            continue
        
        # Keep track of numbers of entries added so far
        n_added_entries = 0

        # print(name)
        
        # We now reduce to the limit. Rather than taking arbitrary members, we try to sample as uniformly
        # as possible among the members of the subcategory level
        # Iterate over sub_category groups - sorted by length - smallest groups first
        for i, group_inner_pair in enumerate(sorted(sub_category_groups, key=lambda k:len(k[1]))):
            name_inner, group_inner = group_inner_pair

            # Calculate how much we are allowed to include. This is simply a matter of spreading
            # what remains evenly over the remaining iterations
            inclusion_size = int(math.ceil((n_entries - n_added_entries) / (len(sub_category_groups) - i)))

            included_entries = group_inner.sort_values([11])[:inclusion_size] 
            group_reduced.append(included_entries)
            n_added_entries += len(included_entries)

            # print(name_inner, inclusion_size, len(group_inner.sort_values([11])))
            # print("\t", name_inner)
            # print("\t", len(included_entries), n_added_entries, inclusion_size, len(group_inner))
            
        # Finally, append added entries to form new dataframe
        groups_reduced.append(pd.concat(group_reduced))

        # print(n_added_entries, n_entries, minimum_group_len)
        assert(n_added_entries == n_entries)

if args.print_group_sizes_only:
    sys.exit()
    
# Merge list of groups into single data frame
data_reduced = pd.concat(groups_reduced)

# Reset index
data_reduced = data_reduced.reset_index(drop=True)


# Create splits by iterating over all sub-categories in sorted order (largest
# first), and adding each sub-category to the split which currently has fewest
# elements. In addition, we keep track of balancing the corresponding main
# categories, so that a split can only receive a sub-category if it does not
# already contain another sub-category from the same category. This counter is
# reset every time all splits have sub-categories from all categories.
# This is the reason that the total number of elements in each split is not the
# same over all splits

# Create splits
splits = [[[],{}] for i in range(args.n_splits)]

# Collect sub_categories across all categories
sub_categories = []
n_categories = len(list(data_reduced.groupby(list(range(1,extract_at_level+1)))))
for name, group in data_reduced.groupby(list(range(1,extract_at_level+1))):
    for name_inner, group_inner in group.groupby(sub_category_col_range):
        sub_categories.append((name, group_inner))
# Sort them by number of elements
sub_categories.sort(key=lambda k: len(k[1]), reverse=True)

print("n_categories: ", n_categories)

# Iterate over sorted sub_categories list
for j, pair in enumerate(sub_categories):
    category_id, sub_category = pair

    # Iterate over splits until split is found in which the
    # category corresponding to this sub-category is not yet present
    # if no suitable entry is found, use the first (smallest)
    split_index = 0
    for i, split_pair in enumerate(splits):
        split, split_category_ids = split_pair
        if category_id not in split_category_ids:
            split_index = i
            break
    
    # Add indices as first element to the splits array
    split, split_category_ids = splits[split_index]
    split += sub_category.index.values.tolist()

    # Register category in split
    split_category_ids[category_id] = True

    # Reset if all splits have seen all categories
    splits_all_complete = True
    for split, split_category_ids in splits:
        splits_all_complete = splits_all_complete and len(split_category_ids) == n_categories
    if splits_all_complete:
        for i, _ in enumerate(splits):
            splits[i][1] = {}
    
    # Sort array so that smallest entrt appears first
    splits.sort(key=lambda k:len(k[0]))

    # print([len(v[0]) for v in splits])
    # print([(len(v[0]),list(v[1].keys())) for v in splits])

# Throw out category id from splits
splits = [v[0] for v in splits]
       
        
print("Split distribution: ", [len(split) for split in splits])
        
# Reorder entries in data frame so that they follow the split division
# (this allows us to keep track of each split using only a start index)
data_reduced_reordered = []
split_start_indices = []
for i, split in enumerate(splits):
    data_reduced_reordered.append(data_reduced.iloc[split])
    if i==0:
        split_start_indices.append(0)
    else:
        split_start_indices.append(split_start_indices[-1]+len(splits[i-1]))
data_reduced = pd.concat(data_reduced_reordered)
data_reduced = data_reduced.reset_index(drop=True)

# Check that members of a sub_category always end in the same split and that all splits contain all main categories
split_ids = {}
print("Split start indices: ", split_start_indices)
for split_index, split_pair in enumerate(zip(split_start_indices[0:], split_start_indices[1:]+[None])):

    sub_category_ids = {}

    # Check that sub_categories always end up in the same split
    split_start, split_end = split_pair
    for index, row in data_reduced.iloc[split_start:split_end].iterrows():
        cath_id = row[1], row[2], row[3], row[4]
        if cath_id not in split_ids:
            split_ids[cath_id] = split_index
        else:
            assert split_ids[cath_id] == split_index

    # Check that each split contains all main categories
    assert len(data_reduced.iloc[split_start:split_end].groupby(list(range(1,extract_at_level+1)))) == n_categories
                    
print ("Processing PDB files...")

positions = []
n_atoms = []
atom_types = []
res_indices = []
labels = []
atom_selector_regexp = re.compile(args.atom_selector_regexp)
for index, row in data_reduced.iterrows():

    # Extract CATH classification and PDB ID
    cath_id = row[0]
    cath_classification = row[1], row[2], row[3], row[4]
    pdb_filename = os.path.join(args.pdb_dir, cath_id)

    print(index, pdb_filename)

    # Parse structure
    structure = pdb_parser.get_structure(pdb_filename, pdb_filename)

    positions_all_atoms = []
    positions_tmp = []
    atom_types_tmp = []
    res_indices_tmp = []
    # Retrieve all atom coordinates
    for atom in structure.get_atoms():
        # The all-atom selection is used to calculate distance from center of mass
        positions_all_atoms.append(atom.get_coord())

        # filter with atom selection
        match = atom_selector_regexp.match(atom.id)
        if match and len(match.group(0)) == len(atom.id):
            positions_tmp.append(atom.get_coord())
            atom_types_tmp.append(match.group(0))
            # assert(match.group(0) == "CA")
            res_indices_tmp.append(int(atom.get_parent().id[1]))

    # Translate to center of mass
    positions_tmp = np.array(positions_tmp)
    positions_tmp = positions_tmp - np.mean(positions_tmp, axis=0)
    positions_all_atoms = np.array(positions_all_atoms)
    positions_all_atoms = positions_all_atoms - np.mean(positions_all_atoms, axis=0)

    # Check that all PDBs have max distance to center of mass within limit
    assert np.max(np.linalg.norm(positions_all_atoms, axis=1)) < args.max_distance
    
    positions.append(positions_tmp)
    n_atoms.append(len(positions_tmp))
    atom_types.append(atom_types_tmp)
    res_indices.append(res_indices_tmp)
    labels.append(cath_classification[:extract_at_level])

# Translate positions to numpy array, by finding maximum number of elements
max_n_atoms = max([len(pos) for pos in positions])
positions_array = np.zeros([len(positions), max_n_atoms, 3])
atom_types_array = np.empty([len(positions), max_n_atoms], dtype='S10')
atom_types_array[:] = ""
res_indices_array = np.zeros([len(positions), max_n_atoms])-1
for i, pos in enumerate(positions):
    positions_array[i, :n_atoms[i]] = pos
    atom_types_array[i, :n_atoms[i]] = atom_types[i]
    res_indices_array[i, :n_atoms[i]] = res_indices[i]

# Save features
print(len(set(labels)))
np.savez_compressed("cath_%d%s_%s"%(len(set(labels)), args.extract_at_level, args.atom_selector_regexp.replace("|","")),
                    n_atoms=np.array(n_atoms),
                    atom_types=atom_types_array,
                    res_indices=res_indices_array,
                    positions=positions_array,
                    labels=np.array(labels),
                    split_start_indices=np.array(split_start_indices))



    




# print(len(data_reduced))

