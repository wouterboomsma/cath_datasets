# CATH protein structure classification data sets 

This repository contains various classification data sets built from the CATH Protein Structure Classification Database (www.cathdb.info). We have processed the data to be useful for the purpose of benchmarking methods for classifying shapes in 3D space. For more information about the hierarchy, see http://cathdb.info/browse/tree.

## cath_3class_ca
The `cath_3class.npz` dataset is the simplest set. It considers the "class" level of the CATH hierarchy, which in the CATH database consists of "Mainly Alpha", "Mainly Beta", "Alpha Beta" and "Few secondary structures". Since the latter category is small, and structurally heterogeneous, we omit it from our set. The three remaining categories are each reduced to have the same number of members (see filtering below), thus creating a balanced set. The three classes differ mainly in the relative quantities of alpha helices and beta strands (protein secondary structure). The main task in this set is thus to detect protein secondary structure in any orientation, and quantify the total amount of the different secondary structure elements in the entire image. The dataset contains only Carbon-alpha positions for each protein (i.e. only a single atom for each amino acid).

## cath_10arch_ca
The `cath_10arch.npz` dataset considers the "architecture" level in the CATH hierarchy. We limit ourselves to architectures with at least 700 members, which leaves us with 10 categories. We then reduce all of these to have the same size (see filtering below), thus creating a balanced set. The differences between the 10 architectures is more subtle than in the 3class case, with considerations not only about the number of secondary structure elements, but also their relative orientations. For an untrained eye, it can be difficult to visually distinguish the architectures. See http://cathdb.info/browse/tree for image examples of the structural variety occurring within architectures. The dataset contains only Carbon-alpha positions for each protein (i.e. only a single atom for each amino acid).

## cath_20topo_ca
The `cath_20_topo.npz` dataset considers the "topology/fold" level in the CATH hierarchy. We limit ourselves to topologies with at least 200 members, which leaves us with 20 categories. We then reduce all of these to have the same size (see filtering below), thus creating a balanced set. A topology/fold uses information not only about the secondary structure and its orientation, but also how the secondary structure elements are connected. See http://cathdb.info/browse/tree for image examples of the structural variety occurring within topologies/folds. The dataset contains only Carbon-alpha positions for each protein (i.e. only a single atom for each amino acid).


## Filtering
The sets are based on a 40% homology-reduced set of PDB structures downloaded from the CATH server. We then filter to obtain balanced categories (equal number of members) at the hierarchy level of interest. For instance, for the `3class` set, there are about 2500 data points for each class. In reducing the size, we use the structures with highest resolution (best experimental quality), with the additional constraint that all included structures with within a 50Å sphere centered around the center of mass of the protein. This last constraint was introduced to allow us to represent all proteins within a well-defined grid size. The constraint is only violated by a small fraction of the original data set. 

## Splits
We provide a 10-fold split of the data, for purposes of separation into train/validation/test sets or cross validation. In addition to the 40% sequence identity cut-off between any entry in the dataset, any two members from different splits are guaranteed to originate from different categories at the "superfamily" level in the CATH hierarchy. In addition, all splits are guaranteed to have members from all categories at the level you are classifying with respect to (i.e., the `class3` set has all 3 classes present in all splits). Note that due to the requirements of non-overlaps at the topology level, the splits are not always entirely of equal size.


## Format and usage
The datasets are in the numpy npz format, and can be read using `np.load()`:

```python
>>> import numpy as np
>>> data = np.load('cath_3class.npz')
```
This object contains a number of keys:
```python
>>> data.keys()
['n_atoms', 'atom_types', 'res_indices', 'positions', 'labels', 'split_start_indices']
```
The primary data is contained in `positions`. For convenience, these are provided as a single numpy array, although the number of atoms generally differs between proteins. The `n_atoms` array is used to extract the relevant positions. For instance, the positions for the first protein are accessible using:
```python
>>> data['positions'][0][:data['n_atoms'][0]]
```
The `data['labels']` entry contains the CATH labels which will generally be used as the classification targets. The arrays have been ordered to that the 10 splits appear sequentially, which means that the start index for each split is sufficient - provided in `data['split_start_indices]`. The `atom_types` key is associated with the atom type for each position and 'res_indices' with the sequential position in the protein chain.

## Licence
This work is derived from the CATH dataset (http://cathdb.info), which was based on work at http://cath.biochem.ucl.ac.uk. It is licenced under the same license as CATH, namely the Creative Commons Attribution 4.0 International License. The code in the scripts directory is licenced under an MIT licence.