data:
ifeq ("$(wildcard cath_data)","")
	mkdir cath_data
	cd cath_data; wget https://sid.erda.dk/share_redirect/Adk8hxLMKZ/cath-domain-list-v4_2_0.txt
	cd cath_data; wget https://sid.erda.dk/share_redirect/Adk8hxLMKZ/cath-dataset-nonredundant-S40-v4_2_0.pdb.tgz
	cd cath_data; tar xvf cath-dataset-nonredundant-S40-v4_2_0.pdb.tgz
endif

cath_3class_ca: data
	python construct_dataset.py --cath-filename cath_data/cath-domain-list-v4_2_0.txt --pdb-dir cath_data/dompdb --min-size 2500 --extract-at-level class --sub-category-level superfam --n-splits 10 --atom-selector-regexp "CA" --min-resolution 3.5

cath_10arch_ca: data
	python construct_dataset.py --cath-filename cath_data/cath-domain-list-v4_2_0.txt --pdb-dir cath_data/dompdb --min-size 700 --extract-at-level arch --sub-category-level superfam --n-splits 10 --atom-selector-regexp "CA" --min-resolution 3.5

cath_20topo_ca: data
	python construct_dataset.py --cath-filename cath_data/cath-domain-list-v4_2_0.txt --pdb-dir cath_data/dompdb --min-size 200 --extract-at-level topo --sub-category-level superfam --n-splits 10 --atom-selector-regexp "CA" --min-resolution 3.5

