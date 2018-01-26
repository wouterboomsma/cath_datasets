data:
ifeq ("$(wildcard cath_data)","")
	mkdir cath_data
	cd cath_data; wget https://sid.erda.dk/share_redirect/Adk8hxLMKZ/cath-domain-list-v4_2_0.txt
	cd cath_data; wget https://sid.erda.dk/share_redirect/Adk8hxLMKZ/cath-dataset-nonredundant-S40-v4_2_0.pdb.tgz
	cd cath_data; tar xvf cath-dataset-nonredundant-S40-v4_2_0.pdb.tgz
endif

cath_3class: data
	python construct_dataset.py --cath-filename cath_data/cath-domain-list-v4_2_0.txt --pdb-dir cath_data/dompdb --min-size 2500 --extract-at-level class --sub-category-level topo --n-splits 10

cath_10arch: data
	python construct_dataset.py --cath-filename cath_data/cath-domain-list-v4_2_0.txt --pdb-dir cath_data/dompdb --min-size 500 --extract-at-level arch --sub-category-level topo --n-splits 10
