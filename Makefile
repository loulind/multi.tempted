.PHONY: clean dirs

# cleans project
clean:
	rm -rf derived_data
	rm -rf figures
	rm -rf interactives
	rm -f report.html

# creates folders for project
dirs:
	mkdir -p derived_data
	mkdir -p figures
	mkdir -p interactives

# creates embeddings vector and vector of paragraphs
derived_data/embeddings.csv derived_data/paragraphs.csv:\
 raw_data/sound_and_fury.txt embeddings.R | dirs
	Rscript embeddings.R

# dimensionality-reduced data
tsne.csv umap.csv: dim_reduce.R\
 derived_data/embeddings.csv\
 derived_data/paragraphs.csv
	Rscript dim_reduce.R

# figures and interactive visuals
figures/tsne_narr_order.png\
 figures/umap_narr_order.png\
 figures/tsne_par_order.png\
 figures/umap_par_order.png\
 interactives/tsne3d.html\
 interactives/umap3d.html\
 interactives/animated.html: tsne.csv umap.csv figures.R
	Rscript figures.R

# report
report.html: figures/tsne_narr_order.png\
 figures/tsne_par_order.png\
 figures/umap_narr_order.png\
 figures/umap_par_order.png\
 report.Rmd | dirs
	R -e "rmarkdown::render('report.Rmd', output_file='report.html')"
