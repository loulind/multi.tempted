.PHONY: clean dirs

# cleans project
clean:
	rm -rf output

# creates folders for project
dirs:
	mkdir -p output

# Creates outputs
output/output.html: data/ipop.rda README.Rmd | dirs
	Rscript -e 'library(devtools);\
	load_all();\
	rmarkdown::render("README.Rmd", output_file="output/output.html", output_format="html_document")'
