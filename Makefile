.PHONY: clean dirs

# cleans project
clean:
	rm -rf output

# creates folders for project
dirs:
	mkdir -p output

# Locate pandoc for rmarkdown: use one on PATH if present, else fall back to the
# copy bundled with RStudio so `make` also works from a plain terminal (where
# RSTUDIO_PANDOC is not set). arm64 Macs name the tools dir "aarch64".
PANDOC_ARCH := $(shell uname -m | sed 's/arm64/aarch64/')
RSTUDIO_PANDOC ?= $(shell command -v pandoc >/dev/null 2>&1 \
	&& dirname "`command -v pandoc`" \
	|| echo /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/$(PANDOC_ARCH))
export RSTUDIO_PANDOC

# Renders the README the way a user would run it: install the package, then
# render README.Rmd (which attaches it with library(multi.tempted)).
output/output.html: data/ipop.rda README.Rmd | dirs
	Rscript -e 'devtools::install(quick = TRUE, upgrade = FALSE, quiet = TRUE);\
	rmarkdown::render("README.Rmd", output_file="output/output.html", output_format="html_document")'
