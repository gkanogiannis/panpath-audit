LATEX = lualatex
BIBER = biber
BUILD = build
JOB = panpath-audit-paper

.PHONY: all clean results

all: $(BUILD)/$(JOB).pdf

results:
	python3 scripts/summarize.py

$(BUILD)/$(JOB).pdf: manuscript.tex references.bib figures/workflow.tex results/derived/results_table.tex results/derived/macros.tex
	mkdir -p $(BUILD)
	$(LATEX) --interaction=nonstopmode --halt-on-error --output-directory=$(BUILD) --jobname=$(JOB) manuscript.tex
	$(BIBER) --input-directory=$(BUILD) --output-directory=$(BUILD) $(JOB)
	$(LATEX) --interaction=nonstopmode --halt-on-error --output-directory=$(BUILD) --jobname=$(JOB) manuscript.tex
	$(LATEX) --interaction=nonstopmode --halt-on-error --output-directory=$(BUILD) --jobname=$(JOB) manuscript.tex

clean:
	rm -f $(BUILD)/$(JOB).aux $(BUILD)/$(JOB).bbl $(BUILD)/$(JOB).bcf \
		$(BUILD)/$(JOB).blg $(BUILD)/$(JOB).log $(BUILD)/$(JOB).out \
		$(BUILD)/$(JOB).run.xml $(BUILD)/$(JOB).toc $(BUILD)/$(JOB).pdf

