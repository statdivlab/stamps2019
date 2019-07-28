# stamps2019

Welcome! This repository contains the statdivlab's teaching materials for STAMPS @ MBL in 2019.

This repository will continue to be updated until the end of STAMPS. All materials are drafts until August 1.


# Who?

The [Stat Div Lab](http://statisticaldiversitylab.com/) is a group of badass statisticians. Your Stat Div Lab representatives at STAMPS in 2019 are:

- Amy Willis (PI), Assistant Professor, UW Biostatistics [@AmyDWillis](https://twitter.com/AmyDWillis)
- Bryan Martin, PhD Candidate, UW Statistics, [@BryanDMartin_](https://twitter.com/BryanDMartin_)
- Pauline Trinh, PhD Candidate, UW Occupational and Environmental Health [@paulinetrinh](https://twitter.com/paulinetrinh)

The Statistical Diversity Lab is the research group of Amy Willis. The Stat Div Lab develops rigorous statistical methods for the analysis of microbiome data. We are excited to be at STAMPS teaching `Statistics Bootcamp`, and `Statistical Estimation`. We are typically very opinionated, though, so feel free to chat to us about anything!

# Tutorials
- Open your cloud RStudio server as we have done before by finding your name on [this page](https://hackmd.io/@astrobiomike/stamps2019) and clicking the link for RStudio.
- When you're in RStudio, copy and paste the following code block to run these commands in order to get all the tutorials and data needed in a folder called `statdivlab`
```
setwd("~/")
dir.create("~/statdivlab/")
setwd("~/statdivlab/")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncob_tutorial/corncob_tutorial.html", destfile = "corncob_tutorial.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncobDESeq2/corncobDESeq2.html", destfile = "corncobDESeq2.html")
download.file("https://github.com/statdivlab/stamps2019/blob/master/labs/corncob_metagenome/corncob_metagenome.html", destfile = "corncob_metagenome.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/diversity-lab.html", destfile = "diversity-lab.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncob_metagenome/count_data.csv", destfile = "count_data.csv")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncob_metagenome/sample_data.csv", destfile = "sample_data.csv")
```

## To begin the `corncob` tutorial
```
file.show("~/statdivlab/corncob_tutorial.html")
```
## To begin the `corncob` and metagenomics tutorial
```
file.show("~/statdivlab/corncob_metagenome.html")
```
## To begin the `corncob` vs. `DESeq2` tutorial
```
file.show("~/statdivlab/corncobDESeq2.html")
```
## To begin the diversity estimation tutorial
```
file.show("~/statdivlab/diversity-lab.html")
```

# Questions

We would love to meet you -- please come and sit with us at mealtimes and tell us about your microbiome data! We best in person.

If you can't catch us, feel free to (1) contact us on Twitter, (2) file an issue on one of our github repositories, or (3) email us. Bryan and Pauline are significantly better at responding to e-mail than Amy, so consider copying one of them in.

# Citations

If you find our teaching materials helpful, please consider citing the following paper

- Willis, A.D. (2019).[Rigorous Statistical Methods for Rigorous Microbiome Science](https://msystems.asm.org/content/4/3/e00117-19). mSystems. doi: 10.1128/mSystems.00117-19

If you use our methods, please cite the relevant paper or preprint. Let us know if you're not sure which article corresponds which method!

# Methods

- Estimating richness: breakaway
- Estimating evenness: DivNet
- Estimating relative abundance: corncob
- Estimating real relative abundance: metacal
- Finding co-abundant genes: CAGs
- Pangenomics: anvi'o
