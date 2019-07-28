# Estimation @ stamps2019

This folder contains the [statdivlab](http://statisticaldiversitylab.com)'s teaching materials for Sunday's lecture on estimating relative abundance and diversity.

# Details

- Lecture: see `estimation-lecture.pdf` above
- Lab:
  - Open your cloud RStudio server as we have done before by finding your name on [this page](https://hackmd.io/@astrobiomike/stamps2019) and clicking the link for RStudio.
  - When you're in RStudio, copy and paste the following code block to run these commands in order to get all the tutorials and data needed in a folder called `statdivlab`
```
setwd("~/")
dir.create("~/statdivlab/")
setwd("~/statdivlab/")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncob_tutorial/corncob_tutorial.html", destfile = "corncob_tutorial.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/corncobDESeq2/corncobDESeq2.html", destfile = "corncobDESeq2.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/diversity-lab.html", destfile = "diversity-lab.html")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/data/count_data.csv", destfile = "count_data.csv")
download.file("https://raw.githubusercontent.com/statdivlab/stamps2019/master/labs/data/sample_data.csv", destfile = "sample_data.csv")
```

## To begin the `corncob` tutorial
```
file.show("corncob_tutorial.html")
```
## To begin the `corncob` vs. `DESeq2` tutorial
```
file.show("corncobDESeq2.html")
```
## To begin the `corncob` and metagenomes tutorial
```
file.show("corncobmetagenome.html")
```
## To begin the diversity estimation tutorial
```
file.show("diversity-lab.html")
```

# Resources

- `corncob`: [paper](https://www.e-publications.org/ims/submission/AOAS/user/submissionFile/39562?confirm=b2fb2331) and [software](https://github.com/bryandmartin/corncob/)
- `breakaway`: [paper 1](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12332) and [paper 2](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssc.12206) and [software](https://github.com/adw96/breakaway)
- `DivNet`: [paper](https://www.biorxiv.org/content/10.1101/305045v1) and [software](https://github.com/adw96/DivNet)
