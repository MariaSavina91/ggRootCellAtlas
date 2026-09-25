# The e2 cross-section map was originally saved with the e1 labels
# (ROI.name, Sections and Atlas all ended in "e1"), so the late elongation
# zone showed e1 expression on the composite plot. Geometry and all other
# annotation were already correct. This script relabels the map and checks
# it against ggPm.At.root.crosssection.e2.txt, which has the right labels
# (but predates the Atlas_reduced column, so cannot be used directly).
#
# Run from the package root: Rscript data_raw/fix_e2_section_labels.R

load("data/ggPm.At.root.crosssection.e2.rda")
m <- ggPm.At.root.crosssection.e2

m$ROI.name <- sub("e1$", "e2", m$ROI.name)
m$Sections <- sub("^e1$", "e2", m$Sections)
m$Atlas <- sub("_e1$", "_e2", m$Atlas)
# Atlas_reduced uses "_e" for both elongation sections and stays unchanged

raw <- utils::read.delim("data_raw/ggPm.At.root.crosssection.e2.txt",
                         stringsAsFactors = FALSE)
stopifnot(isTRUE(all.equal(m[names(raw)], raw, check.attributes = FALSE)))

ggPm.At.root.crosssection.e2 <- m
save(ggPm.At.root.crosssection.e2,
     file = "data/ggPm.At.root.crosssection.e2.rda",
     compress = "gzip", version = 3)
