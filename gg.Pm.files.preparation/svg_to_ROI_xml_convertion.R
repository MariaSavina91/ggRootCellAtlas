library(xml2)
library(ggPlantmap)
library(RCurl)
library(ggplot2)

ra.integrated <- readRDS('C:/Work Folder/R/integrated_7_in_1_2022_09/ra_integrated_final_GO_MS_Rzones2.rds')
ra.integrated.meta<-ra.integrated@meta.data
rm(ra.integrated)

setwd('C:/Work Folder/R/integrated_7_in_1_2022_09/ggPlant_map_for_RootCellAtlas/gg.Pm.files.preporation/')

create_new_ggPm <- function(filename1, filename2, filename3, ra.integrated.meta, filename4) {
  
## Creating new ggPlantmap longitudinal
new.ggPlantmap = data.frame(matrix( 
  vector(), 0, 12, dimnames=list(c(), c("ROI.name","SubCellTypes", "CellTypes", "TissueSubTypes", "TissueTypes", "Zones", "Sections", "Atlas", "ROI.id","point","x","y"))), 
  stringsAsFactors=F) 

# Read the SVG file
svg_file <- readLines(paste0('C:/Work Folder/R/integrated_7_in_1_2022_09/ggPlant_map_for_RootCellAtlas/gg.Pm.files.preporation/',filename1), warn = FALSE)#"long root with dieing LRC_only 1 layer separate paths polygons.svg", warn = FALSE)
svg_content <- paste(svg_file, collapse = "\n")


# Parse the SVG content
svg <- read_xml(svg_content)


# Find all path elements in the SVG
paths <- xml_find_all(svg, "//svg:path")

# Loop through each path and extract information
for (path in paths) {
  
  # Add classname, id, name, and other static elements
   
  path_id<-xml_attr(path, "id")
  if (filename1=='long root with dieing LRC_only 1 layer separate paths polygons.svg') {
    path_id<-stringr::str_extract(string = path_id, pattern = "[0-9][0-9][0-9][0-9]")
  }
  else {
    path_id<-stringr::str_extract(string = path_id, pattern = "[0-9][0-9][0-9][0-9][0-9]")
  }
  
  
  # Extract and add each point
  path_d <- xml_attr(path, "d")
  coordinates <- strsplit(gsub("z", "", path_d)," |,")[[1]]
  
  # Initialize coordinates with the first point
  x <- 0
  y <- 0
  l<-1
  k<-1
  point<-0
  for (i in seq(1, length(coordinates), by = 1)) {
    # Update coordinates using relative values
    
    if ((coordinates[i])=='h') {
      h<-1
      l<-0
      v<-0
      k<-1
    }
    else if ((coordinates[i])=='v') {
      h<-0
      l<-0
      v<-1
      k<-1
    }
    
    else if ((coordinates[i])=='m' | (coordinates[i])=='l') {
      h<-0
      l<-1
      v<-0
      k<-1
    }
    
    else if ((coordinates[i])=='M' | (coordinates[i])=='L') {
      h<-0
      l<-2
      v<-0
      k<-1
    }
    
    else if ((coordinates[i])=='H') {
      h<-2
      l<-0
      v<-0
      k<-1
    }
    
    else if ((coordinates[i])=='V') {
      h<-0
      l<-0
      v<-2
      k<-1
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & l==1 & k==1) {
      x <- x + as.numeric(coordinates[i])
      y <- y + as.numeric(coordinates[i + 1])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
      k<-2
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & l==2 & k==1) {
      x <- as.numeric(coordinates[i])
      y <- as.numeric(coordinates[i + 1])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
      k<-2
    }
    
    else if (k==2) {
      k<-1
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & h==1) {
      x <- x + as.numeric(coordinates[i])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & h==2) {
      x <- as.numeric(coordinates[i])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & v==1) {
      y <- y + as.numeric(coordinates[i])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
    }
    
    else if ((coordinates[i])!='l'&(coordinates[i])!='h'&(coordinates[i])!='v' & v==2) {
      y <- as.numeric(coordinates[i])
      point<-point+1
      new.ggPlantmap[nrow(new.ggPlantmap)+1,]<-c(NA,NA,NA,NA,NA,NA,NA,NA,path_id,point,x,y)
    }
    
  }
}


new.ggPlantmap$ROI.name = as.character(new.ggPlantmap$ROI.name) #character()
new.ggPlantmap$SubCellTypes = as.character(new.ggPlantmap$SubCellTypes) #character()
new.ggPlantmap$CellTypes = as.character(new.ggPlantmap$CellTypes) #character()
new.ggPlantmap$TissueSubTypes = as.character(new.ggPlantmap$TissueSubTypes) #character()
new.ggPlantmap$TissueTypes = as.character(new.ggPlantmap$TissueTypes) #character()
new.ggPlantmap$Zones = as.character(new.ggPlantmap$Zones) #character()
new.ggPlantmap$Sections = as.character(new.ggPlantmap$Sections) #character()
new.ggPlantmap$Atlas = as.character(new.ggPlantmap$Atlas) #character()
new.ggPlantmap$ROI.id = as.numeric(new.ggPlantmap$ROI.id) #integer()
new.ggPlantmap$point = as.numeric(new.ggPlantmap$point) #integer()
new.ggPlantmap$x = as.numeric(new.ggPlantmap$x) #double()
new.ggPlantmap$y = as.numeric(new.ggPlantmap$y) #double()

new.ggPlantmap1<-new.ggPlantmap
new.ggPlantmap1$y<-new.ggPlantmap1$y+2*(mean(max(new.ggPlantmap$y),min(new.ggPlantmap$y))-new.ggPlantmap1$y)



ggPlantmap.plot(new.ggPlantmap,ROI.id,show.legend = FALSE)
ggPlantmap.plot(new.ggPlantmap1,ROI.id,show.legend = FALSE)
labeled<-read.csv(paste0('C:/Work Folder/R/integrated_7_in_1_2022_09/ggPlant_map_for_RootCellAtlas/gg.Pm.files.preporation/',filename2))#'labeled_long_reduced.csv')

new_id<-unique(new.ggPlantmap$ROI.id)
new_id<-as.data.frame(new_id)
Mapping<-new_id
for (i in 1:nrow(Mapping)){
  Mapping$x[i]<-mean(new.ggPlantmap$x[new.ggPlantmap$ROI.id==Mapping$new_id[i]])
  Mapping$y[i]<-mean(new.ggPlantmap$y[new.ggPlantmap$ROI.id==Mapping$new_id[i]])
}
Mapping$x_new<-ceiling(ncol(labeled)*(Mapping$x-min(new.ggPlantmap$x))/(max(new.ggPlantmap$x)-min(new.ggPlantmap$x)))
Mapping$y_new<-ceiling(nrow(labeled)*(Mapping$y-min(new.ggPlantmap$y))/(max(new.ggPlantmap$y)-min(new.ggPlantmap$y)))
Mapping$x_new[Mapping$x_new==0]<-1
Mapping$y_new[Mapping$y_new==0]<-1
a<-vector()
for (i in 1:nrow(Mapping)){
  Mapping$old_id[i]<-(labeled[Mapping$y_new[i],Mapping$x_new[i]])
}
# anyDuplicated(Mapping$old_id) == 0
# print(anyDuplicated(Mapping$old_id) == 0)

# Mapping$old_id[Mapping$old_id==522]

# print(Mapping$old_id[Mapping$old_id==0])
# print(Mapping$new_id[Mapping$old_id==0])
# 
# print(Mapping$old_id[Mapping$old_id==1])
# print(Mapping$new_id[Mapping$old_id==1])
if (filename1=='long root with dieing LRC_only 1 layer separate paths polygons.svg'){
  Mapping$old_id[Mapping$old_id==0]<-522
} 
if (filename1=='root_cross_trich1_only 1 layer separate paths polygons.svg') {
  Mapping$old_id[Mapping$old_id==0][1]<-1
}




for (i in 1:nrow(Mapping)){
  new.ggPlantmap1$ROI.id1[new.ggPlantmap1$ROI.id==Mapping$new_id[i]]<-Mapping$old_id[i]
}
new.ggPlantmap1$ROI.id<-new.ggPlantmap1$ROI.id1
new.ggPlantmap1<-new.ggPlantmap1[,-ncol(new.ggPlantmap1)]

# install.packages("R.matlab")
# library(R.matlab)
# 
# # Specify the path to your mat file
# mat_file_path <- "List_of_subtypes_long_with_dieingLRC_s.mat"
# 
# # Read the mat file
# mat_data <- readMat(mat_file_path)
# 
# df <- data.frame(CellType = character(), CellNumber = numeric(), stringsAsFactors = FALSE)
# 
# # Iterate through the mat_data list and append to the data frame
# for (cell_type in names(mat_data)) {
#   if (length(mat_data[[cell_type]]) > 0) {
#     cell_numbers <- mat_data[[cell_type]][1, ]
#     cell_numbers_df <- data.frame(CellType = rep(sub("\\.cells", "", cell_type), length(cell_numbers)),
#                                   CellNumber = as.numeric(cell_numbers),
#                                   stringsAsFactors = FALSE)
#     df <- rbind(df, cell_numbers_df)
#   }
# }

annotation<-read.csv(paste0('C:/Work Folder/R/integrated_7_in_1_2022_09/ggPlant_map_for_RootCellAtlas/gg.Pm.files.preporation/',filename3))#'cluster_labeled_map_longsec.csv')
new_df <- data.frame(cluster_id = character(), labeled_id = character(), labeled_values = character(), stringsAsFactors = FALSE)

# Iterate through the annotation data and append to the new data frame
for (i in 1:nrow(annotation)) {
  if (annotation$cluster_id[i]=='Young_LRC' | annotation$cluster_id[i]=='Dying_LRC'){
    cluster_id <- gsub("_", " ", annotation$cluster_id[i])
  }
  else{
    cluster_id <- annotation$cluster_id[i]
  }
  labeled_id <- annotation$labeled_id[i]
  labeled_values <- annotation$labeled_values[i]
  
  # Split labeled_values into a numeric vector
  cell_numbers <- as.numeric(unlist(strsplit(labeled_values, "\\s+")))
  
  # Create a data frame with cluster_id, labeled_id, and each cell number
  cell_numbers_df <- data.frame(cluster_id = rep(cluster_id, length(cell_numbers)),
                                labeled_id = rep(labeled_id, length(cell_numbers)),
                                labeled_values = cell_numbers,
                                stringsAsFactors = FALSE)
  
  # Append to the new data frame
  new_df <- rbind(new_df, cell_numbers_df)
}

for (i in 1:nrow(new_df)){
  new.ggPlantmap1$Atlas[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-new_df$cluster_id[i]
  new.ggPlantmap1$ROI.name[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-gsub("_| ", ".", new_df$cluster_id[i])
  new.ggPlantmap1$SubCellTypes[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$SubCellTypes_MS_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
  new.ggPlantmap1$CellTypes[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$CellTypes_MS_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
  new.ggPlantmap1$TissueSubTypes[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$TissueSubTypes_MS_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
  new.ggPlantmap1$TissueTypes[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$TissueTypes_MS_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
  new.ggPlantmap1$Zones[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$Zones_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
  new.ggPlantmap1$Sections[new.ggPlantmap1$ROI.id==new_df$labeled_values[i]]<-
    as.character(ra.integrated.meta$Sections_new[ra.integrated.meta$Atlas_new==new_df$cluster_id[i]][1])
}

write.table(new.ggPlantmap1, file = paste0('C:/Work Folder/R/integrated_7_in_1_2022_09/ggPlant_map_for_RootCellAtlas/gg.Pm.files.preporation/',filename4), sep = "\t", quote = FALSE, row.names = FALSE)#"ggPm.At.longroot.longitudinal.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# View(new.ggPlantmap1)
}

create_new_ggPm('long root with dieing LRC_only 1 layer separate paths polygons.svg', 'labeled_long_reduced.csv', 
                'cluster_labeled_map_longsec.csv', ra.integrated.meta, 'ggPm.At.longroot.longitudinal.txt')

create_new_ggPm('root_cross_only 1 layer separate paths polygons.svg', 'labeled_cross_reduced.csv', 
                'cluster_labeled_map_meristem1.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.m1.txt')

create_new_ggPm('root_cross_only 1 layer separate paths polygons.svg', 'labeled_cross_reduced.csv', 
                'cluster_labeled_map_meristem2.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.m2.txt')

create_new_ggPm('root_cross_only 1 layer separate paths polygons.svg', 'labeled_cross_reduced.csv', 
                'cluster_labeled_map_transition.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.t.txt')

create_new_ggPm('root_cross_only 1 layer separate paths polygons.svg', 'labeled_cross_reduced.csv', 
                'cluster_labeled_map_elongation1.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.e1.txt')

create_new_ggPm('root_cross_only 1 layer separate paths polygons.svg', 'labeled_cross_reduced.csv', 
                'cluster_labeled_map_elongation2.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.e2.txt')

create_new_ggPm('root_cross_trich1_only 1 layer separate paths polygons.svg', 'labeled_cross_diff_reduced.csv', 
                'cluster_labeled_map_differentiation.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.d.txt')

create_new_ggPm('root_cross_trich1_only 1 layer separate paths polygons _withoutPSE.svg', 'labeled_cross_diff_reduced.csv', 
                'cluster_labeled_map_differentiation.csv', ra.integrated.meta, 'ggPm.At.root.crosssection.d.noPSE.txt')
