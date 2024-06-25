
###### PT 1 LIBRARIES & SETUP###########

list.of.packages <- c("readxl", "openxlsx", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}

library(readxl)
library(openxlsx)
library(tidyverse)

#Set working directory for output
print("Please select the folder/directory to save output in.") 
setwd(choose.dir())

#determine whether to start from beginning and generate interactome data, or skip ahead for specific analysis
STEP <- readline("::Do you already have interactome data generated (Y for yes, N for No):: \n")


if(STEP == "N" | STEP == "n"){
  
  ###### PT 2 DATASETS ###########
  
  print("Please select the excel file with ligand data.") 
  lig <- read_excel(file.choose())
  lig_name <- readline("::What tissue or cell type is this:: \n")
  print("Here is a small preview of the table.")
  print(head(lig))
  lig_gene_col <- as.numeric(readline("::Enter the column number of the column with gene symbols:: \n"))
  colnames(lig)[lig_gene_col] <- "lig_gene_symbol"
  lig$lig_gene_symbol <- toupper(lig$lig_gene_symbol)
  lig_col_total <- ncol(lig)
  
  print("Please select the excel file with receptor data.")
  rec <- read_excel(file.choose())
  rec_name <- readline("::What tissue or cell type is this:: \n")
  print("Here is a small preview of the table.")
  print(head(rec))
  rec_gene_col <- as.numeric(readline("::Enter the column number of the column with gene symbols:: \n"))
  colnames(rec)[rec_gene_col] <- "rec_gene_symbol"
  rec$rec_gene_symbol <- toupper(rec$rec_gene_symbol)
  rec_col_total <- ncol(rec)
  
  
  
  ###### PT 3 INTERACTOME LIST (v2.3) ###########
  
  print("Please select the excel file with interactome list.")
  int_list <- read_excel(file.choose())
  print("Here is a small preview of the interactome list.")
  print(head(int_list))
  
  
  
  ###### PT 4 IMAKE INTERACTOME ###########
  make_int <- function(lig_table, lig_cols, rec_table, rec_cols){
    interactome <- inner_join(int_list, lig_table, by = c("Ligand" = "lig_gene_symbol")) #filter in only those genes in the ligand data table
    interactome <- inner_join(interactome, rec_table, by = c("Receptor" = "rec_gene_symbol")) #filter in only those genes in the receptor data table
    interactome <- interactome[, c(1:12, 24:(22+lig_cols), 13:23, (23+lig_cols):(21+lig_cols+rec_cols))] #rearrange
    return(interactome)
  }
  
  interactome1 <- make_int( lig, lig_col_total, rec, rec_col_total)

  ###### PT 5 EXPORT #######
  write.xlsx( interactome1, file = paste0( lig_name, "_to_", rec_name, "_ALLdata.xlsx"))
  #This excel file will be loaded for PT6-10, which are done for each unique interactome analysis derived from this data.
  
}


###### ** PTS 6-11 CAN BE REPEATED FOR DIFFERENT INTERACTOME ANALYSES FROM THE SAME DATASET ** ######
STEP <- "Y"



if(STEP == "Y" | STEP == "y"){
  
  ###### PT 6 LOAD FULL INTERACTOME DATA & CHOOSE ANALYSIS (LIG/REC) ###########
  
  print("Please select the excel file with all of the interactome data.")
  interactome <- read_excel(file.choose())
  print("Here is a small preview of the interactome data.")
  print(head(interactome))
  
  lig_name <- readline("::For this specific analysis, which tissue or cell type are the LIGANDS from:: \n")
  rec_name <- readline("::For this specific analysis, which tissue or cell type are the RECEPTORS from:: \n")
  
  int_cols <- data.frame(colnames(interactome))
  int_cols$number <- as.numeric(rownames(int_cols))
  print("These are the columns of the interactome table.")
  print.data.frame(int_cols)
  
  

  ###### PT 7 REMOVE UNUSED COLUMNS  ###########
  
  print("Enter the numbers of the columns that should be present in the analysis (separated by a space), then click ENTER twice.")
  keep_col <- scan()
  removed_columns <- int_cols$colnames.interactome.[setdiff(int_cols$number, keep_col)]
  interactome <- interactome[,keep_col]
  
  print("Here is a small preview of the interactome analysis table.")
  print(head(interactome))
  
  int_cols <- int_cols[keep_col,]
  rownames(int_cols) <- NULL
  int_cols$number <- as.numeric(rownames(int_cols))
  
  print("These are the columns of the interactome analysis table.")
  print.data.frame(int_cols)
  
  log <- list( Analyzed_Columns = colnames(interactome), Removed_Columns = removed_columns)
  
  
  
  ###### PT 8 FILTER ROWS/INTERACTIONS  ###########
  
  filtering <- as.numeric(readline("::Enter the number of the column from which count values will be filtered.\n
  (Enter 0 if not filtering in any columns):: \n"))
  
  while (filtering > 0){
    
    filt_thresh <- as.double(readline(":: Enter filtering value. (In the next prompt you can choose the filtering direction) :: \n"))
    filt_dir <- readline(":: If you would like to remove all values below this number type 'b'. 
                         If you would like to remove all values above this number type 'a' :: \n")
    
    if(filt_dir == "B" | filt_dir == "b"){
      interactome <- interactome[ which(interactome[ filtering ] >= filt_thresh ) ,]
    } else{
      interactome <- interactome[ which(interactome[ filtering ] <= filt_thresh ) ,]
    }
    
    filter_summary <- data.frame( Filter_Column = c( colnames(interactome)[filtering] ), Filter_Value = c( filt_thresh))
    log <- append( log, list(filter_summary ))
    
    filtering <- as.numeric(readline("::Enter the number of another column from which count values will be filtered.\n
      (Enter 0 if not filtering in any other columns):: \n"))
    if(filtering == 0){ break }
    
  }
  
  
  
  ###### PT 9 RANK INTERACTIONS & EXPORT ###########
  
  print("These are the columns of the interactome table.")
  print.data.frame(int_cols)
  print("Enter the numbers of the columns that should be used to rank the interactions (separated by a space), then click ENTER twice.")
  rank_col <- scan()
  
  # will rank all interactions and assume all given columns have the same weight; will add coded in the future to adjust weights
  rank_table <- interactome[, rank_col]
  rank_table <- data.frame(lapply(rank_table, function(x) x/sum(x)))
  rank_table$sums <- rowSums(rank_table)
  rank_table$rank <- rank(rank_table$sums)
  interactome$rank <- rank_table$rank
  interactome <- interactome[order(interactome$rank,decreasing=TRUE),]
  
  
  write.xlsx( interactome, file = paste0( lig_name, "_to_", rec_name, "_interactome_analysis.xlsx"))
  
  
  
  ###### PT 10 SELECT TOP INTERACTIONS  ###########
  
  top_num <- as.numeric(readline("::How many interactions would you like to show in the figures:: \n"))
  interactome <- interactome[1:top_num,]
  
  
  
  ###### PT 11 PREPARE TXT FILE FOR SANKEYMATIC ###########
  
  gene_types <- data.frame(Types=c('calcium-binding', 'cell adhesion','cell junction','chaperone',
                                   'chromatin-associated','cytokine','cytoskeletal','defense/immunity',
                                   'DNA metabolism','extracellular matrix','gene-specific transcriptional regulator',
                                   'G-protein coupled receptor','growth factor','heparin binding','hormone',
                                   'ion channel','membrane traffic','metabolite interconversion enzyme','neuropeptide',
                                   'other intercellular signal','other transmembrane signal receptor','other/unclassified',
                                   'protein modifying enzyme', 'protein-binding activity modulator','RNA metabolism',
                                   'scaffold/adaptor', 'serine/threonine protein kinase receptor', 'structural',
                                   'transfer/carrier', 'translational', 'transporter', 'tyrosine kinase receptor'))
  
  gene_types$hexcodes <- c('#c62828','#f44336','#e57373','#AD1457','#E91E63','#F48FB1','#6A1B9A','#8E24AA','#CE93D8',
                           '#283593','#3F51B5','#7986CB','#1565C0','#2196F3','#64B5F6','#00838F','#00BCD4','#80DEEA',
                           '#2E7D32','#4CAF50','#A5D6A7','#9E9D24','#CDDC39','#E6EE9C','#EF6C00','#F9A825','#FFEB3B',
                           '#FFFFFF','#3E2723','#795548','#BCAAA4','#263238')
  
  #More Colors: https://hexcolorspicker.com/tools/materialui/colors.html
  
  interactome <- left_join(interactome, gene_types, join_by("lig_final_category" == "Types"))
  interactome <- left_join(interactome, gene_types, join_by("rec_final_category" == "Types"))
  interactome$ligCOLOR <- paste0(":",interactome$Ligand," ",interactome$hexcodes.x)
  interactome$recCOLOR <- paste0(":",interactome$Receptor,". ",interactome$hexcodes.y)
  interactome$FLOW <- paste0(interactome$Ligand," [1] ",interactome$Receptor,".")
  sankey_text <- c(interactome$ligCOLOR, interactome$recCOLOR, interactome$FLOW )
  
  write.table( sankey_text, file = paste0( lig_name, "_to_", rec_name, "_sankey.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Load the text file at https://www.sankeymatic.com/build/
  
  
  
  ###### PT 12 BAR GRAPHS ###########
  
  print("These are the columns of the interactome table.")
  print.data.frame(int_cols)
  
  plotting <- as.numeric(readline("::Enter the number of the column with numeric values to plot.\n
  (Enter 0 if you do not want bar plots of any columns) :: \n"))
  
  while (plotting > 0){
    error_bars_col <- as.numeric(readline("::Enter the number of the column with error bar values. Type 0 if there is not any. :: \n"))
    x_axis_genes <- as.numeric(readline("::Enter the number of the column with gene symbols associated with the values :: \n"))
    color_groups <- as.numeric(readline("::Enter the number of the column with gene categories. Type 0 if there is not any. :: \n"))
    
    y_axis_name <- readline("::Please type y-axis label:: \n")
    y_axis_max <- as.numeric(readline("::Please type y-axis maximum value:: \n"))
    plot_name <- readline("::Please type title for plot:: \n")
    
    numeric_col_name <- colnames(interactome)[plotting]
    gene_col_name <- colnames(interactome)[x_axis_genes]
    
    colnames(interactome)[plotting] <- "VALUE_COL" #doing this so bars in order of largest to smallest, order fxn in next step need col name
    interactome <- interactome[order(interactome$VALUE_COL, decreasing = TRUE),] #doing this so bars in order of largest to smallest
    colnames(interactome)[plotting] <- numeric_col_name #changing it back
    
    if(color_groups>0){
      
      color_col_name <- colnames(interactome)[color_groups] 
      interactome <- left_join(interactome, gene_types, join_by(!!sym(color_col_name) == "Types")) #doing this so last column will have right hexcodes
      hex_colors <- interactome[,c(color_groups,length(interactome))]
      colnames(hex_colors)[1] <- "CATEGORY"
      hex_colors <- hex_colors[order(hex_colors$CATEGORY),]
      hex_colors <- unique(pull(hex_colors,2 ))
      
      plot <- ggplot( data = interactome[!duplicated(interactome[x_axis_genes]),], 
                      aes( x = reorder(!!sym(gene_col_name), -!!sym(numeric_col_name)), y = !!sym(numeric_col_name), fill = !!sym(color_col_name))) 
      bar_colors <- scale_fill_manual(values=hex_colors)
      
    } else{
      plot <- ggplot( data = interactome[!duplicated(interactome[x_axis_genes]),], 
                      aes( x = reorder(!!sym(gene_col_name), -!!sym(numeric_col_name)), y = !!sym(numeric_col_name))) 
      bar_colors <- scale_fill_grey()
    }
    
    bar <- geom_bar(stat="identity", color = 'black')
    plot_labels <- labs(title = plot_name, x = "Genes", y = y_axis_name)
    font <- theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(face = "italic"))
    
    if(error_bars_col>0){
      
      error_col_name <- colnames(interactome)[error_bars_col]
      error_bar <- geom_errorbar(aes( ymin = !!sym(numeric_col_name) - !!sym(error_col_name), ymax = !!sym(numeric_col_name) + !!sym(error_col_name)), width=0.2)
      bar_graph <- plot + bar + bar_colors + error_bar + plot_labels + font + ylim(0,y_axis_max)
      
    } else{
      bar_graph <- plot + bar + bar_colors + plot_labels + font + ylim(0,y_axis_max)
    }
    
    ggsave(width = 20, height = 6, filename = paste0( plot_name, "_sankey.svg"), plot = bar_graph)
   
    plotting <- as.numeric(readline("::Enter the number of another column with numeric values to plot.\n
      (Enter 0 if you do not want any more bar plots) :: \n"))
    if(plotting == 0){ break }
    
    
    
  }
  
  
}




