#libs
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(devtools)
library(arcdiagram)
#funs
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#johnson
local <- read.csv("/Users/.../.csv", row.names = 1)
TARGET <- read.csv("/Users/.../.csv", row.names = 1)

#colnames(ads_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(lietal_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(ped_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(target_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")

# remove single rules
local <- local[grep(",", local$features),]
TARGET <- TARGET[grep(",", TARGET$features),]


#local <- local[which(p.adjust(local$pValue, method = "fdr")<0.2),]
#TARGET <- TARGET[which(p.adjust(TARGET$pValue, method = "fdr")<0.2),]


#genetic
#ads_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Adults_Results2020-09-28/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-28.txt")
#lietal_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Lietal_Results2020-09-28/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-28.txt")
#ped_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Pediatric_Results2020-09-28/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-28.txt")
#target_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_TARGET_Results2020-09-28/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-28.txt")
#colnames(ads_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(lietal_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(ped_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
#colnames(target_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")

#ads_g <- ads_g[grep(",", ads_g$features),]
#lietal_g <- lietal_g[grep(",", lietal_g$features),]
#ped_g <- ped_g[grep(",", ped_g$features),]
#target_g <- target_g[grep(",", target_g$features),]

#ads_g <- ads_g[which(quantile(ads_g$coverageRHS, probs = c(0.05, 0.95))[2] < ads_g$coverageRHS),]
#lietal_g <- lietal_g[which(quantile(lietal_g$coverageRHS, probs = c(0.3, 0.7))[2] < lietal_g$coverageRHS),]
#ped_g <- ped_g[which(quantile(ped_g$coverageRHS, probs = c(0.2, 0.8))[2] < ped_g$coverageRHS),]
#target_g  <- target_g[which(quantile(target_g$coverageRHS, probs = c(0.05, 0.95))[2]<target_g$coverageRHS),]


## aggregate the same
#local <- aggregate(.~features+decision, FUN=mean, data = local, na.action = na.pass)
#local$supportRHS <- round(local$supportRHS)
#li_target <- aggregate(.~features+decision, FUN=mean, data = li_target, na.action = na.pass)
#li_target$supportRHS <- round(li_target$supportRHS)

#nams <- c("Local adults", "Li et al.")
nams <- c("Local Pediatric", "TARGET")

local$decision <- gsub("1","",local$decision)
TARGET$decision <- gsub("1","",TARGET$decision)

source("R/network_gen.R")
vn1 <- network_gen(local, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.5, type = 'RDF', fract=1, topN = 50)
vn2 <- network_gen(local, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.5, type = 'RDF', fract=1, topN = 50)

vn3 <- network_gen(TARGET, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.5, type = 'RDF', fract=1, topN = 50)
vn4 <- network_gen(TARGET, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.5, type = 'RDF', fract=1, topN = 50)

# genetic
#vn1 <- network_gen(local, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
#vn2 <- network_gen(local, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
#vn3 <- network_gen(li_target, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
#vn4 <- network_gen(li_target, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)

### check top nodes
tn1 <- vn1$Diagnosis$nodes$NodeConnection
names(tn1) <- paste0(as.character(vn1$Diagnosis$nodes$label),"(",gsub("3","3_HIGH",gsub("2","2_MEDIUM",gsub("1", "1_LOW", vn1$Diagnosis$nodes$DiscState))),")")
tn2 <- vn2$Relapse$nodes$NodeConnection
names(tn2) <- paste0(as.character(vn2$Relapse$nodes$label),"(",gsub("3","3_HIGH",gsub("2","2_MEDIUM",gsub("1", "1_LOW", vn2$Relapse$nodes$DiscState))),")")

tn3 <- vn3$Diagnosis$nodes$NodeConnection
names(tn3) <- paste0(as.character(vn3$Diagnosis$nodes$label),"(",gsub("3","3_HIGH",gsub("2","2_MEDIUM",gsub("1", "1_LOW", vn3$Diagnosis$nodes$DiscState))),")")
tn4 <- vn4$Relapse$nodes$NodeConnection
names(tn4) <- paste0(as.character(vn4$Relapse$nodes$label),"(",gsub("3","3_HIGH",gsub("2","2_MEDIUM",gsub("1", "1_LOW", vn4$Relapse$nodes$DiscState))),")")

#if(length(which(grepl(",",vn1$Diagnosis$nodes$group)))!=0){
#  tn11 <- vn1$all$nodes$NodeConnection[which(grepl(",",vn1$all$nodes$group))]
#  names(tn11) <- paste0(as.character(vn1$all$nodes$label),"(",gsub("3","HIGH",gsub("2","MEDIUM",gsub("1", "LOW", vn1$all$nodes$DiscState))),")")[which(grepl(",",vn1$all$nodes$group))]
#  tn1 <- c(tn1, tn11)
#  tn2 <- c(tn2, tn11)
#}
#if(length(which(grepl(",",vn2$all$nodes$group)))!=0){
#  tn22 <- vn2$all$nodes$NodeConnection[which(grepl(",",vn2$all$nodes$group))]
#  names(tn22) <- paste0(as.character(vn2$all$nodes$label),"(",gsub("3","HIGH",gsub("2","MEDIUM",gsub("1", "LOW", vn2$all$nodes$DiscState))),")")[which(grepl(",",vn2$all$nodes$group))]
#  tn3 <- c(tn3,tn22)
#  tn4 <- c(tn4,tn22)
#}

allF <- unique(c(names(tn1), names(tn2), names(tn3), names(tn4)))
discS <- gsub("[\\(\\)]", "", regmatches(allF, gregexpr("\\(.*?\\)", allF)))
discS <- tolower(discS)

new_tn1 <- c()
new_tn2 <- c()
new_tn3 <- c()
new_tn4 <- c()
new_tn5 <- c()
new_tn6 <- c()

for(i in 1:length(allF)){
  tn1n <- which(names(tn1)==allF[i])
  if(length(tn1n)==0){
    new_tn1[i] <- 0
  }else{
    new_tn1[i] <- unname(tn1[tn1n])
  }
  tn2n <- which(names(tn2)==allF[i])
  if(length(tn2n)==0){
    new_tn2[i] <- 0
  }else{
    new_tn2[i] <- unname(tn2[tn2n])
  }
  tn3n <- which(names(tn3)==allF[i])
  if(length(tn3n)==0){
    new_tn3[i] <- 0
  }else{
    new_tn3[i] <- unname(tn3[tn3n])
  }
  tn4n <- which(names(tn4)==allF[i])
  if(length(tn4n)==0){
    new_tn4[i] <- 0
  }else{
    new_tn4[i] <- unname(tn4[tn4n])
  }
}

#dftns <- as.data.frame(scale(data.frame(new_tn1, new_tn2, new_tn3, new_tn4)))
dftns <- as.data.frame(data.frame(range01(new_tn1), range01(new_tn2), range01(new_tn3), range01(new_tn4)))
discNames <- c("Local Pediatric Diagnosis", "Local Pediatric Relapse", "TARGET Diagnosis", "TARGET Relapse")
colnames(dftns) <- discNames
rownames(dftns) <- allF
#n1 <- 1
dl2 <- discS #[rowSums(dftns == 0) <= n1]
#dftns <- as.data.frame(scale(dftns[rowSums(dftns == 0) <= n1, ]))

#heatmap(as.matrix(dftns))
#dftns_f <- dftns[-which(is.na(unname(rowMeans(dftns)))),]
dftns_f <- dftns

#shapiro.test(dftns_f$Local_diagnosis) # <0.05 means that it is significantly different from normal dist.
#shapiro.test(dftns_f$Local_relapse)
# we use kendall then


names(dl2) <- "dl"
#as.character(vis_out$`ASPERGER'S DISORDER`$nodes$DiscState)[match(rownames(dftns_f),names(tn2))]
#as.character(vis_out$AUTISM$nodes$DiscState)[match(rownames(dftns_f),names(tn3))]
#as.character(vis_out$CONTROL$nodes$DiscState)[match(rownames(dftns_f),names(tn4))]
M <- as.matrix(dftns_f)
rownames(M) <- gsub("\\s*\\([^\\)]+\\)","",as.character(rownames(M)))
cn <- colnames(M)
Heatmap(M, name="node\nconnection", clustering_distance_columns = "kendall",
        show_column_names = FALSE, cluster_rows = T,
        clustering_distance_rows = "kendall", na_col = "black", row_names_gp = gpar(fontsize = 10),
        border = TRUE,  col = c("seashell","brown1","darkred"), 
        right_annotation = rowAnnotation(dl = dl2, name="levels", col=list(dl=c("1_low"="#56B4E9","2_medium"="#999999","3_high"="#E69F00")), annotation_legend_param = list(dl = list(title = "discretization\nlevel", labels = c("low", "medium", "high")))),
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn, offset = unit(1, "npc")))
)
vec <- seq(0,1,by=1/dim(M)[1])
n <- dim(M)[1]
decorate_heatmap_body("node\nconnection", {
  grid.polyline(x=rep(c(0, 1), n+1),
                y=rep(vec, each = 2),
                id.lengths=rep(2,n+1),
                gp=gpar(lwd=0.5)
  )
  
})

#### arch plot
arcNet <- function(net, icol, dec, mainTitle){
  gene <- rownames(dftns[which.max(dftns[,icol]),])
  discLev <- substr(gsub("[\\(\\)]", "", regmatches(gene, gregexpr("\\(.*?\\)", gene))[[1]]), 1, 1)
  gene2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(gene))
  edges <- net[[dec]]$edges[,1:4]
  edges_from <- as.character(edges$from)
  edges_to <- as.character(edges$to) ## gsub("\\=.*","",as.character(net$short$edges$to))
  indx <- unique(c(which(edges_from == paste0(gene2,"=",discLev)),which(edges_to == paste0(gene2,"=",discLev))))
  edges2 <- edges[indx,]
  M <- as.matrix(edges2[,1:2])
  
  #colors
  valsCols <- round(edges2$connNorm, digits = 1)*10
  edgesCon <- as.numeric(as.factor(as.character(valsCols)))
  colors <- colorRampPalette(c("gainsboro","lavender","darkorchid3"))(length(unique(valsCols)))[edgesCon]
  labelsNodes <- unique(c(t(M)))
  colNodes <- c("#56B4E9","#999999","#E69F00")[as.numeric(sub(".*=","",labelsNodes))]
  
  #nodes values
  nodesDlev <- substr(gsub("[\\(\\)]", "", regmatches(rownames(dftns), gregexpr("\\(.*?\\)", rownames(dftns)))),1,1)
  nodesNams <- gsub("\\s*\\([^\\)]+\\)","",as.character(rownames(dftns)))
  nodes2 <- paste0(nodesNams,"=",nodesDlev)
  nodeSize <- round(dftns[match(labelsNodes, nodes2),icol],digits = 1)*10
  
  ordV <- c(1,order(edgesCon, decreasing = T)+1)
  
  colsLabs <- rep("black",length(labelsNodes))
  colsLabs[which(labelsNodes == paste0(gene2,"=",discLev))] <- "red"
  arcplot(M, lwd.arcs=edgesCon, col.arcs = colors, col.nodes = colNodes, labels=sub("=.*","",labelsNodes),
          ordering=ordV, col.labels=colsLabs, cex.labels=1, font=1, lwd.nodes = nodeSize, main=mainTitle)
}

par(mfcol=c(2,2))
arcNet(vn1, 1, "Diagnosis", discNames[1])
arcNet(vn2, 2, "Relapse", discNames[2])
arcNet(vn3, 3, "Diagnosis", discNames[3])
arcNet(vn4, 4, "Relapse", discNames[4])