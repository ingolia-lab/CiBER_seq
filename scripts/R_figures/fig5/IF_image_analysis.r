#IF image analysis

library(dplyr)

dapi <- read.csv("~/CiBER_seq_package/scripts/R_figures/fig5/merged_dapi.csv", 
                 header=TRUE)
gcn4 <- read.csv("~/CiBER_seq_package/scripts/R_figures/fig5/merged_gcn4.csv", 
                 header=TRUE)

dapi_final <- filter(dapi, dapi$Area != "Area")
dapi_final <- filter(dapi_final, dapi_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv != "test.csv")
gcn4_final <- filter(gcn4, gcn4$Area != "Area")
gcn4_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv <- gcn4_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei_GCN4.csv
gcn4_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv <- gsub("nuclei_GCN4.csv", "nuclei.csv", gcn4_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv)

merge_final <- merge(dapi_final, gcn4_final, by=c("X","Area", "SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv"))

nup133 <- filter(merge_final, grepl("Nup133", merge_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
nup145 <- filter(merge_final, grepl("Nup145", merge_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
ubc9 <- filter(merge_final, grepl("UBC9", merge_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
ulp1 <- filter(merge_final, grepl("ULP1", merge_final$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))

nup133_un <- filter(nup133, grepl("Nup133_UN", nup133$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
nup133_tet <- filter(nup133, grepl("Nup133_TET", nup133$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))

nup133_un$gcn4_norm <- as.numeric(nup133_un$Mean.y)/as.numeric(nup133_un$Mean.x)
nup133_tet$gcn4_norm <- as.numeric(nup133_tet$Mean.y)/as.numeric(nup133_tet$Mean.x)
mean(nup133_un$gcn4_norm)
median(nup133_un$gcn4_norm)
mean(nup133_tet$gcn4_norm)
median(nup133_tet$gcn4_norm)


nup145_un <- filter(nup145, grepl("Nup145_UN", nup145$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
nup145_tet <- filter(nup145, grepl("Nup145_TET", nup145$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))

nup145_un$gcn4_norm <- as.numeric(nup145_un$Mean.y)/as.numeric(nup145_un$Mean.x)
nup145_tet$gcn4_norm <- as.numeric(nup145_tet$Mean.y)/as.numeric(nup145_tet$Mean.x)
mean(nup145_un$gcn4_norm)
median(nup145_un$gcn4_norm)
mean(nup145_tet$gcn4_norm)
median(nup145_tet$gcn4_norm)

ubc9_un <- filter(ubc9, grepl("UBC9_UN", ubc9$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
ubc9_tet <- filter(ubc9, grepl("UBC9_TET", ubc9$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))

ubc9_un$gcn4_norm <- as.numeric(ubc9_un$Mean.y)/as.numeric(ubc9_un$Mean.x)
ubc9_tet$gcn4_norm <- as.numeric(ubc9_tet$Mean.y)/as.numeric(ubc9_tet$Mean.x)
mean(ubc9_un$gcn4_norm)
median(ubc9_un$gcn4_norm)
mean(ubc9_tet$gcn4_norm)
median(ubc9_tet$gcn4_norm)


ulp1_un <- filter(ulp1, grepl("ULP1_UN", ulp1$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))
ulp1_tet <- filter(ulp1, grepl("ULP1_TET", ulp1$SUM_Nup133_TET_01_R3D_D3D.tifnuclei.csv))

ulp1_un$gcn4_norm <- as.numeric(ulp1_un$Mean.y)/as.numeric(ulp1_un$Mean.x)
ulp1_tet$gcn4_norm <- as.numeric(ulp1_tet$Mean.y)/as.numeric(ulp1_tet$Mean.x)
mean(ulp1_un$gcn4_norm)
median(ulp1_un$gcn4_norm)
mean(ulp1_tet$gcn4_norm)
median(ulp1_tet$gcn4_norm)

# Violin with blox plot inside for estradiol titration
nup133_un$cond <- "a.nup133_un"
nup133_tet$cond <- "b.nup133_tet"
nup145_un$cond <- "c.nup145_un"
nup145_tet$cond <- "d.nup145_tet"
ubc9_un$cond <- "e.ubc9_un"
ubc9_tet$cond <- "f.ubc9_tet"
ulp1_un$cond <- "g.ulp1_un"
ulp1_tet$cond <- "h.ulp1_tet"

cond <- append(nup133_un$cond, nup133_tet$cond)
cond <- append(cond, nup145_un$cond)
cond <- append(cond, nup145_tet$cond)
cond <- append(cond, ubc9_un$cond)
cond <- append(cond, ubc9_tet$cond)
cond <- append(cond, ulp1_un$cond)
cond <- append(cond, ulp1_tet$cond)

gcn4_norm_val <- append(nup133_un$gcn4_norm, nup133_tet$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, nup145_un$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, nup145_tet$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, ubc9_un$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, ubc9_tet$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, ulp1_un$gcn4_norm)
gcn4_norm_val <- append(gcn4_norm_val, ulp1_tet$gcn4_norm)

dapi_gcn4 <- data.frame(cond, gcn4_norm_val)

dap_gcn <- ggplot(dapi_gcn4, aes(x=cond, y=gcn4_norm_val, fill=cond)) +
  theme() +  
  geom_violin(trim = TRUE, scale = "area", width=1.1)+
  geom_boxplot(width=0.09, fill="white", outlier.shape = NA)+
  labs(title="Plot of GCN4 signal normalzed to DAPI",x="guide cond", y = "Nuclear GCN4") + 
  scale_y_continuous(breaks=seq(0,1.5,0.1))
  
dap_gcn
  
hist(as.numeric(ulp1_un$Area), breaks = 100)
hist(as.numeric(nup133_un$Area), breaks = 100)

median(as.numeric(nup133_un$Area))
median(as.numeric(nup133_tet$Area))
median(as.numeric(nup145_un$Area))
median(as.numeric(nup145_tet$Area))
median(as.numeric(ubc9_un$Area))
median(as.numeric(ubc9_tet$Area))
median(as.numeric(ulp1_un$Area))
median(as.numeric(ulp1_tet$Area))


