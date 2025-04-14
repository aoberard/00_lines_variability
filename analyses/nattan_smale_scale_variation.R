# Script parameters ----

## Library  ----
library("dplyr")
library("ggplot2")


# Data import ----
lines_modalities <- readxl::read_excel( here::here("data/", "raw-data/", "20241203_lines_variables.xlsx") )

quadrat <- readxl::read_excel(here::here("data", "raw-data", "20250411_quadrat_flora_bocage_2023.xlsx"))

transect <- readxl::read_excel(here::here("data", "raw-data", "20250411_transect_flora_bocage_2023.xlsx"))

strate <- readxl::read_excel(here::here("data", "raw-data", "20250411_strates_flora_bocage_2023.xlsx"))

splist <- readxl::read_excel(here::here("data", "raw-data", "20250411_species_list_flora_bocage_2023.xlsx"))


# Nattan's bocage small-scale data ----

## Trapping lines small-scale variables PCA ----
lines_local_data <- lines_modalities %>% select( c("category","num","type","treatment","hedgerow","broadleaved_class",
                                                   "length","bl_buff100","pine_height_mean","open_width","T_moy",
                                                   "T_max","T_min","H_moy","H_max","H_min","herb_cover",
                                                   "shrub_cover","tree_cover","total_cover","tree_H"))  %>%
  filter(category != "broadleaved_forest") %>% #because small scale data were never collected for those sites
  tibble::column_to_rownames(var = "num") 

#Selection of local scale variables
lines_local_globalvar <- c( "length","bl_buff100","pine_height_mean","open_width","T_moy",
                            "T_max","T_min","H_moy","H_max","H_min","herb_cover",
                            "shrub_cover","tree_cover","total_cover","tree_H")

lines_local_climvar <- c("T_moy", "T_max","T_min","H_moy","H_max","H_min")

lines_local_covervar <- c( "herb_cover", "shrub_cover","tree_cover", "total_cover")

#Choose variable for pca
var_pca <- lines_local_globalvar

#Transform as numeric and delete rows with NA 
lines_local_data <- lines_local_data %>%
  mutate(across(all_of(var_pca), as.numeric))

lines_local_data <- lines_local_data %>%
  filter(if_all(all_of(var_pca), ~ !is.na(.)))

#REMOVING SITE NUMBER 8 BECAUSE OF EFFECT  BEWAREBEWAREBEWAREBEWAREBEWAREBEWAREBEWAREBEWAREBEWAREBEWARE
lines_local_data <- lines_local_data %>%
  .[rownames(lines_local_data) != "8", ]

#Looking for variables correlations  
cor_matrix <- cor(lines_local_data[, var_pca],
                  use = "pairwise.complete.obs", method = "spearman")
print(cor_matrix)
corrplot::corrplot(cor_matrix, method = "color", 
                   type = "upper", 
                   order="hclust",
                   tl.col = "black",
                   tl.srt = 45, 
                   diag = FALSE)

PerformanceAnalytics::chart.Correlation(lines_local_data[, var_pca],
                                        method = "spearman",
                                        histogram = TRUE, 
                                        pch=19)

#Generate PCA
PCA_lines_local <- FactoMineR::PCA(X = lines_local_data[, var_pca],
                                   graph = F ,
                                   scale.unit = T)  

factoextra::fviz_screeplot(PCA_lines_local, addlabels = TRUE)
PCA_lines_local$eig 

factoextra::fviz_pca_var(PCA_lines_local, axe=c(1,2),
                         col.var = "cos2",                                    
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel=TRUE)


var_groupe <- lines_local_data$category 
factoextra::fviz_pca_ind(PCA_lines_local, 
                         axe = c(1,2),
                         label = "ind",
                         habillage = factor(var_groupe),
                         labelsize = 2,
                         palette = (viridis::viridis_pal(option = "viridis", begin = 0.1, end = 0.9)(length(unique(var_groupe)))),
                         addEllipses = F,
                         ggtheme = ggplot2::theme_minimal())


#Extract PCA coordinates values in a new dataframe used to extract multivariate axis results
multivariate_results <- lines_local_data %>%
  tibble::rownames_to_column(var = "num") %>%
  left_join(
    data.frame(PCA_lines_local_axis1 = PCA_lines_local$ind$coord[,1],
               PCA_lines_local_axis2 = PCA_lines_local$ind$coord[,2], 
               line = as.character(rownames(PCA_lines_local$ind$coord))),
    by = c("num" = "line"))


## Test LDA : what variables (used in PCA) segregates the best our a priori groups -----
LDA <- MASS::lda(lines_local_data[, var_pca], lines_local_data$category)

lda_pred <- predict(LDA)
ggplot(data = as.data.frame(lda_pred$x), aes(x = LD1, y = 1, color = lines_local_data$category)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "LDA projection", x = "LD1", y = "LD2", color = "Treatment")

# Inspect the coefficients
LDA$scaling

scaling_df <- as.data.frame(LDA$scaling)
scaling_df$variable <- rownames(scaling_df)

ggplot(scaling_df, aes(x = reorder(variable, LD1), y = LD1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable contributions to LD1", x = "", y = "Coefficient")






# Nattan's flora 2023 data ----

#splist contient plus d'especes différentes que transect et quadrat et toutes les especes qu'ils contiennent (après correction brut data)
setdiff(quadrat$code_esp, splist$code)
setdiff(transect$genre_esp, splist$code)
setdiff(splist$code, quadrat$code_esp)
setdiff(splist$code, transect$genre_esp)
#transect contient plus d'especes que quadrat, toutes les especes de quadrat sont contenues dans transect
setdiff(quadrat$code_esp, transect$genre_esp)
setdiff(transect$genre_esp, quadrat$code_esp)



## Quadrat data : plant species abundances ----

#Add number per quadrat (alternative would have been using mean)
quadrat_summed <- quadrat %>%
  group_by(stand, code_esp) %>%
  summarise(
    across(
      all_of(c("tree", "shrub", "herbaceous")),
      \(x) sum(x, na.rm = TRUE)),
    .groups = "drop")

#Join to species list to obtain species name
quadrat_summed <- quadrat_summed %>%
  left_join(splist,
            by = c("code_esp" = "code"))

#Tidy dataframe by pivoting wide to long
quadrat_summed <- quadrat_summed %>%
  tidyr::pivot_longer(
    cols = c("tree", "shrub", "herbaceous"),
    names_to = "plant_type",
    values_to = "occurence"
  ) %>%
  filter(occurence > 0)

#Generate dataframe with species richness per plant type
quadrat_summed_rich_type <- quadrat_summed %>%
  group_by(stand, plant_type) %>%
  summarise(richness = length(unique(code_esp)),
            .groups = "drop")

quadrat_summed_rich_type_wide <- quadrat_summed_rich_type %>%
  tidyr::pivot_wider(
    names_from = plant_type,
    values_from = richness,
    values_fill = 0  
  )

#Ignore plant_type and group by species only
quadrat_summed_sp <- quadrat_summed %>%
  group_by(stand, code_esp, species) %>%
  summarise(
    across(
      all_of("occurence"),
      \(x) sum(x, na.rm = TRUE)),
    .groups = "drop")

#Generate total species richness and join dataframe with graph containing plant_type richness
quadrat_summed_rich_all <- quadrat_summed_sp %>%
  group_by(stand) %>%
  summarise(plant_richness = length(unique(code_esp)),
            .groups = "drop") %>%
  left_join(quadrat_summed_rich_type_wide)


### Multivariate exploration ----

#Pivot species only (no plant type) dataframe for use in multivariate analysis
commu_quadrat_summed <- quadrat_summed_sp %>%
  select(!species) %>%
  tidyr::pivot_wider(
    names_from = code_esp,
    values_from = occurence
  ) %>%
  mutate(across(all_of(colnames(.)), ~ if_else(is.na(.x), 0, .x)))

#Extract plant species name
plant_species <- unique(quadrat_summed_sp$code_esp)

#Join lines category to dataframe
commu_quadrat_summed <- commu_quadrat_summed %>%
  left_join(lines_modalities %>%
              select(num, category, treatment, broadleaved_class),
            by = c( "stand" = "num")
            )

#PCA

#Looking for variables correlations  
cor_matrix <- cor(commu_quadrat_summed[, plant_species],
                  use = "pairwise.complete.obs", method = "spearman")
print(cor_matrix)
corrplot::corrplot(cor_matrix, method = "color", 
                   type = "upper", 
                   order="hclust",
                   tl.col = "black",
                   tl.srt = 45, 
                   diag = FALSE)

#Generate PCA
PCA_plant_species <- FactoMineR::PCA(X = commu_quadrat_summed[, plant_species],
                                   graph = F ,
                                   scale.unit = T)  
var_groupe <- commu_quadrat_summed$category 



factoextra::fviz_screeplot(PCA_plant_species, addlabels = TRUE)
PCA_plant_species$eig 

factoextra::fviz_pca_var(PCA_plant_species, axe=c(1,2),
                         col.var = "cos2",                                    
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel=TRUE)

factoextra::fviz_pca_ind(PCA_plant_species, 
                         axe = c(1,2),
                         label = "ind",
                         habillage = factor(var_groupe),
                         labelsize = 2,
                         palette = (viridis::viridis_pal(option = "viridis", begin = 0.1, end = 0.9)(length(unique(var_groupe)))),
                         addEllipses = T,
                         ggtheme = ggplot2::theme_minimal())


#Extract PCA coordinates values and add to multivariate dataframe
multivariate_results <- multivariate_results %>%
  left_join(
    data.frame(PCA_plant_species_axis1 = PCA_plant_species$ind$coord[,1],
               PCA_plant_species_axis2 = PCA_plant_species$ind$coord[,2], 
               line = as.character(rownames(PCA_plant_species$ind$coord))),
    by = c("num" = "line"))



#nMDS
nmds_bray <- vegan::metaMDS(
  comm = commu_quadrat_summed[, plant_species],
  distance = "bray",
  k = 2,
  autotransform = FALSE
)
cat(paste0("Final NMDS has a stress of ", round(nmds_bray$stress, 3), "\n"))
vegan::stressplot(nmds_bray)
# Extract scores
nmdspoint <- vegan::scores(nmds_bray)$sites %>%
  as_tibble(rownames = "stand")
nmdsvariable <- vegan::scores(nmds_bray)$species %>%
  as_tibble(rownames = "code_esp")

commu_quadrat_summed$stand <- as.character(commu_quadrat_summed$stand)

nmdspoint %>% 
  left_join(commu_quadrat_summed) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(treatment))) +
  geom_point() +
  stat_ellipse(show.legend = FALSE) +
  geom_text(data = nmdsvariable, aes(x = NMDS1, y = NMDS2, label = code_esp), colour = "grey20") +
  theme_minimal()


## Transect data : plant species presence/absence ----

### Multivariate exploration ----

#Generate presence/absence community data
commu_transect <- transect %>%
  distinct(stand, genre_esp) %>%
  mutate(presence = 1) %>%
  tidyr::pivot_wider(
    names_from = genre_esp,
    values_from = presence,
    values_fill = 0                      
  )

#Extract plant species name
plant_species2 <- unique(transect$genre_esp)

#Join lines category to dataframe
commu_transect <- commu_transect %>%
  left_join(lines_modalities %>%
              select(num, category, treatment, broadleaved_class),
            by = c( "stand" = "num")
  )



#AFC
CA_plant_species = FactoMineR::CA(commu_transect[, plant_species2],
                                  graph = F)

var_groupe <- commu_quadrat_summed$category 



factoextra::fviz_screeplot(CA_plant_species, addlabels = TRUE)
CA_plant_species$eig 

factoextra::fviz_ca_biplot(CA_plant_species,
                           label ="row")

factoextra::fviz_ca_col(CA_plant_species, axe=c(1,2),
                         col.var = "cos2",                                    
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel=TRUE)

factoextra::fviz_ca_row(CA_plant_species, 
                         axe = c(1,2),
                         label = "row",
                         habillage = factor(var_groupe),
                         labelsize = 2,
                         palette = (viridis::viridis_pal(option = "viridis", begin = 0.1, end = 0.9)(length(unique(var_groupe)))),
                         addEllipses = T,
                         ggtheme = ggplot2::theme_minimal())

#Extract PCA coordinates values and add to multivariate dataframe
multivariate_results <- multivariate_results %>%
  left_join(
    data.frame(CA_plant_species_axis1 = CA_plant_species$ind$coord[,1],
               CA_plant_speciess_axis2 = CA_plant_species$ind$coord[,2], 
               line = as.character(rownames(CA_plant_species$ind$coord))),
    by = c("num" = "line"))


#nMDS (Jaccard)
nmds_jaccard <- vegan::metaMDS(
  comm = commu_transect[, plant_species2],
  distance = "jaccard",
  k = 2,
  autotransform = FALSE
)
cat(paste0("Final NMDS has a stress of ", round(nmds_jaccard$stress, 3), "\n"))
vegan::stressplot(nmds_jaccard)
# Extract scores
nmdspoint <- vegan::scores(nmds_jaccard)$sites %>%
  as_tibble(rownames = "stand")
nmdsvariable <- vegan::scores(nmds_jaccard)$species %>%
  as_tibble(rownames = "code_esp")

commu_transect$stand <- as.character(commu_transect$stand)

nmdspoint %>% 
  left_join(commu_transect) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = as.factor(category))) +
  geom_point() +
  stat_ellipse(show.legend = FALSE) +
  geom_text(data = nmdsvariable, aes(x = NMDS1, y = NMDS2, label = code_esp), colour = "grey20") +
  theme_minimal()





# ??? Analyse fonctionnelles ----



## Cover conversion script (from Nattan) ----
#prep
strate$num <- as.character(strate$stand)
#filter strata
strate <- strate %>% filter(strata %in% c("tree", "shrub", "herbaceous"))

#conversion Domin scale into mean cover percentage
strate$cover <- with(strate, 
                     ifelse(
                       domin_scale == 0, 0,
                       ifelse(
                         domin_scale == 1, 0.1,
                         ifelse(
                           domin_scale == 2, 0.5,
                           ifelse(
                             domin_scale == 3, 1,
                             ifelse(
                               domin_scale == 4, 7,
                               ifelse(
                                 domin_scale == 5, 13,
                                 ifelse(
                                   domin_scale == 6, 29.5,
                                   ifelse(
                                     domin_scale == 7, 42,
                                     ifelse(
                                       domin_scale == 8, 63,
                                       ifelse(
                                         domin_scale == 9, 83,
                                         ifelse(
                                           domin_scale == 10, 95.5, NA
                                         )
                                       )
                                     )
                                   )
                                 )
                               )
                             )
                           )
                         )
                       )))


#mean per stand
strate <- strate %>%
  group_by(num, strata) %>%
  summarise(m_cover = mean(cover, na.rm = TRUE))

#pivot
strate <- strate %>%
  tidyr::pivot_wider(names_from = strata, values_from = m_cover)

#cover density
strate$total_cover <- rowSums(strate[,c(2:4)])
#strate$SDHI_str<-diversity(strate[,c(2:4)])
strate <- strate %>% rename(herb_cover = herbaceous)
strate <- strate %>% rename(shrub_cover = shrub)
strate <- strate %>% rename(tree_cover = tree)


# Links between multivariate axis ----

#Correlations






