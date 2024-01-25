## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(neonPlantEcology)
library(vegan)
library(sf)
library(ggpubr)

## ----getdata------------------------------------------------------------------
# D14 <- npe_site_ids(domain = "D14") |>
#   npe_download_plant_div(sites = .)

data("D14")


## ----wrangle------------------------------------------------------------------
comm <- npe_community_matrix(D14)
comm[1:4, 1:4]

data("plot_centroids")
plot_centroids <- sf::st_set_geometry(plot_centroids, NULL)
metadata <- npe_cm_metadata(comm) |>
  dplyr::left_join(plot_centroids)



## ----sac----------------------------------------------------------------------

# checking the metadata file to see where JORN stops and SRER begins
which(metadata$site == "JORN") |> range()
which(metadata$site == "SRER") |> range()

# performing and plotting the species accumulation curve analysis
sp_jorn <- vegan::specaccum(comm[1:152,])
sp_srer <- vegan::specaccum(comm[152:347,])
plot(sp_srer, col = "blue")
plot(sp_jorn, col = "gold", add=T)

## ----nmds, fig.width=7.5, fig.height=6----------------------------------------
nmds <- metaMDS(comm,trace = F)

nmds_sites <- nmds$points |>
  as_tibble(rownames = "rowname") |>
  left_join(metadata)

ggplot(nmds_sites, aes(x=MDS1, y=MDS2, color = eventID, size = elevation, shape = site)) +
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette = "Set1") +
  theme(panel.background = element_rect(fill=NA, color= "black"))


## ----nspp---------------------------------------------------------------------
di <- npe_summary(D14, scale = "site", timescale = "all")

dj<- di |>
  dplyr::select(site, eventID, starts_with("nspp")) |>
  tidyr::pivot_longer(cols = starts_with("nspp")) |>
  dplyr::filter(!name %in% c("nspp_notexotic", "nspp_total")) |>
  dplyr::transmute(origin = str_remove_all(name, "nspp_"),
            nspp = value,
            site=site)

p1<-ggplot(dj, aes(x=nspp, y=site, fill = origin)) +
  geom_bar(stat = "identity", color = "black", lwd = .2) +
  xlab("Species Richness") +
  ylab("Site")


## ----rc-----------------------------------------------------------------------
 dk<- di |>
  dplyr::select(site, eventID, starts_with("rel_cover")) |>
  pivot_longer(cols = starts_with("rel_cover")) |>
  filter(!name %in% c("rel_cover_notexotic", "rel_cover_total")) |>
  transmute(origin = str_remove_all(name, "rel_cover_"),
            relative_cover = value,
            site=site)

p2<- ggplot(dk, aes(x=relative_cover, y=site, fill = origin)) +
  geom_bar(stat = "identity", color = "black", lwd = .2) +
  xlab("Relative Cover")

## ----alphabeta----------------------------------------------------------------
dl <- di|>
  dplyr::select(site, eventID, starts_with("shannon")) |>
  pivot_longer(cols = starts_with("shannon")) |>
  filter(!name %in% c("shannon_notexotic", "shannon_total", "shannon_family")) |>
  transmute(origin = str_remove_all(name, "shannon_"),
            shannon = value,
            site=site)

p3<- ggplot(dl, aes(x = shannon, y = site, fill = origin)) +
  geom_bar(stat = 'identity', color = "black", lwd = .2) +
  xlab("Shannon-Weaver Diveristy")


## ----multipanel---------------------------------------------------------------
ggpubr::ggarrange(p1, 
                  p2 + theme(axis.title.y = element_blank(), 
                                 axis.text.y = element_blank(), 
                                 axis.ticks.y = element_blank()), 
                  p3 + theme(axis.title.y = element_blank(), 
                                 axis.text.y = element_blank(), 
                                 axis.ticks.y = element_blank()), 
                  nrow=1, ncol=3, 
                  common.legend = TRUE, 
                  widths = c(1.35,1,1))



## ----family-------------------------------------------------------------------
data("D14")
lf <- npe_longform(D14, scale = "plot", timescale = "annual")

lf_f <- lf |>
  mutate(family = ifelse(!family %in% c("Poaceae", "Fabaceae"), "Other", family)) |>
  group_by(site, plotID, eventID, family) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()
  
ggplot(lf_f, aes(x=eventID, y=cover, fill = family)) +
  geom_boxplot(position="dodge") +
  facet_wrap(~site, scales = "free_x") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", name = "Family") +
  ylab("Cover") +
  xlab("eventID (Annual Timestep)")


