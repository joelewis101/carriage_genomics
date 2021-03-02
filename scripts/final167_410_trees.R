# final st167/410 trees

# st167 tree

library(ggtree)
library(ape)
library(phytools)
library(tidyverse)
library(blantyreESBL)
library(lubridate)
library(ggnewscale)
library(countrycode)
library(viridis)

st167_metadata <- read_tsv("data_raw/e_coli/st167.tsv")



st167_metadata %>%
  select(
    Uberstrain,
    Name,
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    `Source Niche`,
    `Source Details`,
    Country,
    `Collection Year`,
    ST
  ) %>%  separate_rows(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    sep = ","
  ) %>%
  separate(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    into = c(
      "accession",
      "platform",
      "library",
      "insert_size",
      "experiment"
    ),
    sep = ";"
  ) ->
  st167_metadata


read.tree("data_processed/esco/st167/iqtree/clean.full.filtered_polymorphic_sites.ref_removed.snpsites.fasta.treefile") -> st167_tree

midpoint.root(st167_tree) -> st167_tree

st167_tree$tip.label <- gsub("_filtered","", st167_tree$tip.label)

st167_metadata %>% 
  filter(platform == "ILLUMINA") %>% 
  filter(accession %in% st167_tree$tip.label) %>% 
  bind_rows(
    btESBL_sequence_sample_metadata %>% 
      filter(lane %in% st167_tree$tip.label) %>% 
      transmute(accession = lane,
                `Collection Year` = year(data_date),
                `Source Niche` = "Human",
                Country = "Malawi")
  ) %>% 
  rename(Year = `Collection Year`) %>% 
  as.data.frame() -> st167_metadata

st167_metadata$malawi <- st167_metadata$Country == "Malawi"

st167_metadata$Continent <- 
  countrycode(
    sourcevar = st167_metadata$Country,
    origin = "country.name",
    destination = "un.region.name"
  )

rownames(st167_metadata) <- st167_metadata$accession

st167_metadata %>% 
  select(accession,Country) %>% 
  pivot_wider(id_cols = accession, names_from = Country, values_from = Country,
              values_fn = length, values_fill = 0) %>% 
  mutate(across(everything(), as.character)) %>% 
  as.data.frame() -> country.ohehot

rownames(country.ohehot) <- country.ohehot$accession

amr.ariba <- read_csv("data_processed/esco/st167/ariba/st167_ariba_srst2_summary.csv")

amr.ariba %>% 
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered", "", name),
         name = gsub("#","_", name)) %>% 
  pivot_longer(-name, 
               names_to= c( "cluster", ".value"), 
               names_sep = "\\.") %>% 
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>% 
  filter(match == "yes") %>% 
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>% 
  select(name, gene) %>% 
  pivot_wider(names_from = gene, 
              values_from = gene,
              values_fn = length,
              values_fill = 0) %>%
  mutate(across(everything(), as.character)) %>% 
  as.data.frame() ->
  amr.ariba

rownames(amr.ariba) <- amr.ariba$name


bls <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab"
)


bls %>% 
  rename_with(~ tolower(gsub(" ", "_", .x))) %>% 
  rename_with(~ gsub("#", "", .x)) %>% 
  mutate(allele_name = gsub("-", "_", allele_name),
         class = case_when(
           grepl("carbapenem-hydrolyzing", 
                 curated_gene_product_name) ~ "Carbapenemase",
           grepl("metallo-beta-lactamase", 
                 curated_gene_product_name) ~ "Carbapenemase",
           grepl("extended-spectrum beta-lactamase", 
                 curated_gene_product_name) ~ "ESBL",
           grepl("class C", 
                 curated_gene_product_name) ~ "AmpC",
           grepl("beta-lactamase", 
                 curated_gene_product_name) ~ "Penicillinase"
         )) -> 
  bls

amr.ariba %>% 
  rename(sample = name) %>% 
  pivot_longer(-sample) %>% 
  filter(value == 1) %>% 
  left_join(
    select(bls, allele_name, class),
    by = c("name" = "allele_name")
  ) %>% 
  select(sample, class) %>% 
  filter(!is.na(class)) %>% 
  pivot_wider(id_cols = sample,
              names_from = class,
              values_from = class,
              values_fn = length,
              values_fill = 0) %>% 
  mutate(across(where(is.numeric), function(x) if_else(x > 0,1,0))) %>% 
  mutate(across(where(is.numeric), as.character)) %>% 
  as.data.frame() -> bl.onehot

rownames(bl.onehot) <- bl.onehot$sample  

amr.ariba %>% 
  rename(sample = name) %>% 
  pivot_longer(-sample) %>% 
  filter(value == 1) %>% 
  left_join(
    select(bls, allele_name, class),
    by = c("name" = "allele_name")
  ) %>% 
  filter(class %in% c("ESBL", "Carbapenemase") & value == 1) %>%
  mutate(name = gsub("_", "-", name)) %>% 
  arrange(class, fct_infreq(name)) %>% 
  pivot_wider(id_cols = "sample",
              names_from = name,
              values_from = class) %>% 
  as.data.frame() %>% 
  select(-`TEM-15`) ->    # not in df?
  esbl_cpe.ariba

rownames(esbl_cpe.ariba) <- esbl_cpe.ariba$sample

(
(ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %<+% (select(st167_metadata, accession, Country) %>% mutate(Country = if_else(Country == "Malawi", "Malawi",NA_character_)))  %>%
    gheatmap(
      select(st167_metadata, `Year`),
      font.size = 4,
      width = 0.03,
      colnames_position = "top",
      color = NA,
      colnames_offset_y = 5,
      colnames_angle = 90,
      hjust = 0
    )  + 
    scale_fill_viridis(name = "Year") +
    new_scale_fill()) %>% 
  gheatmap(
    select(st167_metadata, Continent),
    font.size = 4,
    width = 0.03,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00003,
    hjust = 0
  ) + 
  scale_fill_viridis_d(option = "plasma", name = "Continent") + 
  geom_tippoint(aes(color = Country), na.rm = TRUE, size = 0.5) + 
  scale_color_manual(values  = "red", na.translate = FALSE) +
    new_scale_fill() 
) %>% 
  gheatmap(
  select(esbl_cpe.ariba, -sample),
  font.size = 3,
  width = 0.5,
  colnames_position = "top",
  color = "darkgrey",
  colnames_offset_y = 5,
  colnames_angle = 90,
  offset = 0.00006,
  hjust = 0
) + scale_fill_manual(values = viridis_pal(option = "magma")(5)[c(2,3)],
                      na.translate = FALSE,
                      name = "Beta\nLactamase") +
  ylim(NA, 350) +
  geom_treescale(x = 0.0002, y = 250, offset = 2)


# st410 tree ------------



st410_metadata <- read_tsv("data_raw/e_coli/st410.tsv")



st410_metadata %>%
  select(
    Uberstrain,
    Name,
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    `Source Niche`,
    `Source Details`,
    Country,
    `Collection Year`,
    ST
  ) %>%  separate_rows(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    sep = ","
  ) %>%
  separate(
    `Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Status)`,
    into = c(
      "accession",
      "platform",
      "library",
      "insert_size",
      "experiment"
    ),
    sep = ";"
  ) ->
  st410_metadata


read.tree("data_processed/esco/st410/iqtree/clean_full.filtered_pollymorphic_sites.ref_removed.snpsites.fasta.treefile") -> st410_tree

midpoint.root(st410_tree) -> st410_tree

#ggtree(treeio::tree_subset (st410_tree, 1083, levels_back = 0)) + geom_treescale()

st410_tree$tip.label <- gsub("_filtered","", st410_tree$tip.label)

st410_metadata %>% 
  filter(platform == "ILLUMINA") %>% 
  filter(accession %in% st410_tree$tip.label) %>% 
  bind_rows(
    btESBL_sequence_sample_metadata %>% 
      filter(lane %in% st410_tree$tip.label) %>% 
      transmute(accession = lane,
                `Collection Year` = year(data_date),
                `Source Niche` = "Human",
                Country = "Malawi")
  ) %>% 
  rename(
    Year = `Collection Year`) %>% 
  as.data.frame() -> st410_metadata

st410_metadata$malawi <- st410_metadata$Country == "Malawi"


st410_metadata$Continent <- 
  countrycode(
    sourcevar = st410_metadata$Country,
    origin = "country.name",
    destination = "un.region.name"
  )

rownames(st410_metadata) <- st410_metadata$accession

st410_metadata %>% 
  select(accession,Country) %>% 
  pivot_wider(id_cols = accession, names_from = Country, values_from = Country,
              values_fn = length, values_fill = 0) %>% 
  mutate(across(everything(), as.character)) %>% 
  as.data.frame() -> country.ohehot410

rownames(country.ohehot410) <- country.ohehot410$accession

amr.ariba410 <- read_csv("data_processed/esco/st410/ariba/st410_ariba_srst2_summary.csv")

amr.ariba410 %>% 
  mutate(name = gsub("\\./", "", name),
         name = gsub("/report.tsv", "", name),
         name = gsub("_filtered", "", name),
         name = gsub("#","_", name)) %>% 
  pivot_longer(-name, 
               names_to= c( "cluster", ".value"), 
               names_sep = "\\.") %>% 
  mutate(gene = sapply(str_split(ref_seq, "__"), function(x) x[3])) %>% 
  filter(match == "yes") %>% 
  mutate(gene = case_when(
    gene == "TEM_95" ~ "TEM_1",
    TRUE ~ gene
  )) %>% 
  select(name, gene) %>% 
  pivot_wider(names_from = gene, 
              values_from = gene,
              values_fn = length,
              values_fill = 0) %>%
  mutate(across(everything(), as.character)) %>% 
  as.data.frame() ->
  amr.ariba410

rownames(amr.ariba410) <- amr.ariba410$name


amr.ariba410 %>% 
  rename(sample = name) %>% 
  pivot_longer(-sample) %>% 
  filter(value == 1) %>% 
  left_join(
    select(bls, allele_name, class),
    by = c("name" = "allele_name")
  ) %>% 
  select(sample, class) %>% 
  filter(!is.na(class)) %>% 
  pivot_wider(id_cols = sample,
              names_from = class,
              values_from = class,
              values_fn = length,
              values_fill = 0) %>% 
  mutate(across(where(is.numeric), function(x) if_else(x > 0,1,0))) %>% 
  mutate(across(where(is.numeric), as.character)) %>% 
  as.data.frame() -> bl410.onehot

rownames(bl410.onehot) <- bl410.onehot$sample  


amr.ariba410 %>% 
  rename(sample = name) %>% 
  pivot_longer(-sample) %>% 
  filter(value == 1) %>% 
  left_join(
    select(bls, allele_name, class),
    by = c("name" = "allele_name")
  ) %>% 
  filter(class %in% c("ESBL", "Carbapenemase") & value == 1) %>%
  mutate(name = gsub("_", "-", name)) %>% 
  arrange(class, fct_infreq(name)) %>% 
  pivot_wider(id_cols = "sample",
              names_from = name,
              values_from = class) %>% 
  as.data.frame() -> 
#  select(-`TEM-15`) ->    # not in df?
  esbl_cpe.ariba410

rownames(esbl_cpe.ariba410) <- esbl_cpe.ariba410$sample

(  
  (ggtree(treeio::tree_subset (st410_tree, 1083, levels_back = 0)) %<+% (select(st410_metadata, accession, Country) %>% mutate(Country = if_else(Country == "Malawi", "Malawi",NA_character_)))  %>%
     gheatmap(
       select(st410_metadata, `Year`),
       font.size = 4,
       width = 0.03,
       colnames_position = "top",
       color = NA,
       colnames_offset_y = 5,
       colnames_angle = 90,
       hjust = 0
     )  + 
     scale_fill_viridis(name = "Year") +
     new_scale_fill()) %>% 
  gheatmap(
    select(st410_metadata, Continent),
    font.size = 4,
    width = 0.03,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00015,
    hjust = 0
  ) + 
  scale_fill_viridis_d(option = "plasma", name = "Continent") + 
  geom_tippoint(aes(color = Country), na.rm = TRUE, size = 0.5) + 
  scale_color_manual(values  = "red", na.translate = FALSE) +
  new_scale_fill() 
) %>% 
  gheatmap(
    select(esbl_cpe.ariba410, -sample),
    font.size = 3,
    width = 0.5,
    colnames_position = "top",
    color = "darkgrey",
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.0003,
    hjust = 0
  ) + scale_fill_manual(values = viridis_pal(option = "magma")(5)[c(2,3)],
                        na.translate = FALSE,
                        name = "Beta\nLactamase") +
  ylim(NA, 600) +
  geom_treescale(x = 0.0002, y = 250, offset = 2) +
   geom_highlight(aes(node = 568) ,alpha = 0.3) +
  geom_cladelabel(
    node = 568,
    label = "B4/H24RxC",
    angle = 90,
    offset = -0.0006,
    barsize = 0,
    hjust = 0.5,
    fontsize = 3,
    color = "steelblue"
  )

treeio::tree_subset (st410_tree, 564, levels_back = 0) -> st410_subtree

esbl_cpe.ariba410 %>%  
  subset(sample %in% st410_subtree$tip.label) %>% 
  as.data.frame() ->
  esbl_cpe.ariba410_subset
  
  esbl_cpe.ariba410_subset[
    c(TRUE, apply(esbl_cpe.ariba410_subset[-1], 2,function(x) !all(is.na(x))))] ->
    esbl_cpe.ariba410_subset

  (  
  (ggtree(treeio::tree_subset (st410_tree, 564, levels_back = 0)) %<+% (select(st410_metadata, accession, Country) %>% mutate(Country = if_else(Country == "Malawi", "Malawi",NA_character_)))  %>%
     gheatmap(
       select(st410_metadata, `Year`),
       font.size = 4,
       width = 0.05,
       colnames_position = "top",
       color = NA,
       colnames_offset_y = 5,
       colnames_angle = 90,
       hjust = 0
     )  + 
     scale_fill_viridis(name = "Year", limits = c(1975,2020)) +
     new_scale_fill()) %>% 
    gheatmap(
      select(st410_metadata, Continent),
      font.size = 4,
      width = 0.05,
      colnames_position = "top",
      color = NA,
      colnames_offset_y = 5,
      colnames_angle = 90,
      offset = 0.00003,
      hjust = 0
    ) + 
    scale_fill_viridis_d(option = "plasma", name = "Continent") + 
    geom_tippoint(aes(color = Country), na.rm = TRUE, size = 0.5) + 
    scale_color_manual(values  = "red", na.translate = FALSE) +
    new_scale_fill() 
) %>% 
  gheatmap(
    select(esbl_cpe.ariba410_subset, -sample),
    font.size = 3,
    width = 0.5,
    colnames_position = "top",
    color = "darkgrey",
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00006,
    hjust = 0
  ) + scale_fill_manual(values = viridis_pal(option = "magma")(5)[c(2,3)],
                        na.translate = FALSE,
                        name = "Beta\nLactamase") +
  ylim(NA, 200) +
  geom_treescale(x = 0.0002, y = 250, offset = 2) +
    geom_highlight(aes(node = 181) ,alpha = 0.3) +
    geom_cladelabel(node = 181,label = "B4/H24RxC", angle = 90, offset = -0.00028, barsize = 0,hjust = 0.5, fontsize = 4, color = "steelblue")
  