 # st167 tree

library(ggtree)
library(ape)
library(phytools)
library(tidyverse)
library(blantyreESBL)
library(lubridate)
library(ggnewscale)

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
  as.data.frame() -> st167_metadata

st167_metadata$malawi <- st167_metadata$Country == "Malawi"

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
                 curated_gene_product_name) ~ "CPE",
           grepl("metallo-beta-lactamase", 
                 curated_gene_product_name) ~ "CPE",
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
  filter(class %in% c("ESBL", "CPE") & value == 1) %>% 
  arrange(class, fct_infreq(name)) %>% 
  pivot_wider(id_cols = "sample",
              names_from = name,
              values_from = class) %>% 
  as.data.frame() ->
  esbl_cpe.ariba

rownames(esbl_cpe.ariba) <- esbl_cpe.ariba$sample

ggtree(treeio::tree_subset (st167_tree, 321, levels_back = 0)) %>%
  gheatmap(
    select(country.ohehot, -accession),
    font.size = 3,
    width = 2,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00005
  ) %>% 
  gheatmap(
    select(bl.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
  ) + geom_treescale(x = 0.0002, y = 50, offset = 2) +
  ylim(NA, 80)

  
ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
  select(st167_metadata, malawi),
  font.size = 3,
  width = 0.05,
  colnames_position = "top",
  color = NA,
  colnames_offset_y = 5,
  colnames_angle = 90
) %>% 
  gheatmap(
    select(bl.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00005
  )  + geom_treescale()
  
ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
    select(esbl_cpe.ariba, -sample),
    font.size = 3,
    width = 0.3,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00002
  ) %>%
  gheatmap(
    select(st167_metadata, malawi),
    font.size = 3,
    width = 0.02,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0
  )   %>% 
  gheatmap(
    select(st167_metadata, plam),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00011
  )  +
  ylim(NA, 340) + geom_treescale()


ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
    select(amr.ariba, -name),
    font.size = 3,
    width = 2,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00003
  ) %>%
  gheatmap(
    select(st167_metadata, malawi),
    font.size = 3,
    width = 0.02,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0
  ) +
  ylim(NA, 340) + geom_treescale()


read_csv("data_processed/esco/st167/snpdists/st167_snpdists.csv") -> st167_snpdist

st167_snpdist %>%  
  filter(`28099_1#102` > 10000) %>% 
  pull(`snp-dists 0.6.2`)

st167_snpdist %>% 
  filter(`snp-dists 0.6.2` != "SRR11038969_filtered") %>% 
  select(-SRR11038969_filtered) %>% 
  rename_with(function(x) gsub("_filtered","",x)) %>% 
  rename_with(function(x) gsub("#","_",x)) %>% 
  mutate(`snp-dists 0.6.2` = gsub("_filtered","",`snp-dists 0.6.2`),
         `snp-dists 0.6.2` = gsub("#","_",`snp-dists 0.6.2`)) %>% 
  as.data.frame() -> st167_snpdist
 

rownames(st167_snpdist) <- st167_snpdist$`snp-dists 0.6.2`

treeio::tree_subset (st167_tree, 321, levels_back = 0) -> subtree

keeperz <- c("snp-dists 0.6.2", subtree$tip.label)

st167_snpdist %>% 
  dplyr::select(!!!keeperz) %>% 
  filter(`snp-dists 0.6.2` %in% subtree$tip.label) %>% 
  as.data.frame() ->
  subtree.snpdists

rownames(subtree.snpdists) <- subtree.snpdists$`snp-dists 0.6.2`
 
pheatmap(subtree.snpdists[-1],color = viridis(100), 
         display_numbers = TRUE,
         number_format = "%i")

# plasmids

dep_plasm <- read_csv("data_processed/esco/st167/plasmid_mapping/depth_and_cov.csv")
dep_plasm$mean_dep <- dep_plasm$depth_sum /dep_plasm$ref_bases
dep_plasm$cov <- dep_plasm$bases_mapped/ dep_plasm$ref_bases

left_join(
  st167_metadata,
  dep_plasm %>% 
    mutate(lane = gsub("\\./", "", lane),
           lane = gsub("/report.tsv", "", lane),
           lane = gsub("_filtered", "", lane),
           lane = gsub("#","_", lane),
           lane = gsub("/snps.bam", "", lane)) %>% 
    select(lane, cov, mean_dep),
  by = c("accession" = "lane")
) -> st167_metadata

st167_metadata$plam <- st167_metadata$cov > 0.99
 
st167_metadata <- as.data.frame(st167_metadata)
st167_metadata$accession -> rownames(st167_metadata)


ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
    select(esbl_cpe.ariba, -sample),
    font.size = 3,
    width = 0.3,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00004
  ) %>%
  gheatmap(
    select(st167_metadata, malawi),
    font.size = 3,
    width = 0.02,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0
  )   %>% 
  gheatmap(
    select(st167_metadata, plam),
    font.size = 3,
    width = 0.02,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00002
  )  +
  ylim(NA, 340) + geom_treescale()


ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
    select(st167_metadata, malawi),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90
  ) %>% 
  gheatmap(
    select(bl.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00006
  )    %>% 
  gheatmap(
    select(st167_metadata, plam),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00011
  ) %>%   gheatmap(
    select(country.ohehot, -accession),
    font.size = 3,
    width = 2,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00015
  ) +
  ylim(NA, 330)
 
##### ST410 -------------------


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

ggtree(treeio::tree_subset (st410_tree, 1083, levels_back = 0)) + geom_treescale()

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
  as.data.frame() -> st410_metadata

st410_metadata$malawi <- st410_metadata$Country == "Malawi"

rownames(st410_metadata) <- st410_metadata$accession

st410_metadata %>% 
  select(accession,Country) %>% 
  pivot_wider(id_cols = accession, names_from = Country, values_from = Country,
              values_fn = length, values_fill = 0) %>% 
  mutate(across(everything(), as.character)) %>% 
  as.data.frame() -> country.ohehot

rownames(country.ohehot) <- country.ohehot$accession

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


bls <- read_tsv(
  "https://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab"
)


bls %>% 
  rename_with(~ tolower(gsub(" ", "_", .x))) %>% 
  rename_with(~ gsub("#", "", .x)) %>% 
  mutate(allele_name = gsub("-", "_", allele_name),
         class = case_when(
           grepl("carbapenem-hydrolyzing", 
                 curated_gene_product_name) ~ "CPE",
           grepl("metallo-beta-lactamase", 
                 curated_gene_product_name) ~ "CPE",
           grepl("extended-spectrum beta-lactamase", 
                 curated_gene_product_name) ~ "ESBL",
           grepl("class C", 
                 curated_gene_product_name) ~ "AmpC",
           grepl("beta-lactamase", 
                 curated_gene_product_name) ~ "Penicillinase"
         )) -> 
  bls

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


ggtree(treeio::tree_subset (st410_tree, 1083, levels_back = 0)) %>%
  gheatmap(
    select(country.ohehot, -accession),
    font.size = 3,
    width = 2,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.0002
  ) %>% 
  gheatmap(
    select(bl410.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
  )  +
  
  geom_treescale(x = 0.0002, y = 50, offset = 2) +
  ylim(NA, 600)


ggtree(treeio::tree_subset (st167_tree, 320, levels_back = 0)) %>%
  gheatmap(
    select(st167_metadata, malawi),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90
  ) %>% 
  gheatmap(
    select(bl.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.00005
  ) 


gheatmap(
  select(amr.ariba, -name),
  font.size = 3,
  width = 2,
  colnames_position = "top",
  color = NA,
  colnames_offset_y = 5,
  colnames_angle = 90,
  offset = 0.2
) +
  ylim(NA, 340)


data.frame(tipname = st410_tree$tip.label) -> tipz
rownames(tipz) = tipz$tipname
tipz$tipselect <- tipz$tipname == "SRR5714046"

ggtree(treeio::tree_subset (st410_tree, 1083, levels_back = 0)) %>%
  gheatmap(
    select(tipz, tipselect)
  ) 

# B4/H24RxC MDR has lost its carbapenemases in Malawi

ggtree(treeio::tree_subset (st410_tree, 564, levels_back = 0)) %>%
  gheatmap(
    select(country.ohehot, -accession),
    font.size = 3,
    width = 2,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
    offset = 0.0002
  ) %>% 
  gheatmap(
    select(bl410.onehot, ESBL, CPE),
    font.size = 3,
    width = 0.05,
    colnames_position = "top",
    color = NA,
    colnames_offset_y = 5,
    colnames_angle = 90,
  )  + geom_treescale(x = 0.0002, y = 50, offset = 2) +
  ylim(NA, 200) + geom_tiplab(size = 2)

ggtree(treeio::tree_subset (st410_tree, 564, levels_back = 0)) %>%
gheatmap(
  select(amr.ariba410, -name),
  font.size = 3,
  width = 2,
  colnames_position = "top",
  color = NA,
  colnames_offset_y = 5,
  colnames_angle = 90,
) +
  ylim(NA, 200) + geom_treescale()
