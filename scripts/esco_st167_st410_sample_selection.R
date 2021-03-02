
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

# choose reference
  
st167_metadata %>% 
  filter(platform != "ILLUMINA") %>% View()

# US blood culture isolate 2014 - part of US FDARGOS curated reference
# collection
# GCF_002634895

st167_metadata %>% 
  filter(Name == "FDAARGOS_434") %>% View()



# get all illumina read with valid year and country

st167_metadata %>% 
  filter(!is.na(`Collection Year`) &
           !is.na(Country) &
           platform == "ILLUMINA" &
           !grepl("traces", accession)) -> st167_illumina

         
write_csv(st167_illumina, here("data_processed/esco/st167/st167_illumina_included_metadata.csv"))
write_lines(st167_illumina$accession, here("data_processed/esco/st167/st167_illumina_included_lanes.txt"))
  
st167_illumina %>% 
  select(`Collection Year`, Country) %>% 
  ggplot(aes(fct_infreq(Country))) + 
  geom_bar() +
  facet_wrap(~ `Collection Year`, scales = "free") +
  coord_flip()



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

# choose reference

st410_metadata %>% 
  filter(platform != "ILLUMINA") %>% View()

# Canadian rectal swab GCF_002589795, FDAARGOS_433


# get all illumina read with valid year and country

st410_metadata %>% 
  filter(!is.na(`Collection Year`) &
           !is.na(Country) &
           platform == "ILLUMINA" &
           !is.na(`Source Niche`) &
           !grepl("traces", accession)) -> st410_illumina

st410_illumina %>% 
  select(`Collection Year`, Country) %>% 
  ggplot(aes(fct_infreq(Country))) + 
  geom_bar() +
  facet_wrap(~ `Collection Year`, scales = "free") +
  coord_flip()

# too many - sample with equal weighting on year, n = 300
# oh fuck it lets just do all

write_lines(st410_illumina$accession, here("data_processed/esco/st410/st410_illumina_included_lanes.txt"))


write_csv(st167_illumina, here("data_processed/esco/st167/st167_illumina_included_metadata.csv"))
write_lines(st167_illumina$accession, here("data_processed/esco/st167/st167_illumina_included_lanes.txt"))

DASSIM_mlst <- read_csv("~/Documents/PhD/Thesis/bookdown/chapter_7/phylogroup_and_mlst/mlst.csv")

write_lines(
DASSIM_mlst %>%  
  filter(ST == "167") %>% 
  mutate(lane = gsub("_1_", "_1#", lane),
         lane = gsub("_2_", "_2#", lane)) %>% 
  pull(lane),
"data_processed/esco/st167/dassim_st167.txt")

  DASSIM_mlst %>%  
    filter(ST == "167") %>% 
    mutate(lane = gsub("_1_", "_1#", lane),
           lane = gsub("_2_", "_2#", lane)) -> st167 
  
  st167 %>% 
    mutate(lane = gsub("_1#", "_1_", lane),
           lane = gsub("_2#", "_2_", lane)) %>% 
    filter(lane %in% tree$tip.label) 
  
  tree <- read.tree("~/Documents/PhD/Thesis/bookdown/chapter_7/core_genome_tree/core_alnD2ESCO_snp_sites.aln.treefile")
  midpoint.root(tree) -> tree

  
  
  
  ## st410
  write_lines(
  DASSIM_mlst %>% filter(ST == 410) %>% 
    filter(lane %in% tree$tip.label) %>%
    mutate(lane = gsub("1_","1#", lane)) %>% 
    pull(lane),
  "data_processed/esco/st410/dassim_st410_to_include.txt"
  )

