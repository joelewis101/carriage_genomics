# qc on snippy assemblies
# ditch any depth < 20

dep <- read_csv("data_processed/esco/st167/depth_and_cov.csv")
dep$mean_dep <- dep$depth_sum /dep$ref_bases

ggplot(dep, aes(lane, mean_dep, color = mean_dep < 20)) + geom_point()

dep %>% 
  filter(!is.na(mean_dep) & mean_dep > 20) %>% 
  mutate(lane = gsub("\\./", "", lane),
         lane = gsub("/snps\\.bam", "", lane)) %>% 
           pull(lane) ->
st167_mappings_for_gubbins

write_lines(st167_mappings_for_gubbins, "data_processed/esco/st167/st167_mappings_for_gubbins_following_qc.txt")

###


dep <- read_csv("data_processed/esco/st410/snippy/depth_and_cov.csv")
dep$mean_dep <- dep$depth_sum /dep$ref_bases

ggplot(dep, aes(lane, mean_dep, color = mean_dep < 20)) + geom_point()

dep %>% 
  filter(!is.na(mean_dep) & mean_dep > 20) %>% 
  mutate(lane = gsub("\\./", "", lane),
         lane = gsub("/snps\\.bam", "", lane)) %>% 
  pull(lane) ->
  st410_mappings_for_gubbins

write_lines(st410_mappings_for_gubbins, "data_processed/esco/st410/st410_mappings_for_gubbins_following_qc.txt")


