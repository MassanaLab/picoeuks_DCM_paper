library(speedyseq) |> suppressMessages()
library(ggtree) |> suppressMessages()
library(ape) |> suppressMessages()
library(picante) |> suppressMessages()
library(FSA) |> suppressMessages()
library(patchwork) |> suppressMessages()
library(tidyverse) |> suppressMessages()


# General -----------------------------------------------------------------



map <- purrr::map

remove_empty_taxa <- function(physeq){
  
  prune_taxa(x = physeq, taxa = taxa_sums(physeq) > 0)
  
}

percentage <- function(x){100*x/sum(x)}

theme_general <- function() {
  theme(
    axis.text = element_text(color = "black", size = 8),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 8),
    legend.key.size = unit(.4, "cm"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.box.spacing = unit(0.1, "cm"),
    # legend.justification = c(.8, 0),
    legend.key.spacing.y = unit(0.05, "cm")
  )
}

all_colors <- c(
  # --- Set 1: "All" Specific Groups ---
  "MALV-I" = "#FDBF6F", # Light Orange
  "MALV-II" = "#FF7F00", # Dark Orange
  "Ciliophora" = "#5A1846", # Amber/Warm Yellow (Added)
  "Polycystinea" = "#6A3D9A", # Deep Purple
  "Dinoflagellata" = "#B15928", # Brown
  "Acantharea" = "#CAB2D6", # Light Purple
  "MALV-III" = "#1F78B4", # Dark Blue
  "RAD-B" = "#A6CEE3", # Light Blue
  
  # --- Shared Groups (All & Photo) ---
  "Prymnesiophyceae" = "#E31A1C", # Red
  "Pelagophyceae" = "#FB9A99", # Pink
  "Mamiellophyceae" = "#33A02C", # Dark Green
  
  # --- Shared Groups (All & Phago) ---
  "MAST-3" = "#B2DF8A", # Light Green
  
  # --- Set 2: "Photo" Specific Groups ---
  "Dictyochophyceae" = "#FFFF99", # Light Yellow
  "Chrysophyceae" = "#FFD700", # Gold
  "Cryptomonadales" = "#8DD3C7", # Teal/Cyan
  "MOCH-2" = "#4E79A7", # Steel Blue
  "Chloropicophyceae" = "#A6D854", # Lime Green
  "Chlorarachnea" = "#E7298A", # Magenta
  "Prasino-Clade-IX" = "#BEBADA", # Lavender
  "MOCH-5" = "#7570B3", # Slate Blue, change from Pyramimonadales
  
  # --- Set 3: "Phago" Specific Groups ---
  "Picozoa" = "#E6AB02", # Dark Yellow/Mustard
  "Telonemia" = "#A6761D", # Ochre/Tan
  "MAST-1" = "#9400D3", # Dark Violet
  "MAST-4" = "#00F5FF", # Vivid Turquoise
  "MAST-9" = "#F46D43", # Coral Red
  "Centroplasthelida" = "#4575B4", # Medium Blue
  "MAST-8" = "#FDAE61", # Peach
  "Kathablepharidacea" = "#D53E4F", # Crimson
  "MAST-7" = "#228B22", # Forest Green
  "Choanoflagellata" = "#8C510A", # Dark Brown, change from Choano
  
  # --- Set 5: Fungi / New Additions ---
  "Basidiomycota" = "#BC8F8F", # Rosy/Rosy Brown (Added)
  
  # --- Set 4: Other ---
  "Other" = "gray80"
)


# Fuzzy clustering --------------------------------------------------------


plot_memberships <- function(clusters_long, metadata_df, ncol = 7) {
  clusters_membership_level <-
    clusters_long |>
    mutate(
      type =
        case_when(
          cluster == "low_abundance" ~ cluster,
          cluster == "low_membership" ~ "<70%",
          TRUE ~ ">70%"
        ),
      type = factor(type, levels = c(">70%", "<70%", "low_abundance"))
    ) |>
    group_by(sample, type) |>
    summarise(perc = sum(perc), .groups = "drop") |>
    tidyr::complete(
      sample,
      type = factor(c(">70%", "<70%", "low_abundance"), levels = c(">70%", "<70%", "low_abundance")),
      fill = list(perc = 0)
    ) |>
    left_join(metadata_df)
  
  colors <-
    c(
      ">70%" = "deeppink",
      "<70%" = "lightpink",
      "low_abundance" = "gray70"
    )
  
  plot_membership <-
    ggplot(
      data = clusters_membership_level,
    ) +
    geom_line(
      aes(y = depth, x = perc, color = type),
      linewidth = 1,
      orientation = "y"
    ) +
    geom_point(
      aes(y = depth, x = perc, fill = type, color = type),
      pch = 21, size = 1
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Clusters",
      color = "Clusters"
    ) +
    scale_y_reverse() +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25)
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(
        linewidth = .5, fill = "transparent",
      ),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    facet_wrap(~profile_abb, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = colors, guide = guide_legend(direction = "horizontal", nrow = 1), labels = c(">70%", "<70%", "Low abundant")) +
    scale_fill_manual(values = colors, guide = guide_legend(direction = "horizontal", nrow = 1), labels = c(">70%", "<70%", "Low abundant"))
}

make_clusters_long_v2 <- function(clusters) {
  clusters |>
    map(
      ~ .x |>
        dplyr::rename(cluster = cluster_final) |>
        select(-cluster_max) |>
        pivot_longer(names_to = "sample", values_to = "perc", -c(asv_id, Group, Supergroup, cluster, membership))
    ) |>
    bind_rows(.id = "profile")
}

fuzzy_clustering_perc_occur <- function(phy_table, profile, n_clusters, min_mean_abun = 0.01, min_occur = 2, membership_threshold = 0.7, normalize_perc = T) {
  require(e1071)
  
  keep_samples <-
    phy_table |>
    sample_data() |>
    as_tibble() |>
    filter(St == profile) |>
    pull(.sample)
  
  phy_table_profile <-
    phy_table |>
    prune_samples(samples = keep_samples) |>
    remove_empty_taxa()
  
  phy_table_profile_perc <-
    phy_table_profile |>
    transform_sample_counts(percentage)
  
  keep_asvs <-
    tibble(
      asv_id = taxa_names(phy_table_profile_perc),
      mean_perc = rowMeans(otu_table(phy_table_profile_perc)),
      occur = rowSums(otu_table(phy_table_profile_perc) > 0)
    ) |> 
    filter(
      mean_perc >= min_mean_abun,
      occur >= min_occur
    ) |>
    pull(asv_id)
  
  keep_samples <-
    sample_names(phy_table_profile)
  
  asv_table_perc <-
    phy_table_profile_perc |>
    otu_table() |>
    as_tibble() |>
    filter(
      .otu %in% keep_asvs,
      .sample %in% keep_samples
    )
  
  if (normalize_perc == T) {
    input_fuzzy <-
      asv_table_perc |>
      group_by(.otu) |>
      mutate(norm_abun = .abundance * 100 / max(.abundance)) |>
      select(.otu, .sample, .abundance = norm_abun) |>
      pivot_wider(names_from = .sample, values_from = .abundance) |>
      column_to_rownames(".otu") |>
      as.matrix()
  } else {
    input_fuzzy <-
      asv_table_perc |>
      pivot_wider(names_from = .sample, values_from = .abundance) |>
      column_to_rownames(".otu") |>
      as.matrix()
  }
  
  mestimate <- function(df) {
    N <- dim(df)[[1]] # num ASVs
    D <- dim(df)[[2]] # num samples
    m.sj <- 1 + (1418 / N + 22.05) * D^(-2) + (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
    return(m.sj)
  }
  
  m <- mestimate(input_fuzzy)
  k <- n_clusters[[profile]]
  
  set.seed(1789)
  clustering <-
    cmeans(input_fuzzy, centers = k, m = m)
  
  clusters <-
    clustering$membership |>
    as_tibble(rownames = "asv_id") |>
    pivot_longer(-asv_id) |>
    group_by(asv_id) |>
    slice_max(value, n = 1) |>
    mutate(cluster = if_else(value >= membership_threshold, name, "low_membership"))
  
  otu_tab_df <-
    phy_table_profile_perc |>
    otu_table() |>
    as.data.frame() |>
    as_tibble(rownames = "asv_id") |>
    left_join(clusters |> select(asv_id, cluster_max = name, cluster_final = cluster, membership = value)) |>
    left_join(phy_table_profile |> tax_table() |> as_tibble() |> select(asv_id = .otu, Group, Supergroup)) |>
    replace_na(list(
      cluster_final = "low_abundance",
      cluster_max = "low_abundance"
    )) |>
    select(asv_id, Group, Supergroup, contains("cluster"), membership, everything())
  
  if (normalize_perc) {
    norm_df <-
      input_fuzzy |>
      as_tibble(rownames = "asv_id") |>
      rename_with(.cols = -asv_id, ~ str_c("norm_", .x))
    
    
    otu_tab_df <-
      otu_tab_df |>
      left_join(norm_df)
  }
  
  return(otu_tab_df)
}

clusterer_and_plotter_occur <- function(phy_table, n_clusters, min_mean_abun = 0.01, min_occur = 1, membership_threshold = 0.7, normalize_perc = T, get_clusters_all = F, get_clusters_long = F, plot_membership = T, plot_only_memberships = F, ncol = 5, cluster_names, metadata_df) {
  profiles <- names(n_clusters)
  
  clusters <-
    profiles |>
    set_names() |>
    map(
      ~ fuzzy_clustering_perc_occur(
        phy_table = phy_table,
        n_clusters = n_clusters,
        profile = .x,
        normalize_perc = normalize_perc,
        min_occur = min_occur
      )
    )
  
  if (get_clusters_all) {
    clusters_with_names <-
      clusters |>
      imap(
        ~ .x |>
          left_join(cluster_names |> filter(profile == .y) |> select(-contains('profile'), cluster_name, cluster_max = cluster)) |>
          select(asv_id, Group, Supergroup, cluster_max, cluster_name, membership, cluster_final, everything())
      )
    
    return(clusters_with_names)
  }
  
  clusters_long <-
    clusters |>
    make_clusters_long_v2() |>
    filter(!str_detect(sample, "norm_")) |>
    left_join(cluster_names)
  
  if (get_clusters_long) {
    return(clusters_long)
  }
  
  clusters_group_level <-
    clusters_long |>
    filter(!cluster %in% c("low_membership", "low_abundance")) |>
    left_join(sample_data(phy_table) |> as_tibble() |> select(sample = .sample, depth = Depth)) |>
    group_by(cluster_name, profile_abb, depth) |>
    summarise(perc = sum(perc)) |>
    mutate(
      cluster_name_orig = cluster_name,
      cluster_name = case_when(
        str_detect(cluster_name, "^Above") ~ "Above",
        TRUE ~ cluster_name
      ),
      cluster_name = factor(cluster_name, levels = c("Surface", "Above", "DCM", "Below", "Deep"))
    )
  
  colors <-
    c(
      Surface = "tomato",
      Above = "gold",
      Above1 = "gold",
      Above2 = "gold",
      DCM = "darkolivegreen2",
      Below = "plum",
      Deep = "royalblue",
      Deep1 = "royalblue",
      Deep2 = "royalblue"
    )
  
  scale_y_limits <- 
    ggh4x::facetted_pos_scales(
      y = 
        list(
          str_detect(profile_abb, '^R') ~ scale_y_continuous(limits = c(0, 200), breaks = c(0,50,100,150,200)),
          str_detect(profile_abb, '^M') ~ scale_y_continuous(limits = c(0, 150), breaks = c(0,50,100,150), labels = c(0,50,100,150)),
          profile_abb == 'P11' ~ scale_y_continuous(limits = c(0, 280)),
          profile_abb %in% c('P04','P08','P10','P18','P26') ~ scale_y_continuous(limits = c(0, 200), breaks = c(0,50,100,150,200))
        )
    )
  
  clustering_plot <-
    base_plot_fluorescence +
    geom_line(
      data = clusters_group_level,
      aes(y = depth, x = perc, color = cluster_name, group = cluster_name_orig),
      linewidth = 1,
      orientation = "y"
    ) +
    geom_point(
      data = clusters_group_level,
      aes(y = depth, x = perc, fill = cluster_name, color = cluster_name, group = cluster_name_orig),
      pch = 21, size = 1
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Clusters",
      color = "Clusters"
    ) +
    scale_y_limits +
    scale_y_reverse() +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25)
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(
        linewidth = .5, fill = "transparent",
      ),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    facet_wrap(~profile_abb, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = colors, guide = guide_legend(direction = "horizontal", nrow = 1)) +
    scale_fill_manual(values = colors, guide = guide_legend(direction = "horizontal", nrow = 1))
  
  if (plot_membership == F){
    return(clustering_plot)  
  }
  
  membership_plot <- 
    clusters_long |> 
    plot_memberships(metadata_df = metadata_df, ncol = ncol)
  
  if (plot_only_memberships == T){
    return(membership_plot)
  }
  
  merged_plot <- 
    membership_plot + clustering_plot
  
  return(merged_plot)
}

distr_clusters_plotter <- function(group = c("pico", "phago", "photo"), cluster_summary, keep_intermediate_layers = T) {
  group <- match.arg(group)
  
  groups <- cluster_summaries_to_export[[group]]$Group
  
  if (group == "pico") {
    groups <- groups[1:14]
  } else if (group %in% c("photo", "phago")) {
    groups <- groups[1:11]
  } 
  
  if (keep_intermediate_layers == F){
    cluster_summary <- 
      cluster_summary |> 
      filter(cluster_name %in% c('Surface', 'DCM', 'Deep'))
  }
  
  df_processed <-
    cluster_summary |>
    mutate(
      group2 = if_else(Group %in% groups, Group, "Other"),
      group2 = factor(group2, levels = c(groups, "Other"))
    ) |>
    group_by(profile_abb, layer, group2) |>
    summarise(perc = sum(perc)) |>
    mutate(perc = 100 * perc / sum(perc))
  
  plot <-
    df_processed |>
    ggplot(aes(x = perc, y = "")) +
    geom_col(aes(fill = fct_rev(group2)), color = "gray20", linewidth = .2) +
    ggh4x::facet_grid2(profile_abb ~ layer, scales = "free", space = "free") +
    labs(y = NULL, x = "Relative group abundance (%)", ) +
    guides(fill = guide_legend(position = "bottom", reverse = T, title = NULL)) +
    scale_fill_manual(values = all_colors) +
    theme_bw() +
    theme(axis.ticks.y = element_blank()) +
    theme_general()
  
  return(plot)
}

fig5_plotter <- function(df, ngroups = 11) {
  
  top_groups <-
    df |>
    group_by(Group) |>
    summarise(total = sum(mean_perc)) |>
    arrange(-total) |>
    slice_head(n = ngroups) |>
    pull(Group)
  
  df_processed <-
    df |>
    mutate(
      Group = if_else(Group %in% top_groups, Group, "Other"),
      Group = factor(Group, levels = c(top_groups, "Other")),
    ) |>
    group_by(layer, Group) |>
    summarise(perc = sum(mean_perc))
  
  plot <-
    df_processed |>
    ggplot(aes(x = perc, y = fct_rev(layer))) +
    geom_col(aes(fill = fct_rev(Group)),
             color = "gray20",
             linewidth = .2
    ) +
    theme_bw() +
    theme_general() +
    theme(
      legend.key.spacing.x = unit(.2, 'cm')
    ) +
    guides(fill = guide_legend(position = "bottom", reverse = T, title = NULL)) +
    labs(y = NULL, x = "Relative group abundance (%)") +
    scale_fill_manual(
      values = all_colors,
      labels = function(x) stringr::str_pad(x, width = 19, side = "right")
    ) +
    scale_x_continuous(expand = c(0,0))
  
  return(plot)
}

run_kruskal_dunn <- function(df, padj_method = "bh") {
  
  require(broom)
  require(FSA)
  
  kw_res <-
    df |> 
    group_by(Group) |> 
    reframe(broom::tidy(kruskal.test(perc ~ layer))) |> 
    mutate(p.adj = p.adjust(p.value, method = toupper(padj_method)))
  
  sig_groups <- kw_res |> 
    filter(p.adj < 0.05) |> 
    pull(Group)
  
  if (length(sig_groups) == 0) {
    dunn_res <- tibble(
      Group = character(), 
      Comparison = character(), 
      Z = numeric(), 
      P.unadj = numeric(), 
      P.adj = numeric()
    )
  } else {
    dunn_res <- df |> 
      filter(Group %in% sig_groups) |> 
      group_by(Group) |> 
      reframe({
        d_out <- FSA::dunnTest(perc ~ layer, data = pick(everything()), method = padj_method)
        d_out$res
      }) |> 
      rename(P.adj = P.adj)
  }
  
  return(
    list(
      kruskal = kw_res,
      dunn = dunn_res)
  )
}

# Trees -------------------------------------------------------------------

remove_outgroup <- function(tree){
  
  require(phytools)
  
  group_tips <- 
    tree$tip.label[!str_detect(tree$tip.label, '_outgroup')]
  
  rooted_tree <- 
    root_tree_at_outgroup(tree)
  
  tree_node <- 
    MRCA(rooted_tree, group_tips)
  
  tree_without_outgroup <- 
    extract.clade(rooted_tree, tree_node)
  
  return(tree_without_outgroup)
}

add_outgroup <- function(group, n = 5){
  
  outg <- 
    outgroups |> 
    filter(Group == group) |> 
    pull(outgroup)
  
  df_out <- 
    asv_mean_abundance_and_occur |> 
    filter(Group == outg) |> 
    distinct(asv_id, ASV) |> 
    slice_head(n = n) |> 
    mutate(name = str_c(asv_id, '_outgroup')) |>
    select(name, ASV)
  
  return(df_out)
}

root_tree_at_outgroup <- function(tree){
  
  require(phytools)
  
  outgroup_tips <- 
    tree$tip.label[str_detect(tree$tip.label, '_outgroup')]
  
  outgroup_node <- 
    MRCA(tree, outgroup_tips)
  
  rooted_tree <- 
    reroot(tree, outgroup_node)
  
  return(rooted_tree)
  
}

create_base_tree <- function(tree, linewidth = 0.2){
  
  p_tree_base <- 
    tree |> 
    ggtree(
      layout = 'circular',
      color = NA
    ) +
    geom_tree(
      linewidth = linewidth
    )
  
  return(p_tree_base)
  
}

add_geom_col_biome <- function(tree_plot, metadata_biomes, offset = 0.1){
  
  tree_cols <- 
    tree_plot +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      geom = geom_col,
      data = metadata_biomes,
      mapping = aes(x = n, y = label, fill = fct_rev(biome)), 
      color = 'gray20', 
      linewidth = .2,
      axis.params = list(
        axis = "x",
        text.size = 2.5,
      ),
      grid.params=list(),
      offset = offset
    ) +
    scale_fill_manual(values = c(Atlantic = 'burlywood2', Mediterranean = 'lavender', `Not clustered` = 'gray70'),
                      guide = guide_legend(reverse = T, title = NULL, position = 'top', direction = 'horizontal'))
  
  return(tree_cols)
}

add_geom_tile_occurrence <- function(tree_plot, metadata_occurrence){
  
  colors <- viridis::mako(n = 100, direction = -1)
  colors[1] <- "white"
  
  metadata_df <- 
    metadata_occurrence |> 
    complete(
      label, 
      cluster2, 
      fill = list(n = 0)
    )
  
  pwidth = 0.2
  
  if (str_detect(metadata_df$label[1], 'MAST-3')){
    metadata_df <- 
      metadata_df |> 
      filter(cluster2 != 'DCM') |> 
      mutate(cluster2 = factor(cluster2, levels = c('Surface', 'Deep')))
    
    pwidth = .1
  }
  
  tree_tiles <- 
    tree_plot +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = metadata_df,
      mapping = aes(x = cluster2, y = label, fill = n),
      geom = geom_tile,
      axis.params=
        list(
          axis = 'x',
          text.angle = 90,
          text.size = 2.5,
          hjust = 1
        ),
      color = 'gray70',
      linewidth = 0.2,
      pwidth=pwidth
    ) +
    scale_fill_gradientn(
      name = "Num.\nprofiles", 
      colours = colors, 
      limits = c(0, 14), 
      breaks = c(0, 7, 14),
      guide = guide_colorbar(position = 'top', direction = 'horizontal'),
      na.value = 'white'
    )
  
  return(tree_tiles)
  
}

plot_heatmap_layer <- function(tree, perc_col = 'n'){
  
  order_labels = get_taxa_name(ggtree(tree))
  
  group <- str_match(order_labels[1], 'ASV_.*?_(.*?)_')[,2] |> str_remove('[A-I]$')
  if (str_detect(group, 'PIC-')){
    group <- 'Picozoa'
  }
  
  trophic <- groups_to_trophic |> filter(Group == group) |> pull(trophic)
  
  df <- 
    asv_data_for_trees |> 
    filter(label %in% order_labels) |> 
    mutate(label = factor(label, levels = order_labels)) |> 
    select(asv_id, label, cluster2, perc = all_of(perc_col)) |> 
    group_by(asv_id, label) |> 
    tidyr::complete(cluster2)
  
  if (trophic == 'phago'){
    df <- 
      df |> 
      filter(cluster2 != 'DCM') |> 
      mutate(cluster2 = factor(cluster2, levels = c('Surface','Deep')))
  }
  
  p_heatmap <- 
    df |>
    ggplot(aes(x = cluster2, y = fct_rev(label))) + 
    geom_tile(aes(fill= perc), color = 'gray20', linewidth = .2) +
    scale_fill_viridis_c(option ='mako', direction = -1, na.value = 'white', name = 'Num.\nprofiles', 
                         guide = guide_colorbar(position = 'top', direction = 'horizontal'),
                         limits = c(0,14), breaks = c(0,7,14)) +
    scale_y_discrete(drop = F) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  return(p_heatmap)
  
}

plot_bar_nprofiles <- function(tree){
  
  order_labels = get_taxa_name(ggtree(tree))
  
  df <- 
    nprofiles_per_biome |> 
    filter(label %in% order_labels) |> 
    mutate(label = factor(label, levels = order_labels))
  
  p_bar <- 
    df |>
    ggplot(aes(x = n, y = fct_rev(label))) + 
    geom_col(aes(fill= fct_rev(biome)),color = 'gray20', linewidth = .2) +
    scale_y_discrete(drop = F) +
    labs(x = 'N. profiles', y = NULL) +
    theme_minimal() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank()
    )  +
    scale_fill_manual(values = c(Atlantic = 'burlywood2', Mediterranean = 'lavender', `Not clustered` = 'gray70'),
                      guide = guide_legend(reverse = T, title = NULL, position = 'top', direction = 'horizontal'))
}

plot_tree_all <- function(tree, perc_col = 'n', guides = T, add_title = F){
  
  order_labels = get_taxa_name(ggtree(tree))
  
  group <- str_match(order_labels[1], 'ASV_.*?_(.*?)_')[,2] |> str_remove('[A-I]$')
  if (str_detect(group, 'PIC-')){
    group <- 'Picozoa'
  }
  
  p_tree <-
    tree |> 
    ggtree(linewidth = 0.2) + 
    geom_tiplab(aes(label = str_match(label, '(ASV_.*?)_')[,2])) +
    xlim(trees_parameters[[group]][['xlim1']],trees_parameters[[group]][['xlim2']]) +
    coord_cartesian(clip = 'off')
  
  if (add_title != F){
    p_tree <- 
      p_tree + labs(title = add_title)
  }
  
  p_heatmap <- 
    plot_heatmap_layer(tree, perc_col = perc_col)
  
  if (guides == F){
    p_heatmap <- 
      p_heatmap + guides(fill = 'none')
  }
  
  p_bar <- 
    plot_bar_nprofiles(tree)
  
  if (guides == F){
    p_bar <- 
      p_bar + guides(fill = 'none')
  }
  
  p_final <- 
    wrap_plots(
      p_tree,
      p_heatmap, 
      p_bar,
      widths = c(.7,.1,.2),
      guides = 'collect'
    ) & 
    theme(
      legend.position = 'top'
    )
  
  return(p_final)
  
}

get_tree_height <- function(num_seqs, base_height = 2.5, scale_factor = 0.2){
  return(base_height + (scale_factor*num_seqs))
}

order_trees <- 
  c("Prymnesiophyceae", "Pelagophyceae", "Mamiellophyceae", "Dictyochophyceae", "Chrysophyceae", "Cryptomonadales", "MOCH-2", "Chloropicophyceae", "Chlorarachnea", "Prasino-Clade-IX", "MOCH-1", "MAST-3", "Picozoa", "Telonemia", "MAST-1", "MAST-4", "MAST-9", "Centroplasthelida", "Kathablepharidacea", "Cercozoa", "MAST-8", "Choanoflagellata", "MAST-7")

# Environmental variables -------------------------------------------------

plot_cell_counts <- function(profile, fluor_color = 'gray70', neuk_color = 'gold2', peuk_color = 'darkorchid4'){
  
  neuk_df <- neuk_values[[profile]] |> filter(depth < 220) |> arrange(depth)
  peuk_df <- peuk_values[[profile]] |> filter(depth < 220) |> arrange(depth)
  fluor_df <- fluorescence_values[[profile]] |> filter(depth < 220) |> arrange(depth)
  
  axis.titles.cex = 0.8
  axis.numbers.cex = 1.3
  
  # c(bottom, left, top, right)
  par(mar = c(4.5, 3.5, 3.5, 1.5) + 0.1)
  
  max_depth = 220
  min_depth = 0
  
  fluor <- fluor_df$fluorescence
  depth <- fluor_df$depth
  
  plot(fluor, depth, type = "n", 
       ylim = c(max_depth, min_depth), 
       xlim = c(0, max(fluor, na.rm = TRUE) * 1.1),
       xlab = "", ylab = "", axes = FALSE)
  
  polygon(c(0, fluor, 0), c(min_depth, depth, max_depth), col = fluor_color, border = NA)
  
  axis(2, cex.axis = axis.numbers.cex)
  if (profile %in% c('RD1','MM','P10')){
    mtext("Depth (m)", side = 2, line = 2.2, cex = axis.titles.cex)
  }
  box()
  
  neuk <- neuk_df$neuk 
  neuk_depth <- neuk_df$depth  
  
  par(new = TRUE)
  plot(neuk, neuk_depth, type = "n",
       ylim = c(max_depth, min_depth), 
       xlim = c(0, max(neuk, na.rm = TRUE) * 1.1),
       axes = FALSE, xlab = "", ylab = "")
  
  lines(neuk, neuk_depth, lwd = 1.5, col = neuk_color)
  points(neuk, neuk_depth, pch = 19, col = neuk_color, cex = 1.2)
  
  axis(3, col = "black", col.axis = neuk_color, cex.axis = axis.numbers.cex)
  if (profile %in% c('RD1','RD2','RD3','RD4','MD')){
    mtext(expression(Nanoeukaryotes ~ (10^2 ~ cells ~ mL^{-1})), 
          side = 3, line = 2, col = neuk_color, cex = axis.titles.cex)
  }
  
  peuk <- peuk_df$peuk 
  peuk_depth <- peuk_df$depth
  
  par(new = TRUE)
  plot(peuk, peuk_depth, type = "n",
       ylim = c(max_depth, min_depth), 
       xlim = c(0, max(peuk, na.rm = TRUE) * 1.1),
       axes = FALSE, xlab = "", ylab = "")
  
  lines(peuk, peuk_depth, lwd = 1.5, col = peuk_color)
  points(peuk, peuk_depth, pch = 19, col = peuk_color, cex = 1.2)
  
  axis(1, col = "black", col.axis = peuk_color, cex.axis = axis.numbers.cex)
  if (profile %in% c('P10','P11','P18','P26', 'P08')){
    mtext(expression(Picoeukaryotes ~ (10^2 ~ cells ~ mL^{-1})), 
          side = 1, line = 3, col = peuk_color, cex = axis.titles.cex)
  }
  
  text(x = par("usr")[2] * 0.95, y = 195, 
       labels = profile, font = 2, cex = 1.2, adj = c(1, 1))
}

plot_env_variables <- function(profile, fluor_color = 'gray70', temp_color = 'tomato', nitr_color = 'dodgerblue3'){
  
  if (str_detect(profile, '^M')){
    max_depth = 150
  } else {
    max_depth = 220
  }
  min_depth = 0
  
  temp_df <- temperature_values[[profile]] |> filter(depth < max_depth) |> arrange(depth)
  nitr_df <- nitrate_values[[profile]] |> filter(depth < max_depth) |> arrange(depth) |> filter(!is.na(nitrate))
  fluor_df <- fluorescence_values[[profile]] |> filter(depth < max_depth) |> arrange(depth)
  samp_df <- sampling_values[[profile]] |> filter(depth < max_depth) |> arrange(depth)
  
  axis.titles.cex = 0.75
  
  # c(bottom, left, top, right)
  par(mar = c(7, 4, 4, 2) + 0.1)
  
  fluor <- fluor_df$fluorescence
  depth <- fluor_df$depth
  
  plot(fluor, depth, type = "n", 
       ylim = c(max_depth, min_depth), 
       xlim = c(0, max(fluor, na.rm = TRUE) * 1.1),
       xlab = "", ylab = "", axes = FALSE)
  
  polygon(c(0, fluor, 0), c(min_depth, depth, max_depth), col = fluor_color, border = NA)
  
  
  samp <- samp_df$fluorescence
  depth <- samp_df$depth
  
  points(samp, depth, pch = 19, cex = 1.2, col = "black")
  
  axis(1, col.axis = fluor_color)
  
  if (profile == 'P18'){
    mtext("Fluorescence", side = 1, line = 2, col = fluor_color, cex = axis.titles.cex)
  }
  
  axis(2)
  
  if (profile == 'MM'){
    mtext("Depth (m)", side = 2, line = 2.2, cex = axis.titles.cex)
  }
  box()
  
  
  
  temp <- temp_df$temperature
  depth <- temp_df$depth  
  
  par(new = TRUE)
  plot(temp, depth, type = "l", lty = 1, lwd = 1.5,
       ylim = c(max_depth, min_depth), 
       xlim = c(13,30),
       axes = FALSE, xlab = "", ylab = "",
       col = temp_color)
  
  if (profile %in% c('RD1','RD2','RD3','RD4','MD')){
    axis(3, col.axis= temp_color)
  }
  
  if (profile == 'RD3'){
    mtext("Temperature (ºC)", side = 3, line = 2, col = temp_color, cex = axis.titles.cex)
  }
  
  nitr <- nitr_df$nitrate
  nitr_depth <- nitr_df$depth
  
  par(new = TRUE)
  plot(nitr, nitr_depth, type = "n",
       ylim = c(max_depth, min_depth), 
       xlim = c(min(nitr, na.rm = T), ceiling(max(nitr, na.rm = T))),
       axes = FALSE, xlab = "", ylab = "")
  
  if (profile != 'P11'){
    lines(nitr, nitr_depth, lwd = 1.5, col = nitr_color)
    points(nitr, nitr_depth, pch = 21, bg = "white", col = nitr_color, cex = 1.2)
  }
  
  axis(1, line = 3.5, col.axis = nitr_color)
  
  if (profile == 'P18'){
    mtext(expression(NO[3]^"-" ~ "(µM)"), side = 1, line = 6, col = nitr_color, cex = axis.titles.cex)
  }
  
  if (profile == 'P10'){
    title.y = max_depth * .80
  } else {
    title.y = max_depth * .90
  }
  text(x = par("usr")[2]*.95, y = title.y, 
       labels = profile, font = 2, cex = 1.2, adj = c(1, 1))
}