#### Spectral diversity functions
# Etienne Laliberté, February 2019
# Updated using the terra package by Bryce Currey, 2023


### Brightness normalization -----
bright_norm <- function(x) {
  x_norm <- x / sqrt(sum(x^2))
  return(x_norm) 
}


### Brightness normalization for data frame ----
bright_norm_df <- function(x) {
  x$refl_norm <- x$refl / sqrt(sum(x$refl^2))
  return(x)  
}


### Short PCA for matrix-----
pca_mat <- function(x, scaling = c(1, 2), p = 0.99) {
  Y <- scale(x, center = TRUE, scale = FALSE)
  n <- nrow(Y)
  Y.svd = svd(Y)
  values = (1/(n - 1)) * Y.svd$d^2
  epsilon = sqrt(.Machine$double.eps)
  k <- sum(values > epsilon)
  values <- values[1:k]
  prop <- values / sum(values)
  cumprop = cumsum(prop)
  if (p < cumprop[1]) which.values <- c(1, 2) else which.values <- which(cumprop < p)
  values.sel <- values[c(which.values, length(which.values) + 1)]
  n.pcs <- length(values.sel)
  U <- as.matrix(Y.svd$v[, 1:n.pcs])
  if (scaling == 1) {
    obj <- Y %*% U
    descript <- U
  }
  else {
    obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
    descript <- U %*% diag(values.sel^(0.5))
  }
  colnames(obj) <- paste0('PC', 1:n.pcs)
  colnames(descript) <- paste0('PC', 1:n.pcs)
  rownames(descript) <- colnames(x)
  prop.sel <- prop[1:n.pcs]; names(prop.sel) <- colnames(descript)
  cumprop.sel <- cumprop[1:n.pcs]; names(cumprop.sel) <- colnames(descript)
  out <- list(obj = obj, descript = descript, prop = prop.sel, cumprop = cumprop.sel)
  return(out)
}



### PCA for image cube -----
# updated for terra
pca <- function(cube, scaling = c(1, 2), p = 0.99) {
  require(tidyverse)
  require(terra)
  cube.df <- as.data.frame(cube, xy = T)
  xy <- cube.df[, 1:2]
  pixels <- cube.df[, 3:ncol(cube.df)]
  Y <- scale(pixels, center = TRUE, scale = FALSE)
  n <- nrow(Y)
  Y[is.na(Y)] <- 0 # shouldn't be any but could be a few.
  Y.svd = svd(Y) #singular value decomposition
  values = (1/(n - 1)) * Y.svd$d^2
  epsilon = sqrt(.Machine$double.eps)
  k <- sum(values > epsilon)
  values <- values[1:k]
  prop <- values / sum(values)
  cumprop = cumsum(prop)
  if (p < cumprop[1]) which.values <- c(1, 2) else which.values <- which(cumprop < p)
  values.sel <- values[c(which.values, length(which.values) + 1)]
  n.pcs <- length(values.sel)
  U <- as.matrix(Y.svd$v[, 1:n.pcs])
  if (scaling == 1) {
    obj <- Y %*% U
    descript <- U
  } else {
    obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
    descript <- U %*% diag(values.sel^(0.5))
  }
  colnames(obj) <- paste0('PC', 1:n.pcs)
  colnames(descript) <- paste0('PC', 1:n.pcs)
  rownames(descript) <- colnames(pixels)
  points_df <- cbind.data.frame(xy, obj)
  prop.sel <- prop[1:n.pcs];names(prop.sel) <- colnames(descript) #, proj4string = crs(cube)
  cumprop.sel <- cumprop[1:n.pcs]; names(cumprop.sel) <- colnames(descript)
  cube_pc <- rast(points_df, type = 'xyz', crs = crs(cube))
  out <- list(cube_pc = cube_pc, band_contrib = descript, prop = prop.sel, cumprop = cumprop.sel)
  return(out)
}


### Make PC plots ----
make_plot <- function(df) {
  x <- ggplot(df, aes(x = x, y = y) ) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colors = rainbow(20)) +
    ggtitle(label = unique(df$PC)) +
    theme_void() +
    theme(legend.position = 'none',
          plot.title = element_text(face = 'bold', hjust = 0.5, size = 15) ) +
    coord_equal()
  return(x)
}


### Build the PC plots ----
build_pc_plot <- function(cube) {
  require(tidyverse)
  points <- rasterToPoints(cube, spatial = F) %>% 
    as_tibble() %>% 
    tidyr::gather(key = PC, value = value, -x, -y) %>% 
    dplyr::mutate(PC = factor(PC, levels = paste0('PC', 1:n_distinct(PC)) ) ) 
  cube_plots <- points %>%
    group_by(PC) %>% 
    dplyr::do(plots = make_plot(df = .))
return(cube_plots)
}


### Show the PC plots ----
show_pc_plot <- function(x, ...) {
  require(ggplot2)
  require(gridExtra)
  grid.arrange(grobs = x$plots, ...) 
}


# Count non-NA rows
count_noNA <- function(x) {
  n <- nrow(na.omit(x))
  n <- as.data.frame(n)
}

### Count masked/missing pixels -----
# updated
count_pixels <- function(cube, fact = 40) {
  require(terra)
  require(tidyverse)
  # Get plots (or communities)
  new_res <- res(cube) * fact
  cube_plots <- rast(crs = crs(cube))
  ext(cube_plots) <- ext(cube)
  res(cube_plots) <- new_res
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  plot_xy <- as.data.frame(cube_plots, xy=T) %>% dplyr::select(x,y)
  cube_pixels <- disagg(cube_plots, fact = fact)
  # Convert to points
  plot_points <- as.data.frame(cube_pixels, xy=T) %>%
    dplyr::rename_with(~c('x','y','group'))
  cube_points <- as.data.frame(cube, xy=T) %>% 
   right_join(plot_points, by = c('x', 'y'))
  n_pixels <- cube_points %>% 
    group_by(group) %>% 
    do(count_noNA(.)) %>% 
    mutate(n_total = fact^2,
           prop = n / n_total) %>% 
    ungroup() %>%
    dplyr::select(n:prop)
  cube_count <- rast(cbind.data.frame(plot_xy, n_pixels), 'xyz', crs = crs(cube))
  return(cube_count)
}


### SS gamma and alpha ----
# updated to include na.rm = T
sum_squares <- function(Y) {
  n <- nrow(Y)
  Y.cent <- scale(Y, center = T, scale = F)
  sij <- Y.cent^2
  SS.total <- sum(sij, na.rm = T)
  SS.row <- rowSums(sij, na.rm = T)
  SS.col <- colSums(sij, na.rm = T)
  fcsd <- SS.col / SS.total
  lcsd <- SS.row / SS.total
  sdiv <- SS.total / (n - 1)
  out <- list(ss = SS.total, sdiv = sdiv, lcsd = lcsd, fcsd = fcsd)
  return(out)
}


### SS beta ----
# updated to include na.rm = T
sum_squares_beta <- function(Y, m) {
  n <- nrow(Y)
  Y.cent <- bind_cols(dplyr::select(Y, group), as.data.frame(scale(dplyr::select(Y, -group), scale = F)))
  mskj <- Y.cent %>% 
    group_by(group) %>% 
    mutate_at(vars(-group_cols()), function(x) (mean(x))^2) %>% 
    summarise_at(vars(-group_cols()), sum) %>% 
    ungroup() %>% 
    dplyr::select(-group)
  SSbk <- rowSums(mskj, na.rm = T)
  SSbj <- colSums(mskj, na.rm = T)
  SSb <- sum(SSbk, na.rm = T)
  sdiv <- SSb / (n - 1)
  fcsd <- SSbj / SSb
  lcsd <- SSbk / SSb
  out <- list(ss = SSb, sdiv = sdiv, lcss = SSbk, lcsd = lcsd, fcsd = fcsd)
  return(out)
}


### Partitioning spectral diversity ----
# updated to work with terra
specdiv <- function(cube, fact = 40, prop = 0.5, n = 1) {
  # require(raster)
  require(terra)
  require(tidyverse)
  
  # Find minimum number for resampling
  n_pixels <- count_pixels(cube, fact = fact)
  
  # Remove plots with n pixels less than prop
  plot_mask <- subset(n_pixels, 'prop') < prop
  n_pixels_mask <- terra::mask(n_pixels, plot_mask, maskvalue = 1)
  min_pixels <- minmax(n_pixels_mask)[1,1]
  
  # Get cube with plots
  nlayers <- dim(cube)[3]
  cube_plots <- rast(crs = crs(n_pixels_mask))
  ext(cube_plots) <- ext(n_pixels_mask)
  res(cube_plots) <- res(n_pixels_mask)
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  names(cube_plots) <- 'group'
  cube_plots_masked <- terra::mask(cube_plots, plot_mask, maskvalue = 1)
  cube_pixels <- disagg(cube_plots_masked, fact = fact)
  
  # Convert to points
  plots_points <- as.data.frame(cube_plots_masked, xy = T) 
  pixels_points <- as.data.frame(cube_pixels, xy = T)
  
  # Objects to store results
  gamma_ss <- double()
  gamma_sdiv <- double()
  gamma_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(1:n, names(cube)))
  alpha_sdiv <- tibble()
  alpha_fcsd <- tibble()
  alpha_ss <- tibble()
  beta_ss <- double()
  beta_sdiv <- double()
  beta_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(1:n, names(cube)))
  beta_lcsd <- tibble()
  beta_lcss <- tibble()
  
  # Loop to randomly sample min pixels
  for (i in 1:n) {
    cube_points <- as.data.frame(cube, xy = T) %>% 
      inner_join(pixels_points, by = c('x', 'y')) %>% 
      group_by(group) %>% 
      sample_n(size = min_pixels) %>% 
      ungroup()
    xy_all <- cube_points %>% dplyr::select(x,y)
    xy_plots <- plots_points %>% dplyr::select(x,y)
    
    # Gamma diversity
    cube_points_sel_gamma <- cube_points %>%
      dplyr::select(-x, -y, -group)
    sdiv_gamma <- sum_squares(cube_points_sel_gamma)
    gamma_ss[i] <- sdiv_gamma$ss
    gamma_sdiv[i] <- sdiv_gamma$sdiv
    gamma_fcsd[i, ] <- sdiv_gamma$fcsd
    
    # Alpha diversity
    cube_points_sel_alpha <- cube_points %>%
      dplyr::select(-x, -y)
    sdiv_alpha <- cube_points_sel_alpha %>% 
      group_by(group) %>% 
      do(res = sum_squares(Y = dplyr::select(., -group)) )
    # Get sdiv for each community
    alpha_sdiv_tmp <- tibble(rep = i, group = sdiv_alpha$group, sdiv = sapply(sdiv_alpha$res, function(x) x$sdiv))
    # store
    alpha_sdiv <- bind_rows(alpha_sdiv, alpha_sdiv_tmp)
    
    # Get fcsd for each community
    alpha_fcsd_tmp <- bind_cols(rep = rep(i, length(sdiv_alpha$group)), group = sdiv_alpha$group, as.data.frame(t(sapply(sdiv_alpha$res, function(x) x$fcsd))))
    alpha_fcsd <- bind_rows(alpha_fcsd, alpha_fcsd_tmp)
    
    # Get ss for each community
    alpha_ss_tmp <- tibble(rep = i, group = sdiv_alpha$group, ss = sapply(sdiv_alpha$res, function(x) x$ss))
    alpha_ss <- bind_rows(alpha_ss, alpha_ss_tmp)
    
    # Beta diversity
    cube_points_sel_beta <- cube_points %>%
      dplyr::select(-x, -y)
    sdiv_beta <- sum_squares_beta(cube_points_sel_beta, m = min_pixels)
    beta_ss[i] <- sdiv_beta$ss
    beta_sdiv[i] <- sdiv_beta$sdiv
    beta_fcsd[i, ] <- sdiv_beta$fcsd
    beta_lcsd_tmp <- tibble(rep = i, group = 1:length(sdiv_beta$lcsd), lcsd = sdiv_beta$lcsd)
    beta_lcsd <- bind_rows(beta_lcsd, beta_lcsd_tmp)
    beta_lcss_tmp <- tibble(rep = i, group = 1:length(sdiv_beta$lcss), lcss = sdiv_beta$lcss)
    beta_lcss <- bind_rows(beta_lcss, beta_lcss_tmp)
  } 
  
  # Get results together
  
  # LCSD beta
  lcsd_beta_values <- beta_lcsd %>% 
    dplyr::rename(lcsd_beta = lcsd) %>% 
    dplyr::select(-rep) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    dplyr::select(-group)
  lcsd_beta_rast <- rast(cbind.data.frame(xy_plots, lcsd_beta_values), 'xyz', crs = crs(cube))

  
  # LCSS beta
  lcss_beta_values <- beta_lcss %>% 
    dplyr::rename(lcss_beta = lcss) %>% 
    dplyr::select(-rep) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    dplyr::select(-group)
  lcss_beta_rast <- rast(cbind.data.frame(xy_plots, lcss_beta_values), 'xyz', crs = crs(cube))
  
  # FCSD
  fcsd_beta <- colMeans(beta_fcsd, na.rm = T)
  fcsd_gamma <- colMeans(gamma_fcsd, na.rm = T)
  fcsd_alpha_values <- alpha_fcsd %>% 
    dplyr::select(-rep) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    dplyr::select(-group) 
  fcsd_alpha_rast <- rast(cbind.data.frame(xy_plots, fcsd_alpha_values), 'xyz', crs = crs(cube))
  fcsd_alpha_mean <- colMeans(fcsd_alpha_values, na.rm = T)
  
  # SS
  ss_beta <- mean(beta_ss, na.rm = T)
  ss_gamma <- mean(gamma_ss, na.rm = T)
  ss_alpha_sum <- alpha_ss %>%
    dplyr::select(-group) %>% 
    group_by(rep) %>% 
    summarise_all(sum) %>% 
    dplyr::select(ss) %>% 
    summarise_all(mean) %>% 
    as.double()
  
  # SDiv
  sdiv_beta <- mean(beta_sdiv, na.rm = T)
  sdiv_gamma <- mean(gamma_sdiv, na.rm = T)
  sdiv_alpha_values <- alpha_sdiv %>% 
    dplyr::rename(sdiv_alpha = sdiv) %>% 
    dplyr::select(-rep) %>% 
    group_by(group) %>% 
    summarise_all(mean) %>% 
    ungroup() %>% 
    dplyr::select(-group)
  sdiv_alpha_rast <- rast(cbind.data.frame(xy_plots, sdiv_alpha_values), 'xyz', crs = crs(cube))
  sdiv_alpha_mean <- mean(sdiv_alpha_values$sdiv_alpha)
  
  # Prepare outputs
  ss <- tibble(source = c('alpha', 'beta', 'gamma'), sum_squares = c(ss_alpha_sum, ss_beta, ss_gamma) ) %>% 
    mutate(prop_gamma = sum_squares / ss_gamma)
  sdiv <- c(sdiv_alpha_mean, sdiv_beta, sdiv_gamma) ; names(sdiv) <- c('mean_alpha', 'beta', 'gamma')
  fcsd <- bind_cols(source = c('mean_alpha', 'beta', 'gamma'), bind_rows(fcsd_alpha_mean, fcsd_beta, fcsd_gamma))
  rasts <- list(beta_lcsd = lcsd_beta_rast, beta_lcss = lcss_beta_rast, alpha_sdiv = sdiv_alpha_rast, alpha_fcsd = fcsd_alpha_rast)
  out <- list(ss = ss, sdiv = sdiv, fcsd = fcsd, rasters = rasts)
  return(out)
}


