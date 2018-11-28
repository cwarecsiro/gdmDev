#' @title Create site pairs balanced across regions
#' @description Balances the number of pairs formed within and between regions 
#' @param sites (data.frame) Table containing site coordinates (long, lat) and a region identifier. 
#' @param target_pairs (int, optional). Number of requested site pairs. 
#' @param target_ratio (float, default 0.1). Ratio of within to between regional matches.
#' @param generate_statistics_only (bool, default FALSE). Optionally, simply calcualate parameter space.
#' @return A list with two objects: balanced site pairs (data.frame) and a log (data.frame)
#' @examples
#' sites_per_region = c(0, 3, 5, 7, 9, 12)
#' td = gen_testdata(sites_per_region)
#' # visual
#' plot(td$map)
#' points(td$sites, col = 'red', pch = 19)
#' # check parameter space only
#' sitepair_sample_regions(td$sites, target_pairs = 250, generate_statistics_only = TRUE)
#' @export
#' @note Presently requires numeric region 'names.' This dependency needs to be removed. 
sitepair_sample_regions = function(sites, 
                                   target_pairs = NULL, 
                                   target_ratio = 0.1, 
                                   generate_statistics_only = FALSE, 
                                   write_sitepairs = NULL,
                                   write_logs = NULL){
  
  # for logging
  args_supplied = as.list(match.call()[-1])
  input_table = deparse(substitute(sites))
  
  bounds = gen_stats(sites, target_pairs = target_pairs, 
                     target_ratio = target_ratio,
                     generate_statistics_only = generate_statistics_only)
  
  if(generate_statistics_only){
    
    return(bounds)
    
  }
  
  # set up targets for within and between
  targets = make_targets(bounds$output_mat, 
                         target_pairs = target_pairs, 
                         target_ratio = target_ratio, 
                         n_regions = bounds$n_regions)
  
  # fill within
  step1 = fill_diagonal(sites, bounds$output_mat, bounds$reg_col, 
                        bounds$reg_names, targets$target_within, 
                        bounds$regional_sites)
  
  # fill between
  step2 = fill_between(sites, step1$filled_mat, bounds$reg_col, 
                       bounds$reg_names, targets$target_between, 
                       bounds$regional_sites, 
                       step1$within_a, step1$within_b)
  
  # index sites and form pairs
  s1 = sites[step2$within_a, ]
  s2 = sites[step2$within_b, ]
  names(s1) = paste0(names(s1), '_1')
  names(s2) = paste0(names(s2), '_2')
  sitepairs = cbind(s1, s2)
  row.names(sitepairs) = NULL
  
  if(!is.null(write_sitepairs)){
    
    ext = strsplit(basename(write_sitepairs), '\\.')[[1]][2]
    valid_ext = c('feather', 'csv', 'RData')
    if(!ext %in% valid_ext){
      stop(sprintf('%s is not a valid output type. Must be one of %s', 
                   write_sitepairs, paste(valid_ext, collapse = ' ')))
    }
      
    if(ext == 'feather'){
      write_feather(sitepairs, write_sitepairs)
    }
    
    if(ext == 'csv'){
      write.csv(sitepairs, write_sitepairs, row.names = FALSE)
    }
    
    if(ext == 'RData'){
      save(sitepairs, file = write_sitepairs)
    }

  }
  
  if(!is.null(write_logs)){
    
    dst = sprintf('%s/%s_%s', write_logs, input_table, 
                  args_supplied$target_pairs)
    
    n = names(args_supplied)
    collect_args = NULL
    for(i in seq_along(args_supplied)){
      collect_args = c(collect_args, 
                       sprintf('%s=%s', n[i], args_supplied[i]))
    } 
    collect_args = paste(collect_args, collapse = '\n')
    
    perc = calc_proportion(step2$filled_mat)
    
    nmat = data.matrix(step2$filled_mat)
    check_within = diag(nmat)
    check_within = ifelse(check_within < targets$target_within, FALSE, TRUE)
    
    diag(nmat) = NA
    check_between = na.omit(as.vector(nmat))
    check_between = ifelse(check_between < targets$target_between, FALSE, TRUE)
    
    msg = list(m0 = '------------------------------------------------------\n', 
               m1 = sprintf('Sample set with: \n\t%s sites', nrow(sites)),
               m2 = sprintf('\n\t%s possible pairs', N_pairs(nrow(sites))), 
               m3 = sprintf('\n\t%s regions', bounds$n_regions), 
               m4 = sprintf('\n\nUsing a target of %s pairs and a %s within/between ratio', 
                            target_pairs, target_ratio),
               m5 = sprintf('\n\taiming for %s/%s within/between pairs', 
                            targets$target_within, targets$target_between),
               m6 = sprintf('\n\t%s regions met their within targets (%s)', 
                            sum(check_within), targets$target_within),
               m7 = sprintf('\n\t%s between region combinations met their between targets (%s)', 
                            sum(check_between), targets$target_between),
               m8 = sprintf('\n\tFinal percentage of within pairs: %s%%', perc),
               m8 = '\n------------------------------------------------------\n')
    
    sink(sprintf('%s_callsummary.txt', dst))
    cat('USER SUPPLIED ARGUMENTS\n')
    cat('-----------------------\n')
    cat(collect_args)
    cat('\n\nCALL SUMMARY \n')
    for(i in msg) cat(i)
    sink()
    
    
    write.csv(bounds$output_mat, sprintf('%s_MATRIX_ALLPAIRS.csv', dst), 
              row.names = TRUE)
    
    write.csv(step2$filled_mat, sprintf('%s_MATRIX_SAMPLEDPAIRS.csv', dst), 
              row.names = TRUE)
    
  }
  
  return(list(sitepairs = sitepairs,
              matrix_allpairs = bounds$output_mat,
              matrix_sampledpairs = step2$filled_mat, 
              target_within = targets$target_within, 
              target_between = targets$target_between))
  
}


FormPairs_ = function(ns, indexes_only = TRUE, as_df = FALSE, 
                      sites = NULL){
  
  if(!hasArg(ns)){
    stop('Need number of sites')
  }
  
  if(indexes_only & as_df){
    stop('Cannot return indexes as both a list and a data.frame')
  }
  
  pairs = FormPairs(ns)
  
  if(indexes_only){
    
    return(pairs)
    
  } else if (as_df){
    
    return(data.frame(pairs))
    
  } else {
    
    s1 = sites[pairs$pair1, ]
    s2 = sites[pairs$pair2, ]
    names(s1) = paste0(names(s1), '_1')
    names(s2) = paste0(names(s2), '_2')
    
    return(cbind2(s1, s2))
    
  }
  
}

gen_stats = function(sites, 
                     target_pairs = NULL, 
                     target_ratio = 0.1,
                     generate_statistics_only = FALSE){
  
  reg_col = names(sites)[3]
  tbl = table(sites[reg_col])
  reg_names = names(tbl)
  regional_sites = as.vector(tbl)
  n_regions = length(reg_names)
  
  output_mat = data.frame(matrix(nrow = n_regions, ncol = n_regions))
  names(output_mat) = reg_names
  row.names(output_mat) = reg_names
  
  # within
  pairs_within = sapply(regional_sites, N_pairs)
  for (i in seq_along(reg_names)){
    output_mat[i, i] = pairs_within[i]
  }
  
  # between
  pairs_between = list()
  for (i in 1:(n_regions-1)){
    n_site_i = regional_sites[i]
    reg_i = NULL
    for (j in (i+1):n_regions){
      reg_i = c(reg_i, n_site_i*regional_sites[j])
    }
    pairs_between[[i]] = reg_i
  }
  
  for (i in 1:(n_regions-1)){
    output_mat[(i+1):n_regions, i] = unlist(pairs_between[[i]])
  }
  
  output_mat[is.na(output_mat)] = ''
  
  if(generate_statistics_only){
    
    all_n = N_pairs(nrow(sites))
    
    if(!is.null(target_pairs)){
      tg50 = target_pairs
    } else {
      tg50 = ceiling(all_n/2)
    }
    tr10 = ceiling(target_ratio * tg50)
    trw = ceiling(tr10/n_regions)
    
    sampled_mat = output_mat
    
    tg_within = NULL
    for(i in 1:ncol(output_mat)){
      ii = as.numeric(output_mat[i,i])
      tg_within = c(tg_within, ii >= trw)
      
      if(ii < trw){
        sampled_mat[i,i] = ii
      } else{
        sampled_mat[i,i] = trw
      }
      
    }
    
    trb = ceiling((tg50-tr10)/N_pairs(n_regions))
    
    # rounding might make these number not add up...
    tg50 = trb * N_pairs(n_regions) + n_regions * trw
    
    tg_between = NULL
    for(i in 1:(ncol(output_mat)-1)){
      
      for(j in (i+1):nrow(output_mat)){
        
        ij = as.numeric(output_mat[j,i])
        tg_between = c(tg_between, ij >= trb)
        
        if(ij < trb){
          sampled_mat[j,i] = ij
        } else{
          sampled_mat[j,i] = trb
        }
          
      }
      
    }
    
    perc = calc_proportion(sampled_mat)
    
    msg = list(m0 = '------------------------------------------------------\n', 
               m1 = sprintf('Sample set with: \n\t%s sites', nrow(sites)),
               m2 = sprintf('\n\t%s possible pairs', all_n), 
               m3 = sprintf('\n\t%s regions', n_regions), 
               m4 = sprintf('\n\nUsing a target of %s pairs and a %s within/between ratio', tg50, target_ratio),
               m5 = sprintf('\n\taiming for %s pairs: %s/%s', tg50, n_regions*trw, N_pairs(n_regions)*trb),
               m6 = sprintf('\n\t%s regions would meet their within targets (%s)', 
                            sum(tg_within), trw),
               m7 = sprintf('\n\t%s between region combinations would meet their between targets (%s)', 
                            sum(tg_between), trb),
               m8 = sprintf('\n\tFinal percentage of within pairs: %s%%', perc),
               m8 = '\n------------------------------------------------------\n')
    
    for(i in msg) cat(i)
    
    return(list(matrix_allpairs = output_mat))
    
  } else {
    
    return(list(output_mat = output_mat, 
                reg_names = reg_names,
                n_regions = n_regions,
                regional_sites = regional_sites,
                reg_col = reg_col))
  }
  
  
}

N_pairs = function(n) ((n^2)-n)/2

make_targets = function(output_mat, target_pairs, target_ratio, n_regions){
  
  target_within = target_pairs * target_ratio
  target_region = ceiling(target_within/n_regions)
  target_between = ceiling((target_pairs-target_within) / N_pairs(n_regions))
  
  return(list(target_within = target_region,
              target_between = target_between))
  
}

fill_diagonal = function(sites, output_mat, reg_col, reg_names, target_within, 
                         regional_sites){

  # output_mat = bounds$output_mat
  # reg_names = bounds$reg_names
  # target_within = targets$target_within
  # regional_sites = bounds$regional_sites
  
  
  filled_mat = output_mat
  
  within_a = NULL
  within_b = NULL
  
  for(i in 1:nrow(output_mat)){
    
    # i = 4
    # reg_col = 'region'
    
    reg_i = reg_names[i]
    
    possible = as.numeric(output_mat[i,i])
    
    if(possible != 0){
      
      rows_i = which(sites[reg_col] == reg_i)
      
      ns = regional_sites[i]
      pairs_i = FormPairs_(ns, indexes_only = TRUE)
      
      if(possible <= target_within){
        
        # take all
        idx = rows_i[pairs_i$pair1]
        within_a = c(within_a, idx)
        idx = rows_i[pairs_i$pair2]
        within_b = c(within_b, idx)
        
        filled = length(idx)
        
        
      } else {
        
        # sample
        sample_idx = sample.int(as.numeric(possible), target_within)
        idx = rows_i[pairs_i$pair1[sample_idx]]
        within_a = c(within_a, idx)
        idx = rows_i[pairs_i$pair2[sample_idx]]
        within_b = c(within_b, idx)
        
        filled = length(idx)
        
      }
      
    } else {
      
      filled = 0
      
    }
    
    filled_mat[i, i] = filled
    
  }  
  
  return(list(filled_mat = filled_mat, 
              within_a = within_a, 
              within_b = within_b))
  
}

fill_between = function(sites, filled_mat, reg_col, reg_names, target_between, 
                         regional_sites, within_a, within_b){
  
  # filled_mat = step1$filled_mat
  # reg_names = bounds$reg_names
  # target_between = 5
  # regional_sites = bounds$regional_sites
  
  for(i in 1:(ncol(filled_mat)-1)){
    
    # i = 2
    reg_i = reg_names[i]
    
    for(j in (i + 1):nrow(filled_mat)){
      
      # j = i + 1
      possible = as.numeric(filled_mat[j,i])
      reg_j = reg_names[j]
      
      rows_ij = which(sites[reg_col] == reg_i | 
                        sites[reg_col] == reg_j)
      
      regions_ij = sites[rows_ij, reg_col]
      
      pairs_ij = FormPairs2(regions_ij, possible)
      
      if(possible <= target_between){
        
        # take all
        idx = rows_ij[pairs_ij$pair1]
        within_a = c(within_a, idx)
        idx = rows_ij[pairs_ij$pair2]
        within_b = c(within_b, idx)
        
        filled = length(idx)
        
        
      } else {
        
        # sample
        sample_idx = sample.int(as.numeric(possible), target_between)
        idx = rows_ij[pairs_ij$pair1[sample_idx]]
        within_a = c(within_a, idx)
        idx = rows_ij[pairs_ij$pair2[sample_idx]]
        within_b = c(within_b, idx)
        
        filled = length(idx)
        
      }
      
      filled_mat[j, i] = filled
      
    }
    
  }  
  
  return(list(filled_mat = filled_mat, 
              within_a = within_a, 
              within_b = within_b))
  
}

gen_testdata = function(sites_per_region){
  
  ll = length(sites_per_region)
  
  N = ceiling(sqrt(ll))
  rec = ll %% N
  while(!rec == 0){
    N = N+1
    rec = ll %% N
  }
  nr = ll/N
  rast = raster(nrow = nr, ncol = N)
  
  rast[] = 1:ll
  
  polys = rasterToPolygons(rast)
  bbs = lapply(polys@polygons, bbox)
  points = data.frame(matrix(nrow = 0, ncol = 2))
  names(points) = c('x', 'y')
  
  for (i in seq_along(sites_per_region)){
    
    x = runif(sites_per_region[i], bbs[[i]][1,1], bbs[[i]][1,2])
    y = runif(sites_per_region[i],  bbs[[i]][2,1], bbs[[i]][2,2])
    
    points = rbind(points, data.frame(x, y))  
    
  }
  
  points$region = extract(rast, points)
  
  plot(rast, legend = FALSE)
  points(points, pch = 19, col = 'red')
  
  return(list(sites = points, map = rast))
  
}

add_region = function(sites, regions, force_values = TRUE){
  
  regions_class = class(regions)
  if (length(grep('SpatialPolygonsDataFrame', regions_class))){
    # coerce to sf and rasterize
    regions = as(regions, 'sf')
    regions = fasterize(regions, mask, field = field_name, fun = 'first')
    
  } else if(length(grep('sf', regions_class))){
    regions = fasterize(regions, mask, field = field_name, fun = 'first')
  }
  # do extract
  reg = extract(regions, SpatialPoints(sites))
  
  # coord shifter...
  if(force_values){
    
    if(anyNA(reg)){
      
      na_idx = which(is.na(reg))
      na_coords = SpatialPoints(sites[na_idx, ])
      
      buff_width = res(regions)[1]
      if(length(grep('+proj=longlat', crs(regions)))){
        buff_width = buff_width * 100000
      }
      
      valid = NULL
      suppressWarnings(
        
        for (i in seq_along(na_coords)){
          #i = 1
          check = FALSE
          while(!check){
            
            fish = extract(regions, buffer(na_coords[i], buff_width))[[1]]
            if(any(!is.na(fish))){
              valid = c(valid, na.omit(fish)[1])
              check = TRUE
            
            } else {
              
              buff_width  = buff_width * 2
                
            }
            
          }
        
        }
      
      )
        
    }
    
    reg[na_idx] = valid
    sites$region = reg
  
  } else {
    
    sites$region = reg
    p1 = nrow(sites)
    sites = na.omit(sites)
    p2 = nrow(sites)
    if(p1 > p2) warning(sprintf('Dropped %s NA rows', p1-p2))
    
  }
  
  return(sites)
  
}

calc_proportion = function(df){
  
  xmat = data.matrix(df)
  sum_within = 0
  for(i in 1:ncol(xmat)){
    sum_within = sum_within + xmat[i,i]
    xmat[i,i] = 0
  }
  
  sum_between = sum(xmat, na.rm = TRUE)
  
  return(sum_between / sum_within)
  
}