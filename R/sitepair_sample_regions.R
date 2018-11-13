#' @title Create site pairs balanced across regions
#' @description Balances the number of pairs formed within and between regions 
#' @param tab (data.frame) Table containing (minimally) site coordinates. 
#' @param region_name (char, int, vector) Use this to specify the name(s) of the column(s) containing the region IDs if it exists. Default NULL.
#' @param regions (raster, vector) Filepath or object with which create a regionalisation. If it's in memory, it should be of class raste, sp or sf. If on file, it needs to be able to be opened by raster, rgdal, or sf.
#' @param field_name (char) Default NULL. Can be supplied to denote field of vector feature relating to the region to be used. 
#' @param mask (raster) Default NULL. Must be supplied if regions is a vector feature. 
#' @param include_zeros (bool) Relates to regions which are only represented by one site (i.e there can only be zero 'within' region pairs). If TRUE, include a single cross region pair. 
#' @param nreps (int) number of times to attempt try and balance pairs 
#' @return A list with two objects: balanced site pairs (data.frame) and a log (data.frame)
#' @examples
#' @export
#' 
#' @note Written presently only considering that case where all site pairs are to be considered. 
#' @note But it could be applied on top of other site pair selections (random, geo-weighted). Should build this in. 

sitepair_sample_regions = function(tab,
                                   percent_match = 10,
                                   gdm_format = TRUE,
                                   region_name = NULL,
                                   regions = NULL,
                                   field_name = NULL,
                                   mask = NULL,
                                   include_zeros = TRUE,
                                   nreps = 25,
                                   threshold = 10,
                                   verbose = TRUE){
  
  # FORMAT SECTION
  
  EXTRACT = FALSE
  
  if (gdm_format) {
    # gdm format - expecting dist, weight, long, lat, long, lat
    
    # region specified?
    if (!is.null(region_name)){
      # check there are site 1 and site 2
      stopifnot(length(region_name) == 2)
      
      r_idx = sort(unlist(lapply(region_name, grep, names(tab))))
      
      pairs = tab[, c(3:6, r_idx)]
      names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y', 'r1', 'r2')
      
      # go to site-pair balancing
      
    } else {
      # need regions for site pairs - go to extract
      pairs = tab[, -c(1:2)]
      names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y') 
      EXTRACT = TRUE
      
    }
    
  } else {
    # what does tab look like?
    
    cols = ncol(tab)
    if (cols == 2) {
      # sites only - expect this to be long, lat
      
      # ditch duplicates...
      tab = tab[!duplicated(tab, by = names(tab)),]
      
      # all pairs
      n = 1:nrow(tab)
      pairs = all_pairs(n, n)
      s1 = tab[pairs[, 1], ]
      s2 = tab[pairs[, 2], ]
      pairs = cbind(s1, s2)
      names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y')
      
      EXTRACT = TRUE
      
    } else if (cols == 3 & !is.null(region_name)) {
      # sites with regions - expect this to be long, lat, region
      
      # no dupes
      tab = tab[!duplicated(tab, by = names(tab)[1:2]),]
      
      # all pairs
      n = 1:nrow(tab)
      pairs = all_pairs(n, n)
      s1 = tab[pairs[, 1], ]
      s2 = tab[pairs[, 2], ]
      pairs = cbind(s1, s2)
      pairs = pairs[, c(1:2, 4:5, 3, 6)]
      names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y', 'r1', 'r2')
      
      # go to site-pair balancing
      
    } else {
      
      if (!is.null(region_name)) {
        # take first two columns, and region_name - might fail...
        if(class(region_name) == 'character'){
          idx = grep(region_name, names(tab))
        }
        tab = tab[, c(1, 2, idx)]
        tab = tab[!duplicated(tab, by = names(tab)[1:2]),]
        
        # all pairs
        n = 1:nrow(tab)
        pairs = all_pairs(n, n)
        s1 = tab[pairs[, 1], ]
        s2 = tab[pairs[, 2], ]
        pairs = cbind(s1, s2)
        pairs = pairs[, c(1:2, 4:5, 3, 6)]
        names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y', 'r1', 'r2')
        
        # go to site-pair balancing
        
      } else {
        # just assume to use first two columns
        
        # no duplicates...
        tab = tab[!duplicated(tab, by = names(tab)[1:2]), 1:2]
        
        # all pairs
        n = 1:nrow(tab)
        pairs = all_pairs(n, n)
        s1 = tab[pairs[, 1], ]
        s2 = tab[pairs[, 2], ]
        pairs = cbind(s1, s2)
        names(pairs) = c('s1.x', 's1.y', 's2.x', 's2.y')
        
        EXTRACT = TRUE
        
      }

    }
  
  }
  
  # EXTRACT SECTION
  
  if (EXTRACT){
    
    # determine what regions is
    regions_class = class(regions)
    if (length(grep('SpatialPolygonsDataFrame', regions_class))){
      # coerce to sf and rasterize
      regions = as(regions, 'sf')
      regions = fasterize(regions, mask, field = field_name, fun = 'first')
      
    } else if(length(grep('sf', regions_class))){
        regions = fasterize(regions, mask, field = field_name, fun = 'first')
    }
    # do extract
    r1 = extract(regions, SpatialPoints(cbind(pairs$s1.x, pairs$s1.y)))
    r2 = extract(regions, SpatialPoints(cbind(pairs$s2.x, pairs$s2.y)))
    pairs = cbind(pairs, 'r1' = r1, 'r2' = r2)
    
    # NA check: TODO bolt on coordinate shifter 
    # ditch NA rows for now
    p1 = nrow(pairs)
    pairs = na.omit(pairs)
    p2 = nrow(pairs)
    if(p1 > p2) warning(sprintf('Dropped %s NA rows', p1-p2))
    
  }
  
  # BALANCE SECTION
  
  # pairs formed with regions - check regional balance
  mm = pairs$r1 - pairs$r2
  mtch = rep(TRUE, length(mm))
  mtch[which(mm != 0)] = FALSE
  pairs$match = mtch
  
  # regions now becomes a vector unique region ids
  regions = unique(c(pairs$r1, pairs$r2))
  
  # containers for outputs of first loop
  percent_list = NULL
  n_matches = NULL
  matched = data.frame(matrix(nrow = 0, ncol = ncol(pairs)))
  names(matched) = names(pairs)
  starting_matches = NULL
  starting_mismatches = NULL
  
  # loop over pairs and create matches
  for(r in seq_along(regions)){
    
    r_idx = which(pairs$r1 == regions[r] | pairs$r2 == regions[r]) 
    
    if(length(r_idx)){
      reg_df = pairs[r_idx, ]
      match_df = reg_df[reg_df$match, ]
      nomatch_df = reg_df[!reg_df$match, ]
      n_match = nrow(match_df)
      n_diff = nrow(nomatch_df) 
      percent = pc(n_match, n_diff) 
      starting_matches = c(starting_matches, n_match)
      starting_mismatches = c(starting_mismatches, n_diff)
      
      if(percent > percent_match){
        # more matches than target percent - reduce these
        reduce_n = reduce_matches(n_match, n_diff)
        to_remove = sample.int(n_match, reduce_n)
        match_df = match_df[-to_remove, ]
        percent = 10
      }
      
      n_matches = c(n_matches, nrow(match_df))
      matched = rbind(matched, match_df)
      percent_list = c(percent_list, percent)
    }
    if (verbose)  cat('\r', sprintf('region %s done', regions[r])) 
  }
  
  if (include_zeros){
    
    zeros = which(n_matches == 0)
    if(length(zeros)){
      n_matches[zeros] = 1
    }
    
  }
  
  # mismatches
  mismatch = pairs[pairs$match == FALSE, ]
  
  # randomize 
  random_idx = sample.int(nrow(mismatch), nrow(mismatch))
  mismatch = mismatch[random_idx, ]
  
  # set up inputs for next loop
  r1 = as.character(paste(mismatch$r1))
  r2 = as.character(paste(mismatch$r1))
  targets = n_matches * 10
  names(targets) = regions
  tally = rep(0, length(targets))
  names(tally) = regions
  
  # second loop - returns row ids to keep
  selected = SelectMismatches(r1, r2, tally, targets)
  
  calc_percentages = pc(n_matches, selected$tally)
  
  # if (!all(unname(calc_percentages) == percent_match)) {
  # 
  #   # how bad? Good to have some criteria here, but for now, just continue.
  # 
  #   # reduce > 75% quartile of n_matches by 10%
  # 
  #   upper_quartile = unname(quantile(n_matches)['75%'])
  #   #upper_quartile = unname(quantile(final_matches)['75%'])
  #   reduce_regions = regions[which(n_matches >= upper_quartile)]
  #   for(r in seq_along(reduce_regions)){
  # 
  #     r_idx = which(matched$r1 == reduce_regions[r] | matched$r2 == reduce_regions[r])
  # 
  #     if(length(r_idx)){
  #       fact = floor(0.1 * length(r_idx))
  #       to_rm = sample.int(length(r_idx), fact)
  #       matched = matched[-c(to_rm), ]
  # 
  #       final_matches[as.character(paste(reduce_regions[r]))] =
  #         final_matches[as.character(paste(reduce_regions[r]))] - fact
  # 
  #     }
  # 
  #     if (verbose)  cat('\r', sprintf('region %s re-calculated', reduce_regions[r]))
  # 
  #   }
  # 
  # }
  
  if (!all(unname(calc_percentages) == percent_match)) {

    # how bad? Good to have some criteria here, but for now, just reduce matches.

    # Actually, just reduce to a given threshold.
    # I think in the test case, this would be 10%

    reduce_regions = regions[which(calc_percentages > threshold)]
    for(r in seq_along(reduce_regions)){
      # r = 3
      r_idx = which(matched$r1 == reduce_regions[r] | 
                      matched$r2 == reduce_regions[r])

      ll = length(r_idx)
      tl = selected$tally[as.character(paste(reduce_regions[r]))]
      pc_now = pc(ll, tl)
      
      if(ll > 0 & pc_now > threshold){
        reduce_by = reduce_matches(ll, tl, by = threshold)

        # try and reduce only from bloated regions
        # r_idx = which(matched$r1 == reduce_regions[r] | 
        #                 matched$r2 == reduce_regions[r])
        pool_r1 = unlist(lapply(reduce_regions, grep, matched$r1))
        pool_r2 = unlist(lapply(reduce_regions, grep, matched$r2))
        pool_idx = unique(c(pool_r1, pool_r2))
        
        to_rm = pool_idx[sample.int(length(pool_idx), reduce_by)]
        matched = matched[-c(to_rm), ]

      }

      if (verbose)  cat('\r', sprintf('region %s matches reduced', reduce_regions[r]))

    }

  }
  
  # recombine
  to_keep = as.logical(selected$selected)
  # create more objects while testing...
  combined = mismatch[to_keep, ]
  pairs = rbind(matched, combined)

  # re-tally
  final_matches = NULL
  final_mismatches = NULL
  final_percent = NULL
  
  for(r in seq_along(regions)){
    
    r_idx = which(pairs$r1 == regions[r] | pairs$r2 == regions[r]) 
    
    if(length(r_idx)){
      reg_df = pairs[r_idx, ]
      match_df = reg_df[reg_df$match, ]
      nomatch_df = reg_df[!reg_df$match, ]
      n_match = nrow(match_df)
      n_diff = nrow(nomatch_df) 
      percent = pc(n_match, n_diff) 
      final_matches = c(final_matches, n_match)
      final_mismatches = c(final_mismatches, n_diff)
      
      final_percent = c(final_percent, percent)
      
    }
    
  }
  
  log = data.frame(regions = regions, 
                   start_matches = starting_matches, 
                   start_mismatches = starting_mismatches,
                   start_percent = pc(starting_matches, starting_mismatches),
                   final_matches = final_matches,
                   final_mismatches = final_mismatches,
                   final_percent = final_percent)
  
  return(list(pairs = pairs, log = log))
  
}

# helpers
pc = function(x, y) x/y * 100

reduce_mismatches = function(percent, diff){
  fact = ceiling(10/percent)
  reduce_to = floor(diff/fact)
  return(diff - reduce_to)
}

reduce_matches = function(n_match, diff, by = 10){
  reduce_to = ceiling(diff * (by/100))
  return(n_match - reduce_to)
}

all_pairs <- function(a, b){
  allpairs <- data.table(expand.grid(a, b))
  allpairs <- data.frame(allpairs[allpairs[, .I[1], by = list(pmin(Var1, Var2), 
                                                              pmax(Var1, Var2))]$V1])
  del <- which(allpairs$Var1 == allpairs$Var2)
  unique_pairs <- allpairs[-del,]
  return(unique_pairs)
}
