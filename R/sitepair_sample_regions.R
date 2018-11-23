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
#' ap = all_pairs(1:21, 1:21)
#' df = data.frame(LETTERS[ap$Var1], LETTERS[ap$Var2])
#' names(df) = c('a', 'b')
#' tg = rep(3, 21)
#' names(tg) = LETTERS[1:21]
#' tl = rep(0, 21)
#' names(tl) = LETTERS[1:21]
#' op = SelectMismatches(df$a, df$b, tl, tg)
#' df[as.logical(op$selected), ]
#' @export
#' @note Written presently only considering that case where all site pairs are to be considered. 
#' @note But it could be applied on top of other site pair selections (random, geo-weighted). Should build this in. 


# tab = diss
# gen_stats = FALSE
# percent_match = 10
# npairs = 750000
# force_npairs = TRUE
# cap_within_pairs= 1000
# gdm_format = TRUE
# region_name = NULL
# regions = ibra
# field_name = 'REG_CODE_7'
# mask = Aus.domain.mask
# include_zeros = TRUE
# nreps = 25
# tolerance = 5
# verbose = TRUE
  

sitepair_sample_regions = function(tab,
                                   gen_stats = TRUE,
                                   percent_match = 10,
                                   npairs = NULL,
                                   force_npairs = TRUE,
                                   cap_within_pairs= NULL,
                                   gdm_format = TRUE,
                                   region_name = NULL,
                                   regions = NULL,
                                   field_name = NULL,
                                   mask = NULL,
                                   include_zeros = TRUE,
                                   nreps = 25,
                                   tolerance = 5,
                                   write_logs = NULL,
                                   verbose = TRUE){
  
  # for loggin
  args_supplied = as.list(match.call()[-1])
  
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
    
    # r = 32
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
      
      if (!gen_stats){
        if (!is.null(cap_within_pairs)){
          if(n_match > cap_within_pairs){
            to_keep = sample.int(n_match, cap_within_pairs)
            match_df = match_df[to_keep, ]
            percent = pc(cap_within_pairs, n_diff) 
          }
        }
        
        if(percent > percent_match){
          # more matches than target percent - reduce these
          n_match = nrow(match_df)
          reduce_n = reduce_matches(n_match, n_diff)
          to_remove = sample.int(n_match, reduce_n)
          match_df = match_df[-to_remove, ]
          percent = 10
        }
        
      }
        
      n_matches = c(n_matches, nrow(match_df))
      matched = rbind(matched, match_df)
      percent_list = c(percent_list, percent)

    }
    if (verbose)  cat('\r', sprintf('region %s done', regions[r])) 
  }
  
  # build log in stages
  log = data.frame(regions = regions, 
                   starting_within_regions = starting_matches,
                   starting_between_regions = starting_mismatches,
                   starting_percent_within = pc(starting_matches, 
                                             starting_mismatches))
                   
  row.names(log) = NULL
  
  # if a summary only is called
  if(gen_stats){

    cat('gen_stats is TRUE:', 
        '\nReturning summary log of starting matches and mismatches only', 
        sep = '\n')
    return(log)
    
  }
  
  # add first filer to log
  log$first_within_sample = n_matches
  log$percent_within_II = pc(n_matches, log$starting_between_regions)
  
  # targets for regions with single sites?
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
  r2 = as.character(paste(mismatch$r2))
  targets = n_matches * 10
  log$target_between = targets
  names(targets) = regions
  tally = rep(0, length(targets))
  names(tally) = regions
  
  # second loop - returns row ids to keep
  selected = SelectMismatches(r1, r2, tally, targets)
  
  cat('\n')
  
  calc_percentages = pc(n_matches, selected$tally)
  
  # if tally has not been added to (possible), calc_percentages will contain 
  # div 0 Inf 
  calc_percentages = ifelse(is.finite(calc_percentages), calc_percentages, 0)

  log$first_between_sample = unname(selected$tally)
  log$percent_within_III = unname(calc_percentages)
    
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
  
  # I think omit this section... until I work out when it should be called. 
  # if (!all(unname(calc_percentages) == percent_match)) {
  # 
  #   # how bad? Good to have some criteria here, but for now, just reduce matches.
  # 
  #   # Actually, just reduce to a given threshold.
  #   # I think in the test case, this would be +- 5%
  #   
  #   ltol = percent_match - tolerance
  #   utol = percent_match + tolerance
  #   if(mean(calc_percentages) < ltol | mean(calc_percentages) > utol){
  #     
  #     reduce_regions = regions[which(calc_percentages > threshold)]
  #     for(r in seq_along(reduce_regions)){
  # 
  #       r_idx = which(matched$r1 == reduce_regions[r] | 
  #                       matched$r2 == reduce_regions[r])
  # 
  #       ll = length(r_idx)
  #       tl = selected$tally[as.character(paste(reduce_regions[r]))]
  #       pc_now = pc(ll, tl)
  #       
  #       if(ll > 0 & pc_now > threshold){
  #         reduce_by = reduce_matches(ll, tl, by = threshold)
  # 
  #         # try and reduce only from bloated regions
  #         # r_idx = which(matched$r1 == reduce_regions[r] | 
  #         #                 matched$r2 == reduce_regions[r])
  #         pool_r1 = unlist(lapply(reduce_regions, grep, matched$r1))
  #         pool_r2 = unlist(lapply(reduce_regions, grep, matched$r2))
  #         pool_idx = unique(c(pool_r1, pool_r2))
  #         
  #         to_rm = pool_idx[sample.int(length(pool_idx), reduce_by)]
  #         matched = matched[-c(to_rm), ]
  # 
  #       }
  # 
  #       if (verbose)  cat('\r', sprintf('region %s matches reduced', reduce_regions[r]))
  # 
  #     }
  # 
  #   }
  #   
  # }
  
  # select
  combined = mismatch[as.logical(selected$selected), ]
  
  # recombine
  pairs = rbind(matched, combined)
  
  # this will only drop pairs and only when necessary
  if(!is.null(npairs) & !force_npairs){
    
    # Drop pairs to meet npairs target 
    
    if (nrow(pairs) > npairs){
      
      to_drop = sample.int(nrow(pairs), nrow(pairs) - npairs)
      pairs = pairs[-c(to_drop), ]
      
    } 
    
  }
  
  # add more pairs if required
  if(!is.null(npairs) & force_npairs){
    
    # check pair number
    if(nrow(pairs) < npairs){
      
      # reselect from !selected$selected
      reselect = which(selected$selected == 0)
      
      # sample difference from reselect
      diff = npairs - nrow(pairs)
      new_idx = sample(reselect, diff)
      additional = mismatch [new_idx, ]
      
      # recombine
      pairs = rbind(pairs, additional)
      
    }
    
  } 
  
  # re-tally
  final_matches = NULL
  final_mismatches = NULL
  final_percent = NULL
  
  # update this to work on combined table
  for(r in regions){

    m_idx = which(pairs$r1 == r & pairs$match)
    n_idx = which(pairs$r1 == r | pairs$r2 == r & !pairs$match)
                    
    final_matches = c(final_matches, length(m_idx))
    final_mismatches = c(final_mismatches, length(n_idx))
    final_percent = c(final_percent, pc(length(m_idx), length(n_idx)))
    
  }
  
  # recombine
  # pairs = rbind(matched, combined)
  
  # if gdm format, can reshape
  if(gdm_format){
    names(pairs)[1:4] = names(tab)[3:6]
    pairs = join(pairs[, 1:4], tab, by = names(tab)[3:6], match = 'first')
    pairs = pairs[, names(tab)]
  }
  
  log$final_within = final_matches
  log$final_between = final_mismatches
  log$final_percent = final_percent
  log$delta_percent = log$final_percent - log$starting_percent
  
  if(!is.null(write_logs)){
    tag = paste(strsplit(as.character(paste(Sys.time())), ' ')[[1]], 
               collapse = '_')
    tag = gsub(':', '-', tag)
    
    
    
    log_table = sprintf('%s/log_table_%s.csv', write_logs, tag)
    write.csv(log, log_table, row.names = FALSE)
    log_msg = sprintf('%s/log_table_%s.txt', write_logs, tag)
    
    # format this. TODO - add steps
    n = names(args_supplied)
    collect_args = NULL
    for(i in seq_along(args_supplied)){
      collect_args = c(collect_args, 
                       sprintf('%s=%s', n[i], args_supplied[i]))
    } 
    collect_args = paste(collect_args, collapse = '\n')
    sink(log_msg)
    cat('USER SUPPLIED ARGUMENTS\n')
    cat('-----------------------')
    cat(collect_args)
    sink()
  }
  
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
