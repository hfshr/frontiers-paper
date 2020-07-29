splitr <- function(x){cut(x, breaks = c(min(x), median(x), max(x)), labels = c("Low", "High"), 
                          include.lowest = TRUE, ordered_result = TRUE)}


## function to prepare data for network
networkprepare <- function(networkdata){
  filter(networkdata, time == 1) %>% select(-time) -> t1
  filter(networkdata, time == 2) %>% select(-time) -> t2
  filter(networkdata, time == 3) %>% select(-time) -> t3
  explans <- networkdata %>% 
    select(c(id, pi:hours)) %>% 
    names()
  merge(t1, t2, by = explans, suff = c("_1", "_2")) %>% 
    select(-id) -> t12
  merge(t2, t3, by = explans, suff = c("_1", "_2")) %>% 
    select(-id) -> t23
  tallc = rbind(t12, t23)
} # generates tallc files

networkmaker <- function(names){
  networkdata <- datanet %>% 
    select(names) %>% 
    mutate_at(vars(BIS:FFFS, nlebase), splitr)
  tallc <- networkprepare(networkdata)
  
  g1 = tallc %>% select(ends_with("_1")) %>% names()
  g2 = tallc %>% select(ends_with("_2")) %>% names()
  
  explans <- networkdata %>% 
    select(pi:hours, nlebase) %>% 
    names()
  
  bl <-  expand.grid(from = g2, to = g1, stringsAsFactors = FALSE)
  bl <- rbind(bl, expand.grid(from = c(g1, g2), 
                              to = explans, stringsAsFactors = FALSE)) 
  bl <- rbind(bl, data.frame(from = c("clevel", "nlebase"),
                             to = c("gender", "ind_team"), stringsAsFactors = FALSE)) 
  bl <- rbind(bl, data.frame(from = c(g1[2:length(g1)]),
                             to = c(g2[2:length(g2)]), stringsAsFactors = FALSE))
  
  
  wl = data.frame(from = c(paste(vars_select(names(tallc), starts_with("n"))[2]), 
                           paste(vars_select(names(tallc), starts_with("n"))[3])),
                  to = c("injured_1", "injured_2"))
  
  cl = makeCluster(4, type = "SOCK")
  clusterSetRNGStream(cl, 40)
  
  set.seed(155)
  str = boot.strength(tallc, algorithm = "tabu",
                      algorithm.args = list(blacklist = bl, 
                                            #whitelist = wl,
                                            tabu = 50),
                      cluster = cl)
  stopCluster(cl)
  
  return(list(g1 = g1,
              g2 = g2,
              explans = explans,
              tallc = tallc,
              str = str))
  
}

nlesplitr <- function(timenum){
  datanet %>% 
    select(time, nlec) %>%
    filter(time == timenum) %>%
    group_by(time) %>% 
    mutate(lower = log(nlec +1)) %>% 
    mutate_at(vars(lower), 
              function(x) cut(x, breaks = c(min(x), median(x), max(x)), 
                              include.lowest = TRUE, ordered_result = TRUE)) %>% 
    map(~as.numeric(sub('.(.+),.+', '\\1', levels(.x)))[2]) %>% 
    as_tibble() %>% 
    select(-nlec)  %>% 
    mutate(var = paste("nlelg_", timenum, sep = "")) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(nlec) %>% 
            mutate(nlec = log(nlec +1)) %>% 
            map(~min(.)) %>% 
            as_tibble() %>% 
            gather(varr, min) %>% 
            select(-varr)) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(nlec) %>% 
            mutate(nlec = log(nlec +1)) %>% 
            map(~max(.)) %>% 
            as_tibble() %>% 
            gather(varr, max) %>% 
            select(-varr)) %>% 
    mutate(max = round(max, digits = 2)) %>% 
    mutate(Low = paste(min, lower, sep = "-"),
           High = paste(">", lower, sep = ""),
           High = paste(High, max, sep = "-")) %>% 
    mutate(Definition = "") %>% 
    select(var, Definition, Low, High) %>% 
    rename(state_1 = "Low", 
           state_2 = "High")
  
}

tlesplitr <-  function(timenum){
  datanet %>% 
    select(time, tlec) %>%
    filter(time == timenum) %>%
    group_by(time) %>% 
    mutate(lower = log(tlec)) %>% 
    mutate_at(vars(lower), 
              function(x) cut(x, breaks = c(min(x), median(x), max(x)), 
                              include.lowest = TRUE, ordered_result = TRUE)) %>% 
    map(~as.numeric(sub('.(.+),.+', '\\1', levels(.x)))[2]) %>% 
    as_tibble() %>% 
    select(-tlec)  %>% 
    mutate(var = paste("tlelg_", timenum, sep = "")) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(tlec) %>% 
            mutate(tlec = log(tlec)) %>% 
            map(~min(.)) %>% 
            as_tibble() %>% 
            gather(varr, min) %>% 
            select(-varr)) %>% 
    cbind(., datanet %>% 
            filter(time == timenum) %>% 
            select(tlec) %>% 
            mutate(tlec = log(tlec)) %>% 
            map(~max(.)) %>% 
            as_tibble() %>% 
            gather(varr, max) %>% 
            select(-varr)) %>% 
    mutate(max = round(max, digits = 2),
           min = round(min, digits = 2)) %>% 
    mutate(Low = paste(min, lower, sep = "-"),
           High = paste(">", lower, sep = ""),
           High = paste(High, max, sep = "-")) %>% 
    mutate(Definition = "") %>% 
    select(var, Definition, Low, High) %>% 
    rename(state_1 = "Low", 
           state_2 = "High") 
} 

toptail <- function(data, n){
  data %>% 
    top_n(n, prob) %>% 
    rbind(., data %>% 
            top_n(n, -prob))
} # function top and bottom rows

# function for programtically performing CPQs on a given set of vars.
# can combine with map function to loop over a set of vars instead combining them.
newprobtable <- function(tallc, vars, outcome, state, model, repeats = 500000) {
  all.levels = if(any(length(vars) > 1)) { 
    lapply(tallc[, (names(tallc) %in% vars)], levels) 
  } else {
    all.levels <- tallc %>% 
      select(all_of(vars)) %>% 
      sapply(levels) %>% 
      as_tibble() 
  } # makes the code work for when only one variable is used as evidence
  combos <- do.call("expand.grid", c(all.levels, list(stringsAsFactors = FALSE)))  # al combiations
  
  str1 <- "" 
  for (i in seq(nrow(combos))) {
    str1[i] = paste(combos %>% names(), " = '",
                    combos[i, ] %>% sapply(as.character), "'",
                    sep = "", collapse = ", ")
  } # generate character strings for all combinations
  str1 <- rep(str1, times = length(outcome)) # repeat the string for more than one outcome
  str1 <- paste("list(", str1, ")", sep = "")
  
  all.levels.outcome = if(any(length(outcome) > 1)) {
    lapply(tallc[, (names(tallc) %in% outcome)], levels)
  } else {
    all.levels <- tallc %>% 
      select(outcome) %>% 
      sapply(levels) %>% 
      as_tibble()
  } # repeat loop for outcome variables (can have more than one outcome)
  
  combos.outcome <- do.call("expand.grid", c(all.levels.outcome, list(stringsAsFactors = FALSE)))
  
  str3 = rep(paste("(", outcome, " == '", state, "')", sep = ""), each = length(str1)/length(outcome))  # repeat each outcome for the length of combos
  
  fitted <-  bn.fit(avg30, tallc, method = "bayes", iss = 1) # fit the model with bayes method
  

  
  cmd = paste("cpquery(fitted, ", str3, ", ", str1, ", method = 'lw', n = ", repeats, ")", sep = "") # join all elements of string together
 
  
   prob <-  rep(0, length(str1)) # empty vector for probabilities 
  for (i in seq(length(cmd))){
    prob[i] <- eval(parse(text = cmd[i]))
  } # for each combination of strings, what is the probability of outcome
  test <- cbind(combos, prob) %>% 
    mutate(outcome = str3)
  
  
  return(test)
  # output
  # predict <- lm(prob ~ ., data = combos)
  # return_list <- list("combos" = combos, "result" = test, "model" = predict,
  #                     "cmd" = cmd)
  # return(return_list)
}  

## copy from qwraps 2 to have square brackets  
n_perc_x <- function(x, 
                     digits = getOption("qwraps2_frmt_digits", 2), 
                     na_rm = FALSE, 
                     show_denom = "ifNA", 
                     show_symbol = TRUE,
                     markup = getOption("qwraps2_markup", "latex")) { 
  d <- sum(!is.na(x))
  n <- sum(x, na.rm = na_rm)
  p <- frmt(100 * n/d, digits)
  
  if (show_denom == "never") { 
    rtn <- paste0(frmt(as.integer(n)), " [", p, "%]")
  } else { 
    if (show_denom =="always" | any(is.na(x))) { 
      rtn <- paste0(frmt(as.integer(n)), "/", frmt(as.integer(d)), " [", p, "%]")
    } else { 
      rtn <- paste0(frmt(as.integer(n)), " [", p, "%]")
    }
  }
  
  if (!show_symbol) { 
    rtn <- gsub("%", "", rtn)
  }
  
  
  if (markup == "latex") { 
    rtn <- gsub("%", "\\\\%", rtn)
  } 
  
  return(rtn)
}
