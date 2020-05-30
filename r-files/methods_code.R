## Methods
pacman::p_load(tidyverse, knitr, kableExtra, qwraps2, brms, broom.mixed, bnlearn, tidyselect)
library(readxl)
library(janitor)
library(irr)
library(lavaan)
library(caret)
library(parallel)
library(zoo)
library(sjmisc)
## Table 1

x <- read_csv("data/new_cleandata_inj.csv")

id_remove <- x %>% 
  group_by(id) %>% 
  filter(sum(time) == 1) %>% 
  pull(id)

x %>% 
  filter(!id %in% id_remove) %>% 
  group_by(time, id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(time) %>% 
  count()

parchar <- read_csv("data/participant_chars.csv") %>% 
  filter(!id %in% id_remove)

clevel <- parchar %>%
  group_by(gender) %>%
  summarise(
    Recreational = n_perc0(clevel == "Recreational"),
    University = n_perc0(clevel == "University"),
    "National/International" = n_perc0(clevel == "National/International")
  ) %>%
  gather(name, value, -gender) %>%
  spread(gender, value) %>%
  ungroup() %>%
  arrange(match(name, c("Recreational", "University", "National/International")))


npars <- parchar %>%
  group_by(gender) %>%
  summarise(
    Age = mean_sd(age, denote_sd = "paren", digits = 1),
    Height = mean_sd(height, denote_sd = "paren", digits = 1),
    Weight = mean_sd(weight, denote_sd = "paren", digits = 1),
    "Hours per week" = mean_sd(hours, denote_sd = "paren", digits = 1),
    "Student" = n_perc0(student == "y"),
    "Non-student" = n_perc0(student == "n"),
    "Injured" = n_perc0(pic == "1"),
    "Healthy" = n_perc0(pic == "0")
  ) %>%
  gather(name, value, -gender) %>% 
  arrange(gender)%>%
  spread(gender, value) %>%
  rbind(., clevel) %>%
  rename(
    "Male (n = 231) " = "male",
    "Female (n = 120)" = "female"
  ) %>%
  ungroup()

rm(clevel)

## RST reliability

d <- read_csv("data/rst/51itemdata.csv", col_names = T) %>%
  drop_na() %>%
  select(-id) %>%
  select(c(2:52, 1))

colnames(d)[45:51] <- str_replace(names(d[45:51]), "I", "IMP")

# to get character strings
fffs <- paste(names(d[grep("FFFS", colnames(d))]), collapse = "+")
bis <- paste(names(d[grep("BIS", colnames(d))]), collapse = "+")
ri <- paste(names(d[grep("RI", colnames(d))]), collapse = "+")
gdp <- paste(names(d[grep("GDP", colnames(d))]), collapse = "+")
rr <- paste(names(d[grep("RR", colnames(d))]), collapse = "+")
i <- paste(names(d[grep("IMP", colnames(d))]), collapse = "+")


model_fffs <- paste("fffs", fffs, sep = "=~")
model_bis <- paste("bis", bis, sep = "=~")

# full model with all subscales
fullmodel <- "fffs=~V9_FFFS+V39_FFFS+V45_FFFS+V46_FFFS+V48_FFFS+V52_FFFS+V58_FFFS+V59_FFFS 
bis=~V1_BIS+V2_BIS+V6_BIS+V10_BIS+V17_BIS+V18_BIS+V21_BIS+V29_BIS+V33_BIS+V43_BIS+V47_BIS+V50_BIS+V56_BIS+V57_BIS+V60_BIS+V61_BIS+V63_BIS
ri=~V14_RI+V15_RI+V32_RI+V35_RI
gdp=~V5_GDP+V12_GDP+V20_GDP+V31_GDP+V41_GDP+V54_GDP+V65_GDP
rr=~V4_RR+V8_RR+V16_RR+V23_RR+V25_RR+V30_RR+V36_RR+V37_RR
i=~V22_IMP+V27_IMP+V28_IMP+V38_IMP+V40_IMP+V44_IMP+V51_IMP"

fullmodel <- str_remove_all(fullmodel, "_")

# example of cfa model
fit <- cfa(fullmodel, data = d)
sl <- standardizedSolution(fit)
sl <- sl$est.std[sl$op == "=~"]
names(sl) <- names(d[1:51])

rex <- function(x) {
  re <- 1 - x^2
  re
}

rst_rel <- enframe(sl) %>%
  rownames_to_column("id") %>%
  mutate(re = map(.x = value, ~ rex(.x))) %>%
  mutate(re = as.numeric(re)) %>%
  mutate(nnames = str_extract_all(name, "[A-Z]+", simplify = T)[, 2]) %>%
  group_by(nnames) %>%
  summarise(final = sum(value)^2 / (sum(value)^2 + sum(re)))

rm(fit, bis, fffs, fullmodel, gdp, i, model_bis, model_fffs, ri, rr, sl, d, rex)

#### balance icc

dat <- read_excel("data/balancecom/balancecom.xlsx") %>%
  clean_names()

d <- dat %>%
  select(participant, grand_total, set) %>%
  filter(set == 1) %>%
  cbind(
    .,
    dat %>%
      select(participant, grand_total, set) %>%
      filter(set == 2) %>%
      rename(grand_total_2 = "grand_total") %>%
      select(grand_total_2)
  ) %>%
  select(-set, -participant)

icc <- irr::icc(d, model = "twoway", type = "agreement")

##



data <- read_csv("data/fulldatajoined.csv") %>%
  filter(time != 4) # remove last time point as only questionnaire data is available

####################### Cleaning up variables + impute #################################

# all possible explanatory
explan <- data %>%
  select(c(1, 2, 7, 9:15, 17, 19, 37:47))

# basic set of explanatory
basic.explan <- explan %>%
  select(id, time, variable, gender, sportg, clevel, hours, pic) %>%
  mutate(
    pic = factor(pic, labels = c("no injury", "injury")),
    clevel = ifelse(clevel == "notplaying", "club_university_county", clevel)
  ) %>%
  group_by(id) %>%
  mutate(hours = na.locf(hours)) %>%
  ungroup()

# impute missing data

# select variables of interest from full data set and join with explanatory variables
voi <- data %>%
  select(
    id,
    time,
    variable,
    injured,
    nle,
    ple,
    tle,
    st,
    nst,
    gt,
    ndt,
    dt,
    stiffness,
    frequency,
    decrement,
    rmssd,
    sdnn,
    meanHR,
    bis,
    bas,
    fffs
  ) %>%
  left_join(basic.explan, ., by = c("id", "time", "variable")) %>%
  mutate(
    gt = ifelse(gt == 0, NA, gt),
    bal_asym = abs(ndt - dt)
  )

# additional variables that wont be imputed
ids <- data %>%
  select(id, time, variable, change, delta, pre, post, days_missed)

## check how many missing

missingcount <- tibble(
  var = c("hrv", "myoton"),
  missing = c(map(voi, ~ sum(is.na(.)))$rmssd, map(voi, ~ sum(is.na(.)))$stiffness),
  total = c(668, 668)
) %>%
  mutate(percdiff = missing / total * 100)

missingtot <- voi %>%
  select(id, rmssd, stiffness) %>%
  mutate(
    rmssd = tidyr::replace_na(rmssd, 0),
    stiffness = tidyr::replace_na(stiffness, 0)
  )

## impute missing values using bag impute

prepro <- preProcess(voi[, c(4:length(voi))], "bagImpute")
datanew <- predict(prepro, voi[, c(4:length(voi))]) %>%
  cbind(ids, .)

lescacounts <- read_excel("data/lesca/lescacounts.xlsx", sheet = "Sheet3")

datanew <- datanew %>%
  inner_join(., lescacounts, by = c("id", "time")) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    negsev = nle / totneg,
    totsev = tle / sum(totneg, totpos)
  ) %>%
  mutate(
    negsev = ifelse(negsev == "NaN", 0, negsev),
    totsev = ifelse(totsev == "NaN", 0, totsev)
  ) %>%
  ungroup()

rm(voi, ids, prepro)

### Prepare data for network

datanet <- datanew %>%
  select(-bis, -bas, -fffs, -delta, -post, -pre, -change) %>%
  left_join(., read_csv("data/rst/finalrst.csv"), by = c("id", "time")) %>% # updated 51 item rst
  filter(variable == 1) %>%
  rename(
    "balance" = gt,
    "pi" = pic
  ) %>%
  mutate(
    injured = factor(ifelse(injured == 0, "healthy", "injured")),
    clevel = factor(
      clevel,
      levels = c(
        "club_university_county",
        "national",
        "international"
      )
    ),
    clevel = fct_collapse(clevel, "national_international" = c("national", "international")),
    pi = factor(pi),
    gender = factor(gender),
    sportg = factor(sportg),
    hours = factor(ifelse(hours < median(hours), "Low", "High"), levels = c("Low", "High")),
    sportg = factor(sportg),
    ind_team = fct_collapse(
      sportg,
      "individual" = c("athletics", "gym", "other"),
      "team" = c(
        "football",
        "hockey",
        "netball",
        "cricket",
        "rugby",
        "basketball"
      )
    )
  ) %>%
  group_by(id) %>%
  mutate(
    nlec = cumsum(nle),
    tlec = cumsum(tle),
    nlebase = ifelse(time == 1, nle, NA),
    nlebase = na.locf(nlebase),
    injnum = ifelse(injured == "injured", 1, 0),
    injcount = sum(injnum),
    oneinj = ifelse(injcount >= 1, 1, 0),
    oneinj = factor(oneinj, levels = c(0, 1), labels = c("healty", "injured"))
  ) %>%
  ungroup() %>%
  group_by(time) %>%
  mutate(
    nlel = log(nlec + 1),
    tlel = log(tlec)
  ) %>%
  ungroup() %>%
  mutate(negsev = ifelse(negsev == "NaN", 0, negsev)) %>%
  ungroup() %>%
  mutate(rmssd = round(log(rmssd), digits = 2)) %>%
  group_by(time) %>%
  mutate(
    totneg_d = splitr(totneg),
    "nlelg" = splitr(log(nlec + 1)),
    "tlelg" = splitr(log(tlec)),
    nlel = splitr(log(nle + 1)),
    tlel = splitr(log(tle + 1))
  ) %>%
  ungroup()



nlename <- list(
  nlelg = datanet %>%
    select(
      id,
      time,
      injured,
      nlelg,
      pi,
      nlebase,
      gender,
      ind_team,
      clevel,
      hours,
      BIS,
      stiffness,
      RR,
      I,
      GDP,
      RI,
      sdnn,
      rmssd,
      balance,
      FFFS
    ) %>% names(),
  tlelg = datanet %>%
    select(
      id,
      time,
      injured,
      tlelg,
      pi,
      nlebase,
      gender,
      ind_team,
      clevel,
      hours,
      BIS,
      stiffness,
      RR,
      I,
      GDP,
      RI,
      sdnn,
      rmssd,
      balance,
      FFFS
    ) %>% names(),
  both = datanet %>%
    select(
      id,
      time,
      injured,
      nlelg,
      tlelg,
      pi,
      nlebase,
      gender,
      ind_team,
      clevel,
      hours,
      BIS,
      stiffness,
      RR,
      I,
      GDP,
      RI,
      sdnn,
      rmssd,
      balance,
      FFFS
    ) %>% names()
)
rstnames <- list(BAS = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          BAS, rmssd,sdnn, balance, FFFS) %>% names(),
                 RI = datanet %>% 
                   select(id, time, injured, nlelg, tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RI, rmssd,sdnn, balance, FFFS) %>% names(),
                 GDP = datanet %>% 
                   select(id, time, injured,nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          GDP, rmssd,sdnn, balance, FFFS) %>% names(),
                 I = datanet %>% 
                   select(id, time, injured, nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          I, rmssd,sdnn, balance, FFFS) %>% names(),
                 RR = datanet %>% 
                   select(id, time, injured,nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RR, rmssd,sdnn, balance, FFFS) %>% names(),
                 ALL = datanet %>% 
                   select(id, time, injured, nlelg,tlelg, pi, nlebase, gender, ind_team, clevel, hours, BIS, stiffness, 
                          RR, I, GDP, RI, rmssd,sdnn, balance, FFFS) %>% names())

scorer <- function(names) {
  test <- map(names, ~ networkmaker(.x))
  avg <- map(.x = names(test), ~ averaged.network(test[[.x]][[5]], threshold = 0.35)) %>%
    set_names(c(names(test)))
  "Variable" <- names(test)
  x <- map_df(names(test), ~ data.frame(score = bnlearn::score(avg[[.x]], test[[.x]][[4]], type = "bic"))) %>%
    cbind(Variable, .) %>%
    arrange(-score)
}


le <- scorer(nlename)
rstmodel <- scorer(rstnames)


#### Code used for table 3

networkdata <- datanet %>%
  select(
    id,
    time,
    injured,
    nlelg,
    pi,
    nlebase,
    gender,
    ind_team,
    clevel,
    hours,
    stiffness,
    RI,
    BIS,
    rmssd,
    balance,
    FFFS
  ) %>%
  mutate_at(vars(stiffness:FFFS, nlebase), splitr)

## prepare 2 turn network
tallc <- networkprepare(networkdata)

names <- datanet %>%
  select(nlebase, stiffness:FFFS) %>%
  names()

varcutoffs <- tallc %>%
  select(pi:hours) %>%
  select(-nlebase) %>%
  map(~ levels(.)) %>%
  as_tibble() %>%
  gather(var, state) %>%
  mutate(ind = rep(c("state_1", "state_2"), length.out = n())) %>%
  spread(ind, state) %>%
  mutate(Definition = c(
    "Current competitive level",
    "Gender of the participant",
    "Number of hours spent training per week",
    "Participate in an individual or team based sport",
    "Previous injury - Whether an injury had been sustained in the previous 12 months prior to the study"
  )) %>%
  mutate_at(vars(state_1, state_2), funs(str_to_title(.))) %>%
  select(var, Definition, state_1, state_2) %>%
  rbind(
    .,
    datanet %>%
      select(
        id,
        time,
        injured,
        pi,
        nlebase,
        gender,
        ind_team,
        clevel,
        hours,
        stiffness,
        RI,
        RR,
        I,
        GDP,
        BIS,
        rmssd,
        sdnn,
        bal_asym,
        balance,
        FFFS
      ) %>%
      mutate_at(
        vars(stiffness:FFFS, nlebase),
        function(x) {
          cut(
            x,
            breaks = c(min(x), median(x), max(x)),
            include.lowest = TRUE,
            ordered_result = TRUE
          )
        }
      ) %>%
      select(stiffness:FFFS, nlebase) %>%
      map(~ as.numeric(sub(".(.+),.+", "\\1", levels(.x)))[2]) %>%
      as_tibble() %>%
      gather(var, lower) %>%
      {
        . ->> namess
      } %>%
      cbind(
        .,
        datanet %>%
          select(namess$var) %>%
          map(~ min(.)) %>%
          as_tibble() %>%
          gather(varr, min) %>%
          select(-varr)
      ) %>%
      cbind(
        .,
        datanet %>%
          select(namess$var) %>%
          map(~ max(.)) %>%
          as_tibble() %>%
          gather(varr, max) %>%
          select(-varr)
      ) %>%
      mutate(
        Low = paste(min, lower, sep = "-"),
        High = paste(">", lower, sep = ""),
        High = paste(High, max, sep = "-")
      ) %>%
      select(var, Low, High) %>%
      mutate(var = factor(
        var,
        levels = c(
          "nlebase",
          "FFFS",
          "BIS",
          "RI",
          "RR",
          "I",
          "GDP",
          "stiffness",
          "rmssd",
          "sdnn",
          "bal_asym",
          "balance"
        )
      )) %>%
      arrange(var) %>%
      mutate(Definition = c(
        "Untransformed NLE at TP 1",
        "Fight-Flight-Freeze System",
        "Behavioural Inhibition System",
        "Reward Interest",
        "Reward reactivity",
        "Impulsivity",
        "Goal drive persistence",
        "Sum of all stiffness locations",
        "Root mean squared difference of successive RR intervals",
        "Standard deviation of RR series",
        "Percentage difference between left and right leg balance score",
        "Total balance score"
      )) %>%
      rename(
        state_1 = "Low",
        state_2 = "High"
      ) %>%
      select(var, Definition, state_1, state_2)
  ) %>%
  rbind(., map_df(1:3, nlesplitr)) %>%
  rbind(., map_df(1:3, tlesplitr)) %>%
  rowid_to_column("tempid") %>%
  mutate(
    state_1 = ifelse(tempid > 5, paste(state_1, "(Low)"), state_1),
    state_2 = ifelse(tempid > 5, paste(state_2, "(High)"), state_2),
    state_1 = ifelse(tempid == 3, "0-9 (Low)", state_1),
    state_2 = ifelse(tempid == 3, ">9-35 (High)", state_2)
  ) %>%
  select(-tempid)

varcutoffs[18, 2] <- "Log NLE at TP 1"
varcutoffs[19, 2] <- "Log NLE at TP 2"
varcutoffs[20, 2] <- "Log NLE at TP 3"
varcutoffs[21, 2] <- "Log TLE at TP 1"
varcutoffs[22, 2] <- "Log TLE at TP 2"
varcutoffs[23, 2] <- "Log TLE at TP 3"


####


