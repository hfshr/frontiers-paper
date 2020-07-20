pacman::p_load(tidyverse, knitr, kableExtra, qwraps2, brms, broom.mixed, bnlearn, tidyselect,
               Rgraphviz, sjmisc)

dinj <- read_csv("data/new_cleandata_inj.csv")

injtable <- dinj %>% 
  select(id, time, variable, gender, injured:typeofinjuryg) %>% 
  filter(time != 4 & injured == 1) 


inj_summary <- injtable %>% 
  mutate(bodylocation = ifelse(bodypart %in% c("upperarm_shoulder",
                                               "head_neck",
                                               "torso",
                                               "elbow_lowerarm_hand"), "Upper body", "Lower body")) %>% 
  group_by(gender, bodylocation) %>% 
  summarise("Joint / ligament" = n_perc0(typeofinjuryg == "joint_ligament"),
            "Muscle / tendon" = n_perc0(typeofinjuryg == "muscle_tendon"),
            "Other" = n_perc0(typeofinjuryg %in% c("other", "brain", "bone", "skin"))) %>% 
  gather(variable, value, -gender, -bodylocation)  %>% 
  spread(bodylocation, value) 


inj_mechanism <- injtable %>% 
  group_by(gender) %>% 
  summarise("Acute" = n_perc_x(acu_over == "acute", digits = 0),
            "Chronic" = n_perc_x(acu_over == "chronic", digits = 0),
            "Contact" = n_perc_x(contact_non == "contact", digits = 0),
            "Non-contact" = n_perc_x(contact_non == "noncontact", digits = 0)) %>% 
  gather(variable, value, -gender)%>% 
  spread(gender, value) %>% 
  ungroup()



tallc <- networkprepare(networkdata)
## extract names for blacklist + plotting

g1 <- tallc %>%
  select(ends_with("_1")) %>%
  names()
g2 <- tallc %>%
  select(ends_with("_2")) %>%
  names()

## get explanatory variable names
explans <- networkdata %>%
  select(pi:hours, nlebase) %>%
  names()

## create blacklist
bl <- expand.grid(from = g2, to = g1, stringsAsFactors = FALSE) %>%
  rbind(., expand.grid(from = c(g1, g2), to = explans, stringsAsFactors = FALSE)) %>%
  rbind(
    .,
    data.frame(
      from = c("clevel", "nlebase", "nlebase"),
      to = c("gender", "ind_team", "gender"),
      stringsAsFactors = FALSE
    )
  ) %>%
  rbind(
    .,
    data.frame(
      from = c(g1[2:length(g1)]),
      to = c(g2[2:length(g2)]),
      stringsAsFactors = FALSE
    )
  )
## create whitelist
wl <- data.frame(
  from = c("nlelg_1", "nlelg_2"),
  to = c("injured_1", "injured_2")
)

## generate model


library(parallel)

cl = makeCluster(8, type = "SOCK")
clusterSetRNGStream(cl, 30)

str <- boot.strength(
  tallc,
  algorithm = "tabu",
  algorithm.args = list(
    blacklist = bl,
    whitelist = wl,
    tabu = 1000
  ),
  cluster = cl
)


stopCluster(cl)

## get average network with 30% cut off value
avg30 <- averaged.network(str, threshold = 0.30)
fitted <- bn.fit(avg30, tallc)

## plot the network
gR <- strength.plot(
  avg30,
  str,
  shape = "rectangle",
  render = FALSE,
  groups = list(g1, g2, explans),
  layout = "dot"
)

nodeRenderInfo(gR)$fill[g1] = "deepskyblue"
nodeRenderInfo(gR)$fill[g2] = "tomato"
nodeRenderInfo(gR)$fill["injured_1"] = "gold"
nodeRenderInfo(gR)$fill["injured_2"] = "gold"
#gR <- renderGraph(gR)

# tiff("main_network.tif", res = 300, width = 15, height = 10, units = "in")
#renderGraph(gR)
# dev.off()

# mb(avg30, "injured_2")
## change network - see notes from first network

networkdatacont <- datanet %>%
  select(
    id,
    time,
    oneinj,
    pi,
    nlebase,
    gender,
    ind_team,
    clevel,
    hours,
    stiffness,
    nlec,
    RI,
    BIS,
    rmssd,
    balance,
    FFFS
  ) %>%
  mutate(nlebase = splitr(nlebase))

tallcc <- networkprepare(networkdatacont)

g1 <- tallcc %>%
  select(ends_with("_1")) %>%
  names()
g2 <- tallcc %>%
  select(ends_with("_2")) %>%
  names()

diffgtable <- data.frame(tallcc[, c("pi", "gender", "clevel", "hours")], tallcc[, g2] - tallcc[, g1]) %>%
  `names<-`(sub("_2", "", names(.))) %>%
  mutate(injured = tallcc$oneinj_1) %>%
  select(-oneinj) %>%
  select(stiffness:FFFS) %>%
  gather(varname, value) %>%
  group_by(varname) %>%
  summarise(
    Mean = mean(value),
    SD = sd(value)
  ) %>%
  dendroTools::round_df(2)



diffg <- data.frame(tallcc[, c("pi", "gender", "clevel", "hours")], tallcc[, g2] - tallcc[, g1]) %>%
  `names<-`(sub("_2", "", names(.))) %>%
  mutate(injured = tallcc$oneinj_1) %>%
  select(-oneinj) %>%
  std(stiffness:FFFS, append = T, suffix = "") %>%
  mutate(nlec = log(nlec + 1))

explans <- diffg %>%
  select(pi:hours) %>%
  names()

indivars <- diffg %>%
  select(stiffness:injured) %>%
  names()

bl <- expand.grid(
  from = setdiff(names(diffg), c("pi", "gender", "clevel")),
  to = c("pi", "gender", "clevel"),
  stringsAsFactors = FALSE
) %>%
  rbind(
    .,
    data.frame(
      from = c("gender", "clevel", "injured"),
      to = c("clevel", "gender", "hours"),
      stringsAsFactors = FALSE
    )
  )


cl = makeCluster(4, type = "SOCK")
clusterSetRNGStream(cl, 44)


strchange <- boot.strength(
  diffg,
  algorithm = "tabu",
  algorithm.args = list(blacklist = bl, 
                        tabu = 1000),
  cluster = cl
)

stopCluster(cl)

avgcont <- averaged.network(strchange, threshold = 0.30)


change <- strength.plot(
  avgcont,
  strchange,
  shape = "rectangle",
  render = F,
  groups = list(explans, indivars),
  layout = "dot"
)
nodeRenderInfo(change)$fill[indivars] = "tomato"
nodeRenderInfo(change)$fill["injured"] = "gold"
# 
# tiff("change_network.tif", res = 300, width = 15, height = 10, units = "in")
# renderGraph(change)
# dev.off()
fittedcont <- bn.fit(avgcont, diffg)


## Bayesian lm

testdataa <- cpdist(fittedcont, nodes = c("FFFS", "rmssd", "BIS"), evidence = TRUE, n = 1000) 

prior1 <- prior(normal(0, 5), class = "b")

bm1 <- brm(FFFS ~ BIS*rmssd, data = testdataa, prior = prior1)



comp <- newprobtable(
  tallc,
  vars = names(tallc[c(1:6, 8:14)]),
  outcome = "injured_1",
  state = "injured",
  avg30,
  repeats = 10000
)

injured1comp <- comp %>%
  top_n(1, prob) %>%
  as_tibble() %>%
  mutate(
    id = "highrisk",
    prob = round(prob, digits = 8)
  ) %>%
  rbind(
    .,
    comp %>%
      top_n(1, -prob) %>%
      as_tibble() %>%
      mutate(
        id = "lowrisk",
        prob = round(prob, digits = 8)
      )
  ) %>%
  mutate(prob = round(prob, digits = 2)) %>%
  gather(variable, value, -id) %>%
  spread(id, value) %>%
  mutate(variable = factor(
    variable,
    levels = c(
      "outcome",
      "prob",
      "pi",
      "clevel",
      "hours",
      "ind_team",
      "nlebase",
      "stiffness_1",
      "balance_1",
      "RI_1",
      "BIS_1",
      "FFFS_1",
      "nlelg_1",
      "rmssd_1",
      "gender"
    )
  )) %>%
  arrange(variable)

## MB 2

