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
set.seed(44)
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
gR <- renderGraph(gR)


jpeg("./figures_doc/Fig1.jpeg", res = 300, width = 20, height = 15, units = "cm")
renderGraph(gR)
dev.off()
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
    hrv,
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

set.seed(55)
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

jpeg("./figures_doc/Fig4.jpeg", res = 300, width = 20, height = 15, units = "cm")
renderGraph(change)
dev.off()
fittedcont <- bn.fit(avgcont, diffg)


## Bayesian lm

testdataa <- cpdist(fittedcont, nodes = c("FFFS", "hrv", "BIS"), evidence = TRUE, n = 1000) 

prior1 <- prior(normal(0, 5), class = "b")
set.seed(66)
bm1 <- brm(FFFS ~ BIS*hrv, data = testdataa, prior = prior1)



# comp <- newprobtable(
#   tallc,
#   vars = names(tallc[c(1:6, 8:14)]),
#   outcome = "injured_1",
#   state = "injured",
#   avg30,
#   repeats = 10000
# )
# 
# injured1comp <- comp %>%
#   top_n(1, prob) %>%
#   as_tibble() %>%
#   mutate(
#     id = "highrisk",
#     prob = round(prob, digits = 8)
#   ) %>%
#   rbind(
#     .,
#     comp %>%
#       top_n(1, -prob) %>%
#       as_tibble() %>%
#       mutate(
#         id = "lowrisk",
#         prob = round(prob, digits = 8)
#       )
#   ) %>%
#   mutate(prob = round(prob, digits = 2)) %>%
#   gather(variable, value, -id) %>%
#   spread(id, value) %>%
#   mutate(variable = factor(
#     variable,
#     levels = c(
#       "outcome",
#       "prob",
#       "pi",
#       "clevel",
#       "hours",
#       "ind_team",
#       "nlebase",
#       "stiffness_1",
#       "balance_1",
#       "RI_1",
#       "BIS_1",
#       "FFFS_1",
#       "nlelg_1",
#       "hrv_1",
#       "gender"
#     )
#   )) %>%
#   arrange(variable)

## MB 2

plotsg1 <- model2network('[hours][stiffness_1][nlelg_1][clevel][injured_1|hours:nlelg_1:stiffness_1][balance_1|clevel:injured_1]')
q1p <- graphviz.plot(plotsg1,
                     shape = "rectangle",
                     render = F,
                     groups = list(c("clevel", "hours"), "injured_1", c("stiffness_1", "nlelg_1", "balance_1")))

strq1 <- str %>%
  filter(from %in% c("clevel", "hours", "nlelg_1", "stiffness_1", "injured_1") &
           to %in% c("balance_1", "injured_1")) %>%
  filter(strength > 0.30) %>%
  mutate(strength = format(round(strength, digits = 2), nsmall = 2))

labels <- as.character(strq1$strength)
labels[4] <- str_pad(labels[4], side = "left", width = 2)

edgenames <- edgeNames(q1p)
names(labels) <- edgenames
nodeRenderInfo(q1p)$fill[c("injured_1", "nlelg_1")] = "gold"
edgeRenderInfo(q1p) <- list(lwd = c("nlelg_1~injured_1" = 6,
                                    "stiffness_1~injured_1" = 3,
                                    "clevel~balance_1" = 3),
                            lty = c("injured_1~balance_1" = "dashed",
                                    "hours~injured_1" = "dashed"))

q1p <- layoutGraph(q1p, edgeAttrs=list(label = labels))
graph.par(list(edges=list(fontsize = 15)))

nodeRenderInfo(q1p)$fill[c("nlelg_1", "stiffness_1", "balance_1")] = "deepskyblue"
nodeRenderInfo(q1p)$fill["injured_1"] = "gold"

jpeg("figures_doc/Fig2.jpeg", res = 300, width = 20, height = 15, units = "cm")
renderGraph(q1p)
dev.off()


plotsg2 <- model2network('[nlelg_2][gender][FFFS_1][injured_2|FFFS_1:nlelg_2][stiffness_2|gender:injured_2][balance_2|injured_2][hrv_2|injured_2]')
q2p <- graphviz.plot(plotsg2,
                     render = F,
                     shape = "rectangle",
                     groups = list(c("gender"), "injured_2", c("stiffness_2", "nlelg_2", "balance_2", "hrv_2")))

strq2 <- str %>%
  filter(from %in% c("FFFS_1", "nlelg_2", "stiffness_1", "injured_2", "gender") &
           to %in% c("balance_2", "injured_2", "nlelg_2", "stiffness_2", "hrv_2")) %>%
  filter(strength > 0.29) %>%
  filter(direction != 0) %>%
  mutate(strength = format(round(strength, digits = 2), nsmall = 2))

labels <- as.character(strq2$strength)
labels[6] <- str_pad(labels[6], side = "left", width = 2)
edgenames <- c("gender~stiffness_2", 
               "FFFS_1~injured_2",
               "injured_2~stiffness_2", 
               "injured_2~balance_2", 
               "injured_2~hrv_2", 
               "nlelg_2~injured_2")
names(labels) <- edgenames
nodeRenderInfo(q2p)$fill[c("injured_2", "nlelg_2")] = "gold"
edgeRenderInfo(q2p) <- list(lwd = c("injured_2~stiffness_2" = 3,
                                    "nlelg_2~injured_2" = 6,
                                    "clevel~balance_1" = 3,
                                    "gender~stiffness_2" = 4),
                            lty = c("FFFS_1~injured_2" = "dashed",
                                    "injured_2~balance_2" = "dashed",
                                    "injured_2~hrv_2" = "dashed"))

q2p <- layoutGraph(q2p, edgeAttrs=list(label = labels))
graph.par(list(edges=list(fontsize = 15)))
nodeRenderInfo(q2p)$fill[c("nlelg_2", "stiffness_2", "balance_2", "FFFS_1", "hrv_2")] = "tomato"
nodeRenderInfo(q2p)$fill["injured_2"] = "gold"
jpeg("figures_doc/Fig3.jpeg", res = 300, width = 20, height = 15, units = "cm")
renderGraph(q2p)
dev.off()


plotsg3 <- model2network('[gender][pi][hours][injured|pi:hours][stiffness|gender:injured][nlec|injured]')
q3p <- graphviz.plot(plotsg3,
                     render = F,
                     shape = "rectangle",
                     groups = list(c("gender", "pi", "hours"), "injured", c("stiffness", "nlec")))
strq3 <- strchange %>%
  filter(from %in% c("pi", "hours", "injured", "gender") &
           to %in% c("stiffness", "nlec", "injured")) %>%
  filter(strength > 0.30) %>%
  mutate(strength = format(round(strength, digits = 2), nsmall = 2))

labels <- as.character(strq3$strength)
edgenames <- c("pi~injured", "gender~stiffness", "hours~injured",
               "injured~stiffness", "injured~nlec")
names(labels) <- edgenames
nodeRenderInfo(q3p)$fill[c("injured")] = "gold"
edgeRenderInfo(q3p) <- list(lwd = c("injured~stiffness" = 4,
                                    "hours~injured" = 5,
                                    "gender~stiffness" = 4),
                            lty = c("pi~injured" = "dashed"))

q3p <- layoutGraph(q3p, edgeAttrs=list(label = labels))
graph.par(list(edges=list(fontsize = 15)))
nodeRenderInfo(q3p)$fill[c("nlec", "stiffness")] = "tomato"
nodeRenderInfo(q3p)$fill["injured"] = "gold"
jpeg("./figures_doc/Fig5.jpeg", res = 300, width = 20, height = 15, units = "cm")
renderGraph(q3p)
dev.off()
rm(plotsg3, strq3, labels)

# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(rsvg)
# grViz(
#   "
#  digraph {
# 
# # graph attributes
# graph [rankdir = LR, 
# overlap = true, 
# fontsize = 12]
# 
# # node attributes
# node [shape = rectangle]
# 
# # edge attributes
# edge [color = black,
# len = 0.5]
# 
# # node statements
# 
# T1 [label = <T1 <br/> <i>n</i>  = 351 <br/> All measures <br/> (no injury reporting)>]
# T2 [label = <T2 <br/> <i>n</i>  = 236 <br/> All measures>]
# T3 [label = <T3 <br/> <i>n</i>  = 157 <br/> All measures>]
# T4 [label = <T4 <br/> <i>n</i>  = 163 <br/> Only injury reporting>]
# 
# d1 [label = 'Sep \n 2016 and 2017', shape = plaintext]
# d2 [label = 'Jan \n 2017 and 2018', shape = plaintext]
# d3 [label = 'May \n 2017 and 2018', shape = plaintext]
# d4 [label = 'Sep \n 2017 and 2018', shape = plaintext]
# 
# d1 -> d2 -> d3 -> d4
# 
# T1 -> T2 -> T3 -> T4 
# 
# }
# "
# ) %>%
#   export_svg %>% charToRaw %>% rsvg_png("./figures_doc/pic.png")
# 
# 
# 
# grViz(" 
#       digraph CFA {
#       a [label = 'Node 1', shape = rectangle, color=CornflowerBlue ]; 
# 
#       node [shape = ellipse, color=CornflowerBlue]
#       T1    [label = <Node 2 <br/> <u>extra detail</u>>]; 
#       T2    [label = 'Node 3']; 
# 
#       {rank = same; a T1 T2}
# 
#       # Connect nodes with edges and labels
#       a -> T1
#       T2 -> a[dir=back]
#        }
# 
#       ")
# 
# sessionprotocol <- grViz("digraph {
# 
# # graph attributes
# graph [
# rankdir = TB, 
# overlap = true, 
# fontsize = 14]
# 
# node [shape = rectangle, 
# width = 2, 
# fixedsize = true]
# 
# start [label = 'Participants separated into groups \n Informed consent (T1 only)',
# fixedsize = false]
# pm [label = 'Computer-based \n measures', shape = oval]
# phys [label = 'Physical measures', shape = oval]
# les [label = 'LESCA \n (15 minutes)']
# rst [label = 'RST-PQ \n (10 minutes)']
# in [label = 'Injury reporting \n (5 minutes)']
# hrv [label = 'HRV \n (10 minutes)']
# myo [label = 'Myoton \n (5 minutes)']
# bal [label = 'Postural stability \n (5 minutes)']
# endpm [label = 'Start computer based \n measures']
# endphys [label = 'Start physical \n measures']
# 
# start -> pm [label = 'Group 1']
# 
# pm -> les -> rst -> in -> endphys
# 
# start -> phys [label = 'Group 2']
# 
# phys -> hrv -> myo -> bal -> endpm
# 
# }")
# sessionprotocol