#' Create histogram
#' PI: Dr Arimori
#' date created: 2020/10/09
#' ---

# subroutines-------

dir.sub <- './src/R/sub'
Bibtex <- FALSE

list.files.dir.sub <- list.files(path = dir.sub)

for(i in 1:length(list.files.dir.sub))
  source(sprintf("%s/%s", dir.sub, list.files.dir.sub[i]))

sink('sessionInfo.txt')
sessionInfo()
sink()

# Load data ---------------------------------------------------------------

load(file = sprintf("%s/%s", dir.ADS, fn.ADS))

# models ------------

vars.1 <- data.frame(colinfo)[colinfo$prop_model.0==1,"col_names"]
vars.2 <- data.frame(colinfo)[colinfo$prop_model.1==1,"col_names"]
vars.3 <- data.frame(colinfo)[colinfo$prop_model.2==1,"col_names"]
vars.4 <- data.frame(colinfo)[colinfo$prop_model.3==1,"col_names"]

vars.smd.1 <- 
  data.frame(colinfo)[
    colinfo$smd==1,
    "col_names"
    ]

vars <- "vars.1"
vars.smd <- "vars.smd.1"

fml.ps_model <- sprintf(
  "%s ~ %s", 
  data.frame(colinfo)[colinfo$exposure==1,"col_names"],
  paste(eval(parse(text = vars)),collapse = "+")
  )


## Construct a table
tabUnmatched <-
  CreateTableOne(
    vars = eval(parse(text=vars.smd)), 
    strata = data.frame(colinfo)[colinfo$exposure==1,"col_names"], 
    data = data, 
    test = FALSE)
## Show table with SMD
sink(
  sprintf("%s/TableOne.txt", dir.output)
  )
print(tabUnmatched, smd = TRUE)
sink()

# Propensity score model ---------------

propensityScoreModel <-
  glm(
    gsub(
      data.frame(colinfo)[colinfo$exposure==1,"col_names"],
      sprintf("factor\\(%s\\)",data.frame(colinfo)[colinfo$exposure==1,"col_names"]),
      fml.ps_model
      ),
    family  = binomial(link = "logit"),
    data = data, na.action = na.exclude
    )

# propensityScoreModel <-
#   brglm::brglm(
#     gsub("treatment","factor\\(treatment\\)",fml.ps_model),
#     family = binomial(link="logit"),
#     data = data %>% data.frame(),
#     na.action = na.exclude
#   )

# propensityScoreModel <- 
#   rpart::rpart(
#     fml.ps_model,
#     data = data
#     )

data.propensityScores <-
  data %>% data.frame()
  
data.propensityScores$propensity_score <- 
  predict(
    propensityScoreModel,
    # data=data, 
    type= "response",
    # type= "prob",
    na.action = na.exclude()
    ) %>% unlist()

data.propensityScores_IPW <-
  IPW_weights(
    treatment = data.propensityScores[,data.frame(colinfo)[colinfo$exposure==1,"col_names"]],
    propensity_score = data.propensityScores$propensity_score,
    dat = data.propensityScores
    )

res.roc.propensity_score <- roc(
  response = 
    as.factor(data.propensityScores_IPW[,data.frame(colinfo)[colinfo$exposure==1,"col_names"]]),
  predictor = 
    data.propensityScores_IPW$propensity_score
  )

# Plot distribution of propensity score and weighted counts. ----------------
#'
ggdata.propensityScores <-
  ggplot(
    data =
      data.propensityScores_IPW,
    aes(
      x = propensity_score
      )
    )

ggdata.propensityScores.weighted_count <-
  ggplot(
    data =
      data.propensityScores_IPW %>%
      pivot_longer(
        cols =
          c(starts_with("w_at")),
        values_to = "weight",
        names_to = "target_pop"
      ) %>%
      dplyr::filter(target_pop=="w_ato"),
    aes(
      x = propensity_score,
      weight = weight
    )
  )

pdf(
  file = 
    sprintf("%s/IPWcount.%s.pdf", dir.output,vars),
  width = 21
  )
plot(
  ggdata.propensityScores + 
    geom_density(
      aes(
        fill = 
          as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"]))
        ),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5,
      position="identity"
      ) +
  geom_point(
    aes(
      y=as.numeric(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"])),
      x=propensity_score,
      color=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"])),
      size=2
      )
    ) +
    theme_bw()
  )
plot(
  ggdata.propensityScores.weighted_count + 
    geom_density(
      aes(fill=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"]))),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5#,
#      position="dodge"
      ) + 
    facet_grid(~target_pop) + 
    theme_bw()
    )
plot(
  res.roc.propensity_score,
  print.thres=TRUE 
  )
legend(
  x = 0.6, y=0.5,cex = 0.7, 
  # lwd = c(2,2,0), lty = 1:2,
  legend = c(
    sprintf(
      "AUC = %s (0.95CI: %s, %s)",
      round(auc(res.roc.propensity_score),3),
      round(ci(auc(res.roc.propensity_score))[1],3),
      round(ci(auc(res.roc.propensity_score))[3],3)
    )
  ),
  bty = "n"
)
dev.off()


# Balancing assessment ----------------------------------------------------

#' The standardized mean differences between the two patients' groups of *treatment* 
#' Reference:
#' https://cran.r-project.org/web/packages/tableone/vignettes/smd.html
#' [accessed:2020/08/31])

res.svydesign.weighted <- 
  survey::svydesign(
    ids = ~ 1, 
    data = data.propensityScores_IPW %>% dplyr::filter(!is.na(propensity_score)),
    weights = ~ w_ato
    )
  

## Construct a table
tabWeighted.weighted <- 
  svyCreateTableOne(
    vars = eval(parse(text = vars.smd)),
    strata = data.frame(colinfo)[colinfo$exposure==1,"col_names"], 
    data = res.svydesign.weighted, 
    test = FALSE
    )
## Show table with SMD

sink(
  sprintf("%s/TableOne._weighted.txt",dir.output)
  )
print(tabWeighted.weighted, smd = TRUE)
sink()


## Construct a data frame containing variable name and SMD from all methods
dataPlot <- data.frame(
  variable  = rownames(ExtractSmd(tabUnmatched)),
  rawdata = as.numeric(ExtractSmd(tabUnmatched)),
  weighted_data = as.numeric(ExtractSmd(tabWeighted.weighted))
  ) %>% 
  left_join(
    colinfo, 
    by = c("variable"="col_names")
  )
dataPlot <- dataPlot[,c("col_labels","rawdata","weighted_data")]

## Create long-format data for ggplot2
dataPlotMelt <-
  melt(
    data          = dataPlot,
    id.vars       = c("col_labels"),
    variable.name = "Method",
    value.name    = "SMD"
    )

## Order variable names by magnitude of SMD
varNames <- unique(
  as.character(dataPlot$col_label)[
    order(dataPlot$rawdata)
    ]
  )


## Order factor levels in the same order
dataPlotMelt$col_label <- 
  factor(
    dataPlotMelt$col_label,
    levels = varNames
    )

## Plot using ggplot2

quartz(
  family = "Arial",type = "pdf",
  file =   sprintf("%s/smd_%s.pdf",dir.output,vars)
  )
ggplot(
  data = dataPlotMelt,
  mapping = aes(x = col_labels, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + 
  theme(legend.key = element_blank())
dev.off()


# Make histogram ---------------------------------------------------------

ggdata.histogram.orig_data <- ggplot(
  data = data.propensityScores_IPW, 
  aes(x=eval(parse(text=data.frame(colinfo)[colinfo$var.y==1,"col_names"])))
  )

ggdata.histogram.weighted_data <- ggplot(
  data = data.propensityScores_IPW, 
  aes(
    x = eval(parse(
      text = 
        data.frame(colinfo)[colinfo$var.y==1,"col_names"])
      ),
    weight = w_ato
    )
  )

quartz(
  family = "Arial",type = "pdf",
  file =   sprintf("%s/histogram_%s.pdf",dir.output,vars)
  )
plot(
  ggdata.histogram.orig_data +
    geom_density(
      aes(fill=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"]))),
      bw="SJ",
      #      binwidth = FD,
      alpha=0.5#,
      #      position="dodge"
      ) +
    labs(
      fill=data.frame(colinfo)[colinfo$exposure==1,"col_names"]
      ) +
    theme_bw()
    )
plot(
  ggdata.histogram.weighted_data +
    geom_density(
      aes(fill=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"]))),
      bw="SJ",
      #      binwidth = FD,
      alpha=0.5#,
      #      position="dodge"
    ) +
    labs(
      fill=data.frame(colinfo)[colinfo$exposure==1,"col_names"]
    ) +
    theme_bw()
)
plot(
  ggdata.histogram.orig_data +
    stat_ecdf(
      aes(color=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"])))
      ) +
    labs(
      color=data.frame(colinfo)[colinfo$exposure==1,"col_names"]
    ) +
    theme_bw()
  )
plot(
  ggdata.histogram.weighted_data +
    stat_ecdf(
      aes(color=as.factor(get(data.frame(colinfo)[colinfo$exposure==1,"col_names"])))
      ) +
    labs(
      color=data.frame(colinfo)[colinfo$exposure==1,"col_names"]
    ) +
    theme_bw()
  )
dev.off()

# End.

#' NOTE: 2020/10/09
#'  The src for the Main analyses will refer the src 
#' in the mrmtshmp/sampleDesign_20.03 repository.