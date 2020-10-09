#' Comparison of changes in disease activities between treatment groups.
#' PI: Dr Endoh
#' date created: 2020/09/28
#' ---

# subroutines-------

dir.sub <- './src/R/sub'
Bibtex <- FALSE

list.files.dir.sub <- list.files(path = dir.sub)

for(i in 1:length(list.files.dir.sub))
  source(sprintf("%s/%s", dir.sub, list.files.dir.sub[i]))

sink('sessionInfo.')
sessionInfo()
sink()

# Load data ---------------------------------------------------------------

load(file = sprintf("%s/%s", dir.ADS, fn.ADS))

colinfo$col_label <- 
  colinfo$col_label.2

data$tratment <- factor(data$treatment)


# models ------------
vars.0 <- data.frame(colinfo)[colinfo$prop_model==1,"col_names"]
vars.1 <- data.frame(colinfo)[colinfo$prop_model.0==1,"col_names"]
vars.2 <- data.frame(colinfo)[colinfo$prop_model.1==1,"col_names"]
vars.3 <- data.frame(colinfo)[colinfo$prop_model.2==1,"col_names"]
vars.4 <- data.frame(colinfo)[colinfo$prop_model.3==1,"col_names"]
vars.5 <- data.frame(colinfo)[colinfo$prop_model.4==1,"col_names"]
vars.6 <- data.frame(colinfo)[colinfo$prop_model.5==1,"col_names"]
vars.7 <- data.frame(colinfo)[colinfo$prop_model.6==1,"col_names"]
vars.8 <- data.frame(colinfo)[colinfo$prop_model.7==1,"col_names"]
vars.9 <- data.frame(colinfo)[colinfo$prop_model.8==1,"col_names"]
vars.10 <- data.frame(colinfo)[colinfo$prop_model.9==1,"col_names"]
# vars.11 <- data.frame(colinfo)[colinfo$prop_model.10==1,"col_names"]
# vars.12 <- data.frame(colinfo)[colinfo$prop_model.11==1,"col_names"]
# vars.13<- data.frame(colinfo)[colinfo$prop_model.12==1,"col_names"]
# vars.14 <- data.frame(colinfo)[colinfo$prop_model.13==1,"col_names"]
# vars.15 <- data.frame(colinfo)[colinfo$prop_model.14==1,"col_names"]
# vars.16 <- data.frame(colinfo)[colinfo$prop_model.15==1,"col_names"]
# vars.17 <- data.frame(colinfo)[colinfo$prop_model.16==1,"col_names"]
# vars.18 <- data.frame(colinfo)[colinfo$prop_model.17==1,"col_names"]
# vars.19 <- data.frame(colinfo)[colinfo$prop_model.18==1,"col_names"]
# vars.20 <- data.frame(colinfo)[colinfo$prop_model.19==1,"col_names"]
# vars.21 <- data.frame(colinfo)[colinfo$prop_model.20==1,"col_names"]
# vars.22 <- data.frame(colinfo)[colinfo$prop_model.21==1,"col_names"]


vars.smd.1 <- 
  data.frame(colinfo)[
    colinfo$smd==1,
    "col_names"
    ]

vars <- "vars.0"
vars.smd <- "vars.smd.1"

fml.ps_model <- sprintf(
  "treatment ~ %s", 
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
    strata = "treatment", 
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
dataPlot <- dataPlot[,c("col_label","rawdata","weighted_data")]

## Create long-format data for ggplot2
dataPlotMelt <-
  melt(
    data          = dataPlot,
    id.vars       = c("col_label"),
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
  mapping = aes(x = col_label, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + 
  theme(legend.key = element_blank())
dev.off()

# Make long data ----------------------------------------------------------

ADS <- data %>%
  pivot_longer(
    cols =ends_with("oc_init")
  ) %>%
  mutate(
    name=gsub("(.+)\\.oc_init","\\1",name)
  ) %>%
  dplyr::rename(
    value.pre=value
  ) %>%
  dplyr::left_join(
    data %>%
      pivot_longer(
        cols =ends_with("oc_post")
      ) %>%
      mutate(
        name=gsub("(.+)\\.oc_post","\\1",name)
      ) %>%
      dplyr::select(ID, name, value),
    by = c("ID","name") 
  ) %>%
  mutate(
    value.change = value.pre - value
  ) %>%
  left_join(
    data.propensityScores_IPW[,grep("^ID$|^w_.+",colnames(data.propensityScores_IPW))]
  )


## Weighted model (use svyglm function)

#' From help(svyglm):
#' Note
#' svyglm always returns 'model-robust' standard errors;
#'  the Horvitz-Thompson-type standard errors used 
#'  everywhere in the survey package are a generalisation
#'   of the model-robust 'sandwich' estimators. 
#'  In particular, a quasi-Poisson svyglm will return
#'   correct standard errors for relative risk regression models.

list.res.glm <-
  dlply(
    ADS,
    .(name),
    function(D){
      res.svydesign.weighted <- 
        survey::svydesign(
          ids = ~ 1, 
          data = D %>% dplyr::filter(!is.na(w_ato)),
          weights = ~ w_ato
        )
      res.glmWeighted <-
        svyglm(
          formula = value.change ~ treatment,
          # family  = quasibinomial(link = "logit"),                  
          design  = res.svydesign.weighted
        )
      }
    )


res.svydesign.weighted.binomial_outcome <- 
  survey::svydesign(
    ids = ~ 1, 
    data = data.propensityScores_IPW %>% dplyr::filter(!is.na(w_ato)),
    weights = ~ w_ato
  )
res.glmWeighted.binomial_outcome <-
  svyglm(
    formula = as.factor(outcome) ~ treatment,
    family  = quasibinomial(link = "logit"),                  
    design  = res.svydesign.weighted.binomial_outcome
    )

res.glmWeighted.ordinal_outcome <-
  svyolr(
    formula = factor(EULAR_response.oc_post) ~ treatment,
    design  = res.svydesign.weighted.binomial_outcome
  )

#' calculate p-value for regression coefficients of
#'  the ord-logistic regression.

summ.res.glmWeighted.ordinal_outcome <-
  tidy(
    res.glmWeighted.ordinal_outcome,
    conf.int = TRUE,
    conf.level = 0.95,
    exponentiate = TRUE,
    p.values = FALSE
  ) %>% data.frame()
ctable <- coef(summary(res.glmWeighted.ordinal_outcome))[, "t value"]
p <- pnorm(abs(ctable), lower.tail = FALSE) * 2
table.summ.res.glmWeighted.ordinal_outcome <- cbind(summ.res.glmWeighted.ordinal_outcome, "p value" = p)

sink(sprintf("%s/res.glm.txt", dir.output))
print("+++++++++++++++++")
print("$outcome")
print("+++++++++++++++++")
print(summary(res.glmWeighted.binomial_outcome))
print(exp(res.glmWeighted.binomial_outcome$coefficients))
print(exp(confint(res.glmWeighted.binomial_outcome)))

print("+++++++++++++++++")
print("$EULAR_response")
print("+++++++++++++++++")
print(table.summ.res.glmWeighted.ordinal_outcome)
lapply(list.res.glm, function(L){list(summary(L),confint(L))})
sink()


#' Investigation for the inconsistency 
#' in the estimates between 'Response' 
#' and other numeric outcomes.

quartz(type = 'pdf',
  file = sprintf('%s/%s', dir.output,'outcome.changeVal_vs_response.pdf')
  )
dlply(
  ADS %>% data.frame(),
  .(name),
  function(D){
    ylab <- unique(D$name)
    ggdata <- ggplot(data = D,
                     aes(x = as.factor(outcome), y=value.change)
    )
    plot(ggdata+geom_boxplot()+labs(y=ylab)+theme_bw())
  }
)
dev.off()