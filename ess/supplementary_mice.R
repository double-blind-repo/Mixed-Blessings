# Running Parralel Analysis. 19 Factors.
n.factors <- fa.parallel(data.mice[wbvars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact
factorvars <- paste0("Factor",1:n.factors, collapse = " + ")

# Factor Names Table Builder
# Variable descriptions source: https://www.europeansocialsurvey.org/docs/round6/survey/ESS6_data_protocol_e01_4.pdf
make.factor.names.table <- function(model, no.general = FALSE) {
  start.factor <- ifelse(no.general==TRUE,1,2)
  variable.descriptions = read.csv("variable_descriptions.csv")
  names(variable.descriptions) <- c("variable","description")
  factor.loadings <- parameterestimates(model)
  
  factor.names.table <- data.frame(factor=character(),
                                   description=character(), 
                                   est=numeric()) 
  for (i in start.factor:n.factors) {
    table.temp <- factor.loadings[factor.loadings$lhs==paste0("F",i),]
    table.temp <- table.temp[table.temp$op=="=~",]
    table.temp <- table.temp[table.temp$est>=0.2,]
    table.temp <- table.temp[c("rhs","est")]
    table.temp <- merge(table.temp,variable.descriptions,by.x="rhs",by.y="variable")
    table.temp$factor <- paste0("Factor ", i)
    table.temp <- table.temp[c("factor","description","est")]
    table.temp <- table.temp[order(table.temp$est, decreasing = TRUE),]
    factor.names.table <- rbind(factor.names.table,table.temp)
  }
  return(factor.names.table)
}


# Orthogonal Bifactor Rotation ----
# ESEM
efa.mice.orthogonal <- fa(data.mice[wbvars], nfact = n.factors, rotate = "bifactor", fm = "minres", maxit=10000, weight = data.mice$newweight)
loadmat.mice.orthogonal <- zapsmall(matrix(round(efa.mice.orthogonal$loadings, 2), nrow = 47, ncol = n.factors))
rownames(loadmat.mice.orthogonal) <- colnames(data.mice[wbvars])
terms <- vector()
for (i in 1:n.factors) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.orthogonal[,i]), "*", names(loadmat.mice.orthogonal[,1]), collapse = "+"))
}
esem.formula.mice.orthogonal <- paste(terms, collapse = "\n")
esem.mice.orthogonal <- lavaan::cfa(esem.formula.mice.orthogonal, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing="listwise")
esem.fit.mice.orthogonal <- fitmeasures(esem.mice.orthogonal, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.orthogonal <- make.factor.names.table(esem.mice.orthogonal)
write.csv(table.mice.orthogonal,"supplementary/bifactor orthogonal/ess_imputed_mice_orthogonal.csv", row.names = FALSE)
# Factors are labelled qualitatively.
factor.names <- c("General","Engagement","Trust","Sadness","Social Closeness", "Social Optimism",
                  "Help and Support", "Positive Affect", "Social Engagement", "Evaluation", "Purpose",
                  "Self-Esteem", "Learning", "Resilience", "Effort", "Energy", "Calmness",
                  "Self-Efficacy", "Freedom")

for (i in n.factors:1) {
  table.mice.orthogonal$factor <- str_replace(table.mice.orthogonal$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(table.mice.orthogonal,"supplementary/bifactor orthogonal/ess_labelled_imputed_mice_orthogonal.csv", row.names = FALSE)

# Building Full Table
full.table.mice.orthogonal <- 1:n.factors
y <- parameterEstimates(esem.mice.orthogonal, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in wbvars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.factors)
  full.table.mice.orthogonal <- rbind(full.table.mice.orthogonal,line)
}
full.table.mice.orthogonal <- as.data.frame(full.table.mice.orthogonal)
names(full.table.mice.orthogonal) <- paste0("Factor",1:n.factors)
rownames(full.table.mice.orthogonal) <- c("Label", wbvars)
full.table.mice.orthogonal[1,] <- factor.names
write.csv(full.table.mice.orthogonal,"supplementary/bifactor orthogonal/ess_full_factor_loadings_mice_orthogonal.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice.orthogonal <- esem.mice.orthogonal %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice.orthogonal <- esem.mice.orthogonal %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice.orthogonal <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice.orthogonal, variance.table.mice.orthogonal) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice.orthogonal,"supplementary/bifactor orthogonal/ess_mean_variance_mice_orthogonal.csv", row.names = FALSE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.factors,"~ hinctnta + ",personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.orthogonal, x)
mimic.1.mice.orthogonal <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.1.fit.mice.orthogonal <- fitmeasures(mimic.1.mice.orthogonal, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.orthogonal)
all_models <- all_models %>% filter(grepl(" ~ hinctnta", term))
all_models$model <- paste0("Factor",1:n.factors)
all_models$term <- "hinctnta"
all_models <- arrange(all_models, estimate)
for (i in n.factors:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.orthogonal <- all_models
mimic.1.graph.mice.orthogonal <- dwplot(all_models)
mimic.1.graph.mice.orthogonal <- mimic.1.graph.mice.orthogonal %>% relabel_predictors(hinctnta = "Income Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Bifactor Orthogonal") + xlim(-0.075,0.075)
ggsave(file="supplementary/bifactor orthogonal/ess_income_decile_controlled_imputed_mice_orthogonal.pdf", mimic.1.graph.mice.orthogonal, width = 276, height = 145, units = "mm", device = pdf())
write.csv(mimic.results.table.mice.orthogonal,"supplementary/bifactor orthogonal/ess_mimic_income_decile_controlled_imputed_mice_orthogonal.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.orthogonal, x)
mimic.2.mice.orthogonal <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice.orthogonal <- fitmeasures(mimic.2.mice.orthogonal, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice.orthogonal)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.mice.orthogonal <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=7, ncol=3)
ggsave(file="supplementary/bifactor orthogonal/ess_income_controls_categorical_controlled_imputed_mice_orthogonal.pdf", mimic.2.graph.mice.orthogonal, width = 210, height = 297, units = "mm", device = pdf())


# Oblimin  Rotation ----
# ESEM
efa.mice.oblimin <- fa(data.mice[wbvars], nfact = n.factors, rotate = "oblimin", fm = "minres", maxit=10000, weight = data.mice$newweight)
loadmat.mice.oblimin <- zapsmall(matrix(round(efa.mice.oblimin$loadings, 2), nrow = 47, ncol = n.factors))
rownames(loadmat.mice.oblimin) <- colnames(data.mice[wbvars])
terms <- vector()
for (i in 1:n.factors) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.oblimin[,i]), "*", names(loadmat.mice.oblimin[,1]), collapse = "+"))
}
esem.formula.mice.oblimin <- paste(terms, collapse = "\n")
esem.mice.oblimin <- lavaan::cfa(esem.formula.mice.oblimin, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing="listwise")
esem.fit.mice.oblimin <- fitmeasures(esem.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.oblimin <- make.factor.names.table(esem.mice.oblimin, no.general = TRUE)
write.csv(table.mice.oblimin,"supplementary/oblimin/ess_imputed_mice_oblimin.csv", row.names = FALSE)
# Factors are labelled qualitatively.
factor.names <- c("Engagement", "Trust", "Evaluation", "Social Closeness",
                  "Sadness", "Self-Esteem","Help and Support", "Positive Affect",
                  "Effort", "Social Optimism", "Social Engagement", "Purpose", "Appreciated",
                  "Resilience", "Learning", "Energy", "Freedom", "Calmness", "Self-Efficacy")
for (i in n.factors:1) {
  table.mice.oblimin$factor <- str_replace(table.mice.oblimin$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(table.mice.oblimin,"supplementary/oblimin/ess_labelled_imputed_mice_oblimin.csv", row.names = FALSE)

# Building Full Table
full.table.mice.oblimin <- 1:n.factors
y <- parameterEstimates(esem.mice.oblimin, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in wbvars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.factors)
  full.table.mice.oblimin <- rbind(full.table.mice.oblimin,line)
}
full.table.mice.oblimin <- as.data.frame(full.table.mice.oblimin)
names(full.table.mice.oblimin) <- paste0("Factor",1:n.factors)
rownames(full.table.mice.oblimin) <- c("Label", wbvars)
full.table.mice.oblimin[1,] <- factor.names
write.csv(full.table.mice.oblimin,"supplementary/oblimin/ess_full_factor_loadings_mice_oblimin.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice.oblimin <- esem.mice.oblimin %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice.oblimin <- esem.mice.oblimin %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice.oblimin <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice.oblimin, variance.table.mice.oblimin) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice.oblimin,"supplementary/oblimin/psid_mean_variance_mice_oblimin.csv", row.names = FALSE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.factors,"~ hinctnta + ",personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.oblimin, x)
mimic.1.mice.oblimin <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.1.fit.mice.oblimin <- fitmeasures(mimic.1.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.oblimin)
all_models <- all_models %>% filter(grepl(" ~ hinctnta", term))
all_models$model <- paste0("Factor",1:n.factors)
all_models$term <- "hinctnta"
all_models <- arrange(all_models, estimate)
for (i in n.factors:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.oblimin <- all_models
mimic.1.graph.mice.oblimin <- dwplot(all_models)
mimic.1.graph.mice.oblimin <- mimic.1.graph.mice.oblimin %>% relabel_predictors(hinctnta = "Income Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Oblimin (Oblique)") + xlim(-0.075,0.075)
ggsave(file="supplementary/oblimin/ess_income_decile_controlled_imputed_mice.pdf", mimic.1.graph.mice.oblimin, width = 276, height = 145, units = "mm", device = pdf())
write.csv(mimic.results.table.mice.oblimin,"supplementary/oblimin/ess_mimic_income_decile_controlled_imputed_mice_oblimin.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.oblimin, x)
mimic.2.mice.oblimin <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice.oblimin <- fitmeasures(mimic.2.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice.oblimin)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.mice.oblimin <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=7, ncol=3)
ggsave(file="supplementary/oblimin/ess_income_controls_categorical_controlled_imputed_mice_oblimin.pdf", mimic.2.graph.mice.oblimin, width = 210, height = 297, units = "mm", device = pdf())

#Income Categories Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.oblimin, x)
mimic.2.mice.oblimin <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice.oblimin <- fitmeasures(mimic.2.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice.oblimin)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.mice.oblimin <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=5, ncol=4)
ggsave(file="supplementary/oblimin/ess_income_controls_categorical_controlled_imputed_mice_oblimin.pdf", mimic.2.graph.mice.oblimin, width = 210, height = 297, units = "mm", device = pdf())

# Varimax Rotation ----
# ESEM
efa.mice.varimax <- fa(data.mice[wbvars], nfact = n.factors, rotate = "varimax", fm = "minres", maxit=10000, weight = data.mice$newweight)
loadmat.mice.varimax <- zapsmall(matrix(round(efa.mice.varimax$loadings, 2), nrow = 47, ncol = n.factors))
rownames(loadmat.mice.varimax) <- colnames(data.mice[wbvars])
terms <- vector()
for (i in 1:n.factors) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.varimax[,i]), "*", names(loadmat.mice.varimax[,1]), collapse = "+"))
}
esem.formula.mice.varimax <- paste(terms, collapse = "\n")
esem.mice.varimax <- lavaan::cfa(esem.formula.mice.varimax, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing="listwise")
esem.fit.mice.varimax <- fitmeasures(esem.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.varimax <- make.factor.names.table(esem.mice.varimax,no.general = TRUE)
write.csv(table.mice.varimax,"supplementary/varimax/ess_imputed_mice_varimax.csv", row.names = FALSE)
# Factors are labelled qualitatively.
factor.names <- c("Sadness", "Engagement", "Purpose", "Trust", "Positive Affect", 
                  "Help and Support", "Social Closeness", "Evaluation", "Social Engagement",
                  "Social Optimism", "Self Esteem", "Resilience", "Appreciated",
                  "Learning", "Energy", "Calmness", "Selflessness", "Self-Efficacy", "Sadness")

for (i in n.factors:1) {
  table.mice.varimax$factor <- str_replace(table.mice.varimax$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(table.mice.varimax,"supplementary/varimax/ess_labelled_imputed_mice_varimax.csv", row.names = FALSE)

# Building Full Table
full.table.mice.varimax <- 1:n.factors
y <- parameterEstimates(esem.mice.varimax, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in wbvars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.factors)
  full.table.mice.varimax <- rbind(full.table.mice.varimax,line)
}
full.table.mice.varimax <- as.data.frame(full.table.mice.varimax)
names(full.table.mice.varimax) <- paste0("Factor",1:n.factors)
rownames(full.table.mice.varimax) <- c("Label", wbvars)
full.table.mice.varimax[1,] <- factor.names
write.csv(full.table.mice.varimax,"supplementary/varimax/ess_full_factor_loadings_mice_varimax.csv", row.names = TRUE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.factors,"~ hinctnta + ",personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.varimax, x)
mimic.1.mice.varimax <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.1.fit.mice.varimax <- fitmeasures(mimic.1.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.varimax)
all_models <- all_models %>% filter(grepl(" ~ hinctnta", term))
all_models$model <- paste0("Factor",1:n.factors)
all_models$term <- "hinctnta"
all_models <- arrange(all_models, estimate)
for (i in n.factors:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.varimax <- all_models
mimic.1.graph.mice.varimax <- dwplot(all_models)
mimic.1.graph.mice.varimax <- mimic.1.graph.mice.varimax %>% relabel_predictors(hinctnta = "Income Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Varimax (Orthogonal)") + xlim(-0.075,0.075)
ggsave(file="supplementary/varimax/ess_income_decile_controlled_imputed_mice.pdf", mimic.1.graph.mice.varimax, width = 276, height = 145, units = "mm", device = pdf())
write.csv(mimic.results.table.mice.varimax,"supplementary/varimax/ess_mimic_income_decile_controlled_imputed_mice_varimax.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.varimax, x)
mimic.2.mice.varimax <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice.varimax <- fitmeasures(mimic.2.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice.varimax)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.mice.varimax <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=5, ncol=4)
ggsave(file="supplementary/varimax/ess_income_controls_categorical_controlled_imputed_mice_varimax.pdf", mimic.2.graph.mice.varimax, width = 210, height = 297, units = "mm", device = pdf())

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice.varimax, x)
mimic.2.mice.varimax <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice.varimax <- fitmeasures(mimic.2.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice.varimax)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.mice.varimax <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=5, ncol=4)
ggsave(file="supplementary/varimax/ess_income_controls_categorical_controlled_imputed_mice_varimax.pdf", mimic.2.graph.mice.varimax, width = 210, height = 297, units = "mm", device = pdf())

save.image("supplementary_session_complete.RData")