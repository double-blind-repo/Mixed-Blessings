# 15 factors.
n.parallel <- fa.parallel(psid15.wb16.data.mice[total.wb.vars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact

# Factor Names Table Function
make.factor.names.table <- function(model, no.general = FALSE) {
  start.factor <- ifelse(no.general==TRUE,1,2)
  variable.descriptions = read.csv("variable_descriptions.csv")
  names(variable.descriptions) <- c("variable","description")
  factor.loadings <- parameterestimates(model)
  
  factor.names.table <- data.frame(factor=character(),
                                   description=character(), 
                                   est=numeric()) 
  for (i in start.factor:n.parallel) {
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
efa.mice.orthogonal <- fa(psid15.wb16.data.mice[total.wb.vars], nfact = n.parallel, rotate = "bifactor", fm = "minres", maxit=10000, weight=psid15.wb16.data.mice$WB16WT)
loadmat.mice.orthogonal <- zapsmall(matrix(round(efa.mice.orthogonal$loadings, 2), nrow = 41, ncol = n.parallel))
rownames(loadmat.mice.orthogonal) = colnames(psid15.wb16.data.mice[total.wb.vars])
terms <- vector()
for (i in 1:n.parallel) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.orthogonal[,i]), "*", names(loadmat.mice.orthogonal[,1]), collapse = "+"))
}
esem.formula.mice.orthogonal <- paste(terms, collapse = "\n")
esem.mice.orthogonal <- lavaan::cfa(esem.formula.mice.orthogonal, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = total.wb.vars, missing="listwise")
esem.fit.mice.orthogonal <- fitmeasures(esem.mice.orthogonal, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.orthogonal <- make.factor.names.table(esem.mice.orthogonal)

# Writing .csv so factors can be labelled.
write.csv(table.mice,"supplementary/bifactor orthogonal/psid_imputed_mice_orthogonal.csv", row.names = FALSE)

# Factors are labelled qualitatively.
factor.names <- c("General","Functioning","Satisfaction","Negative Affect (Month)", "Positive Affect (Month)",
                  "Positive Affect (Day)", "Nervousness","Loneliness",
                  "Pain/Tiredness", "Calmness (Month)", "Respect", "Anxiety / Stress (Day)",
                  "Social Relationships", "Full of Life", "Important Things in Life")

for (i in n.parallel:1) {
  table.mice.orthogonal$factor <- str_replace(table.mice.orthogonal$factor,paste0("Factor ",i),factor.names[i])
}

# Writing .csv with labelled factors.
write.csv(table.mice.orthogonal,"supplementary/bifactor orthogonal/psid_labelled_imputed_mice_orthogonal.csv", row.names = FALSE)

# Building Full Table
full.table.mice.orthogonal <- 1:n.parallel
y <- parameterEstimates(esem.mice.orthogonal, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in total.wb.vars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.parallel)
  full.table.mice.orthogonal <- rbind(full.table.mice.orthogonal,line)
}
full.table.mice.orthogonal <- as.data.frame(full.table.mice.orthogonal)
names(full.table.mice.orthogonal) <- paste0("Factor",1:n.parallel)
rownames(full.table.mice.orthogonal) <- c("Label", total.wb.vars)
full.table.mice.orthogonal[1,] <- factor.names
write.csv(full.table.mice.orthogonal,"supplementary/bifactor orthogonal/psid_full_factor_loadings_mice_orthogonal.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice.orthogonal <- esem.mice.orthogonal %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice.orthogonal <- esem.mice.orthogonal %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice.orthogonal <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice.orthogonal, variance.table.mice.orthogonal) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice.orthogonal,"supplementary/bifactor orthogonal/psid_mean_variance_mice_orthogonal.csv", row.names = FALSE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.parallel,"~", paste("income.decile", collapse = " + "), " + ", "female + age + arechildren +", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice.orthogonal, x)
mimic.1.mice.orthogonal <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = total.wb.vars, missing = "listwise")
mimic.1.fit.mice.orthogonal <- fitmeasures(mimic.1.mice.orthogonal, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.orthogonal)
all_models <- all_models %>% filter(grepl(paste0(" ~ ", "income.decile"), term))
all_models$model <- paste0("Factor",1:n.parallel)
all_models$term <- "income.decile"
all_models <- arrange(all_models, estimate)
for (i in n.parallel:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.orthogonal <- all_models
mimic.1.graph.mice.orthogonal <- dwplot(all_models)
mimic.1.graph.mice.orthogonal <- mimic.1.graph.mice.orthogonal %>% relabel_predictors(income.total = "Income") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Bifactor Orthogonal")
ggsave(file="supplementary/bifactor orthogonal/psid_income_controlled_imputed_mice_orthogonal.pdf", mimic.1.graph.mice.orthogonal, width = 276, height = 145, units = "mm")
write.csv(mimic.results.table.mice.orthogonal,"supplementary/bifactor orthogonal/psid_mimic_income_decile_controlled_imputed_mice_orthogonal.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.parallel,"~", paste(decile.vars, collapse = " + "), " + ", "female + age + arechildren + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice.orthogonal, x)
mimic.2.mice.orthogonal <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = total.wb.vars, missing="listwise")
all_models <- tidy(mimic.2.mice.orthogonal)
inc_models <- all_models %>% filter(grepl(" ~ income.decile_", term))
graphs = list()
for (i in 1:n.parallel){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("income.decile_",2:10)
  finalgraphs <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- finalgraphs %>% relabel_predictors(income.decile_2 = "2nd Decile", income.decile_3 = "3rd Decile", income.decile_4 = "4th Decile", income.decile_5 = "5th Decile", income.decile_6 = "6th Decile", income.decile_7 = "7th Decile", income.decile_8 = "8th Decile", income.decile_9 = "9th Decile", income.decile_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-1, 1.2)
}
mimic.2.graph.mice.orthogonal <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], nrow=5, ncol=3)
ggsave(file="supplementary/bifactor orthogonal/psid_income_categories_controlled_oblique_mice_orthogonal.pdf", mimic.2.graph.mice.orthogonal, width = 210, height = 297, units = "mm")


# Oblimin Rotation ----
efa.mice.oblimin <- fa(psid15.wb16.data.mice[total.wb.vars], nfact = n.parallel, rotate = "oblimin", fm = "minres", maxit=10000, weight=psid15.wb16.data.mice$WB16WT)
loadmat.mice.oblimin <- zapsmall(matrix(round(efa.mice.oblimin$loadings, 2), nrow = 41, ncol = n.parallel))
rownames(loadmat.mice.oblimin) = colnames(psid15.wb16.data.mice[total.wb.vars])
terms <- vector()
for (i in 1:n.parallel) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.oblimin[,i]), "*", names(loadmat.mice.oblimin[,1]), collapse = "+"))
}
esem.formula.mice.oblimin <- paste(terms, collapse = "\n")
esem.mice.oblimin <- lavaan::cfa(esem.formula.mice.oblimin, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
esem.fit.mice.oblimin  <- fitmeasures(esem.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.oblimin  <- make.factor.names.table(esem.mice.oblimin, no.general = TRUE)

# Writing .csv so factors can be labelled.
write.csv(table.mice.oblimin,"supplementary/oblimin/psid_imputed_mice_oblimin.csv", row.names = FALSE)

# Factors are labelled qualitatively.
factor.names <- c("Positive Affect (Day)", "Negative Affect (Month)","Satisfaction","Positive Affect (Month)","Frustration",
                  "Functioning","Nervousness","Evaluation","Respect","Calmness (Month)",
                  "Pain/Tiredness","Loneliness","Anxiety / Stress (Day)","Capability","Full of Life")

for (i in n.parallel:1) {
  table.mice.oblimin$factor <- str_replace(table.mice.oblimin$factor,paste0("Factor ",i),factor.names[i])
}

# Writing .csv with labelled factors.
write.csv(table.mice.oblimin,"supplementary/oblimin/psid_labelled_imputed_mice_oblimin.csv", row.names = FALSE)

# Building Full Table
full.table.mice.oblimin <- 1:n.parallel
y <- parameterEstimates(esem.mice.oblimin, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in total.wb.vars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.parallel)
  full.table.mice.oblimin <- rbind(full.table.mice.oblimin,line)
}
full.table.mice.oblimin <- as.data.frame(full.table.mice.oblimin)
names(full.table.mice.oblimin) <- paste0("Factor",1:n.parallel)
rownames(full.table.mice.oblimin) <- c("Label", total.wb.vars)
full.table.mice.oblimin[1,] <- factor.names
write.csv(full.table.mice.oblimin,"supplementary/oblimin/psid_full_factor_loadings_mice_oblimin.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice.oblimin <- esem.mice.oblimin %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice.oblimin <- esem.mice.oblimin %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice.oblimin <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice.oblimin, variance.table.mice.oblimin) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice.oblimin,"supplementary/oblimin/psid_mean_variance_mice_oblimin.csv", row.names = FALSE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.parallel,"~", paste("income.decile", collapse = " + "), " + ", "female + age + arechildren +", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice.oblimin, x)
mimic.1.mice.oblimin <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing = "listwise")
mimic.1.fit.mice.oblimin <- fitmeasures(mimic.1.mice.oblimin, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.oblimin)
all_models <- all_models %>% filter(grepl(paste0(" ~ ", "income.decile"), term))
all_models$model <- paste0("Factor",1:n.parallel)
all_models$term <- "income.decile"
all_models <- arrange(all_models, estimate)
for (i in n.parallel:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.oblimin <- all_models
mimic.1.graph.mice.oblimin <- dwplot(all_models)
mimic.1.graph.mice.oblimin <- mimic.1.graph.mice.oblimin %>% relabel_predictors(income.total = "Income") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Oblimin (Oblique)")
ggsave(file="supplementary/oblimin/psid_income_controlled__imputed_mice_oblimin.pdf", mimic.1.graph.mice.oblimin, width = 276, height = 145, units = "mm")
write.csv(mimic.results.table.mice.oblimin,"supplementary/oblimin/psid_mimic_income_decile_controlled_imputed_mice_oblimin.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.parallel,"~", paste(decile.vars, collapse = " + "), " + ", "female + age + arechildren + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice.oblimin, x)
mimic.2.mice.oblimin <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
all_models <- tidy(mimic.2.mice.oblimin)
inc_models <- all_models %>% filter(grepl(" ~ income.decile_", term))
graphs = list()
for (i in 1:n.parallel){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("income.decile_",2:10)
  finalgraphs <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- finalgraphs %>% relabel_predictors(income.decile_2 = "2nd Decile", income.decile_3 = "3rd Decile", income.decile_4 = "4th Decile", income.decile_5 = "5th Decile", income.decile_6 = "6th Decile", income.decile_7 = "7th Decile", income.decile_8 = "8th Decile", income.decile_9 = "9th Decile", income.decile_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-1, 1.2)
}
mimic.2.graph.mice.oblimin <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], nrow=5, ncol=3)
ggsave(file="supplementary/oblimin/psid_income_categories_controlled_oblique_mice_oblimin.pdf", mimic.2.graph.mice.oblimin, width = 210, height = 297, units = "mm")


# Varimax Rotation ----
# Factor Names Table Function
make.factor.names.table <- function(model, no.general = FALSE) {
  start.factor <- ifelse(no.general==TRUE,1,2)
  variable.descriptions = read.csv("variable_descriptions.csv")
  names(variable.descriptions) <- c("variable","description")
  factor.loadings <- parameterestimates(model)
  
  factor.names.table <- data.frame(factor=character(),
                                   description=character(), 
                                   est=numeric()) 
  for (i in start.factor:n.parallel) {
    table.temp <- factor.loadings[factor.loadings$lhs==paste0("F",i),]
    table.temp <- table.temp[table.temp$op=="=~",]
    table.temp <- table.temp[table.temp$est>=0.19,]
    table.temp <- table.temp[c("rhs","est")]
    table.temp <- merge(table.temp,variable.descriptions,by.x="rhs",by.y="variable")
    table.temp$factor <- paste0("Factor ", i)
    table.temp <- table.temp[c("factor","description","est")]
    table.temp <- table.temp[order(table.temp$est, decreasing = TRUE),]
    factor.names.table <- rbind(factor.names.table,table.temp)
  }
  return(factor.names.table)
}

efa.mice.varimax <- fa(psid15.wb16.data.mice[total.wb.vars], nfact = n.parallel, rotate = "varimax", fm = "minres", maxit=10000, weight=psid15.wb16.data.mice$WB16WT)
loadmat.mice.varimax <- zapsmall(matrix(round(efa.mice.varimax$loadings, 2), nrow = 41, ncol = n.parallel))
rownames(loadmat.mice.varimax) = colnames(psid15.wb16.data.mice[total.wb.vars])
terms <- vector()
for (i in 1:n.parallel) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice.varimax[,i]), "*", names(loadmat.mice.varimax[,1]), collapse = "+"))
}
esem.formula.mice.varimax <- paste(terms, collapse = "\n")
esem.mice.varimax <- lavaan::cfa(esem.formula.mice.varimax, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = total.wb.vars, missing="listwise")
esem.fit.mice.varimax  <- fitmeasures(esem.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice.varimax  <- make.factor.names.table(esem.mice.varimax, no.general = TRUE)

# Writing .csv so factors can be labelled.
write.csv(table.mice.varimax,"supplementary/varimax/psid_imputed_mice_varimax.csv", row.names = FALSE)

# Factors are labelled qualitatively.
factor.names <- c("Evaluative","Functioning","Positive Affect (Month)","Positive Affect (Day)","Frustration",
                  "Negative Affect (Month)","Nervousness","Pain/Tiredness","Loneliness",
                  "Anxiety / Stress (Day)","Respect","Calmness (Month)","Social Relationships","Evaluation", "Full of Life")

for (i in n.parallel:1) {
  table.mice.varimax$factor <- str_replace(table.mice.varimax$factor,paste0("Factor ",i),factor.names[i])
}

# Writing .csv with labelled factors.
write.csv(table.mice.varimax,"supplementary/varimax/psid_labelled_imputed_mice_varimax.csv", row.names = FALSE)

# Building Full Table
full.table.mice.varimax <- 1:n.parallel
y <- parameterEstimates(esem.mice.varimax, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in total.wb.vars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.parallel)
  full.table.mice.varimax <- rbind(full.table.mice.varimax,line)
}
full.table.mice.varimax <- as.data.frame(full.table.mice.varimax)
names(full.table.mice.varimax) <- paste0("Factor",1:n.parallel)
rownames(full.table.mice.varimax) <- c("Label", total.wb.vars)
full.table.mice.varimax[1,] <- factor.names
write.csv(full.table.mice.varimax,"supplementary/varimax/psid_full_factor_loadings_mice_varimax.csv", row.names = TRUE)

# MIMIC
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.parallel,"~", paste("income.decile", collapse = " + "), " + ", "female + age + arechildren +", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice.varimax, x)
mimic.1.mice.varimax <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = TRUE, ordered = total.wb.vars, missing = "listwise")
mimic.1.fit.mice.varimax <- fitmeasures(mimic.1.mice.varimax, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice.varimax)
all_models <- all_models %>% filter(grepl(paste0(" ~ ", "income.decile"), term))
all_models$model <- paste0("Factor",1:n.parallel)
all_models$term <- "income.decile"
all_models <- arrange(all_models, estimate)
for (i in n.parallel:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice.varimax <- all_models
mimic.1.graph.mice.varimax <- dwplot(all_models)
mimic.1.graph.mice.varimax <- mimic.1.graph.mice.varimax %>% relabel_predictors(income.total = "Income") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Varimax (Orthogonal)")
ggsave(file="supplementary/varimax/psid_income_controlled_imputed_mice_varimax.pdf", mimic.1.graph.mice.varimax, width = 276, height = 145, units = "mm")
write.csv(mimic.results.table.mice.varimax,"supplementary/varimax/psid_mimic_income_decile_controlled_imputed_mice_varimax.csv")
