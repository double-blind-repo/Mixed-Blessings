# Imputation ----
# Setting wellbeing variable types to factor (categorical)
factor.vars <- c(total.wb.vars,"marital","empstatus","health","sex")
psid15.wb16.data.for.imputation[factor.vars] <- lapply(psid15.wb16.data.for.imputation[factor.vars],as.factor)
psid15.wb16.data.mice <- psid15.wb16.data.for.imputation

# Setting up Parrallel Computing
n.cores <- detectCores()
n.cores <- ifelse(n.cores>5, 5, n.cores - 1)

# Imputing Data
system.time(
  imputed.datasets <- parlmice(as.data.frame(psid15.wb16.data.for.imputation), m = 5, maxit = 5, method="cart", cluster.seed = "34", n.core=n.cores, n.imp.core = 1, printFlag = TRUE)
)

# Saving workspace in case of crash.
save.image(file='imputation_session_complete_children.RData')

# Combining imputed datasets.
imputed_complete <- mice::complete(imputed.datasets)

# Filling NA wellbeing with imputed values.
psid15.wb16.data.mice[total.wb.vars] <- imputed_complete[total.wb.vars]

# Filling NA personality with imputed values.
psid15.wb16.data.mice[personality.recoded.vars] <- imputed_complete[personality.recoded.vars]

# Saving .csv in case of crash.
write.csv(psid15.wb16.data.mice,"psid_data_imputed_mice_children.csv", row.names = FALSE)



# Post imputation wrangling ----
psid15.wb16.data.mice <- cbind(psid15.wb16.data.mice,dummy_cols(psid15.wb16.data.mice["marital"]))
psid15.wb16.data.mice <- cbind(psid15.wb16.data.mice,dummy_cols(psid15.wb16.data.mice["empstatus"]))
psid15.wb16.data.mice <- cbind(psid15.wb16.data.mice,dummy_cols(psid15.wb16.data.mice["health"]))
psid15.wb16.data.mice["female"] <- ifelse(psid15.wb16.data.mice$sex==2,1,0)
psid15.wb16.data.mice$age2 <- psid15.wb16.data.mice$age^2
psid15.wb16.data.mice$arechildren <- ifelse(psid15.wb16.data.mice$children>0,1,0)
psid15.wb16.data.mice$arechildren <- as.factor(psid15.wb16.data.mice$arechildren)

psid15.wb16.data.mice$no.education <- ifelse(psid15.wb16.data.mice$education==0,1,0)
psid15.wb16.data.mice$less.than.high.school <- ifelse(psid15.wb16.data.mice$education<12&psid15.wb16.data.mice$education>0,1,0)
psid15.wb16.data.mice$high.school <- ifelse(psid15.wb16.data.mice$education==12,1,0)
psid15.wb16.data.mice$some.college <- ifelse(psid15.wb16.data.mice$education>12&psid15.wb16.data.mice$education<16,1,0)
psid15.wb16.data.mice$college <- ifelse(psid15.wb16.data.mice$education==16,1,0)
psid15.wb16.data.mice$some.postgraduate <- ifelse(psid15.wb16.data.mice$education==17,1,0)
edu.vars <- c("less.than.high.school","high.school",
              "some.college", "college","some.postgraduate")

names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'marital_1'] <- 'married'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'marital_2'] <- 'never.married'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'marital_3'] <- 'widowed'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'marital_4'] <- 'divorced.annulled'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'marital_5'] <- 'separated'
marital.vars <- c("married", "widowed", "divorced.annulled", "separated")

names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_1'] <- 'employed'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_2'] <- 'temp.unemployed.sick.maternity'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_3'] <- 'unemployed'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_4'] <- 'retired'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_5'] <- 'disabled'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_6'] <- 'keeping.house'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'empstatus_7'] <- 'student'
employment.vars <- c("employed","retired","disabled","keeping.house","student")

names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'health_1'] <- 'health.excellent'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'health_2'] <- 'health.very.good'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'health_3'] <- 'health.good'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'health_4'] <- 'health.fair'
names(psid15.wb16.data.mice)[names(psid15.wb16.data.mice) == 'health_5'] <- 'health.poor'
health.vars <- c("health.excellent","health.very.good","health.good","health.fair")

# Recoding Income Decile
psid15.wb16.data.mice$income.decile <- ntile(psid15.wb16.data.mice$income.total,10)
#psid15.wb16.data.mice$income.decile[is.na(psid15.wb16.data.mice$income.total)] <- 99
a <- fastDummies::dummy_cols(psid15.wb16.data.mice["income.decile"], ignore_na = TRUE)
psid15.wb16.data.mice <- cbind(psid15.wb16.data.mice,a)
decile.vars <- c(paste0("income.decile_",2:10))
#psid15.wb16.data.mice$income.decile[is.na(psid15.wb16.data.mice$income.total)] <- NA

# Log Income
psid15.wb16.data.mice$log.income.total = log(psid15.wb16.data.mice$income.total + 1)

# Recoding Income Decile
psid15.wb16.data.mice$income.decile <- ntile(psid15.wb16.data.mice$income.total,10)
#psid15.wb16.data$income.decile[is.na(psid15.wb16.data$income.total)] <- 99
a <- fastDummies::dummy_cols(psid15.wb16.data.mice["income.decile"], ignore_na = TRUE)
psid15.wb16.data.mice <- cbind(psid15.wb16.data.mice,a)
decile.vars <- c(paste0("income.decile_",2:10))



# ESEM ----
# Converting type back to numeric for FA.
psid15.wb16.data.mice[total.wb.vars] <- lapply(psid15.wb16.data.mice[total.wb.vars],as.numeric)

# Determining whether number of factors has changed.
# Remains 15.
n.parallel <- fa.parallel(psid15.wb16.data.mice[total.wb.vars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact

# Factor Names Table Function
make.factor.names.table <- function(model) {
  variable.descriptions = read.csv("variable_descriptions.csv")
  names(variable.descriptions) <- c("variable","description")
  factor.loadings <- parameterestimates(model)
  
  factor.names.table <- data.frame(factor=character(),
                                   description=character(), 
                                   est=numeric()) 
  for (i in 2:n.parallel) {
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

# Oblique Rotation Using All Items And Sample Weights
efa.mice <- fa(psid15.wb16.data.mice[total.wb.vars], nfact = n.parallel, rotate = "biquartimin", fm = "minres", maxit=10000, weight=psid15.wb16.data.mice$WB16WT)
loadmat.mice <- zapsmall(matrix(round(efa.mice$loadings, 2), nrow = 41, ncol = n.parallel))
rownames(loadmat.mice) = colnames(psid15.wb16.data.mice[total.wb.vars])
terms <- vector()
for (i in 1:n.parallel) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice[,i]), "*", names(loadmat.mice[,1]), collapse = "+"))
}
esem.formula.mice <- paste(terms, collapse = "\n")
esem.mice <- lavaan::cfa(esem.formula.mice, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
oblique.fit.mice <- fitmeasures(esem.mice, c("cfi.scaled","tli.scaled","rmsea.scaled"))
oblique.table.mice <- make.factor.names.table(esem.mice)

# Writing .csv so factors can be labelled.
write.csv(oblique.table.mice,"psid_oblique_imputed_mice.csv", row.names = FALSE)

# Factors are labelled qualitatively.

factor.names <- c("General","Negative Affect (Month)","Satisfaction","Positive Affect (Day)",
                  "Nervousness","Functioning","Positive Affect (Month)", "Evaluation",
                  "Respect", "Pain/Tiredness","Loneliness", "Calmness (Month)",
                  "Capability","Anxiety (Day)","Life Fullness")

for (i in n.parallel:1) {
  oblique.table.mice$factor <- str_replace(oblique.table.mice$factor,paste0("Factor ",i),factor.names[i])
}

# Writing .csv with labelled factors.
write.csv(oblique.table.mice,"psid_oblique_labled_imputed_mice.csv", row.names = FALSE)


# Building Full Table
full.table.mice <- 1:n.parallel
y <- parameterEstimates(esem.mice, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in total.wb.vars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.parallel)
  full.table.mice <- rbind(full.table.mice,line)
}
full.table.mice <- as.data.frame(full.table.mice)
names(full.table.mice) <- paste0("Factor",1:n.parallel)
rownames(full.table.mice) <- c("Label", total.wb.vars)
full.table.mice[1,] <- factor.names
write.csv(full.table.mice,"psid_full_factor_loadings_mice.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice <- esem.mice %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice <- esem.mice %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice, variance.table.mice) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice,"psid_mean_variance_mice.csv", row.names = FALSE)

# Income Decile Regression in SEM (With Controls) ----
x <- paste0("\nF",1:n.parallel,"~", paste("income.decile", collapse = " + "), " + ", "female + age + arechildren + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice, x)
mimic.1.mice <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing = "listwise")
mimic.1.fit.mice <- fitmeasures(mimic.1.mice, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice)
all_models <- all_models %>% filter(grepl(paste0(" ~ ", "income.decile"), term))
all_models$model <- paste0("Factor",1:n.parallel)
all_models$term <- "income.decile"
all_models <- arrange(all_models, estimate)
for (i in n.parallel:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice <- all_models
mimic.1.graph.mice <- dwplot(all_models)
mimic.1.graph.mice <- mimic.1.graph.mice %>% relabel_predictors(income.total = "Income") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("PSID: MIMIC Coefficients (With Controls)")
ggsave(file="psid_income_controlled_oblique_mice.pdf", mimic.1.graph.mice, width = 276, height = 145, units = "mm")
write.csv(mimic.results.table.mice,"psid_mimic_income_decile_controlled_mice.csv")

x <- paste0("\nF",1:n.parallel,"~", paste(decile.vars, collapse = " + "), " + ", "female + age + arechildren + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.mice, x)
mimic.2.mice <- lavaan::cfa(y, data = psid15.wb16.data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
all_models <- tidy(mimic.2.mice)
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
mimic.2.graph.mice <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], nrow=5, ncol=3)
ggsave(file="psid_income_categories_controlled_oblique_mice.pdf", mimic.2.graph.mice, width = 210, height = 297, units = "mm")