#Relabelling
psid15.wb16.data <- cbind(psid15.wb16.data,dummy_cols(psid15.wb16.data["marital"]))
psid15.wb16.data <- cbind(psid15.wb16.data,dummy_cols(psid15.wb16.data["empstatus"]))
psid15.wb16.data <- cbind(psid15.wb16.data,dummy_cols(psid15.wb16.data["health"]))
psid15.wb16.data["female"] <- psid15.wb16.data["sex"]-1
psid15.wb16.data$arechildren <- ifelse(psid15.wb16.data$children>0,TRUE,FALSE)
psid15.wb16.data$arechildren <- as.factor(psid15.wb16.data$arechildren)

psid15.wb16.data$no.education <- ifelse(psid15.wb16.data$education==0,1,0)
psid15.wb16.data$less.than.high.school <- ifelse(psid15.wb16.data$education<12&psid15.wb16.data$education>0,1,0)
psid15.wb16.data$high.school <- ifelse(psid15.wb16.data$education==12,1,0)
psid15.wb16.data$some.college <- ifelse(psid15.wb16.data$education>12&psid15.wb16.data$education<16,1,0)
psid15.wb16.data$college <- ifelse(psid15.wb16.data$education==16,1,0)
psid15.wb16.data$some.postgraduate <- ifelse(psid15.wb16.data$education==17,1,0)
edu.vars <- c("less.than.high.school","high.school",
              "some.college", "college","some.postgraduate")

names(psid15.wb16.data)[names(psid15.wb16.data) == 'marital_1'] <- 'married'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'marital_2'] <- 'never.married'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'marital_3'] <- 'widowed'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'marital_4'] <- 'divorced.annulled'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'marital_5'] <- 'separated'
marital.vars <- c("married", "widowed", "divorced.annulled", "separated")

names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_1'] <- 'employed'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_2'] <- 'temp.unemployed.sick.maternity'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_3'] <- 'unemployed'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_4'] <- 'retired'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_5'] <- 'disabled'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_6'] <- 'keeping.house'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'empstatus_7'] <- 'student'
employment.vars <- c("employed","retired","disabled","keeping.house","student")

names(psid15.wb16.data)[names(psid15.wb16.data) == 'health_1'] <- 'health.excellent'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'health_2'] <- 'health.very.good'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'health_3'] <- 'health.good'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'health_4'] <- 'health.fair'
names(psid15.wb16.data)[names(psid15.wb16.data) == 'health_5'] <- 'health.poor'
health.vars <- c("health.excellent","health.very.good","health.good","health.fair")

# Recoding Income Decile
psid15.wb16.data$income.decile <- ntile(psid15.wb16.data$income.total,10)
#psid15.wb16.data$income.decile[is.na(psid15.wb16.data$income.total)] <- 99
a <- fastDummies::dummy_cols(psid15.wb16.data["income.decile"], ignore_na = TRUE)
psid15.wb16.data <- cbind(psid15.wb16.data,a)
decile.vars <- c(paste0("income.decile_",2:10))
#psid15.wb16.data$income.decile[is.na(psid15.wb16.data$income.total)] <- NA

# Log Income
psid15.wb16.data$log.income.total = log(psid15.wb16.data$income.total + 1)

# ESEM ----
#15 Factors
n.parallel <- fa.parallel(psid15.wb16.data[total.wb.vars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact

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
b5.efa = fa(psid15.wb16.data[total.wb.vars], nfact = n.parallel, rotate = "biquartimin", fm = "minres", maxit=10000, weight=psid15.wb16.data$WB16WT)
b5.loadmat = zapsmall(matrix(round(b5.efa$loadings, 2), nrow = 41, ncol = n.parallel))
rownames(b5.loadmat) = colnames(psid15.wb16.data[total.wb.vars])
terms = vector()
for (i in 1:n.parallel) {
  terms[i] =
    paste0("F",i,"=~ ", paste0(c(b5.loadmat[,i]), "*", names(b5.loadmat[,1]), collapse = "+"))
}
esem.formula.listwise <- paste(terms, collapse = "\n")
esem.listwise <- lavaan::cfa(esem.formula.listwise, data = psid15.wb16.data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
oblique.fit.listwise <- fitmeasures(esem.listwise, c("cfi.scaled","tli.scaled","rmsea.scaled"))
oblique.table.listwise <- make.factor.names.table(esem.listwise)
write.csv(oblique.table.listwise,"psid_oblique_listwise.csv", row.names = FALSE)
factor.names <- c("General","Negative Affect (Month)","Positive Affect (Day)","Satisfaction","Nervousness",
                  "Frustration","Functioning","Positive Affect (Month)","Respected","Pain/Tiredness",
                  "Loneliness", "Evaluation", "Capability", "Calmness (Month)", "Full of Life (Month)")
for (i in n.parallel:1) {
  oblique.table.listwise$factor <- str_replace(oblique.table.listwise$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(oblique.table.listwise,"psid_oblique_labled_listwise.csv", row.names = FALSE)

# Building Full Table
full.table.listwise <- 1:n.parallel
y <- parameterEstimates(esem.listwise, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in total.wb.vars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.parallel)
  full.table.listwise <- rbind(full.table.listwise,line)
}
full.table.listwise <- as.data.frame(full.table.listwise)
names(full.table.listwise) <- paste0("Factor",1:n.parallel)
rownames(full.table.listwise) <- c("Label", total.wb.vars)
full.table.listwise[1,] <- factor.names
write.csv(full.table.listwise,"psid_full_factor_loadings_listwise.csv", row.names = TRUE)

# Mean Variance Table
variance.table.listwise <- esem.listwise %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.listwise <- esem.listwise %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.listwise <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.listwise, variance.table.listwise) %>% select(factor, mean, variance)
write.csv(mean.variance.table.listwise,"psid_mean_variance_listwise.csv", row.names = FALSE)

# MIMIC ----
# Income Decile Regression in SEM (With Controls) ----
# arechidren control variable omitted as it causes the model to fail to converge.
x <- paste0("\nF",1:n.parallel,"~", paste("income.decile", collapse = " + "), " + ", "female + age + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.listwise, x)
mimic.1.listwise <- lavaan::cfa(y, data = psid15.wb16.data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing = "listwise")
all_models <- tidy(mimic.1.listwise)
all_models <- all_models %>% filter(grepl(paste0(" ~ ", "income.decile"), term))
all_models$model <- paste0("Factor",1:n.parallel)
all_models$term <- "income.decile"
all_models <- arrange(all_models, estimate)
for (i in n.parallel:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.listwise <- all_models
mimic.1.graph.listwise <- dwplot(all_models)
mimic.1.graph.listwise <- mimic.1.graph.listwise %>% relabel_predictors(income.total = "Income") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("PSID: MIMIC Coefficients (With Controls)")
ggsave(file="psid_income_controlled_oblique_listwise.pdf", mimic.1.graph.listwise, width = 210, height = 145, units = "mm")
write.csv(mimic.results.table.listwise,"psid_mimic_income_decile_controlled_listwise.csv")

x <- paste0("\nF",1:n.parallel,"~", paste(decile.vars, collapse = " + "), " + ", "female + age + arechildren + ", paste(health.vars, collapse = " + ")," + ",paste(employment.vars, collapse = " + ")," + ",paste(marital.vars, collapse = " + "), " + ", paste(personality.recoded.vars, collapse = "+"), " + ", paste(edu.vars, collapse = "+"), collapse = "")
y <- paste0(esem.formula.listwise, x)
mimic.2.listwise <- lavaan::cfa(y, data = psid15.wb16.data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = total.wb.vars, missing="listwise")
all_models <- tidy(mimic.2.listwise)
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
mimic.2.graph.listwise <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], nrow=5, ncol=3)
ggsave(file="psid_income_categories_controlled_oblique_listwise.pdf", mimic.2.graph.listwise, width = 210, height = 297, units = "mm")