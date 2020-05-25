# Creating Weight
data["newweight"] <- data["pspwght"]*data["pweight"]

# Reversing Variables so that larger numbers correspond to better outcomes. ----
reverse <- c("dngval", "dclvlf", "accdng", "lotsgot", "optftr", "pstvms", "fltdpr", "flteeff", "slprl", "fltlnl", "fltsd", "cldgng",
             "fltanx", "flclpla","wkvlorg")
rev_fun <- function(variable) {
  max <- max(data[variable], na.rm = TRUE)
  reversed <- (max + 1) - data[variable]
  return(reversed)
}
reversed <- lapply(reverse,rev_fun)
data[reverse] <- as.data.frame(reversed)

# Recoding Variables ----
data$edu <- ifelse(data$edulvlb<5555, round(data$edulvlb/100) , 5555) 
data$agea2 <- data$agea^2
data$female <- ifelse(data$gndr==2,1,0)

data <- cbind(data,dummy_cols(data["hinctnta"]))
data <- cbind(data,dummy_cols(data["health"]))
names(data)[names(data) == 'health_1'] <- 'healthvg'
names(data)[names(data) == 'health_2'] <- 'healthg'
names(data)[names(data) == 'health_3'] <- 'healthf'
names(data)[names(data) == 'health_4'] <- 'healthb'
names(data)[names(data) == 'health_5'] <- 'healthvb'

data <- cbind(data,dummy_cols(data["maritalb"]))
names(data)[names(data) == 'maritalb_1'] <- 'married'
names(data)[names(data) == 'maritalb_2'] <- 'civilu'
names(data)[names(data) == 'maritalb_3'] <- 'separated'
names(data)[names(data) == 'maritalb_4'] <- 'divorced'
names(data)[names(data) == 'maritalb_5'] <- 'widowed'
names(data)[names(data) == 'maritalb_6'] <- 'single'

data["children"] <- as.numeric(data["chldhhe"]==2)

data <- cbind(data, dummy.code(data$cntry))
data <- cbind(data,dummy_cols(data["hinctnta"]))
data <- cbind(data,dummy_cols(data["edu"]))
data <- cbind(data,dummy_cols(data["mnactic"]))
names(data)[names(data) == 'mnactic_1'] <- 'paidwork'
names(data)[names(data) == 'mnactic_2'] <- 'education'
names(data)[names(data) == 'mnactic_3'] <- 'unemployedlooking'
names(data)[names(data) == 'mnactic_4'] <- 'unemployednot'
names(data)[names(data) == 'mnactic_5'] <- 'sickdisabled'
names(data)[names(data) == 'mnactic_6'] <- 'retired'
names(data)[names(data) == 'mnactic_7'] <- 'service'
names(data)[names(data) == 'mnactic_8'] <- 'housework'
names(data)[names(data) == 'mnactic_9'] <- 'other'

data$otheractivity <- ifelse(data$other==1,1,ifelse(data$service==1,1,0))

# Creating Control Variable Strings ----
personalvars <- "agea + female"
familyvars <- "children + married + civilu + separated + divorced + widowed"
healthvars <- "healthvg + healthg + healthf + healthb"
relincvars <- "hinctnta_2 + hinctnta_3 + hinctnta_4 + hinctnta_5 + hinctnta_6 + hinctnta_7 + hinctnta_8 + hinctnta_9 + hinctnta_10"
eduvars <- "edu_1 + edu_2 + edu_3 + edu_4 + edu_5 + edu_6 + edu_7 + edu_8 + edu_5555"
activityvars <- "paidwork + education + unemployednot + sickdisabled + retired + service + housework + other"
countryvars <- "BE + BG + CH + CY + CZ + DE + DK + EE + ES + FI + FR + HU + IE + IL + IS + IT + LT + NL + NO + PL + PT + RU + SE + SI + SK + UA + XK"
factorvars <- paste0("Factor",1:18, collapse = " + ")

# ESEM ----
# Running Parralel Analysis
# 18 Retained
n.factors <- fa.parallel(data[wbvars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact
factorvars <- paste0("Factor",1:n.factors, collapse = " + ")

# Factor Names Table Builder
# Variable descriptions source: https://www.europeansocialsurvey.org/docs/round6/survey/ESS6_data_protocol_e01_4.pdf
make.factor.names.table <- function(model) {
  variable.descriptions = read.csv("variable_descriptions.csv")
  names(variable.descriptions) <- c("variable","description")
  factor.loadings <- parameterestimates(model)
  
  factor.names.table <- data.frame(factor=character(),
                                   description=character(), 
                                   est=numeric()) 
  for (i in 2:n.factors) {
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
efa.listwise <- fa(data[wbvars], nfact = n.factors, rotate = "biquartimin", fm = "minres", maxit=10000, weight = data$newweight)
loadmat.listwise <- zapsmall(matrix(round(efa.listwise$loadings, 2), nrow = 47, ncol = n.factors))
rownames(loadmat.listwise) <- colnames(data[wbvars])
terms <- vector()
for (i in 1:n.factors) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.listwise[,i]), "*", names(loadmat.listwise[,1]), collapse = "+"))
}
esem.formula.listwise <- paste(terms, collapse = "\n")
esem.listwise <- lavaan::cfa(esem.formula.listwise, data = data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing="listwise")
oblique.fit.listwise <- fitmeasures(esem.listwise , c("cfi.scaled","tli.scaled","rmsea.scaled"))
oblique.table <- make.factor.names.table(esem.listwise)
write.csv(oblique.table,"ess_oblique_listwise.csv", row.names = FALSE)
# Factors are labelled qualitatively.
factor.names <- c("General","Engagement","Trust","Social Closeness","Social Optimism", "Help and Support",
                  "Evaluation","Social Engagement","Positive Affect","Self-Esteem", "Sadness",
                  "Purpose","Effort","Appreciated","Learning","Resilience","Self-Efficacy","Selflessness")
for (i in n.factors:1) {
  oblique.table$factor <- str_replace(oblique.table$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(oblique.table,"ess_oblique_labelled_listwise.csv", row.names = FALSE)

# Building Full Table ----
full.table.listwise <- 1:n.factors
y <- parameterEstimates(esem.listwise, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in wbvars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.factors)
  full.table.listwise <- rbind(full.table.listwise,line)
}
full.table.listwise <- as.data.frame(full.table.listwise)
names(full.table.listwise) <- paste0("Factor",1:n.factors)
rownames(full.table.listwise) <- c("Label", wbvars)
full.table.listwise[1,] <- factor.names
write.csv(full.table.listwise,"ess_full_factor_loadings_listwise.csv", row.names = TRUE)

# Mean Variance Table
variance.table.listwise <- esem.listwise %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.listwise <- esem.listwise %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.listwise <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.listwise, variance.table.listwise) %>% select(factor, mean, variance)
write.csv(mean.variance.table.listwise,"ess_mean_variance_listwise.csv", row.names = FALSE)

# MIMIC ----
# Income Decile Regression in SEM (With Controls) ----
x <- paste0("\nF",1:18,"~ hinctnta + ",personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.listwise, x)
mimic.1.listwise <- lavaan::cfa(y, data = data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.1.fit.listwise <- fitmeasures(mimic.1.listwise, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.listwise)
all_models <- all_models %>% filter(grepl(" ~ hinctnta", term))
all_models$model <- paste0("Factor",1:18)
all_models$term <- "hinctnta"
all_models <- arrange(all_models, estimate)
for (i in 18:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.listwise <- all_models
mimic.1.graph.listwise <- dwplot(all_models)
mimic.1.graph.listwise <- mimic.1.graph.listwise %>% relabel_predictors(hinctnta = "Income Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Regression Coefficients (With Controls)") + xlim(-0.075,0.075)
ggsave(file="income_controls_continous_oblique_listwise.pdf", mimic.1.graph.listwise, width = 276, height = 145, units = "mm", device = pdf())
write.csv(mimic.results.table.listwise,"ess_mimic_income_decile_controlled_listwise.csv")

# Income Decile Regression in SEM (With Controls) - Categorical ----
x <- paste0("\nF",1:18,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.listwise, x)
mimic.2.listwise <- lavaan::cfa(y, data = data, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.2.fit.listwise <- fitmeasures(mimic.2.listwise, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.listwise)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:18){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.6, 0.8)
}
mimic.2.graph.listwise <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], nrow=6, ncol=3)
ggsave(file="income_controls_categorical_oblique_listwise.pdf", mimic.2.graph.listwise, width = 210, height = 297, units = "mm", device = pdf())
