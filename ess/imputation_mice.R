# Imputation (Takes over 12 hrs) ----
# Setting wellbeing variable types to factor (categorical)
data.mice <- data.for.imputation
data.mice[wbvars] <- lapply(data.mice[wbvars],as.factor)

# Setting up Parrallel Computing
# Cores - 1 to allow for other work to continue on this PC. Using all cores uses 100% of CPU.
n.cores <- detectCores()
n.cores <- ifelse(n.cores>5, 5, n.cores - 1)

# Imputing Data
system.time(
  imputed.datasets <- parlmice(as.data.frame(data.mice), m = 5, maxit = 5, method="cart", cluster.seed = "34", n.core=n.cores, n.imp.core = 1, printFlag = TRUE)
)

# Saving workspace in case of crash.
save.image(file='imputation_session_complete.RData')

# Combining imputed datasets.
imputed_complete <- mice::complete(imputed.datasets)

# Filling NA wellbeing with imputed values.
data.mice <- imputed_complete

# Saving .csv in case of crash.
write.csv(data.mice,"ess_data_imputed_mice.csv", row.names = FALSE)


# Post imputation wrangling ----
#Creating Weight
data.mice["newweight"] <- data.mice["pspwght"]*data.mice["pweight"]

# Converting type back to numeric for FA.
data.mice[wbvars] <- lapply(data.mice[wbvars],as.numeric)

# Reversing Variables so that larger numbers mean better outcomes.
reverse <- c("dngval", "dclvlf", "accdng", "lotsgot", "optftr", "pstvms", "fltdpr", "flteeff", "slprl", "fltlnl", "fltsd", "cldgng",
             "fltanx", "flclpla","wkvlorg")
rev_fun <- function(variable) {
  max <- max(data.mice[variable], na.rm = TRUE)
  reversed <- (max + 1) - data.mice[variable]
  return(reversed)
}
reversed <- lapply(reverse,rev_fun)
data.mice[reverse] <- as.data.frame(reversed)

# Recoding Variables
data.mice$edu <- ifelse(data.mice$edulvlb<5555, round(data.mice$edulvlb/100) , 5555) 
data.mice$female <- ifelse(data.mice$gndr==2,1,0)
data.mice$agea2 <- data.mice$agea^2

data.mice <- cbind(data.mice,dummy_cols(data.mice["hinctnta"]))
data.mice <- cbind(data.mice,dummy_cols(data.mice["health"]))
names(data.mice)[names(data.mice) == 'health_1'] <- 'healthvg'
names(data.mice)[names(data.mice) == 'health_2'] <- 'healthg'
names(data.mice)[names(data.mice) == 'health_3'] <- 'healthf'
names(data.mice)[names(data.mice) == 'health_4'] <- 'healthb'
names(data.mice)[names(data.mice) == 'health_5'] <- 'healthvb'

data.mice <- cbind(data.mice,dummy_cols(data.mice["maritalb"]))
names(data.mice)[names(data.mice) == 'maritalb_1'] <- 'married'
names(data.mice)[names(data.mice) == 'maritalb_2'] <- 'civilu'
names(data.mice)[names(data.mice) == 'maritalb_3'] <- 'separated'
names(data.mice)[names(data.mice) == 'maritalb_4'] <- 'divorced'
names(data.mice)[names(data.mice) == 'maritalb_5'] <- 'widowed'
names(data.mice)[names(data.mice) == 'maritalb_6'] <- 'single'

data.mice["children"] <- as.numeric(data.mice["chldhhe"]==2)

data.mice <- cbind(data.mice, dummy.code(data.mice$cntry))
data.mice <- cbind(data.mice,dummy_cols(data.mice["hinctnta"]))
data.mice <- cbind(data.mice,dummy_cols(data.mice["edu"]))

data.mice <- cbind(data.mice,dummy_cols(data.mice["mnactic"]))
names(data.mice)[names(data.mice) == 'mnactic_1'] <- 'paidwork'
names(data.mice)[names(data.mice) == 'mnactic_2'] <- 'education'
names(data.mice)[names(data.mice) == 'mnactic_3'] <- 'unemployedlooking'
names(data.mice)[names(data.mice) == 'mnactic_4'] <- 'unemployednot'
names(data.mice)[names(data.mice) == 'mnactic_5'] <- 'sickdisabled'
names(data.mice)[names(data.mice) == 'mnactic_6'] <- 'retired'
names(data.mice)[names(data.mice) == 'mnactic_7'] <- 'service'
names(data.mice)[names(data.mice) == 'mnactic_8'] <- 'housework'
names(data.mice)[names(data.mice) == 'mnactic_9'] <- 'other'
data.mice$otheractivity <- ifelse(data.mice$other==1,1,ifelse(data.mice$service==1,1,0))

# Creating Control Variable Strings
personalvars <- "agea + female"
familyvars <- "children + married + civilu + separated + divorced + widowed"
healthvars <- "healthvg + healthg + healthf + healthb"
relincvars <- "hinctnta_2 + hinctnta_3 + hinctnta_4 + hinctnta_5 + hinctnta_6 + hinctnta_7 + hinctnta_8 + hinctnta_9 + hinctnta_10"
eduvars <- "edu_1 + edu_2 + edu_3 + edu_4 + edu_5 + edu_6 + edu_7 + edu_8 + edu_5555"
activityvars <- "paidwork + education + unemployednot + sickdisabled + retired + service + housework + other"
countryvars <- "BE + BG + CH + CY + CZ + DE + DK + EE + ES + FI + FR + HU + IE + IL + IS + IT + LT + NL + NO + PL + PT + RU + SE + SI + SK + UA + XK"


# ESEM ----
# Running Parralel Analysis. 19 Factors.
n.factors <- fa.parallel(data.mice[wbvars], fm = 'ml', fa = 'fa', n.iter = 50, SMC = TRUE, quant = .95)$nfact
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
efa.mice <- fa(data.mice[wbvars], nfact = n.factors, rotate = "biquartimin", fm = "minres", maxit=10000, weight = data.mice$newweight)
loadmat.mice <- zapsmall(matrix(efa.mice$loadings, nrow = 47, ncol = n.factors))
rownames(loadmat.mice) <- colnames(data.mice[wbvars])
terms <- vector()
for (i in 1:n.factors) {
  terms[i] <-
    paste0("F",i,"=~ ", paste0(c(loadmat.mice[,i]), "*", names(loadmat.mice[,1]), collapse = "+"))
}
esem.formula.mice <- paste(terms, collapse = "\n")
esem.mice <- lavaan::cfa(esem.formula.mice, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing="listwise")
esem.fit.mice <- fitmeasures(esem.mice, c("cfi.scaled","tli.scaled","rmsea.scaled"))
table.mice <- make.factor.names.table(esem.mice)
write.csv(table.mice,"ess_oblique_imputed_mice.csv", row.names = FALSE)
# Factors are labelled qualitatively.
factor.names <- c("General","Engagement","Trust","Social Closeness","Social Optimism", "Help and Support",
                  "Evaluation","Social Engagement","Positive Affect","Effort","Self-esteem",
                  "Sadness","Appreciation","Learning","Resilience","Calmness","Energy","Freedom","Self-Efficacy")
for (i in n.factors:1) {
  table.mice$factor <- str_replace(table.mice$factor,paste0("Factor ",i),factor.names[i])
}
write.csv(table.mice,"ess_oblique_labelled_imputed_mice.csv", row.names = FALSE)


# Building Full Table
full.table.mice <- 1:n.factors
y <- parameterEstimates(esem.mice, se = FALSE, zstat = FALSE, ci = FALSE)
for (value in wbvars) {
  line <- y$est[y$rhs == value]
  line <- head(line, n.factors)
  full.table.mice <- rbind(full.table.mice,line)
}
full.table.mice <- as.data.frame(full.table.mice)
names(full.table.mice) <- paste0("Factor",1:n.factors)
rownames(full.table.mice) <- c("Label", wbvars)
full.table.mice[1,] <- factor.names
write.csv(full.table.mice,"ess_full_factor_loadings_mice.csv", row.names = TRUE)

# Mean Variance Table
variance.table.mice <- esem.mice %>% parameterestimates() %>% filter(op == "~~" & grepl("F", lhs) & lhs == rhs) %>% select(lhs, est) %>% mutate(est = round(est, digits = 3)) %>% rename(variance = est)
mean.table.mice <- esem.mice %>% parameterestimates() %>% filter(op == "~1" & grepl("F", lhs)) %>% select(lhs, est) %>% rename(mean = est)
mean.variance.table.mice <- tibble(factor = factor.names) %>% bind_cols(.,mean.table.mice, variance.table.mice) %>% select(factor, mean, variance)
write.csv(mean.variance.table.mice,"ess_mean_variance_mice.csv", row.names = FALSE)

# MIMIC ----
# Income Decile Regression in SEM (With Controls)
x <- paste0("\nF",1:n.factors,"~ hinctnta + ",personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice, x)
mimic.1.mice <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.1.fit.mice <- fitmeasures(mimic.1.mice, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.1.mice)
all_models <- all_models %>% filter(grepl(" ~ hinctnta", term))
all_models$model <- paste0("Factor",1:n.factors)
all_models$term <- "hinctnta"
all_models <- arrange(all_models, estimate)
for (i in n.factors:1){
  all_models$model <- gsub(paste0("Factor",i), factor.names[i], all_models$model)
}
mimic.results.table.mice <- all_models
mimic.1.graph.mice <- dwplot(all_models)
mimic.1.graph.mice <- mimic.1.graph.mice %>% relabel_predictors(hinctnta = "Income Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank()) + ggtitle("Regression Coefficients (With Controls)") + xlim(-0.075,0.075)
ggsave(file="ess_income_decile_controlled_oblique_imputed_mice.pdf", mimic.1.graph.mice, width = 276, height = 145, units = "mm", device = pdf())
write.csv(mimic.results.table.mice,"ess_mimic_income_decile_controlled_imputed_mice.csv")

# Income Decile Regression in SEM (With Controls) - Categorical
x <- paste0("\nF",1:n.factors,"~", relincvars, " + ", personalvars," + ",familyvars," + ",healthvars, " + ", countryvars, "+", eduvars, "+", activityvars,collapse = "")
y <- paste0(esem.formula.mice, x)
mimic.2.mice <- lavaan::cfa(y, data = data.mice, verbose = F, estimator = "WLSMV", orthogonal = FALSE, ordered = wbvars, missing = "listwise")
mimic.2.fit.mice <- fitmeasures(mimic.2.mice, c("cfi.scaled","tli.scaled","rmsea.scaled"))
all_models <- tidy(mimic.2.mice)
inc_models <- all_models %>% filter(grepl("hinctnta", term))
graphs = list()
for (i in 1:n.factors){
  factor <- paste0("F",i," ~ ")
  factor.name <- factor.names[i]
  temp_models <- inc_models %>% filter(grepl(factor, term)) 
  temp_models["term"] <- paste0("hinctnta_",2:10)
  temp_plot <- dwplot(temp_models, alpha = 0.05, dodge_size = 0.00001)
  graphs[[i]] <- temp_plot %>% relabel_predictors(hinctnta_2 = "2nd Decile", hinctnta_3 = "3rd Decile", hinctnta_4 = "4th Decile", hinctnta_5 = "5th Decile", hinctnta_6 = "6th Decile", hinctnta_7 = "7th Decile", hinctnta_8 = "8th Decile", hinctnta_9 = "9th Decile", hinctnta_10 = "10th Decile") + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + theme_bw()  + guides(colour = guide_legend(reverse = TRUE)) + theme(legend.title = element_blank(), plot.title = element_text(size=8, hjust=0.3), legend.position="none") + ggtitle(paste0(factor.name)) + xlim(-0.4, 0.6)
}
mimic.2.graph.mice <- gridExtra::grid.arrange(graphs[1][[1]], graphs[2][[1]], graphs[3][[1]], graphs[4][[1]], graphs[5][[1]], graphs[6][[1]], graphs[7][[1]], graphs[8][[1]], graphs[9][[1]],graphs[10][[1]], graphs[11][[1]], graphs[12][[1]], graphs[13][[1]], graphs[14][[1]], graphs[15][[1]], graphs[16][[1]], graphs[17][[1]], graphs[18][[1]], graphs[19][[1]], nrow=7, ncol=3)
ggsave(file="ess_income_controls_categorical_controlled_imputed_mice.pdf", mimic.2.graph.mice, width = 210, height = 297, units = "mm", device = pdf())


