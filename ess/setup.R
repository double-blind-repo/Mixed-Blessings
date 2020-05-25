# Downloading ESS Data ----
set_email("doubleblind@email.com")
data <- import_rounds(6)
data <- recode_missings(data)

# These are the wellbeing variables from the Round 6 ESS Module. Excluding those concerning work.
wbvars <- c("accdng","cldgng","dclvlf","deaimpp", "dngval","enjlf",
            "enrglot", "flapppl", "flclpla", "flrms", "fltanx","fltdpr",
            "flteeff", "fltlnl","fltpcfl", "fltsd", "happy", "inprdsc",
            "lchshcp", "lfwrs", "lotsgot", "lrnntlf", "nhpftr","optftr", 
            "physact", "plinsoc", "pplahlp", "pplfair", "pplhlp","ppltrst",
            "prhlppl", "pstvms","rehlppl", "sclact","sclmeet", "sedirlf",
            "slprl", "stflife", "tmabdng", "tmdotwa", "tmendng", "tmimdng",
            "tnapsur", "trtrsp","wkvlorg", "wrbknrm", "wrhpp")

other.vars <- c("pspwght","pweight","inwyys","agea","edulvlb","gndr","hinctnta",
                "health","maritalb","chldhhe","mnactic",
                "cntry")

data <- data[c(wbvars,other.vars)]

# No well-being data for Albania
data <- data[!data$cntry == "AL",]
data$cntry <- as.factor(data$cntry)

# Creating Copy For Imputation
data.for.imputation <- data