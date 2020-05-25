# Setting Working Directory and File Paths
data.directory = paste0(getwd(),"/data/")
main.data.directory = paste0(data.directory,"main/")
wb2016.data.path = paste0(data.directory,"WB2016/wbdata.dta")

# Build PSID ----
build.psid(datadr = main.data.directory, small = TRUE)

# Crosswalk file downloaded from http://psidonline.isr.umich.edu/help/xyr/psid.xlsx
cwf <- read_excel("psid_crosswalk.xlsx")
years.list <- 2015

family.id <- getNamesPSID("ER34301", cwf, years.list)
sequence.n <- getNamesPSID("ER34302", cwf, years.list)

income.total <- getNamesPSID("ER65349", cwf, years.list)
children <- getNamesPSID("ER60021", cwf, years.list)
sex <- getNamesPSID("ER32000", cwf, years.list)
age <- getNamesPSID("ER30004", cwf, years.list)
marital <- getNamesPSID("ER32049", cwf, years.list)
empstatus <- getNamesPSID("ER30293", cwf, years.list)
health <- getNamesPSID("ER33128", cwf, years.list)
education <- getNamesPSID("ER30010", cwf, years.list)

famvars <- data.frame(year=years.list, income.total=income.total, children=children)
indvars <- data.frame(year=years.list, sequence.n=sequence.n, family.id=family.id,
                      age=age,marital=marital,empstatus=empstatus,health=health, sex=sex, education=education)
                     
my.dir <- main.data.directory
data <- build.panel(datadir=my.dir,fam.vars=famvars, ind.vars=indvars, heads.only=FALSE)
data <- as_tibble(data)

# WB2016 Data ----
# Download the 2016 WellBeing and Daily Life supplement from https://simba.isr.umich.edu/Zips/ZipMain.aspx
# Run the STATA.do file and save the output as wbdata.dta in the WB2016 subfolder.
wb.data <- read.dta(wb2016.data.path)
names(wb.data)[names(wb.data) == 'WB16YRID'] <- 'family.id'
names(wb.data)[names(wb.data) == 'WB16SN'] <- 'sequence.n'

names(wb.data)[names(wb.data) == 'WB16A1'] <- 'life.sat.whole'
names(wb.data)[names(wb.data) == 'WB16A2'] <- 'best.life'
names(wb.data)[names(wb.data) == 'WB16A3A'] <- 'ideal.life'
names(wb.data)[names(wb.data) == 'WB16A3B'] <- 'life.conditions'
names(wb.data)[names(wb.data) == 'WB16A3C'] <- 'sat.with.life'
names(wb.data)[names(wb.data) == 'WB16A3D'] <- 'got.important.things'
names(wb.data)[names(wb.data) == 'WB16A3E'] <- 'change.nothing'
global.wb.vars <- c('life.sat.whole','best.life', 'ideal.life','life.conditions','sat.with.life','got.important.things','change.nothing')

names(wb.data)[names(wb.data) == 'WB16A6A'] <- 'purpose.meaning'
names(wb.data)[names(wb.data) == 'WB16A6B'] <- 'supportive.relationships'
names(wb.data)[names(wb.data) == 'WB16A6C'] <- 'engaged'
names(wb.data)[names(wb.data) == 'WB16A6D'] <- 'contribute.happiness'
names(wb.data)[names(wb.data) == 'WB16A6E'] <- 'capable'
names(wb.data)[names(wb.data) == 'WB16A6F'] <- 'good.person.life'
names(wb.data)[names(wb.data) == 'WB16A6G'] <- 'optimistic'
names(wb.data)[names(wb.data) == 'WB16A6H'] <- 'respected'
self.eval.vars <- c('purpose.meaning','supportive.relationships','engaged','contribute.happiness','capable','good.person.life','optimistic','respected')

names(wb.data)[names(wb.data) == 'WB16B1A'] <- 'cheerful'
names(wb.data)[names(wb.data) == 'WB16B1B'] <- 'good.spirits'
names(wb.data)[names(wb.data) == 'WB16B1C'] <- 'extreme.happiness'
names(wb.data)[names(wb.data) == 'WB16B1D'] <- 'calm.peaceful'
names(wb.data)[names(wb.data) == 'WB16B1E'] <- 'felt.satisfied'
names(wb.data)[names(wb.data) == 'WB16B1F'] <- 'full.of.life'
names(wb.data)[names(wb.data) == 'WB16B2A'] <- 'not.cheer.up'
names(wb.data)[names(wb.data) == 'WB16B2B'] <- 'nervous'
names(wb.data)[names(wb.data) == 'WB16B2C'] <- 'restless.figity'
names(wb.data)[names(wb.data) == 'WB16B2D'] <- 'hopeless'
names(wb.data)[names(wb.data) == 'WB16B2E'] <- 'everything.effort'
names(wb.data)[names(wb.data) == 'WB16B2F'] <- 'felt.worthless'
month.vars <- c('cheerful','good.spirits','extreme.happiness','calm.peaceful','felt.satisfied','full.of.life',
                'not.cheer.up','nervous','restless.figity','hopeless','everything.effort','felt.worthless')


names(wb.data)[names(wb.data) == 'WB16C14A'] <- 'day.calm'
names(wb.data)[names(wb.data) == 'WB16C14B'] <- 'day.happy'
names(wb.data)[names(wb.data) == 'WB16C14C'] <- 'day.enthusiastic'
names(wb.data)[names(wb.data) == 'WB16C14D'] <- 'day.content'
names(wb.data)[names(wb.data) == 'WB16C14E'] <- 'day.interested'
names(wb.data)[names(wb.data) == 'WB16C15A'] <- 'day.angry'
names(wb.data)[names(wb.data) == 'WB16C15B'] <- 'day.frustrated'
names(wb.data)[names(wb.data) == 'WB16C15C'] <- 'day.sad'
names(wb.data)[names(wb.data) == 'WB16C15D'] <- 'day.stressed'
names(wb.data)[names(wb.data) == 'WB16C15E'] <- 'day.lonely'
names(wb.data)[names(wb.data) == 'WB16C15F'] <- 'day.worried'
names(wb.data)[names(wb.data) == 'WB16C15G'] <- 'day.bored'
names(wb.data)[names(wb.data) == 'WB16C16A'] <- 'day.tired'
names(wb.data)[names(wb.data) == 'WB16C16B'] <- 'day.pain'

day.vars <- c('day.calm','day.happy','day.enthusiastic','day.content','day.interested','day.angry',
              'day.frustrated','day.sad','day.stressed','day.lonely','day.worried','day.bored','day.tired','day.pain')


total.wb.vars <- c(global.wb.vars, self.eval.vars, month.vars, day.vars)
total.wb.vars.temp <- total.wb.vars[!total.wb.vars %in% "best.life"]

# Replacing Missing Values
wb.data <- as_tibble(wb.data)
wb.data[total.wb.vars.temp]<- replace_with_na_all(data = wb.data[total.wb.vars.temp], condition = ~.x == "9")
wb.data["best.life"] <- replace_with_na_all(data = wb.data["best.life"], condition = ~.x == "99")

# Reversing some indicators such that the largest positive number corresponds with the 'best' outcome.
reverse <- c("life.sat.whole", "ideal.life", "life.conditions", "sat.with.life", "got.important.things", "change.nothing",
             "purpose.meaning", "supportive.relationships", "engaged", "contribute.happiness",
             "capable", "good.person.life", "optimistic", "respected", "cheerful", "good.spirits",
             "extreme.happiness", "calm.peaceful", "felt.satisfied", "full.of.life",
             "day.calm", "day.happy", "day.enthusiastic", "day.content", "day.interested")
rev_fun <- function(variable) {
  max <- max(wb.data[variable], na.rm = TRUE)
  reversed <- (max + 1) - wb.data[variable]
  return(reversed)
}
reversed <- lapply(reverse,rev_fun)
wb.data[reverse] <- as.data.frame(reversed)

# Personality
conscientiousness.vars <- c("WB16D1A", "WB16D1G", "WB16D1K")
extraversion.vars <- c("WB16D1B", "WB16D1H", "WB16D1L")
agreeableness.vars <- c("WB16D1C", "WB16D1F", "WB16D1M")
openness.vars <- c("WB16D1D", "WB16D1I", "WB16D1N")
neuroticism.vars <- c("WB16D1E", "WB16D1J", "WB16D1O")
personality.vars <- c(conscientiousness.vars,extraversion.vars,agreeableness.vars,openness.vars,neuroticism.vars)

wb.data[personality.vars]<- replace_with_na_all(data = wb.data[personality.vars], condition = ~.x == "9")

personality.reverse <- c("WB16D1A", "WB16D1K", "WB16D1B", "WB16D1H", "WB16D1F", "WB16D1M", "WB16D1D", "WB16D1I", "WB16D1N",
                         "WB16D1E", "WB16D1J")
reversed <- lapply(personality.reverse,rev_fun)
wb.data[personality.reverse] <- as.data.frame(reversed)

wb.data['conscientiousness'] <- rowMeans(wb.data[conscientiousness.vars], na.rm = TRUE)
wb.data['extraversion'] <- rowMeans(wb.data[extraversion.vars], na.rm = TRUE)
wb.data['agreeableness'] <- rowMeans(wb.data[agreeableness.vars], na.rm = TRUE)
wb.data['openness'] <- rowMeans(wb.data[openness.vars], na.rm = TRUE)
wb.data['neuroticism'] <- rowMeans(wb.data[neuroticism.vars], na.rm = TRUE)

personality.recoded.vars <- c("conscientiousness", "extraversion", "agreeableness", "openness", "neuroticism")

wb.data <- wb.data[c("family.id", "sequence.n","WB16WT",total.wb.vars,personality.recoded.vars)]

# Merge
psid15.wb16.data <- merge(data,wb.data,by=c("family.id", "sequence.n"))

# Creating Copy For Imputation
psid15.wb16.data.for.imputation <- psid15.wb16.data