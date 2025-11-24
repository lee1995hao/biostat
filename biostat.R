file_list_list <- list.files("/Users/holee/Downloads/EDC_Raw_Datasets_Respiratory_100subjects_20251105_080749_qqocoa")
keio_caonima_1 <- paste0("/Users/holee/Downloads/EDC_Raw_Datasets_Respiratory_100subjects_20251105_080749_qqocoa","/",file_list_list)
ae <- read.csv(keio_caonima_1[which(keio_caonima_1 %like% "ae")]) 
keio_caonima_row <- read.csv(keio_caonima_1[8])
keio_caonima_row %>% mutate(
  STUDYID = study_id,
  VISIT   = visit_name,
  USUBJID = paste0(study_id,"-",subject_id),
  VSDTC = paste0(collection_date,"T",collection_time),
  LBTEST  = test_name,
  LBORRES = result_value,
  LBORRESU = result_unit,
  LBSTNRLO = reference_range_low,
  LBSTNRHI = reference_range_high,
  LBSTNRFL = abnormal_flag,
  LBNRIND  = clinical_significance,
  LBLOC    = laboratory_name,
  LBTECH   = technician_id
  ) %>% left_join(test_map,by = "LBTEST")  # 补充 LBTESTCD

keio_caonima_row["test_name"]


test_map <- data.frame(
  LBTEST = c("Hemoglobin", "Glucose", "Sodium", "Potassium", "Chloride", "Creatinine"),
  LBTESTCD = c("HGB", "GLUC", "NA", "K", "CL", "CREAT"),
  c(1:6)
)



library(dplyr)

df <- tibble(
  group = c("A","A","B","A","B"),
  value = c(10, 5, 8, NA, 3)
)

data.frame(df%>% drop_na())
df %>% group_by(group)%>%
  mutate(ddd = rank(value))


df <- tibble(
  student = c("Alice", "Bob", "Charlie"),
  scores  = list(
    c(90, 95, 100),    # Alice 的三次成绩
    c(82, 88),         # Bob 的两次成绩
    c(70, 75, 80, 85)  # Charlie 的四次成绩
  )
)

library(survival)
library(survminer)   # ggsurvplot() 属于这个包


data(lung)
lung$sex <- factor(lung$sex, labels = c("Male","Female"))

fit_cox <- coxph(
  Surv(time, status == 2) ~ age + sex,
  data = lung
)

newdata <- data.frame(
  age = c(60, 60),
  sex = factor(c("Male", "Female"), levels = c("Male", "Female"))
)
fit_surv <- survfit(fit_cox, newdata = newdata)




ggsurvplot(
  fit_surv,
  data = newdata
)



summary(fit_surv, times = 1000)$surv


set.seed(123)

n <- 120                          # 样本量
df <- data.frame(
  id = 1:n,
  treatment = sample(c("Treatment","Control"), n, replace=TRUE),
  ALT = c(rnorm(n/2, 65, 20), rnorm(n/2, 55, 17)),   # 生物指标
  age = rnorm(n, 60, 8),
  time = round(runif(n, 30, 180)),                   # 随访天数
  event = rbinom(n, 1, 0.4)                          # 1=死亡 / 事件
)
head(df)

library(survival)
library(survminer)
ssss <- Surv(df$time,df$event)
data_s <- survfit(ssss~ treatment,data = df)
ggsurvplot(data_s,data = df,conf.int = T,risk.table = T,pval = T)

ggplot(repeat_df, aes(x = visit_day, y = ALT, group = id, color = treatment)) +
  geom_line(alpha = 0.2)



repeat_df <- df %>%
  slice(rep(1:n(), each = 4)) %>%
  mutate(visit_day = rep(seq(0, 180, length.out = 4), times = n),
         ALT = ALT + rnorm(n * 4, 0, 10))  # 随时间波动



model <- glm(event ~ treatment + age, data = df, family = "binomial")
summary <- tidy(model, exponentiate = TRUE, conf.int = TRUE)

fit <- coxph(Surv(time, event) ~ treatment + age + ALT, data = df)
tb <- tidy(fit)


tb <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)

# 按 HR 远离 1 的程度排序（可按 estimate 或 p 值）
tb <- tb %>% mutate(term_clean = fct_reorder(term_clean, estimate))

p <- ggplot(tb, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = 2) +          # 参考线 HR=1
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_x_log10(                                       # HR 常用对数刻度
    breaks = c(0.5, 0.7, 1, 1.5, 2, 3, 5),
    limits = c(min(0.5, min(tb$conf.low, na.rm=TRUE)), 
               max(5,   max(tb$conf.high, na.rm=TRUE)))
  ) +
  labs(
    title = "Cox PH Forest Plot",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal()

# 例：按 sex 做子组（假设 df$sex 已存在且为 factor）

library(survival)
library(broom)
library(ggplot2)

fit <- coxph(Surv(time, event) ~ treatment + age + ALT, data = df)
zph <- cox.zph(fit)  
tb <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)

ggplot(tb, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "", title = "Forest Plot (Cox Model)") +
  theme_minimal()


check_cox_diagnostics(fit)

basehaz(fit, centered = F)  
cox.zph(fit)  



df <- tibble::tibble(
  group = rep(c("Treatment", "Placebo"), each = 50),
  hba1c = c(rnorm(50, mean = 6.8, sd = 0.4), rnorm(50, mean = 7.2, sd = 0.5))
)
new_df <- df %>% group_by(group) %>% mutate(numbeer = row_number()) 

ggplot(data = new_df, mapping = aes(x = numbeer,y = hba1c, group = group, colour = group,))+
  geom_line()



ggplot(mtcars, aes(x = wt, y = mpg, size = hp, colour = factor(am))) +
  geom_point(aes(alpha = drat,shape = as.factor(gear))) + 
  scale_colour_manual(values = c("1" = "blue", "0" = "green"))



ggplot(mtcars, aes(x = wt, y = mpg, label = rownames(mtcars))) +
  geom_text()



data(lung)
lung$sex <- factor(lung$sex, labels = c("Male", "Female"))

fit <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
summary(fit)


ggforest(fit, data = lung)


data_out <- data.frame(x = c(1:20), y = (c(1:20))*3 + rnorm(20), z = c(rep(1,10),rep(2,10))) 
out_idx <- sample(1:nrow(data_out),3)
data_out[out_idx,2] <- data_out[out_idx,2] + 1000
 
ggplot(data_out, aes(x = z, y = y,group = z)) +
  geom_boxplot(
    outlier.colour = "red"    # 异常值颜色
  )


panduanyichang <- function(x){
  center <- median(x)             # 也可以换成 mean(x)
  spread <- sd(x, na.rm = TRUE)
  outlier_flag <- ifelse(x < center - 1.5 * spread|x > center + 1.5 * spread,1,2)
  return(outlier_flag)            # 返回逻辑向量 TRUE/FALSE
}


data_out <- data_out %>% group_by(z) %>% mutate(outind_x = panduanyichang(y))


ggplot(data_out, aes(x = factor(z), y = y)) +
  geom_violin(outlier.shape = NA) +
  geom_point(mapping = aes(colour = factor(outind_x)))+
  scale_color_manual(values = c("1" = "red", "2" = "black"),
                     labels = c("Outlier", "Normal")) 



data <- data.frame(murder = USArrests$Murder,state = tolower(rownames(USArrests))) 
map <- map_data("state") 
k <- ggplot(data, aes(fill = murder))
k + geom_map(aes(map_id = state), map = map)  + expand_limits(x = map$long, y = map$lat) 






library(ggplot2)

# 模拟数据（3组）
df <- data.frame(
  group = c("A", "B", "C"),
  mean = c(10, 15, 12),
  sd = c(2, 1.5, 1)
)

# 用 geom_col() 显示 bar + geom_errorbar()
ggplot(df, aes(x = group, y = mean)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1)




library(ggplot2)

# 1. 生成网格数据
x <- seq(-10,10, length.out = 100)
y <- seq(-10,10, length.out = 100)
eg <- expand.grid(x = x,y = y)
library(mvtnorm)
mu <- c(0, 0)
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # 协方差矩阵
eg$z <- dmvnorm(eg,mean = mu,sigma= sigma)
# 2. 画出等高线图
ggplot(eg, aes(x = x, y = y, z = z)) +
  geom_contour(color = "blue")

data.frame(eg) %>% slice(c(1,4))
table3 <- tibble(
  name = c("Tom", "Lisa"),
  rate = c("5/6/7", "3/2")
)









# 安装（如未安装）
# install.packages(c("survival", "survminer", "broom", "dplyr"))

library(survival)
library(survminer)  # 若装不上，可看后面的“无 survminer 方案”
library(dplyr)
library(broom)

data(lung)

# 构造：时间 = time；事件 = status==2；分组 = sex
# （在 ADTTE 中通常对应：AVAL, CNSR(0/1), TRTxxA）
df <- lung |>
  mutate(
    event = ifelse(status == 2, 1, 0),
    trt   = factor(sex, labels = c("Male", "Female"))  # 仅作演示
  ) |>
  select(time, event, trt, age, ph.ecog)
df <- na.omit(df)
#缺失

library(mice)
imp <- mice(df, m = 5, method = "pmm", seed = 123)
df_filled <- complete(imp)

#异常
data_3 <- df %>% filter(ph.ecog != 3)
table(df$ph.ecog)
select_trated <- surv_fit(Surv(time,event)~ ph.ecog, data = data_3)

summary(is.na(data_3))%>% fill(trt,age, .direction = "down")

ggsurvplot(select_trated,df,  conf.int = TRUE,
           conf.int.style = "ribbon",   # 带状CI
           pval = TRUE)
select_trated <- coxph(Surv(time,event)~ ph.ecog, data = data_3)
ggforest(select_trated, data = df)



















n <- 100
x <- rnorm(n)
beta0_true <- 2;beta1_true <- 0; sigma_true <- 1
y <- beta0_true + beta1_true * x + rnorm(n, sd = sigma_true)
d <- list(x = x, y = y)

## 先验（与你一致；注意 spike 的“点质量”几乎采不到，见注释）
log_prior <- function(pars){
  beta0 <- pars[1]
  beta1 <- pars[2]
  logmix <- log(0.9 * dnorm(beta1, 0, 0.0001) + 0.1 * dnorm(beta1, 0, 10))
  dnorm(beta0, 0, 10, log = TRUE) + logmix
}

## 似然

log_lik <- function(pars, data){
  beta0 <- pars[1]; beta1 <- pars[2]; sigma <- pars[3]
  if (sigma <= 0) return(-Inf)
  mu <- beta0 + beta1 * data$x
  sum(dnorm(data$y, mu, sigma, log = TRUE))
}

log_posterior <- function(pars, data){
  sigma <- pars[3]
  if (sigma <= 0) return(-Inf)
  lp_sigma <- dnorm(sigma, 0, 10, log = TRUE)  # 建议改成 half-normal 或对数正态，见说明
  log_lik(pars, data) + log_prior(pars[1:2]) + lp_sigma
}
## 初始值

init <- c(0, 0, 1)

## 采样（注意：pars=，data=，par_names=）

set.seed(42)
mcmc <- Metro_Hastings(
  li_func   = log_posterior,
  pars      = init,
  data      = d
)

apply(mcmc$trace,2,mean)
mcmc$acceptance_rate


















# ----- 数据与依赖 -----
set.seed(1)
n <- 100
X <- rnorm(n)
beta <- 0.5
mu <- 2
sigma <- 0.6
logT <- mu + beta * X + sigma * rnorm(n)
T_y <- exp(logT)
status <- rep(1, n)
dat <- data.frame(time = T_y, status = status, X = X)

# 用到的包
library(MHadaptive)   # Metro_Hastings
library(extraDistr)   # dinvgamma

# 统一成 list，用 t/x/status 命名，便于向量化
d <- list(t = dat$time, x = dat$X, status = dat$status)

# ----- 先验：放宽以免过度收缩 -----
log_prior <- function(pars) {
  mu    <- pars[1]
  beta  <- pars[2]
  sigma <- pars[3]
  if (!is.finite(sigma) || sigma <= 0) return(-Inf)
  
  # 宽松正态先验 + 逆Gamma 在 sigma（注意是 sigma 本身，不是方差）
  lp_mu    <- dnorm(mu,   mean = 0, sd = 5,  log = TRUE)
  lp_beta  <- dnorm(beta, mean = 0, sd = 10, log = TRUE)
  lp_sigma <- extraDistr::dinvgamma(sigma, a = 2, b = 1, log = TRUE)
  lp_mu + lp_beta + lp_sigma
}

# ----- Gumbel(AFT) 对数似然：log T = mu + beta x + sigma * epsilon,  epsilon ~ Gumbel(0,1) -----
loglik_gumbel_aft <- function(par, data) {
  mu    <- par[1]
  beta  <- par[2]
  sigma <- par[3]
  if (!is.finite(sigma) || sigma <= 0) return(-Inf)
  
  e <- (log(data$t) - mu - beta * data$x) / sigma
  
  # 事件
  logf <- -e - exp(-e) - log(sigma * data$t)
  # 删失
  logs <- -exp(-e)
  
  sum(ifelse(data$status == 1, logf, logs))
}

log_posterior <- function(pars, data) {
  loglik_gumbel_aft(pars, data) + log_prior(pars)
}

# ----- 初值与建议协方差（很重要，决定接受率） -----
# 也可以用 survreg() 得到更贴近真值的初值
init <- c(mu = 0, beta = 0, sigma = 1)
prop_cov <- diag(c(0.05, 0.05, 0.02))  # 适中步长，后面可按接受率微调

set.seed(42)
mcmc <- Metro_Hastings(
  li_func    = log_posterior,
  pars       = init,
  data       = d,
  iterations = 200000
  
 
)

colMeans(mcmc$trace)          # 后验均值
mcmc$acceptance_rate          # 接受率（目标区间 0.2 ~ 0.4 左右）



















library(MHadaptive)
library(extraDistr)

set.seed(1)
n <- 100
X <- cbind(
  X1 = rnorm(n),
  X2 = rbinom(n, 1, 0.5),
  X3 = runif(n, -1, 1),
  X4 = rpois(n, 2)
)
p <- ncol(X)

beta_true <- c(0.5, -0.8, 0.3, -0.4)  # 4 维 β 真值
mu_true   <- 2
sigma_true<- 0.6

logT <- as.vector(mu_true + X %*% beta_true + sigma_true * rnorm(n))
T_y  <- exp(logT)
status <- rep(1, n)

dat <- data.frame(time = T_y, status = status, X)
d <- list(t = dat$time, X = as.matrix(dat[, colnames(X)]), status = dat$status)

# ---- 对数似然（Gumbel AFT）----
loglik_gumbel_aft <- function(pars, data){
  p <- ncol(data$X)
  mu    <- pars[1]
  beta  <- pars[2:(1+p)]
  sigma <- pars[length(pars)]
  if (!is.finite(sigma) || sigma <= 0) return(-Inf)
  
  e <- (log(data$t) - mu - as.vector(data$X %*% beta)) / sigma
  logf <- -e - exp(-e) - log(sigma * data$t)  # 事件
  logs <- -exp(-e)                             # 删失
  sum(ifelse(data$status == 1, logf, logs))
}

# ---- 先验 ----
log_prior <- function(pars, p){
  mu    <- pars[1]
  beta  <- pars[2:(1+p)]
  sigma <- pars[length(pars)]
  if (!is.finite(sigma) || sigma <= 0) return(-Inf)
  
  lp_mu    <- dnorm(mu, 0, 5,  log = TRUE)
  lp_beta  <- sum(dnorm(beta, 0, 10, log = TRUE))
  lp_sigma <- extraDistr::dinvgamma(sigma, a = 2, b = 1, log = TRUE)
  lp_mu + lp_beta + lp_sigma
}

log_posterior <- function(pars, data){
  p <- ncol(data$X)
  loglik_gumbel_aft(pars, data) + log_prior(pars, p)
}

# ---- 初值：长度 = 1(μ) + p(β) + 1(σ) ----
init <- c(mu = 0, rep(0, p), sigma = 1)

# 建议协方差（步长）：维度 (p+2)×(p+2)
step_mu    <- 0.05
step_beta  <- rep(0.05, p)
step_sigma <- 0.02
prop_cov <- diag(c(step_mu, step_beta, step_sigma))

par_names <- c("mu", paste0("beta", 1:p), "sigma")

set.seed(42)
mcmc <- Metro_Hastings(
  li_func    = log_posterior,
  pars       = init,
  data       = d
)

colMeans(mcmc$trace)
mcmc$acceptance_rate


#######################cdisc


row_DM_data <- read.csv("/Users/holee/Desktop/EDC_Raw_Datasets_Oncology_100subjects_20251113_121419_wv8oc2/dm_raw_20251113_121419.csv",header = TRUE)
head(row_DM_data,5)
data_body <- row_DM_data%>%mutate(
  STUDYID = study_id,
  DOMAIN  = "DM",
  USUBJID = paste0(str_trim(study_id),"-",str_trim(subject_id)),
  subjectid = subject_id) %>% left_join(RFSTDTC_RFENDTC_data,
                                        by = c("subjectid" = "subject_id")) %>% distinct()


row_DM_data_3 <- read.csv("/Users/holee/Desktop/EDC_Raw_Datasets_Oncology_100subjects_20251113_121419_wv8oc2/rawlb3_raw_20251113_121419.csv",header = TRUE)

head(row_DM_data_3)

RFSTDTC_RFENDTC_data <- row_DM_data_3 %>% 
  mutate(reform_data =as.Date(collection_date, "%Y-%m-%d")) %>% 
  group_by(subject_id) %>% mutate(RFSTDTC = min(reform_data,na.rm = TRUE),
                                     RFENDTC = max(reform_data,na.rm = TRUE))%>% select(RFSTDTC,RFENDTC,subject_id)%>% distinct()






row_DM_dead_data <- read.csv("/Users/holee/Desktop/EDC_Raw_Datasets_Oncology_100subjects_20251113_121419_wv8oc2/rawvs2_raw_20251113_121419.csv",header = TRUE)

row_DM_dead_data <- row_DM_dead_data%>% filter(visit_name == "Week 12") %>% select(subject_id,assessment_date) %>% mutate(
  assessment_date = as.Date(assessment_date,"%Y-%m-%d"))



data_process_finish <- data_body %>% left_join(RFSTDTC_RFENDTC_data,by=c("subjectid"="subject_id"))  %>% 
   left_join(row_DM_dead_data, by = c("subjectid" = "subject_id")) %>%
  mutate(sitid = as.character(str_sub(USUBJID,-3,-1))) %>% mutate(brthdct = as.Date(informed_consent_date,"%Y-%m-%d")) %>%
  mutate(age = as.numeric(as.Date("2050-3-07") - as.Date(brthdct)) %/% 365.25)%>% 
  mutate(AGEU = "YEARS") %>%
  mutate(sex = ifelse(gender == "Male","M","F")) %>%
  mutate(arm= ifelse(str_sub(study_id,0,1) == "S","treatment","Control")) %>%
  mutate(armcd = ifelse(str_sub(study_id,0,1) == "S","TRT","PLA")) %>%
  mutate(actamct = arm,actarm = armcd) %>%
  select(get_namme)
  
head(data_body,1)
 
dim(data_process_finish)
get_namme <- str_trim(unlist(str_split("STUDYID,DOMAIN,USUBJID,subjectid, RFSTDTC.x,RFENDTC.x ,
           RFSTDTC.y ,RFENDTC.y,assessment_date,sitid,brthdct,age,AGEU,sex,arm,armcd,actamct,actarm",",")))
# data.frame(model.matrix(~ race_category - 1, data = data_body)) %>% 
#   mutate(race =  case_when(race_categoryAsian == 1 ~ "Asian",
#   race_categoryBlack.or.African.American == 1 ~ "Black.or.African.American",
#   race_categoryHispanic == 1 ~ "race_categoryHispanic",
#   race_categoryOther == 1 ~ "Other",
#   race_categoryWhite == 1 ~ "White"
#                                                                                      ))
length(get_namme)

countrycode("japan","country.name","iso3c")
library(dplyr)
library(purrr)

# excel_data: 包含 Variable Name 与 Variable Label
# df: 你的 SDTM/ADaM 数据集
label_map <- excel_data$`Variable Label`


excel_data <- data.frame(
  Variable_Name = c(
    "STUDYID", "DOMAIN", "USUBJID",
    "BRTHDTC", "SEX", "AGE",
    "AETERM", "AEDECOD", "AESTDTC", "AEENDTC","subjectid"
  ),
  Variable_Label = c(
    "Study Identifier",
    "Domain Abbreviation",
    "Unique Subject Identifier",
    "Date/Time of Birth",
    "Sex",
    "Age",
    "Reported Term for the Adverse Event",
    "Dictionary-Derived Term",
    "Start Date/Time of Adverse Event",
    "End Date/Time of Adverse Event",
    "butongrendeid"
  )
)


connect_name_label <- function(x){
  colname <- cur_column()
  if(colname %in% excel_data$Variable_Name){
    idxxx <- which(excel_data$Variable_Name == colname)
    attr(x,"label") <- excel_data$Variable_Label[idxxx]
  }
  return(x)
}



connect_name_label<- function(coll){
  col_name <- cur_column()
  if(col_name %in% excel_data$Variable_Name){
    idxx <- which(excel_data$Variable_Name == col_name)
    attr(coll,"label") <- excel_data$Variable_Label[idxx]
  }
  return(coll)
}

data_process_finish <- data_process_finish %>% rename(c(AGE =age))

data_labeled <- data_process_finish %>%
  mutate(across(everything(), connect_name_label))


## 查看是否变量已将添加label
attr(data_labeled$AGE, "label")



df2 <- tibble(
  SEX = c("M","M","F","F","M","F"),
  RACE = c("WHITE","BLACK","WHITE","BLACK","WHITE","BLACK"),
  AGE = c(20,30,25,35,40,32)
)


df2 %>% mutate(niange =  str_c(SEX,",",AGE)) %>% 
  mutate(niange_l = str_split(niange,",")) %>% unnest_wider(niange_l,"niange") %>% select(-niange)

df2 %>%
  group_by(SEX, RACE) %>%
  summarise(
    group_key = cur_group()
  )




data.frame(
sentence_cc = c("This is a simple example to illustrate word splitting function now")
) %>% mutate(chaifen = chaifenjuzi(sentence_cc,2)) %>%
  unnest_wider(chaifen,"") %>% select(- sentence_cc)

chaifenjuzi <- function(sentence_cc,n){
  words <- unlist(str_split(sentence_cc, " "))
  
  total_long <- length(words)
  celling <- ceiling(total_long / n)
  
  split_sentences <- c()  # 初始化结果向量
  
  for(i in 0:(celling-1)){
    start_idx <- i * n + 1
    end_idx <- min((i + 1) * n, total_long)
    # 按单词索引分组
    split_sentences <- c(split_sentences, str_c(words[start_idx:end_idx], collapse = " "))
  }
  
  return(list(split_sentences))
}

as.Date("2021a9a10", "%Ya%ma%d")

datetime_str <- "2021-09-01-13:14:09"

start_zifu <- str_locate_all(datetime_str,"-")[[1]][,1][3]
end_zifu <- str_length(datetime_str)
str_sub(datetime_str,start_zifu + 1,end_zifu)


# 数据表1：患者基本信息
patients <- data.frame(
  ID = c("P001", "P002"),
  NAME = c("张三", "李四")
)

# 数据表2：体重数据
weights <- data.frame(
  ID = c("P001", "P002"),
  WEIGHT = c(70, 80)
)


derive_vars_merged(
  patients,
  weights,
  by_vars = exprs(ID),
  new_vars = exprs(WEIGHT)
)





df <- tibble(
  USUBJID = c("01", "02"),
  AESTDTC = c("2020-06-01", "2020-07-01")
)

df2 <- derive_vars_dt(
  df,
  dtc = AESTDTC,
  new_vars_prefix = "AEST"
)



library(admiral)
library(tidyverse)

adae <- tibble(
  USUBJID = c("01", "02"),
  TRTSDT  = as.Date(c("2020-06-01", "2020-06-10")),
  AESTDTC = c("2020-06-05", "2020-06-12"),
  AEENDTC = c("2020-06-08", "2020-06-15")
)

adae %>% select(names(adae)[which(names(adae) %like% "DT")])
adae %>%
  derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AEST") %>%  #AESTDT
  derive_vars_dt(dtc = AEENDTC, new_vars_prefix = "AEEN") %>%  #AEENDT
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(AESTDT, AEENDT)
  )
#  reference_date = TRTSDT 裸变量名：就是用这个列。
# source_vars = exprs(AESTDT, AEENDT) exprs()：把一个或者多个列名打包成一个列的列表以方便进行进一步操作
# 
# A tibble: 2 × 8
# USUBJID TRTSDT     AESTDTC    AEENDTC    AESTDT     AEENDT     AESTDY AEENDY
# <chr>   <date>     <chr>      <chr>      <date>     <date>      <dbl>  <dbl>
#   1 01      2020-06-01 2020-06-05 2020-06-08 2020-06-05 2020-06-08      5      8
# 2 02      2020-06-10 2020-06-12 2020-06-15 2020-06-12 2020-06-15      3      6
# AESTDY、AEENDY 这种变量名不是你自己随便取的；
# 而是 ADaM 标准已经规定好的，你必须按照官方名称生成。








dm <- tibble(
  USUBJID = c("01", "02")
)

vs <- tibble(
  USUBJ_ID   = c("01","01","02"),
  VSDTC     = as.Date(c("2021-01-01","2021-02-01","2021-01-10")),
  VSTESTCD  = "WEIGHT",
  VSSTRESN  = c(60, 62, 55)
)

adsl <- dm %>% derive_vars_merged(
  vs,
  by_vars     = exprs(USUBJID=USUBJ_ID),
  order       = exprs(VSDTC),
  mode        = "last",                    # 每人最后一次
  new_vars    = exprs(VSSTRESN)     # 新变量名 BLWT
)

# by_vars = exprs(
#   USUBJID  = USUBJ_ID,
#   VISITNUM = VISIT,
#   DY       = STUDYDAY
# )
 


ae <- tibble(
  USUBJID = c("01","01"),
  AEDT    =  as.Date(c("2021-02-05","2021-03-10"))
)

ex <- tibble(
  USUBJID = c("01","01","01"),
  EXSTDT  = as.Date(c("2021-01-10","2021-02-01","2021-03-01"))
)

ae2 <- ae %>%
  derive_vars_joined(
    dataset_add = ex,
    by_vars     = exprs(USUBJID),
    join_type   = "before",                 # ⭐ 必须明确写
    order       = exprs(EXSTDT),
    mode        = "last",                 # 取最近一次剂量
    filter_join = EXSTDT <= AEDT          # 条件 join
  )

library(admiral)
library(dplyr)
library(tidyr)

vs <- tibble(
  USUBJID  = c("01","01","01","02","02"),
  VSTESTCD = c("SYSBP","DIABP","PULSE","SYSBP","DIABP"),
  VSSTRESN = c(120, 80, 70, 110, 75)
)

vs_tr <- derive_vars_transposed(
  dataset    = vs,
  by_vars    = exprs(USUBJID),
  key_var    = VSTESTCD,        # pivot_wider 的 names_from
  value_var  = VSSTRESN        # pivot_wider 的 values_from
)

vs_wide <- derive_vars_transposed(
  dataset   = vs_long,
  by_vars   = exprs(USUBJID, VISIT), # by_vars 接受表达式列表，这里正确
  key_var   = VSTESTCD,              # 错误修正：移除 exprs()
  value_var = VSSTRESN               # 错误修正：移除 exprs()
)

library(dplyr)
library(admiral)
library(tibble)

# 1. 你的原始长格式数据
vs_long <- tribble(
  ~USUBJID, ~VISIT,     ~VSTESTCD, ~VSSTRESN,
  "001",    "SCREENING", "SYSBP",   120,
  "001",    "SCREENING", "DIABP",   80,
  "001",    "SCREENING", "PULSE",   70,
  "001",    "WEEK 1",    "SYSBP",   122,
  "001",    "WEEK 1",    "DIABP",   81,
  "001",    "WEEK 1",    "PULSE",   72,
  "002",    "SCREENING", "SYSBP",   130,
  "002",    "SCREENING", "DIABP",   85,
  "002",    "SCREENING", "PULSE",   75
)

# 2. (关键步骤) 创建你的“基础”数据集
#    它只包含你希望保留的唯一“键”。
base_data <- vs_long %>%
  distinct(USUBJID, VISIT)

print("--- 基础数据集 (base_data): ---")
print(base_data)

# 3. 现在正确调用 derive_vars_transposed
vs_wide <- derive_vars_transposed(
  dataset       = base_data,            # 基础 (你希望添加列的数据集)
  dataset_merge = vs_long,              # 源 (你要转置的数据集)
  by_vars       = exprs(USUBJID, VISIT), # 用于合并的键
  key_var       = VSTESTCD,            
  value_var     = VSSTRESN              
)

print("--- 转换后 (宽格式): ---")
print(vs_wide)



library(admiral)
library(tibble)
library(dplyr)
library(rlang) # for exprs()

# 基础数据集 (ADSL) - 每位受试者一行
adsl <- tribble(
  ~USUBJID, ~AGE, ~SEX,
  "001",    65,   "M",
  "002",    70,   "F"
)

# 源数据集 (VS) - 长格式
vs <- tribble(
  ~USUBJID, ~VSTESTCD, ~VSSTRESN,
  "001",    "SYSBP",   120,
  "001",    "PULSE",   70,        # 001 的脉搏
  "001",    "SYSBP",   122,
  "001",    "PULSE",   74,        # 001 的脉搏
  "002",    "SYSBP",   130,
  "002",    "PULSE",   80,        # 002 的脉搏
  "002",    "SYSBP",   135,
  "002",    "PULSE",   82         # 002 的脉搏
)
adsl_with_avgpulse <- derive_vars_merged_summary(
  dataset       = adsl,                   # 1. 基础数据集
  dataset_merge = vs,                     # 2. 源数据集
  by_vars       = exprs(USUBJID),         # 3. 合并的键
  
  new_var       = AVGPULSE,               # 4. 要创建的新列的名称
  
  # --- 汇总的逻辑 ---
  analysis_var  = VSSTRESN,               # 5. 要汇总的数值列
  summary_fun   = mean,                   # 6. 要使用的汇总函数
  
  # 7. (关键) 如何在 "源" 数据集中筛选出我们想要的数据？
  filter_merge  = VSTESTCD == "PULSE"
)
# group_by() + summarise()（dplyr）然后如果需要「所有 visit 都有同一个 summary」，相当于 summarise() 结果再 left_join()




adae <- tibble(
  USUBJID = c("01","01"),
  TRTSDT  = as.Date("2024-01-01"),
  AESTDT  = as.Date(c("2023-12-31","2024-01-10"))
)

adae2 <- adae %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars    = exprs(AESTDT)
  )

adae2   # 自动生成 AESTDY


library(dplyr)
library(admiral)   # 里面有 derive_var_extreme_flag()

# 原始 AE 数据 ---------------------------------------------------------
adae <- tribble(
  ~USUBJID, ~ADY, ~AESEQ, ~AETERM,
  "001",        10,     2, "HEADACHE",
  "001",         5,     1, "COUGH",
  "001",        30,     3, "NAUSEA",
  "002",         1,     1, "FEVER",
  "002",        20,     2, "RASH",
  "002",        20,     3, "DIZZINESS",
  "003",         8,     1, "SORE THROAT"
)

# 按照你说的逻辑：
# (1) 按 USUBJID 分组
# (2) 先按 ADY 升序，再按 AESEQ 升序
# (3) 新建标志变量 AFIRSTFL
# (4) mode = "first" → 选“最早那一条”
adae_first <- adae %>% 
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID),
    order   = exprs(ADY, AESEQ),
    new_var = AFIRSTFL,
    mode    = "first"
  )


df2 <- df %>% 
  group_by(grp) %>%
  summarise(mean_x = mean(x))
df2%>% mutate(flag = row_number())


df2 %>% group_by(grp) %>% mutate(flag = row_number())



dat$TRT <- relevel(factor(dat$TRT), ref = "Placebo")
fit_intercept <- lmer(
  SBP ~ week + TRT + (1 | USUBJID),
  data = dat
)

summary(fit_intercept)

dat$USUBJID

set.seed(123)

n <- 80

d

dat <- tibble(
  USUBJID = 1:n,
  TRT = rep(c("Drug", "Placebo"), each = n/2),
  BASE = c(rnorm(n/2, 150, 10), rnorm(n/2, 148, 10)),
  # treatment effect: Drug improves more
  AVAL = BASE - ifelse(TRT=="Drug", rnorm(n/2, 15, 4), rnorm(n/2, 8, 4))
)


fit <- lm(AVAL ~ BASE + TRT, data = dat)
summary(fit)
emm <- emmeans(fit, ~ TRT)
emm
contrast(emm, list(Placebo_minus_Drug = c(-1, 1)))










which("White Cell Count" %like% str_trim(lb_map$LBTEST))
lb_map


#############LB

library(tibble)

lb_map <- tribble(
  ~LBCAT,           ~LBTEST,          ~LBTESTCD, ~LBORRES,     ~FACTOR,
  "Urine Chemistry","Creatinine",     "CREAT",   "mg/dL",      0.0884,
  "Hematology",     "Red Cell Count", "RBC",     "x10E6/uL",       NA,
  "Chemistry",      "Sodium",         "SODIUM",  "mmol/L",         1,
  "Chemistry",      "Sodium",         "SODIUM",  "mmol/L",         1
)
lb_map$LBCAT <- str_to_upper(lb_map$LBCAT)
LB_data_row <- read.csv("/Users/holee/Downloads/sdtm_lb_20251118_093920.csv")

dim(LB_data_row)
  
dddd <- LB_data_row %>% group_by(USUBJID) %>% mutate(LBSEQ_2 =  rank(LBSTNRHI,ties.method = "first", na.last = "keep")) %>% 
  left_join(lb_map,by = c("LBCAT" = "LBCAT")) %>% 
  mutate(LBSTRLO_1 = LBSTNRLO*FACTOR) %>%
  mutate(LBORRES_1 = as.character(LBORRES.y))%>%
  mutate(if_else(is.na(LBSTRLO_1),"NOT DONE",NA)) %>%
  ungroup() %>% 
  mutate(strat_time = as.Date(LBDTC) - as.Date("2025-01-08")) %>%
  group_by(USUBJID) %>% 
  mutate(LBSEQ_2_rank = rank(strat_time,ties.method = "first", na.last = "keep"))%>%
  filter(LBSEQ_2_rank == 1)%>%
  ungroup() %>% select(USUBJID,LBSTRESC)



cccc <- LB_data_row %>% group_by(USUBJID) %>% mutate(LBSEQ_2 =  rank(LBSTNRHI,ties.method = "first", na.last = "keep")) %>% 
  left_join(lb_map,by = c("LBCAT" = "LBCAT")) %>% 
  mutate(LBSTRLO_1 = LBSTNRLO*FACTOR) %>%
  mutate(LBORRES_1 = as.character(LBORRES.y))%>%
  mutate(if_else(is.na(LBSTRLO_1),"NOT DONE",NA)) %>%
  ungroup() %>%
  left_join(dddd,by = c("USUBJID"="USUBJID")) %>%
  mutate(jizhun = if_else(LBSTRESC.y == LBSTRESC.x,"Y",NA))
  
bbbb <- cccc %>% arrange(LBSTRLO_1,LBSTNRHI) %>%
  group_by(USUBJID) %>%
  mutate(LBSEQ_33 = row_number())


#######################ADam
vss <- read.csv("/Users/holee/Desktop/SDTM/sdtm_vs_20251119_055522.csv")
head(vss,20)


vss %>% derive_vars_merged(
  dataset_add = vss,
  by_vars = exprs(USUBJID,VSTESTCD),
  new_var = exprs(weak_1_time = VSDTC),
  filter_add = VISIT == "WEEK 1"
)


vss_1 <- vss %>% mutate(AVAL = log(VSSTRESN))

vss %>% derive_var_base(
  by_var = exprs(USUBJID,VSTESTCD),
  source_var = VSDTC,
  new_var = WEEK_1_TIME,
  filter = VISIT == "WEEK 1"
)%>%
  derive_vars_dt(dtc = WEEK_1_TIME, new_vars_prefix = "WEEK1") %>%  # WEEK1.DT
  derive_vars_dt(dtc = VSDTC, new_vars_prefix = "VS") %>%           # VS.DT
  derive_vars_duration(
    new_var = DAYS_SINCE_WEEK_1,
    start_date = WEEK1DT,      # 用 Date 类型变量
    end_date   = VSDT,         # 用 Date 类型变量
    out_unit = "years"
  )

vss %>% filter(VISIT == "SCREENING" & VSTEST == "Height") %>% select(USUBJID,VSTEST) %>%
  derive_vars_merged(
    dataset_add =  vss,
    by_vars = exprs(USUBJID),
    order = exprs(VSSTRESN),
    mode = "first",       # 取 VSSTRESN 最大那条
    new_vars = exprs(HEIGHT = VSSTRESN)
  )


vss %>%
  derive_var_base(
    by_vars = exprs(USUBJID,VSTESTCD),
    source_var = VSSTRESN,
    new_var = BASE_VSSTRESN,
    filter = VISIT == "BASELINE"
  ) %>%  
  filter(VSSTRESN > 100) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID,VSTESTCD),
    order  = exprs(VSSTRESN),
    mode = "first",
    new_var = MAX_FLAG
  )

VSSS_ind <- vss %>%
  derive_var_base(
    by_vars = exprs(USUBJID, VSTESTCD),
    source_var = VSSTRESN,
    new_var = BASE_VSSTRESN,
    filter = VISIT == "BASELINE"
  ) %>%  
  filter(VISIT != "BASELINE" & VISIT != "SCREENING") %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, VSTESTCD),
    order  = exprs(VSSTRESN),
    mode = "first",
    new_var = MAX_FLAG
  )



vss %>% derive_param_computed(
  by_vars = exprs(USUBJID,VSTESTCD),
  
)


vss_2 <- vss_1 %>%
  derive_var_base(
    by_vars    = exprs(USUBJID, VSTESTCD), 
    source_var = AVAL,
    new_var    = BASE,
    filter     = VSBLFL == "Y"    # 哪些行算 baseline
  ) %>%
  mutate(
    CHG  = AVAL - BASE
  )
vs_bmi <- vs %>%
  derive_param_computed(
    by_vars = exprs(USUBJID, VISIT, VSDTC),
    parameters = c("HEIGHT", "WEIGHT"),
    analysis_value = VSSTRESN,   # 用 VSSTRESN 计算
    compute_expr = VSSTRESC.WEIGHT / (VSSTRESC.HEIGHT/100)^2,
    set_values_to = exprs(
      VSTESTCD = "BMI",
      VSTEST   = "Body Mass Index",
      VSSTRESN = NA_real_,      # 下面自动赋值
      VSSTRESU = "kg/m^2"
    )
  )

uuuuu1 <- vss %>% mutate(PARAMCD = VSTESTCD) %>%  
  derive_param_computed(
  by_vars = exprs(USUBJID, VISIT, VSDTC),
  parameters = c("HEIGHT", "WEIGHT"),
  set_values_to = exprs(
    VSSTRESN =(log(VSSTRESN.HEIGHT) + 2 * VSSTRESN.WEIGHT) / 2,
    VSTESTCD = "aaaaa"
  ),
  keep_nas = TRUE
) %>% filter(USUBJID == "STUDY001-SUBJ0001")

advs <- tribble(
  ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
  "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
  "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
  "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
  "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
  "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
  "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
  "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
  "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
) %>%
  mutate(
    AVALU = "mmHg",
    ADT = case_when(
      VISIT == "BASELINE" ~ as.Date("2024-01-10"),
      VISIT == "WEEK 2" ~ as.Date("2024-01-24")
    ),
    ADTF = NA_character_
  )



library(tibble)
library(dplyr)

# Example 1: Derive BMIBL
adsl <- tribble(
  ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
  "PILOT01", "01-1302",   61,   "YEARS",
  "PILOT01", "17-1344",   64,   "YEARS"
)

advs <- tribble(
  ~STUDYID,  ~USUBJID,  ~PARAMCD, ~PARAM,        ~VISIT,      ~AVAL, ~AVALU, ~ABLFL,
  "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm",   "Y",
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE",   82.1, "kg",   "Y",
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2",    81.19, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4",    82.56, "kg",   NA,
  "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6",    80.74, "kg",   NA,
  "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm",   "Y",
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE",  58.06, "kg",   "Y",
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2",    58.97, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4",    57.97, "kg",   NA,
  "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6",    58.97, "kg",   NA
)

# 1. 从 advs 筛选需要的参数行（HEIGHT, WEIGHT）
# 2. 将长格式转为宽格式（pivot_wider）
# 3. left_join 到 adsl
# 4. 计算新变量（BMIBL）
adsl%>%derive_vars_computed(
  dataset_add = advs,
  by_vars = exprs(STUDYID, USUBJID), #的连接键（join keys）
  parameters = c("WEIGHT", "HEIGHT"),
  new_vars = exprs(BMIBL = AVAL.HEIGHT - AVAL.WEIGHT),
  filter_add = ABLFL == "Y")

derive_vars_computed(
  dataset = adsl,
  dataset_add = advs,
  by_vars = exprs(STUDYID, USUBJID),
  parameters = c("WEIGHT", "HEIGHT"),
  new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)))

advs %>% derive_param_computed(
  by_vars = exprs(USUBJID,VISIT),
  parameters = c("DIABP","SYSBP"),
  set_values_to = exprs(
    AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 2,
    PARAM = "Mean Arterial Pressure (mmHg)",
    PARAMCD = "MAP",
    AVALU = "mmHg",
    ADTF = NA
  ),
  keep_nas = TRUE
)



head(vss_2)

vss_aval <- vss_2 %>%
  mutate(
    PARAMCD = VSTESTCD
  ) %>%
  select(-CHG)

vss_chg <- vss_2 %>%
  mutate(
    PARAMCD = paste0(VSTESTCD, "_CHG"),
    AVAL    = CHG
  ) %>%
  select(-CHG)

vss_final <- bind_rows(vss_aval, vss_chg) %>%
  arrange(USUBJID)



vss_final %>% select(-BASE) ### delete unused variable


library(admiral)
library(dplyr)
library(lubridate)

# 模拟数据
adae <- tibble(
  USUBJID = c("01", "01", "01", "02", "02"),
  AESTDT  = as.Date(c("2024-06-01", "2024-06-05", "2024-06-10", "2024-06-03", "2024-06-08")),
  TRTSDT  = as.Date(c("2024-06-01", "2024-06-01", "2024-06-01", "2024-06-04", "2024-06-04"))
)



adae %>%
  derive_vars_event(
    by_vars   = exprs(USUBJID),
    condition = AESTDT >= TRTSDT,
    order     = exprs(AESTDT),
    new_var   = TRTEMFL
  )




derive_param_computed(
  advs,
  by_vars = exprs(USUBJID, VISIT),
  parameters = "WEIGHT",
  set_values_to = exprs(
    AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
    PARAMCD = "BMI",
    PARAM = "Body Mass Index (kg/m^2)",
    AVALU = "kg/m^2"
  ),
  constant_parameters = c("HEIGHT"),
  constant_by_vars = exprs(USUBJID)
)



adsl <- tribble(
  ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
  "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
  "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""),
  "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""),
  "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""),
  "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
)

lb <- tribble(
  ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
  "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
  "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
  "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
  "PILOT01",    "LB", "01-1133",    304, "2013-05-01T10:13",
  "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
  "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
  "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
  "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
  "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
  "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
) %>%
  mutate(
    ADT = convert_dtc_to_dt(LBDTC)
  )

adsl %>% derive_vars_extreme_event(
  by_vars = exprs(STUDYID, USUBJID),
  events = list(
    # 事件1：死亡
    event(
      dataset_name = "adsl",
      set_values_to = exprs(LSTALVDT = DTHDT,DTHFL = "Y"),
      condition = !is.na(DTHDT),
          # 最后生存日期 = 死亡日期            # 死亡标志
    )),
  source_datasets = list(adsl = adsl),
  tmp_event_nr_var = event_nr,
  order = exprs(LSTALVDT, event_nr),
  mode = "last",
  new_vars = exprs(LSTALVDT, DTHFL)
)



library(admiral)
library(dplyr)

# 创建示例数据
data <- tibble(
  USUBJID = c("001", "002", "003", "004", "005"),
  AGE = c(25, 45, 60, 70, 35),
  BMI = c(18.5, 24.5, 28.0, 32.5, 22.0)
)

# 使用 derive_vars_cat 创建年龄分类变量
result <- data %>%
  derive_vars_cat(
    definition = exprs(
      AGEGR1 = case_when(
        AGE < 18 ~ "< 18",
        AGE >= 18 & AGE < 65 ~ "18-64",
        AGE >= 65 ~ ">= 65"
      ),
      AGEGR1N = case_when(
        AGE < 18 ~ 1,
        AGE >= 18 & AGE < 65 ~ 2,
        AGE >= 65 ~ 3
      )
    )
  )




advs <- tibble::tribble(
  ~USUBJID,       ~VSTEST,  ~AVAL,
  "01-701-1015", "Height", 147.32,
  "01-701-1015", "Weight",  53.98,
  "01-701-1023", "Height", 162.56,
  "01-701-1023", "Weight",     NA,
  "01-701-1028", "Height",     NA,
  "01-701-1028", "Weight",     NA,
  "01-701-1033", "Height", 175.26,
  "01-701-1033", "Weight",  88.45
)

definition <- exprs(
  ~condition,                        ~AVALCAT1, ~AVALCA1N,  ~NEWCOL,
  VSTEST == "Height" & AVAL > 160,   ">160 cm",         1, "extra1",
  VSTEST == "Height" & AVAL <= 160, "<=160 cm",         2, "extra2"
)
derive_vars_cat(
  dataset = advs,
  definition = definition
)


library(admiral)
library(dplyr)

# 假设实验室数据 adlb
adlb <- tibble::tribble(
  ~USUBJID,     ~PARAM, ~AVAL, ~AVALU,  ~ANRHI,
  "01-701-1015", "ALT",   150,  "U/L",      40,
  "01-701-1023", "ALT",    70,  "U/L",      40,
  "01-701-1036", "ALT",   130,  "U/L",      40,
  "01-701-1048", "ALT",    30,  "U/L",      40,
  "01-701-1015", "AST",    50,  "U/L",      35
)

definition_mcrit <- exprs(
  ~PARAM,                      ~condition,    ~MCRIT1ML, ~MCRIT1MN,
  "ALT",                    AVAL <= 10,    "<=ANRHI",         1,
  "ALT", 10 < AVAL & AVAL <= 39, ">1-3*ANRHI",         2,
  "ALT",                 40 < AVAL,   ">3*ANRHI",         3
)
adlb %>%
  derive_vars_cat(
    definition = definition_mcrit,
    by_vars = exprs(PARAM)
  )





