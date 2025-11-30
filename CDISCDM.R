
dm_data <- read.csv("/Users/holee/Desktop/EDC_Raw_Datasets_Oncology_100subjects_20251113_121419_wv8oc2/dm_raw_20251113_121419.csv")
ae_data <- read.csv("/Users/holee/Desktop/EDC_Raw_Datasets_Oncology_100subjects_20251113_121419_wv8oc2/mh_raw_20251113_121419.csv")
ae_data$examination_date <- as.Date(ae_data$examination_date, format="%Y-%m-%d")
max_min_dataset <- ae_data %>% 
  group_by(subject_id) %>%
  summarise(max_date = max(end_date),
            min_date = min(start_date))

dm_data


dm_data_1 <- dm_data %>% mutate(study_id_new = "STUDY001",
                                subject_id_new = as.character(str_c("SUB", str_pad(1:100,side = "left",width = 4,pad = "0"))))



dm_data_2 <- dm_data_1 %>% left_join(max_min_dataset, by = c("subject_id_new" = "subject_id"))%>%
  rename(rfstdtc = max_date,rfscdtc = min_date)

dead_time <- ae_data %>% 
  group_by(subject_id) %>%
  summarise(max_date = max(as.Date(end_date)))


final_data <- dm_data_2 %>% left_join(dead_time, by = c("subject_id_new" = "subject_id")) %>% 
  mutate(dthfl = if_else(is.na(max_date),"", "Y")) %>%
  mutate(dthdtc =  as.character(max_date)) %>% select(-max_date) %>%
  mutate(site_id = str_sub(study_id,-3,-1)) %>%
  mutate(age_days = difftime(Sys.Date(), as.Date(birth_date), units = "days")) %>%
  mutate(age = as.integer(age_days/365.25)) %>%
  select(-age_days) %>%
  mutate(ageu = "YEARS") %>%
  mutate(sex = case_when(
    gender == "Male" ~ "M",
    gender == "Female" ~ "F",
    TRUE ~ "other")
  ) %>%
  mutate(ethnicity_new = "Not Hispanic or Latino") %>%
  mutate(actam = if_else(tttt == 0,"placedo", "treatment")) %>%
  mutate(armcd = if_else(tttt == 0,str_to_upper("pla"),str_to_upper("tre"))) %>% 
  select(-tttt) %>% 
  select(c(study_id_new,subject_id_new,rfstdtc,rfscdtc,dthfl,dthdtc,site_id,age,ageu,sex,ethnicity_new,actam,armcd))

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



final_data <- final_data %>% mutate(age = as.numeric(age),......)

str(final_data)

output_final_data <- final_data %>% mutate(across(everything(), connect_name_label))

write_xpt(output_final_data,"out_DM.xpt")


### updata from EC2





dataset %>% derive_var_base(
  by_vars = exprs(USUBJID,PARAMCD),
  source_var = AVAL,
  new_var = BASELINE_AVAL
)

vs %>% derive_vars_merged(
  dataset_add =  dm %>% select(STUDYID, USUBJID,AGE,AGEU),
  by_vars = exprs(STUDYID, USUBJID),
  filter_add = AGEU == "YEARS",
  new_vars = exprs(level_age = case_when(
    AGE <= 62 ~ "LOW",
    TRUE ~ "HIGH"
  ))
)

derive_vars_duration(
  strat_date = exprs(TRTEDT),
  end_date   = exprs(LSTALVDT),
  new_var    = exprs(TRTDUR),
  out_unit   = "days"
)

adsl %>% derive_vars_dt(
  dtc = TRTEDT,
  new_vars_prefix = "TRTEDT_T"
) %>% derive_vars_dt(
  dtc = LSTALVDT,
  new_vars_prefix = "LSTALVDT_T"
) %>%
  derive_vars_duration(
    start_date = TRTEDT_TDT,
    end_date   = LSTALVDT_TDT,
    new_var    = TRTDUR,
    out_unit   = "day"
  )




sdtm_ae <- read.csv("/Users/holee/Desktop/SDTM/sdtm_ae_20251119_055449.csv")
sdtm_dm <- read.csv("/Users/holee/Desktop/SDTM/sdtm_dm_20251119_055437.csv")


adsl <- sdtm_dm %>% select(STUDYID,USUBJID,AGE,SEX,RACE)

sdtm_ae$AEACN

add_indix_variable <- sdtm_ae %>% mutate(
  AESER = case_when(
    AESER == "Y" ~ "Yes",
    TRUE ~ "No"),
  AERELFL = case_when(AEREL == "RELATED" ~ "Y",
                      TRUE ~ "N"),
  AESTDISP = case_when(
    AEACN == "DRUG WITHDRAWN" ~ "Y",
    TRUE ~ "N"
  )
  
)

add_indix_exflag_variable <- add_indix_variable %>% derive_vars_dtm(dtc = AESTDTC,new_vars_prefix = "AST")%>%
  derive_vars_dtm(dtc = AEENDTC,new_vars_prefix = "AEN")%>%
  derive_vars_dtm_to_dt (source_vars = exprs (ASTDTM, AENDTM))%>% 
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs (STUDYID, USUBJID))%>%
  derive_var_extreme_flag( ## 给极大值极小值加flag
    by_vars = exprs(STUDYID, USUBJID,AEDECOD),
    order = exprs(ASTDTM),
    mode = "first", 
    new_var = AE_FIRST_FL) 



aaaaadad <- add_indix_exflag_variable %>% 
  derive_vars_merged(
    dataset_add = add_indix_exflag_variable,
    by_vars = exprs(USUBJID,STUDYID,AEDECOD),
    order = exprs(ASTDT),
    mode = "first",
    new_var = exprs(date_porgess = ASTDT)
  ) %>%
  # 第二步：计算时长
  derive_vars_duration(
    new_var = DURATION_days,    # 新变量名，不需要 exprs()
    start_date = date_porgess,        # 开始日期
    end_date = ASTDT,   # 结束日期 (即刚刚生成的第一个日期)
    out_unit = "days"
  )



add_indix_exflag_variable_1 <- add_indix_exflag_variable %>%
  restrict_derivation(
    filter = AEENDY > 200,
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID),
      order = exprs(ASTDT),
      mode = "first",
      new_var = AEENDY_MAXFL,
    )
  )





add_indix_exflag_variable_1 <- add_indix_exflag_variable %>%
  restrict_derivation( #根据 filter 条件创建临时子集。（例如：filter = AVAL > 100）对这个临时子集执行derivation的函数
    #将派生出的新变量结果leftjoin回原始的数据集
    derivation = derive_vars_merged,
    args = params(
      dataset_add = add_indix_exflag_variable,
      by_vars = exprs(USUBJID),
      order = exprs(ASTDT),
      mode = "first",
      new_vars = exprs(DATE_GT_100 = ASTDT)
    )
    ,filter = AEENDY > 200)

add_indix_exflag_variable_2 <- add_indix_exflag_variable %>% group_by(USUBJID) %>% arrange(USUBJID,AENDT) %>%
  mutate(
    PREV_AENDT = lag(AEENDY),
    PRGINDX = case_when(AEENDY > 1.1 * PREV_AENDT ~"PD",TRUE ~ "no PD")
  ) %>%
  ungroup()


pd_event <- event(dataset_name = "add_indix_exflag_variable_2",#event函数生成一个临时表只包含 PRGINDX == "PD" 的行 + 新生成的变量
                  condition = PRGINDX == "PD",
                  set_values_to = exprs(
                    jinzhanshijian = AENDT,
                    jinzhanshijian_lab = "PD"
                  ))



# 2. 合并回主表[对应 derive_vars_extreme_event 的整体行为]
caoleee <- add_indix_exflag_variable_2 %>% 
  derive_vars_extreme_event(
    by_vars = exprs(USUBJID),#Join Key
    events = list(pd_event), 
    source_datasets = list(add_indix_exflag_variable_2 = add_indix_exflag_variable_2), 
    order = exprs(AENDT), 
    mode = "first", 
    new_vars = exprs(
      PD_DATE = jinzhanshijian     # 将 event 里的 jinzhanshijian 拿出来，改名为 PD_DATE
    )
  )


add_indix_exflag_variable_2 %>% 
  derive_extreme_event(
    by_vars = exprs(USUBJID),
    event(dataset_name = "add_indix_exflag_variable_2",
          condition = lag(AESTDY) - AESTDY > 30),
    set_values_to = exprs(
      PD_DATE = AENDT,
      PD_DATE_LAB = "PD"
    ),
    order = exprs(AENDT),
    mode = "first",
    new_vars = exprs(PROG_EVENT, PROG_DATE)
  )


add_indix_exflag_variable_2 %>% group_by()




add_indix_exflag_variable_2






library(dplyr)
library(lubridate)
library(admiral)

# 简化版 ADRS 数据：每次影像评估一条
adrs <- tribble(
  ~USUBJID, ~AVALC, ~ADT,
  "01",     "SD",   ymd("2024-01-10"),
  "01",     "PD",   ymd("2024-03-01"),
  "01",     "PD",   ymd("2024-05-01"),
  
  "02",     "PR",   ymd("2024-01-05"),
  "02",     "SD",   ymd("2024-02-10"),
  
  "03",     "SD",   ymd("2024-02-01"),
  "03",     "PD",   ymd("2024-04-12")
)


add_indix_exflag_variable_3 <- add_indix_exflag_variable_2 %>%
  derive_vars_extreme_event(source_datasets = list(add_indix_exflag_variable_2 = add_indix_exflag_variable_2),
                            #########附表需要的原始表格
                            events = list(
                              event(
                                dataset_name = "add_indix_exflag_variable_2",
                                condition = PRGINDX == "PD",
                                set_values_to = exprs(
                                  EVENT = "PD",
                                  DATE  = AENDTM,
                                  REASON = AESEV
                                )
                              )
                            ),
                            ##################以上建立附表
                            by_vars = exprs(USUBJID),
                            order = exprs(DATE),
                            mode  = "first",
                            ###################取附表中的group排序后取第一个值
                            new_vars = exprs(EVENT, DATE,REASON)
                            #########取附表中的那个变量横着加回去
  )

