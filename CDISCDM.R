
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


