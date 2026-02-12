library(dplyr)

df_1 <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/AALSDXFX.csv"))
meta <- read.csv(paste0(here(), "/proc/meta.csv"))


df_1_proc <- df_1 %>%
  select(-Form_Name, -Visit_Date, -Visit_Name, -alsdx2, -alsdx3, -alsdxdt, -blbelmn, -lleelmn, -lueelmn, -rleelmn,
         -rueelmn, -trnkclmn, -trnkcumn, -trnkelmn, -elescrlr) %>%
  mutate(across(alsdx1:ruecumn, ~case_when(
    .x == 90 ~ NA_integer_,
    .x == 2 ~ 0,
    TRUE ~ .x
  )) )


df_2 <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/AALSHXFX.csv"))

df_2_proc <- df_2 %>%
  select(-Form_Name,-SubjectUID, -Visit_Date, -Visit_Name, -alsdxloc, -diagdt, -hxotsp, -onsetdt)




US_vocab <- read.csv(paste0(here(), "raw/als-metadata/answerals/clinical/RUCA-codes-2020-zipcode.csv"))
df_3 <- read.csv(paste0(here(),"raw/als-metadata/answerals/clinical/Environmental_Questionnaire.csv"))


US_vocab <- US_vocab %>%
  mutate(joinder = paste(POName, State, sep=",") ) %>%
  select(joinder,PrimaryRUCA) %>%
  filter(!duplicated(joinder))


df_3_proc <- df_3 %>%
  filter(!duplicated(SubjectUID)) %>%
  mutate(concussrb = if_else(is.na(concussrb), 2, concussrb ),
         concusstb = as.numeric(if_else(concusstb == "", "0", concusstb)),
         driavgtb = as.numeric(if_else(driavgtb == "","0", driavgtb)),
         drinktb = as.numeric(if_else(drinktb == "", "0", drinktb)),
         concus = concussrb * concusstb) %>%
  select(Participant_ID,smokerb,concus,driavgtb,drinktb,edrb,headrb,milirb,teck1,teck2,teck3, city1, state1)


corrections_cities <- tibble::tribble(
  ~wrong, ~correct,
  "Rockyface", "Rocky Face",
  "Malborough", "Marlborough",
  "St. Charles", "Saint Charles",
  "Gadston", "Gadsden",
  "Ft. Worth", "Fort Worth",
  "Tuscon", "Tucson",
  "Lake St. Louis", "Lake Saint Louis",
  "St. Louis", "Saint Louis",
  "Wilmington,", "Wilmington",
  "WILMINGTON", "Wilmington",
  "Dallax", "Dallas",
  "Dixs", "Dix",
  "St. Charles", "Saint Charles",
  "Queens", "Queens Village",
  "New London Township", "New London",
  "Lagrange", "LaGrange",
  "Washington and New Orleans", "Washington",
  "Seatle", "Seattle",
  "Los Feliz", "Los Angeles",
  "Liberty Township", "Cincinnati",
  "Landsdale", "Lansdale",
  "plano", "Plano",
  "Bellfontaine", "	Bellefontaine",
  "Selingsgrove", "Selinsgrove",
  "Roland Heights", "Rowland Heights",
  "poplar bluff", "Poplar Bluff",
  "Talahasee", "Tallahassee",
  'St. Ann', "Saint Ann",
  "Galipelise", "Gallipolis",
  "JAMESTOWN", "Jamestown",
  "ARLINGTON", "Arlington",
  "BROOKVILLE", "Brookville",
  "Jacksons' Gap", "Jacksons Gap",
  "WILTON", "Wilton",
  "Bluerock", "Blue Rock",
  "St. Augustine", 'Saint Augustine',
  "Hokomis", "Nokomis",
  "Mt. Vernon", "Mount Vernon",
  "Grotton", "Groton",
  "Burchwood", "Birchwood",
  "Charlotsville", "Charlottesville",
  "pittsfield", "Pittsfield",
  "Overland", "Saint Louis",
  "O'Fallon", "Saint Louis",
  "O'fallon", "Saint Louis",
  "Colombus", "Columbus",
  "cottekill", "Cottekill",
  "GERMANTOWN", "Germantown",
  "Shaker Heights", "Cleveland",
  "Mt. Olive", "Mount Olive",
  "VanBuren", "Van Buren",
  "Weldon Spring", "Saint Louis",
  "Nashvile", "Nashville",
  "St. Leonard", "Saint Leonard",
  "St. Peters", "Saint Peters",
  "Creve Couer", "Saint Louis",
  "Prairie Du Sac", "Prairie du Sac",
  "Agora Hills", "Agoura Hills",
  "Downington", "Downingtown",
  "Fredrick", "Frederick",
  "MARTINS FERRY", "Martins Ferry",
  "Mt. Airy", "Mount Airy",
  "winfield", "Winfield",
  "St. Thomas", "Saint Thomas",
  "irving", "Irving",
  "St. James", "Saint James",
  "Skameateles", "Skaneateles",
  "Staceyville", 'Stacyville',
  "Mcdonough", "McDonough",
  "S. Glastonbury", "South Glastonbury",
  "Ruxton", "Baltimore",
  "Oak Leaf", "Dallas",
  "Los Angles", "Los Angeles",
  "University city", "Saint Louis",
  "CAMBRIDGE", "Cambridge",
  "New York City", "New York"
)


corrections_states <- tibble::tribble(
  ~wrong, ~correct,
  "Arizona", "AZ",
  "Arkansas", "AR",
  "Ca", "CA",
  "California", "CA",
  "Connecticut", "CT",
  "DC, LA", "DC",
  "District of Columbia", "DC",
  "Florida", "FL",
  "Il", "IL",
  "Illinois", "IL",
  "Indiana", "IN",
  "Iowa", "IA",
  "Kentucky", "KY",
  "Maine", "ME",
  "Maryland", "MD",
  "Md", "MD",
  "Massachusetts", "MA",
  "Michigan", "MI",
  "Missouri", "MO",
  "Mo", "MO",
  "New Hampshire", "NH",
  "New Hapshire", "NH",
  "New York", "NY",
  "North Carolina", "NC",
  "Oh", "OH",
  "Ohio", "OH",
  "Oklahoma", "OK",
  "Pennsylvania", "PA",
  "Rhode Island", "RI",
  "Texas", "TX",
  "Vermont", "VT",
  "Virginia", "VA",
  "Wyoming", "WY"
)
corr_vec_cities <- setNames(corrections_cities$correct, corrections_cities$wrong)
df_3_proc$city <- dplyr::recode(df_3_proc$city1, !!!corr_vec_cities)

corr_vec_states <- setNames(corrections_states$correct, corrections_states$wrong)
df_3_proc$state <- dplyr::recode(df_3_proc$state1, !!!corr_vec_states)

df_stat_3 <- df_3_proc %>%
  mutate(joinder = paste(city, state, sep=",") ) %>%
  inner_join(US_vocab, by="joinder") %>%
  select(-city, -city1, -state, -state1, -joinder)
  





df_mut <- read.csv(paste0(here(), "/proc/ALS_Gene_Mutations.csv"))



mut_df_proc <- df_mut %>%
  select(-mutotsp) %>%
  mutate(across(ang:vcp, ~ifelse(.x == 2, NA, .x) ),
         general_mut = as.integer(rowSums(across(ang:vcp), na.rm = TRUE) > 0)) %>%
  select(Participant_ID, general_mut)



df_vital <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/Vital_Signs.csv"))


df_vital_proc <- df_vital %>%
  filter(Visit_Date == 0) %>% 
  select(Participant_ID, bmi, bpdias, bpsys) %>% 
  mutate(across(bmi:bpsys, ~as.numeric(.x) ))


df_history <- read.csv(paste0(here(), "/raw/als-metadata/answerals/clinical/Family_History_Log.csv"))

df_history_proc <- df_history %>% 
  select(Participant_ID, fhals, fhalz, fhdem, fhftd, fhpd) %>%
  mutate(history_tot = as.integer(rowSums(across(fhals:fhpd), na.rm = TRUE)>1) ) %>%
  filter(!duplicated(Participant_ID))


df_general <- df_1_proc %>%
  left_join(df_2_proc, by="Participant_ID") %>%
  left_join(df_stat_3, by="Participant_ID") %>%
  left_join(mut_df_proc, by="Participant_ID") %>%
  left_join(df_vital_proc, by="Participant_ID") %>%
  left_join(df_history_proc, by="Participant_ID")

write.csv(df_general, paste0(here(), "/proc/clin_data.csv"), row.names = F)

