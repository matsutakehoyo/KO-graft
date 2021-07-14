library(tidyverse)
library(readxl)
source("take.R")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path))
	return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

files <- list.files(pattern = "mERG_peaks.csv", recursive=TRUE)
data <- list()
for (i in seq_along(files)){
	# i=1
	print(files[i])

	data[[i]] <- read_csv(files[i]) %>%
		add_column(mouse = split_path(files[i])[3], 
		           graft = split_path(files[i])[4]) %>%
		separate(graft, into=c("graft", "after_transplant"), sep="_")
}
data <- bind_rows(data)

data <- data %>% 
	separate(file, into=c("mouse", "stimuli", "condition"), sep="_") %>%
	separate(mouse, into=c("mouse", "eye"), sep=-1)
	


data <- data %>% 
	mutate(condition=toupper(condition)) %>% #force upper case to avoid errors
	mutate(condition=factor(condition, levels = c("PRE", "AMES", "L-AP4", "WASHOUT"))) %>%
	mutate(condition = forcats::fct_recode(condition, 
                                         "Pre" = "PRE", 
                                         "AMES" = "AMES", 
                                         "L-AP4" = "L-AP4",
                                         "Washout" = "WASHOUT")) %>%
	mutate(stimuli=factor(stimuli, 
	       levels =c("allND2mAmERG", "allND10mAmERG"), 
	       # labels=c("10.56 log photons/cm2/s", "12.84 log photons/cm2/s"))) %>%
	       labels=c("weak", "strong"))) %>%
	mutate(graft = factor(graft, levels=c("WT", "B4", "ISL1"), labels=c("wt", "Bhlhb4", "Islet1")),
	       after_transplant = factor(after_transplant, levels=c("5W", "8W", "12W")))

data %>% distinct(condition) %>% pull()
data %>% distinct(stimuli) %>% pull()
data %>% distinct(graft) %>% pull()
data %>% distinct(after_transplant) %>% pull()

data <- data %>% add_column(graft_covered = NA)
ch_file <- list.files(pattern="[Ee]lectrodes")
if (length(ch_file)){
	ch_data <- read_excel(ch_file) %>%
	separate(`Mouse ID`, into = c("mouse", "eye"), sep=-1)

	for (m in seq_along(ch_data$mouse)){
		mouse_id <- ch_data$mouse[m] 
		eye_id <-  ch_data$eye[m]
		ch_id <- ch_data[m,] %>% select(X__2:X__65) %>%  as.integer()
		
		data <- data %>%	
			mutate(graft_covered = ifelse(mouse==mouse_id & eye==eye_id, channel %in% ch_id, graft_covered)) 
	}
}
data

saveRDS(data, "merg_data.rds")

df1 <- data %>%
	filter(b_time < 65) %>%
	# filter(condition=="L-AP4") 
	filter(condition!="Pre") %>%
	# filter(!is.na(b_amplitude)) %>%
	filter(b_amplitude>0) %>%
	filter(after_transplant!="5W") %>%
	filter(stimuli=="12.84 log photons/cm2/s")



p1 <- ggplot(df1, aes(x=graft_covered, y=b_amplitude, color=graft)) +
	geom_violin() +
	geom_jitter(alpha=1/4)	+
	facet_grid(graft~condition) +
	ylim(0, 0.07) +
	labs(y="b-wave amp (mV)", x="graft") +
	scale_x_discrete(labels=c("FALSE" = "off graft", "TRUE" = "on graft")) +
	scale_color_manual(values=cbPalette) +
	basic_style
ggsave("b wave summary1.pdf", width=5, height=5)
ggsave("b wave summary1.png", width=5, height=5)

p2 <- ggplot(df1, aes(x=condition, y=b_amplitude, color=graft)) +
	geom_violin() +
	geom_jitter(alpha=1/4)	+
	facet_grid(after_transplant+graft~.) +
	ylim(0, 0.07) +
	labs(y="b-wave amp (mV)", x="condition") +
	scale_color_manual(values=cbPalette) +
	basic_style + keynote
ggsave("b wave summary2.pdf", width=16, height=12)
ggsave("b wave summary2.png", width=16, height=12)

summary1 <- df1 %>%
	group_by(graft, condition, after_transplant, mouse, eye) %>%
	summarise(b_mean=mean(b_amplitude, na.rm =TRUE), 
	          graft_area=sum(graft_covered))

p3 <- ggplot(summary1, aes(x=graft_area, y=b_mean, color=graft)) +
	geom_point()+
	scale_color_manual(values=cbPalette) +
	facet_grid(condition) +
	labs(y="b-wave amp mean (mV)", x="graft size (n channel)") +
	basic_style + keynote
ggsave("b wave summary3 (array mean).pdf", width=16, height=12)
ggsave("b wave summary3 (array mean).png", width=16, height=12)

p4 <- ggplot(summary1, aes(x=condition, y=b_mean/graft_area, color=graft)) +
	geom_violin() +
	geom_jitter()+
	labs(y="b-wave amp mean / graft area") +
	scale_color_manual(values=cbPalette) +
	facet_grid(graft~.) +
	basic_style + keynote
ggsave("b wave summary4 (array mean).pdf", width=16, height=12)
ggsave("b wave summary4 (array mean).png", width=16, height=12)

df2 <- data %>%
	filter(b_time < 65) %>%
	# filter(condition=="L-AP4") 
	filter(condition!="Pre") %>%
	filter(after_transplant!="5W")
	

p5 <- ggplot(df2, aes(x=graft_covered, y=b_amplitude, color=graft)) +
	geom_violin() +
	geom_jitter(alpha=1/4) +
	facet_grid(graft~stimuli) +
	labs(y="b-wave amp (mV)") +
	scale_x_discrete(labels=c("FALSE" = "off graft", "TRUE" = "on graft")) +
	scale_color_manual(values=cbPalette) +
	basic_style + keynote
ggsave("b wave summary1 (Pre).pdf", width=16, height=12)
ggsave("b wave summary1 (Pre).png", width=16, height=12)


p6 <- ggplot(df2, aes(x=stimuli, y=b_amplitude, color=graft)) +
	geom_violin() +
	geom_jitter(alpha=1/4) +
	facet_grid(graft~.) +
	labs(y="b-wave amp (mV)") +
	scale_color_manual(values=cbPalette) +
	basic_style + keynote
ggsave("b wave summary2 (Pre).pdf", width=16, height=12)
ggsave("b wave summary2 (Pre).png", width=16, height=12)

ggplot(data, aes(x=a_time)) +
	geom_histogram(binwidth=0.001) +
	facet_grid(graft~graft_covered) +
	xlim(-0.02, 0.02)

# data <- data %>%
# 	mutate(a_abs_ongraft = a_bas * graft_covered, 
# 	       a_abs_offgraft = a_bas * (1-graft_covered),
# 	       b_abs_ongraft = b_bas * graft_covered, 
# 	       b_abs_offgraft = b_bas * (1-graft_covered),
# 	       a_time_ongraft = a_time * graft_covered, 
# 	       a_time_offgraft = a_time * (1-graft_covered),
# 	       b_time_ongraft = b_time * graft_covered, 
# 	       b_time_offgraft = b_time * (1-graft_covered))
# 	group_by(graft, stimuli, condition) %>%
# 	summarise(a_abs = mean(a_abs),
# 	          b_abs = mean(b_abs),
# 	          a_time = mean(a_time),
# 	          b_time = mean(b_time))





