library(EBImage)
library(tidyverse)
source('take.R')

split_path <- function(path) {
	if (dirname(path) %in% c(".", path))
		return(basename(path))
	return(c(basename(path), split_path(dirname(path))))
}

extract_img_value <- function(img, color){
	# color = "red"
	df <- imageData(channel(img, color)) %>%
		as.vector %>%
		as.tibble %>%
		mutate(pixel = ifelse(value>otsu(channel(img, color)), "signal", "bg")) %>%
		group_by(pixel) %>%
		summarise(mean = mean(value),
					 sd = sd(value),
					 mode = dens_mode(value), 
					 n = n()) %>%
		unite("mean_sd_mode_n", mean, sd, mode, n) %>%
		spread(key=pixel, value=mean_sd_mode_n) %>%
		separate(bg, c("mean_bg", "sd_bg", "mode_bg", "n_bg"), sep="_", convert = TRUE) %>%
		separate(signal, c("mean_signal", "sd_signal", "mode_signal", "n_signal"), sep="_", convert = TRUE) %>%
		add_column(channel = color)
	
	return(df)
}
# extract_img_value(img, "blue") %>% glimpse

file <- list.files(pattern=".tif$", recursive=TRUE)
meta <- list.files(pattern="_meta.txt$", recursive=TRUE)

pb <- progress_estimated(length(file))
df_red <- list()
df_green <- list()
df_blue <- list()
for (i in seq_along(file)){
	# i <- 100
	pb$tick()$print()
	
	exposure_time <- read_delim(meta[i], delim=":", col_names = FALSE) %>%
		spread(key=X1, value=X2, convert = TRUE) %>% #glimpse()
		select(`Shutter Speed                   `) %>%
		pull() %>%
		parse(text = .) %>% eval()
	
	f <- file[i]
	img = readImage(f)
	# df_img[[i]] <- tibble(file = f, image = list(img), exposure = exposure_time)
	# display(img)
	# display(green)
	# hist(img)
	# hist(green)
	# otsu(green)
	# str(img)
	# green <- channel(img, "green")
	# blue <- channel(img, "blue")
	# red <- channel(img, "red")
	df_red[[i]] <- extract_img_value(img, "red") %>%
		add_column(
			filename = split_path(f)[1],
			line = split_path(f)[3],
			dd = as.numeric(gsub("[^\\d]+", "", split_path(f)[2], perl=TRUE)),
			exposure = exposure_time)
	df_green[[i]] <- extract_img_value(img, "green") %>%
		add_column(
			filename = split_path(f)[1],
			line = split_path(f)[3],
			dd = as.numeric(gsub("[^\\d]+", "", split_path(f)[2], perl=TRUE)),
			exposure = exposure_time)
	df_blue[[i]] <- extract_img_value(img, "blue") %>%
		add_column(
			filename = split_path(f)[1],
			line = split_path(f)[3],
			dd = as.numeric(gsub("[^\\d]+", "", split_path(f)[2], perl=TRUE)),
			exposure = exposure_time)
}

df_red <- bind_rows(df_red)
df_green <- bind_rows(df_green)
df_blue <- bind_rows(df_blue)

save(df_green, df_red, df_blue, file = "channel_data.Rda")

df_green <- df_green %>%
	mutate(bg = df_blue$mean_signal - df_blue$mean_bg)
df_blue <- df_blue %>% 
	mutate(bg = mean_signal - mean_bg)

df <-bind_rows(df_red, df_green, df_blue)
glimpse(df)

ggplot(df, aes(x=mean_bg)) +
	geom_histogram(binwidth = 1/255) +
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=mode_signal)) +
	geom_point(binwidth = 1/255) +
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=(mode_signal)/mode_bg)) +
	geom_point() +
	facet_grid(.~channel) +
	geom_text(aes(label=filename))

ggplot(df, aes(x=dd, y=(mean_signal-mean_bg))) +
	geom_point() +
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=(mode_signal-mode_bg))) +
	geom_point() +
	facet_grid(.~channel) +
	basic_style
	
ggplot(df, aes(x=mean_signal, y=mode_signal)) +
	geom_point() +
	facet_grid(.~line) +
	basic_style

ggplot(df, aes(x=mean_bg, y=mode_bg)) +
	geom_point() +
	facet_grid(.~line) +
	basic_style

ggplot(df_green, aes(x=dd, y=(mean_signal-mean_bg)/bg)) +
	geom_point() +
	basic_style

