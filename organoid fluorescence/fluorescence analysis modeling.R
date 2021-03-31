library(tidyverse)
source('take.R')
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('matsu_stan.R')

gompertz <-function(A, mu, lambda, time){
	# modified Gompertz growth curve
	A * exp(-exp(mu * exp(1) / A * (lambda - time) + 1))
}



load(file = "channel_data.Rda")

df <-bind_rows(df_red, df_green, df_blue)

# df %>% select(line, dd, filename, exposure) %>% print(n=Inf)

# get background intensity from blue channel
ggplot(df_blue, aes(x=dd, y=mean_signal)) +
	geom_point() + 
	basic_style

ggplot(df_blue, aes(x=dd, y=mean_bg)) +
	geom_point() + 
	basic_style

ggplot(df_blue, aes(x=dd, y=mean_signal-mean_bg)) +
	geom_point() + 
	basic_style

df <- df %>% 
	group_by(channel) %>%
	mutate(bg = df_blue$mean_signal - df_blue$mean_bg)

glimpse(df)

ggplot(df, aes(x=dd, y=n_signal)) +
	geom_point() + 
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=mean_signal-mean_bg)) +
	geom_point() + 
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=(mean_signal-mean_bg) / bg -1 )) +
	geom_point() +
	facet_grid(.~channel) +
	basic_style

ggplot(df %>% filter(channel=="green"), aes(x=dd, y=(mean_signal-mean_bg) / bg - 1)) +
	geom_point() +
	facet_grid(.~line) +
	basic_style

# ggplot(df, aes(x=mode_bg)) +
# 	geom_histogram(binwidth = 1/255) +
# 	facet_grid(.~channel) +
# 	basic_style

ggplot(df, aes(x=dd, y=mean_signal-mean_bg)) +
	geom_point(binwidth = 1/255) +
	facet_grid(.~channel) +
	basic_style

ggplot(df, aes(x=dd, y=(mode_signal)/mode_bg)) +
	geom_point() +
	facet_grid(.~channel) +
	geom_text(aes(label=filename))


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


# modeling data
df_model <- df_green %>% 
	filter(dd<41) %>%
	mutate(rel_signal = (mean_signal-mean_bg) / bg,
			 #rel_signal = mean_signal,
			 line = factor(line, levels = c("wt", "Bhlhb4", "Islet1")), 
			 time_id = group_indices(., dd))

p <- ggplot() +
	geom_point(aes(x=dd, y=rel_signal, color=line), data=df_model) +
	facet_grid(line~.) +
	basic_style

df_model %>% group_by(line) %>%
	summarise(n=n())

model_data <- 
	list(n_obs = nrow(df_model), 
		  n_line = df_model$line %>% unique() %>% length(), 
		  x = df_model$dd,
		  y = df_model$rel_signal,
		  obs2line = df_model$line %>% as.numeric())

# nested model
model_name <- "growth_curve_nested.stan"
model <- 
	stan(file=model_name,
		  data = model_data,
		  cores = getOption("mc.cores", 1L),
		  control = list(adapt_delta = 0.9),
		  chains = 8,
		  thin=1, iter=2000)
model_name <- gsub(".stan", "", model_name)

get_elapsed_time(model)
check_div(model)
check_treedepth(model,max_depth=10)
check_energy(model)
check_n_eff(model)
check_rhat(model)


# diagnostics
posterior_summary <- mcmc_summary(model)
saveRDS( posterior_summary, paste0(model_name, " posterior summary.rds") )
# posterior_summary <- readRDS(paste0(model_name, " posterior summary.rds"))


posterior <- tidy_posterior(model)
saveRDS( posterior, paste0(model_name, " posterior tidy.rds") )
# posterior <- readRDS( paste0(model_name, " posterior tidy.rds") )


model_diagnose(posterior, posterior_summary, name = model_name)
mcmc_pairs(model, model_name, parameters=c("a0", "a_l", "m0", "m_l", "l0", "l_l", "rate"))

# posteriors
model_posterior(posterior, name = model_name)

# reviewer requested "conventioal inteval" hdi -> .95
pretty_posterior(
	posterior %>% filter(Parameter=="rate"),
	model_name,
	c("rate"), splatoon, 
	width=5, height=5, hdi1=.95, text_size=13
)

pretty_posterior(
	posterior %>% filter(Parameter=="a0"),
	model_name,
	c("overall"), splatoon, 
	width=5, height=5, hdi1=.95, text_size=13
) 
 

pretty_posterior(
	posterior %>% filter(Parameter=="a_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="a_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,  
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="m0"),
	model_name,
	c("overall"), splatoon, 
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="m_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="m_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,  
	width=5, height=5, hdi1=.95, text_size=13

) 

pretty_posterior(
	posterior %>% filter(Parameter=="l0"),
	model_name,
	c("overall"), splatoon, 
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="l_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=5, height=5, hdi1=.95, text_size=13
) 

pretty_posterior(
	posterior %>% filter(Parameter=="l_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,
	width=5, height=5, hdi1=.95, text_size=13
) 



# posterior predictive check

# prediction for line
post_pred1 <- posterior %>% filter(Parameter %in% c("a_l", "m_l", "l_l")) %>% 
	spread(Parameter, value) %>%
	left_join( posterior %>% filter(Parameter=="rate") %>% 
	          	select(-ind, -ind1, -Parameter) %>% rename( rate=value)  ) %>%
	# filter(Iteration==1) %>%
	mutate( dd = list(seq(0, 40, 1)) ) %>% unnest() %>%
	rowwise() %>%
	mutate( sim = gompertz(a_l, m_l, l_l, dd)) %>%
	unnest() %>%
	group_by(ind, dd) %>%
	summarise(mean = mean(sim), 
	          hdi_lower = hdi(sim, 0.95)[1], 
	          hdi_upper = hdi(sim, 0.95)[2]) %>%
	ungroup() %>%
	mutate(line = factor(ind, labels=c("wt", "Bhlhb4", "Islet1")))

# prediction for observation
post_pred2 <- posterior %>% filter(Parameter %in% c("a_l", "m_l", "l_l")) %>% 
	spread(Parameter, value) %>%
	left_join( posterior %>% filter(Parameter=="rate") %>% 
	          	select(-ind, -ind1, -Parameter) %>% rename( rate=value)  ) %>%
	# filter(Iteration==1) %>%
	mutate( dd = list(seq(0, 40, 1)) ) %>% unnest() %>%
	rowwise() %>%
	mutate( sim = rgamma(100, gompertz(a_l, m_l, l_l, dd)*rate, rate) %>% list() ) %>%
	unnest() %>%
	group_by(ind, dd) %>%
	summarise(mean = mean(sim), 
	          hdi_lower = hdi(sim, 0.95)[1], 
	          hdi_upper = hdi(sim, 0.95)[2]) %>%
	ungroup() %>%
	mutate(line = factor(ind, labels=c("wt", "Bhlhb4", "Islet1")))

p_pred_check  <- p + 
	geom_ribbon( aes(x=dd, ymin=hdi_lower, ymax=hdi_upper, fill=as.factor(line)), 
                 post_pred1, alpha=1/3) +
	geom_ribbon( aes(x=dd, ymin=hdi_lower, ymax=hdi_upper, fill=as.factor(line)), 
                 post_pred2, alpha=1/3) +	
	facet_grid(line~.) +
	scale_fill_manual(values=cbPalette) + 
	scale_color_manual(values=cbPalette) + 
	theme(legend.position="none") +
	scale_y_continuous(breaks=scales::pretty_breaks(n=4), expand = c(0, 0), limits=c(0,19) )+
	scale_x_continuous(breaks=scales::pretty_breaks(n=4), expand = c(0, 0), limits=c(5,38)) +
	labs(y="Fluorescence (AU)", x="DD (day)") +
	theme(text = element_text(size=12)) 
ggsave(paste0(model_name, " posterior/posterior predictive check.pdf"), width=2, height=6)

pred_hdi <- ggplot() + 
	geom_ribbon(aes(x=dd, ymin=hdi_lower, ymax=hdi_upper, fill=as.factor(line)), 
                post_pred1, alpha=1/3) + 
	scale_fill_manual(values=cbPalette) + basic_style

pred_mean <- ggplot() + 
	geom_line(aes(x=dd, y=mean, color=as.factor(line)), post_pred1) +
	scale_color_manual(values=cbPalette) + 
	basic_style	

ggsave(paste0(model_name, " posterior/posterior predictive check faceted.pdf"), 
       gridExtra::grid.arrange(p, pred_hdi, pred_mean, ncol=1),
       width=4, height=8)


# crossed model
model_name <- "growth_curve_crossed.stan"
model <- 
	stan(file=model_name,
		  data = model_data,
		  cores = getOption("mc.cores", 1L),
		  control = list(adapt_delta = 0.99, max_treedepth = 10),
		  chains = 8,
		  iter=2000)
model_name <- gsub(".stan", "", model_name)

get_elapsed_time(model)
check_div(model)
check_treedepth(model,max_depth=10)
check_energy(model)
check_n_eff(model)
check_rhat(model)

posterior_summary <- mcmc_summary(model)
saveRDS( posterior_summary, paste0(model_name, " posterior summary.rds") )
# posterior_summary <- readRDS(paste0(model_name, " posterior summary.rds"))

posterior <- tidy_posterior(model)
saveRDS( posterior, paste0(model_name, " posterior tidy.rds") )
# posterior <- readRDS( paste0(model_name, " posterior tidy.rds") )

# diagnostics
model_diagnose(posterior, posterior_summary, name = model_name)
mcmc_pairs(model, model_name, parameters=c("a0", "a_l", "m0", "m_l", "l0", "l_l", "rate"))
mcmc_pairs(model, model_name, parameters=c("A0", "A_l", "M0", "M_l", "L0", "L_l", "rate"), "sum-to-zero")

# posteriors
model_posterior(posterior, name = model_name)

pretty_posterior(
	posterior %>% filter(Parameter=="A0"),
	model_name,
	c("overall"), splatoon, 
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="A_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="A_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,  width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="M0"),
	model_name,
	c("overall"), splatoon, 
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="M_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="M_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,  width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="L0"),
	model_name,
	c("overall"), splatoon, 
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="L_l"),
	model_name,
	c("wt", "Bhlhb4", "Islet1"),
	width=2.5, height=2.5
) 

pretty_posterior(
	posterior %>% filter(Parameter=="L_l") %>% posterior_diff(),
	model_name,
	c("wt", "Bhlhb4", "Islet1") %>% make_diff_label(),
	wes_palette("Moonrise3"),
	diff=TRUE,  width=2.5, height=2.5
) 


# posterior predictive check
post_pred <- posterior %>% filter(Parameter %in% c("n_a", "n_m", "n_l")) %>% 
	spread(Parameter, value) %>%
	left_join( posterior %>% filter(Parameter=="rate") %>% 
	          	select(-ind, -ind1, -Parameter) %>% rename( rate=value)  ) %>%
	# filter(Iteration==1) %>%
	mutate( n_a=exp(n_a), n_m=exp(n_m), n_l=exp(n_l) ) %>% 
	mutate( dd = list(seq(0, 40, 1)) ) %>% unnest() %>%
	rowwise() %>%
	mutate( sim = rgamma(100, gompertz(n_a, n_m, n_l, dd)*rate, rate) %>% list() ) %>%
	unnest() %>%
	group_by(ind, dd) %>%
	summarise(mean = mean(sim), 
	          hdi_lower = hdi(sim, 0.89)[1], 
	          hdi_upper = hdi(sim, 0.89)[2]) %>%
	ungroup() %>%
	mutate(line = factor(ind, labels=c("wt", "Bhlhb4", "Islet1")))

p_pred_check <- p + geom_ribbon(aes(x=dd, ymin=hdi_lower, ymax=hdi_upper, fill=as.factor(line)), 
                post_pred, alpha=1/3) +
	geom_line(aes(x=dd, y=mean, color=as.factor(line)), post_pred) +
	facet_grid(line~.) +
	scale_fill_manual(values=cbPalette) + 
	scale_color_manual(values=cbPalette) + 
	theme(legend.position="none") +
	scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
	scale_x_continuous(breaks=scales::pretty_breaks(n=4)) +
	xlim(5, 38) + ylim(0, 19)
ggsave(paste0(model_name, " posterior/posterior predictive check.pdf"), width=2, height=3)

pred_hdi <- ggplot() + 
	geom_ribbon(aes(x=dd, ymin=hdi_lower, ymax=hdi_upper, fill=as.factor(line)), 
                post_pred, alpha=1/3) + 
	scale_fill_manual(values=cbPalette) + basic_style

pred_mean <- ggplot() + 
	geom_line(aes(x=dd, y=mean, color=as.factor(line)), post_pred) +
	scale_color_manual(values=cbPalette) + 
	basic_style	

ggsave(paste0(model_name, " posterior/posterior predictive check faceted.pdf"), 
       gridExtra::grid.arrange(p, pred_hdi, pred_mean, ncol=1),
       width=4, height=8)



