library(tidyverse)
source('take.R')
library(gridExtra)
library(scales)
library(GO.db)
options(tibble.print_max = 30, tibble.print_min = 20)
options(tibble.width = 80)
# for foreach parallel processing
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())
getDoParWorkers()
# for map parallel proccessing
library(furrr)
plan(multiprocess)
library(ggrastr)

# find the probes for a particular gene
findgene <- function (keyword){
	# look in gene neme
	genes_name <- grep(keyword, data$GeneName,
							 ignore.case=TRUE, perl=TRUE, value=TRUE) %>% unique()
	
	# look in description keyword in Description
	description <- grep(keyword, data$Description,
							  ignore.case=TRUE, perl=TRUE, value=TRUE) %>% unique()
	# genes that have keyword in description 
	genes_desc <- data %>% 
		filter(Description %in% description) %>% 
		distinct(GeneName) %>% 
		pull() %>% as.character()
	
	# look in probe
	genes_prob <- data %>% 
		filter(ProbeName == keyword) %>%
		distinct(GeneName) %>% 
		pull() %>% as.character()
	
	genes <- c(genes_name, genes_desc, genes_prob)
	
	print(paste(c("Genes:", genes), collapse="  ", sep="\t"))
	df <- data %>% 
		filter(GeneName %in% genes) %>%
		distinct(GeneName, ProbeName, Description)
	if (length(df)){
		return(df)
	} else {
		print(paste("NO Gene Found"))
	}
}

#find gene from its name
# findgene("OPN4")
#find genes usign regex 
# findgene("^opn.*")
#retina related genes
#findgene("retina")
#find by probe name 
#findgene("A_52_P542570")


plot_expression <- function(probe){
	# probe = "A_55_P2061879"
	probe_data <- data %>% filter(ProbeName %in% probe)
	ggplot() +
		geom_bar(aes(x=dd, y=Normalized, fill=line), data=probe_data, 
					stat="identity", position="dodge") +
		scale_fill_manual(values=cbPalette) +
		scale_y_continuous(breaks= pretty_breaks()) +
		theme(axis.title.x=element_blank(), 
				axis.title.y=element_blank()) +
		ggtitle(label=paste(probe_data$GeneName %>% unique()) ) +
		facet_wrap(~GeneName) +
		basic_style + theme(legend.position="none") 
	
}

# the data
load("data.Rda")
data <- as_tibble(data)

length(unique(data$ProbeName))
length(unique(data$GeneName))

# some data wrangling
data <- data %>%
	filter(ControlType==FALSE) %>% #remove controls
	filter(gIsWellAboveBG==TRUE) %>% #remove low signals
	# filter(Row>50) %>%
	#discard columns bigining with gIs (flags for quality control), and Raw values
	# select(grep("Raw", grep("^gIs", colnames(data), value=TRUE, invert=TRUE), 
	#             value=TRUE, invert=TRUE))
	dplyr::rename(line = Cell.Line, dd = dif.day) %>% 
	mutate(line = fct_recode(line, Islet1 = "Isl1"), 
			 GO = as.character(GO)) 

levels(data$line)
levels(data$dd)

# some times the same probe is present multiple times 
data <- data %>% 
	group_by(dd, line, ProbeName, GeneName, GO, Description) %>%
	summarise(
		Normalized = mean(Normalized, na.rm = TRUE), 
		n_probe = n()
	) %>%
	ungroup()

data$GeneName %>% unique %>% length
data$ProbeName %>% unique %>% length

glimpse(data)

plot_expression(findgene("opsin") %>% pull(ProbeName))
plot_expression(findgene("OPN4") %>% pull(ProbeName))
findgene("BDNF") %>% pull(ProbeName) %>% plot_expression()

findgene("spns2") %>% pull(ProbeName) %>% plot_expression()
findgene("s1p")

findgene("S1pr3") %>% pull(ProbeName) %>% plot_expression()
findgene("chat") %>% pull(ProbeName) %>% plot_expression()
findgene("tyrosine hy") %>% pull(ProbeName) %>% plot_expression()
findgene("chat") %>% pull(Description)

#genes of interest Mandai sensei list
# mandai_list <- read.csv("Genes of interest (ProbeName).csv", header=TRUE)$ProbeName
# findgene(as.character(mandai_list))
# mandai_genes <- list()
# for (i in seq_along(mandai_list)){
# 	mandai_genes[[i]] <- findgene(as.character(mandai_list[i]))	
# }
# mandai_genes <- bind_rows(mandai_genes)
# mandai_genes <- mandai_genes %>% distinct(GeneName)
# write.csv(mandai_genes, file = "mandai_genes.csv")

# import gene list
library(readxl)
gene_list <- read_excel("../genes.xlsx") %>%
	mutate(group = factor(group))

# check if all the gene names match
genes <- gene_list$gene_name
for (i in seq_along(genes)){
	test <- data %>% filter(GeneName==genes[i])
	if (!nrow(test))
		print(paste(genes[i], "error !!!!"))
	#else print(paste(genes[i], "OK"))
	
}
levels(gene_list$group)

# gene plots by groups
groups <- gene_list %>% filter(!is.na(group)) %>%
	distinct(group) %>% pull() %>% as.character()


for (g in groups){
	# gnees in group
	gene_group <- gene_list %>% filter(group==g) %>%
		distinct(gene_name) %>%
		pull() %>% as.character()
	
	probe_list <- match(gene_group, data$GeneName) %>% 
		data$ProbeName[.]	%>% 
		as.character() 
	
	plots <- list()
	# subset data 
	for (i in seq_along(probe_list)){
		probe_data <- data %>% filter(ProbeName==probe_list[i])
		plots[[i]] <- ggplot() +
			geom_bar(aes(x=dd, y=Normalized, fill=line), data=probe_data, 
						stat="identity", position="dodge") +
			scale_fill_manual(values=cbPalette) +
			scale_y_continuous(breaks= pretty_breaks()) +
			theme(axis.title.x=element_blank(), 
					axis.title.y=element_blank()) +
			ggtitle(label=paste(probe_data$GeneName %>% unique()), 
					  subtitle = paste(probe_list[i])) +
			basic_style + theme(legend.position="none") 
	}
	layout <- matrix(1:30, nrow=6, ncol=5) %>% t
	gene_plots <- marrangeGrob(grobs = plots, layout_matrix=layout)
	pdf(paste(g, "genes.pdf"), width = 9, height = 8)
	print(gene_plots)
	dev.off()
}

glimpse(data)

DD10 <- data %>% 
	spread(key=dd, value=Normalized) %>% 
	gather(key=dd, value=Normalized, DD16:DD23)
glimpse(DD10)


# use ggrastr to plot points as raster as there are too many points for the pdf


compare_dd10 <- ggplot(DD10, aes(x=DD10, y=Normalized, color=line)) +
	geom_point_rast(alpha=1/5, size=1.5, stroke=0) +
	scale_x_continuous(trans=log10_trans()) +
	scale_y_continuous(trans=log10_trans()) +
	# geom_line(data=data.frame(x=c(0.005, 1000), y=c(0.01, 1000)), aes(x=x, y=y), colour="gray50", size=0.3) +
	geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.01*4, 1000*4)), aes(x=x, y=y), colour="gray30", size=0.3, alpha=1/3) +
	geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.01/4, 1000/4)), aes(x=x, y=y), colour="gray30", size=0.3, alpha=1/3) +
	facet_grid(dd~line) +
	scale_color_manual(values = cbPalette) +
	labs(x="DD10") +
	theme(axis.title.y=element_blank()) +
	basic_style + theme(legend.position = "none")
ggsave("Compare DD10.pdf", width=5, height=3)
# ggsave("Compare DD10.png", width=5, height=3, bg = "transparent")

# make scatter plots comparing to wt
wt <- data %>% 
	spread(key=line, value=Normalized) %>% 
	gather(key=line, value=Normalized, Bhlhb4:Islet1)
compare_wt <- ggplot(wt, aes(x=wt, y=Normalized, color=line)) +
	geom_point_rast(alpha=1/5, size=1.5, stroke=0) +
	scale_x_continuous(trans=log10_trans()) +
	scale_y_continuous(trans=log10_trans()) +
	# geom_line(data=data.frame(x=c(0.005, 1000), y=c(0.01, 1000)), aes(x=x, y=y), colour="gray50", size=0.3) +
	geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.01*4, 1000*4)), aes(x=x, y=y), colour="gray30", size=0.3, alpha=1/3) +
	geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.01/4, 1000/4)), aes(x=x, y=y), colour="gray30", size=0.3, alpha=1/3) +
	facet_grid(line~dd) +
	scale_color_manual(values = c("#E69F00", "#56B4E9")) +
	# labs(x="Normalized mRNA of wt", y="Normalized mRNA") +
	labs(x="wt") +
	theme(axis.title.y=element_blank()) +
	basic_style + theme(legend.position = "none")
ggsave("Compare wt.pdf", width=5, height=3)
# ggsave("Compare wt.png", width=5, height=3, bg = "transparent")

####################### testing #######################

#######################Make Tibble for Gene Ontology Analysis ####################################
# go_text <- readLines("go-slim.txt")
go_text <- readLines("go.obo")

# go term always start with [Term] and at the end there are [Typedef]
entries <- grep("\\[(Term|Typedef)\\]", go_text)
# go_hierlim[entries+1] 

pb <- progress_estimated(length(entries)-1)
go <- list()
for (i in 1:(length(entries)-1)){
	pb$tick()$print()
	id <- gsub("id: ", "", go_text[entries[i]+1], perl=TRUE)
	alt_id <- gsub("alt_id: ", "", 
	               grep("alt_id: ", c(go_text[entries[i]:entries[i+1]]),
	               perl=TRUE, value=TRUE), perl=TRUE) %>% list()
	name <- gsub("name: ", "", 
	             grep("name: ", c(go_text[entries[i]:entries[i+1]]),
	               perl=TRUE, value=TRUE))
	namespace <- gsub("namespace: ", "", 
	                  grep("namespace: ", c(go_text[entries[i]:entries[i+1]]),
	                       perl=TRUE, value=TRUE))
	def <- gsub("def: ", "", 
	            grep("def: ", c(go_text[entries[i]:entries[i+1]]),
	                 perl=TRUE, value=TRUE))
	intersection_of <- gsub("intersection_of: ", "", 
	                        grep("intersection_of: ", c(go_text[entries[i]:entries[i+1]]),
	                        	perl=TRUE, value=TRUE) %>% unique()) %>% list() 
	relationship <- gsub("relationship: ", "", 
	                     grep("relationship: ", c(go_text[entries[i]:entries[i+1]]),
	                          perl=TRUE, value=TRUE) %>% unique()) %>% list()
	subset <- gsub("subset: ", "", 
	                     grep("subset: ", c(go_text[entries[i]:entries[i+1]]),
	                          perl=TRUE, value=TRUE) %>% unique()) %>% list()
	is_a <- gsub("(is_a: )(GO:[0-9]{7})( ! .*)", "\\2", 
	             grep("is_a: GO:.*+", c(go_text[entries[i]:entries[i+1]]), 
	                  perl=TRUE, value=TRUE)) %>% 
	             unique() %>% list()

	go[[i]] <- tibble(id=id,
	                  alt_id=alt_id, 
	                  name=name, 
	                  namespace=namespace, 
	                  def=def, 
	                  intersection_of=intersection_of,
	                  relationship=relationship, 
	                  subset=subset, 
	                  is_a=is_a) %>%
	replace_na(list(alt_id = list(NA))) %>%
	replace_na(list(intersection_of = list(NA))) %>%
	replace_na(list(relationship = list(NA))) %>%
	replace_na(list(subset = list(NA))) %>%
	replace_na(list(is_a = list(NA)))
}
go <- bind_rows(go)

go <- go %>% 
	# extract part_of from relationships and intersection
	mutate(part_of1  = ifelse(grepl("part_of ", relationship),
	                         relationship, NA), 
		   part_of2  = ifelse(grepl("part_of ", intersection_of),
	                         intersection_of, NA)) %>%
	rowwise() %>%
	mutate(part_of1 = grep("part_of", part_of1 %>% unlist(), value=TRUE) %>%
					  gsub("(part_of )(GO:[0-9]{7})( ! .*)", "\\2", .) %>% list(),
		   part_of2 = grep("part_of", part_of2 %>% unlist(), value=TRUE) %>%
					  gsub("(part_of )(GO:[0-9]{7})( ! .*)", "\\2", .) %>% list()) %>%
	mutate(part_of = c(part_of1 %>% unlist() %>% na.omit(), 
	                   part_of2 %>% unlist() %>% na.omit()) %>% as.character() %>% list()) %>%
	replace_na(list(part_of = list(NA))) %>%
	ungroup()

# go <- go %>% 
# 	# extract part_of from relationships and intersection
# 	mutate(is_a  = ifelse(grepl("is_a ", is_a),
# 	                         is_a, NA), 
# 	rowwise() %>%
# 	mutate(is_a = grep("is_a", is_a %>% unlist(), value=TRUE) %>%
# 					  gsub("(is_a )(GO:[0-9]{7})( ! .*)", "\\2", .) %>% list()) %>%
# 	replace_na(list(is_a = list(NA))) %>%
# 	ungroup()

# go[70:80,] %>% 
# identify do not annotate subset gocheck_do_not_annotate|gocheck_do_not_manually_annotate
go <- go %>%
	# mutate(do_not_annotate = ifelse(grepl("gocheck_do_not_annotate|gocheck_do_not_manually_annotate", subset), TRUE, FALSE)) 
	mutate(do_not_annotate = ifelse(grepl("gocheck_do_not_annotate", subset), TRUE, FALSE)) 

do_not_annotate <- go %>% filter(do_not_annotate) %>% pull(id) %>% unique()

go %>% filter(do_not_annotate) %>% pull(name)

go <- go %>% 
	unnest(is_a, .drop=FALSE) %>%
	unnest(part_of, .drop=FALSE) %>%
	unnest(alt_id, .drop=FALSE)

glimpse(go)

# saveRDS(go_all, file = "GO_all_data.rds")

space <- c("biological_process", "molecular_function", "cellular_component")

####################### find parental nodes by namespace #######################
# function to look up parental GO in go_ref by id
find_ref <- function(go_id){
	go_ref %>% filter(GO_id==go_id) %>% pull(relationship) %>% unlist()
}
# find_ref("GO:0045211") %>% unlist()

find_namespace <- function(go_id){
	go %>% filter(id==go_id) %>% pull(namespace) %>% unique()
}
# find_namespace("GO:0000991")

find_name <- function(go_id){
	go %>% filter(id==go_id) %>% pull(name) %>% unique()
}

find_level1 <- function(go_id){
	go_hier %>% filter(level2==go_id) %>% pull(level1) %>% unique()
}


go_hier <- list(); k=1
for (s in space){
	print(paste(s))

	go_hier[[k]] <- go %>% 
		# process each namespace
		filter(namespace==s) %>% 
		# exclude namespce GO that are namespace
		filter(!(id %in% c("GO:0003674", "GO:0008150", "GO:0005575"))) %>%
		# exclude GO that lead to namespace GO
		filter(!(is_a %in% c("GO:0003674", "GO:0008150", "GO:0005575"))) %>%
		filter(!(part_of %in% c("GO:0003674", "GO:0008150", "GO:0005575")))  #%>%
		# exclude GO in teh do_not_annotate list
		# filter(!(id %in% do_not_annotate))

	# this is to reference parental GOs
	go_ref <- go_hier[[k]] %>% 
		dplyr::select(id, is_a, part_of) %>%
		# remove orphan GOs 
		filter(!(is.na(is_a) & is.na(part_of))) %>%
		# exclude GO in that lead to do_not_annotate list
		# filter(!(is_a %in% do_not_annotate)) %>%
		# filter(!(part_of %in% do_not_annotate)) %>%
		# exclude GO taht lead to different namespace
		gather(key=rel, value=GOid, is_a:part_of) %>%
		mutate(namespace=future_map(GOid, ~find_namespace(.x)), .progress=TRUE) %>% 
		unnest(namespace) %>%
		filter(namespace==s) %>%
		# summarise all the parental GOs
		group_by(id) %>%
		summarise(relationship = GOid %>% unique() %>% list()) %>% 
		# renema id column to avoid confusion
		dplyr::rename(., GO_id=id)
	
	# go_ref[1,2] %>% unlist()
	go_hier[[k]] <- go_hier[[k]] %>% dplyr::select(id, name, namespace) %>% distinct()

	# this works but its VERY slow, extremely inneficient
	# system.time({	
	# go_hier[[k]] <- go_hier[[k]] %>% 
	# 	rowwise() %>%
	# 	mutate(parent0 = ifelse(id %in% go_ref$GO_id,
	# 	                        go_ref %>% filter(., id==GO_id) %>% pull(relationship) %>% as.list(),
	# 	                        id %>% as.list())) %>% 
	# 	unnest(parent0, .drop=FALSE) %>% distinct()
	# })

	# this is also very slow
	# system.time({	
	# go_hier[[k]] <- go_hier[[k]] %>% 
	# 	mutate(parent0 = 
	# 	       ifelse(id %in% go_ref$GO_id,
	# 	              id %>% map(~find_ref(.x)),
	# 	              id %>% as.list())) %>% 
	# 	unnest(parent0, .drop=FALSE) %>% distinct()
	# })

	system.time({	
	go_hier[[k]] <- go_hier[[k]] %>% 
		mutate(parent0 = 
		       ifelse(id %in% go_ref$GO_id,
		              id %>% future_map(~find_ref(.x), .progress = TRUE),
		              id %>% as.list())) %>% 
		unnest(parent0, .drop=FALSE) %>% distinct()
	})	
	# system.time({	
	# go_hier[[k]] <- go_hier[[k]] %>% 
	# 	mutate(parent0 = 
	# 	       ifelse(id %in% go_ref$GO_id,
	# 	              go_ref %>% filter(GO_id %in% id) %>% pull(relationship),
	# 	              id %>% as.list())) %>% 
	# 	unnest(parent0, .drop=FALSE) %>% distinct()
	# })
	# print(go_hier[[k]], width=Inf)
	# this is a variable to check is final parent node is reached 
	check <- nrow(go_hier[[k]])

	# find parent node
	p <- 1
	while(check>0){
		print(paste("parent node: ", p))
		start <- Sys.time()
		go_hier[[k]] <- go_hier[[k]] %>% 
			mutate(!!paste0("parent", p):= 
			       ifelse(!!rlang::parse_quosure(paste0("parent", p-1)) %in% go_ref$GO_id,
			              !!rlang::parse_quosure(paste0("parent", p-1)) %>% future_map(~find_ref(.x), .progress=TRUE), 
			              !!rlang::parse_quosure(paste0("parent", p-1)) %>% as.list())) %>%
			unnest(!!rlang::parse_quosure(paste0("parent", p)), .drop=FALSE) %>%
			distinct()

		end <- Sys.time()
		print(paste(Sys.time(), "parent node", p, "check: ", check))
		print(end-start)
		check <- as.logical(go_hier[[k]][[paste0("parent", p-1)]]!=go_hier[[k]][[paste0("parent", p)]]) %>% sum()
		p <- p+1
	}
	
	# find the parent node
	nodes <- gsub("(parent)([0-9])", "\\2", colnames(go_hier[[k]][ncol(go_hier[[k]])])) %>% as.numeric()
	go_hier[[k]] <- go_hier[[k]] %>% 
		mutate(level1=!!rlang::parse_quosure(paste0("parent", nodes)),
		       level2=id)
	for (n in 1:nodes){
		go_hier[[k]] <- go_hier[[k]] %>% 
			mutate(level2 = ifelse(!!rlang::parse_quosure(paste0("parent", n-1))==level1,
			                       level2, 
			                       !!rlang::parse_quosure(paste0("parent", n-1))))
	}
	# go_hier[[k]] %>% print(width=Inf)
	k=k+1
}
go_hier <- bind_rows(go_hier)

go_hier[[3]] %>% print(width=Inf)
#go_hier %>% group_by(namespace) %>% nest()
# saveRDS(go_hier, "go_hier.rds")
####################### find parental nodes by namespace #######################
####################### Make Tibble for Gene Ontology Analysis ####################################
# readRDS("go_hier.rds")

find_parent <- function(df, name_space="molecular_function"){
	df <- df %>%
		# separate GO terms
		mutate(GO=strsplit(GO, split="\\|")) %>%
		unnest() %>%
		separate(GO, into=c("GO_id", "GO_description"), sep=10) %>%
		# filter out namespace GOs
		filter(GO_id!="GO:0008150") %>%	# GO:0008150 biological_process
		filter(GO_id!="GO:0003674") %>%	# GO:0003674 molecular_function
		filter(GO_id!="GO:0005575") # GO:0009295 cellular_component


	go_hier_sp <- go_hier %>% 
		filter(namespace==name_space) %>%
		distinct(id, name, namespace, level1, level2)
	
	# correlate GO information releationship (level2)
	df <- df %>%	
		mutate(level2 = ifelse(GO_id %in% go_hier_sp$id, go_hier_sp$level2, NA), 
		       name = ifelse(GO_id %in% go_hier_sp$id, go_hier_sp$name, NA),
		       namespace = ifelse(GO_id %in% go_hier_sp$id, go_hier_sp$namespace, NA)) 

	# Count GOs by gene
	df <- df %>%
		# filter(namespace=="biological_process") %>%
		filter(namespace==name_space) %>%
		# filter(namespace=="cellular_component") %>%
		filter(!is.na(level2)) %>%
		group_by(GeneName) %>% #select(GO_id, level2, level1) %>% print(n=100)
		summarise(level2=unique(level2) %>% list()) %>%
		# unnest(level1, .drop=FALSE) %>% 
		unnest(level2, .drop=FALSE) %>%
		mutate(n_genes = unique(GeneName) %>% length())
	
	# find level1 parent node	
	df <- df %>%	
		mutate(level1 = ifelse(level2 %in% go_hier_sp$id,
		                       level2 %>% future_map(~find_level1(.)),
		                       NA)) %>%
		unnest(level1)
	

	# Calculate level1 freq
	df <- df %>%
		group_by(level2) %>%
		summarise(n=n(), level1=list(unique(level1)), n_genes=unique(n_genes)) %>% 
		unnest(level1) %>%
		mutate(freq = n / sum(n) * 100) %>% 
		arrange(desc(n)) #%>% print(n=100)

	# find the name of the GO terms 
	df <- df %>%	
		mutate(name_level1 = ifelse(level1 %in% go$id, 
		                            level1 %>% future_map(~find_name(.)), NA), 
		       name_level2 = ifelse(level2 %in% go$id, 
		                            level2 %>% future_map(~find_name(.)), NA)) %>%
		unnest(name_level1) %>%
		unnest(name_level2)
	
	return(df)
}

####################### Gene Ontology Analysis #######################

# function to compare conditions and output a list of genes
compare_genes_dd <- function(ref, target_dd, target_line, change="up"){
	genes <- data %>%
	spread(key=dd, value=Normalized) %>%
	mutate(ref = eval(parse(text=ref))) %>%
	gather(key=dd, value=Normalized, DD10:DD23) %>%
	mutate(log2_exp = log2(Normalized/ref)) %>%
	filter(dd==target_dd, line==target_line) %>%
	arrange(desc(log2_exp))
	if (change=="up"){
		genes <- genes %>%
			filter(log2_exp > 2)
	} else if (change=="down"){
		genes <- genes %>%
			filter(log2_exp < -2)
	} else {
		return(NA)
	}	

	write.table(pull(genes, GeneName) %>% as.character(), 
	            paste0(target_line, "_", target_dd, "(", change, ").txt"),
	            quote = FALSE, row.names = FALSE, col.names = FALSE)
	return(genes)
}

compare_genes_line <- function(ref, target_line, target_dd, change="up"){
	genes <- data %>%
		spread(key=line, value=Normalized) %>%
		mutate(ref = eval(parse(text=ref))) %>%
		gather(key=line, value=Normalized, wt:Islet1) %>%
		mutate(log2_exp = log2(Normalized/ref)) %>%
		filter(dd==target_dd, line==target_line) %>%
		arrange(desc(log2_exp))
	if (change=="up"){
		genes <- genes %>%
			filter(log2_exp > 2)
	} else if (change=="down"){
		genes <- genes %>%
			filter(log2_exp < -2)
	} else {
		return(NA)
	}
	
	write.table(pull(genes, GeneName) %>% as.character(), 
	            paste0(target_dd, "_", target_line, "(", change, ").txt"),
	            quote = FALSE, row.names = FALSE, col.names = FALSE)
	return(genes)
}


system.time({
line_data <- foreach  (s=space, .combine=rbind) %:% 
	foreach (r = c("up", "down"), .combine=rbind) %:% 
		foreach (l = c("Bhlhb4", "Islet1"), .combine=rbind) %:% 
			foreach (d=c("DD10", "DD16", "DD23"), .combine=rbind) %dopar% {
				print(paste(s, l, d, r))
				compare_genes_line("wt", l, d, change=r) %>%
					find_parent(name_space=s) %>%
					mutate(line = l, dd=d, namespace=s, regulation=r, i=i)												
			}
		
})
# saveRDS(line_data, "line_data.rds")
line_data

system.time({
dd_data <- foreach  (s=space, .combine=rbind) %:% 
	foreach (r = c("up", "down"), .combine=rbind) %:% 
		foreach (l = c("wt", "Bhlhb4", "Islet1"), .combine=rbind) %:% 
			foreach (d=c("DD16", "DD23"), .combine=rbind) %dopar% {
				print(paste(s, l, d, r))
				compare_genes_dd("DD10", d, l, change=r) %>%
					find_parent(name_space=s) %>%
					mutate(line = l, dd=d, namespace=s, regulation=r, i=i)												
			}
		
})
# saveRDS(dd_data, "dd_data.rds")

# df <- compare_genes_line("wt", "Bhlhb4", d, change=r) %>% find_parent()
line_data <- line_data %>% 
	mutate(line = factor(line, levels=c("wt", "Bhlhb4", "Islet1")), 
	       dd = factor(dd, levels=c("DD10", "DD16", "DD23")))
line_data %>% distinct(line, dd, namespace, regulation)

dd_data <- dd_data %>% 
	mutate(line = factor(line, levels=c("wt", "Bhlhb4", "Islet1")), 
	       dd = factor(dd, levels=c("DD10", "DD16", "DD23")))
dd_data %>% distinct(line, dd, namespace, regulation)


glimpse(dd_data)

# make the plots for changes acroos lines
line_plots <- list(); line_plots_data <- list(); i=1
for (r in c("up", "down")){	
	for (s in space){
		# # top 5 categories in level1
		groups_level1 <- line_data %>% 
			filter(regulation==r) %>%
			filter(namespace==s) %>%
			group_by(name_level1, name_level2) %>%
			summarise(n=sum(freq)) %>%
			arrange(desc(n)) %>%
			pull(name_level1) %>% unique()
		groups_level1 <- groups_level1[1:5]

		# select top 10 
		groups_level2 <- line_data %>% 
			filter(regulation==r) %>%
			filter(namespace==s) %>%
			group_by(name_level1, name_level2) %>%
			summarise(n=sum(freq)) %>%
			arrange(desc(n)) %>%
			pull(name_level2) %>% unique()
		groups_level2 <- groups_level2[1:10]

		# line_data$name_level1 %in% groups_level1
		# line_plots_data[[i]] %>% distinct(name_level1)
		line_plots_data[[i]] <- line_data %>% 
			filter(namespace==s) %>% 
			filter(regulation==r) %>% 
			mutate(name_level2 = ifelse(name_level2 %in% groups_level2, name_level2, "other"), 
			       name_level1 = ifelse(name_level1 %in% groups_level1, name_level1, "other")) %>% 
			mutate(name_level1 = str_wrap(name_level1, width = 20), 
			       name_level2 = str_wrap(name_level2, width = 40)) 
	
		# reorganise level1 to top 5 rank
		line_plots_data[[i]] <- line_plots_data[[i]] %>% 
			mutate(name_level1=factor(name_level1,
			                          levels=c(str_wrap(groups_level1, width=20), "other")),
					name_level2=factor(name_level2, 
					                   levels=c(str_wrap(groups_level2, width=40), "other")) 
			)	
	
		line_plots[[i]] <- ggplot(line_plots_data[[i]], aes(x=name_level1, y=n)) +
			geom_bar(aes(fill=name_level2), stat="identity") +
			geom_text(aes(label=paste(n_genes, "genes")), 
			          data = line_plots_data[[i]] %>% distinct(dd, line, n_genes), 
			          y=Inf, x=Inf, vjust="inward", hjust="inward", 
			          alpha=1/2) +
			labs(x="Category", y="# of GOs") + 
			ggtitle(s) +
			scale_y_continuous(breaks= pretty_breaks()) +
			scale_fill_manual(values=c(rev(hue_pal()(10)), "grey")) +
			theme(axis.text.x=element_text(angle = 60, hjust = 1)) +
			facet_grid(dd~line) +
			basic_style
		ggsave(paste0(s, "GO entrichment (line ", r, ").pdf"), width=10, height=5)
		i=i+1
	}
}

# make the plots for changes acroos dd
dd_plots <- list(); dd_plots_data <- list(); i=1
for (r in c("up", "down")){	
	for (s in space){
		# # top 5 categories in level1
		groups_level1 <- dd_data %>% 
			filter(regulation==r) %>%
			filter(namespace==s) %>%
			group_by(name_level1, name_level2) %>%
			summarise(n=sum(freq)) %>%
			arrange(desc(n)) %>%
			pull(name_level1) %>% unique()
		groups_level1 <- groups_level1[1:5]

		# select top 10 
		groups_level2 <- dd_data %>% 
			filter(regulation==r) %>%
			filter(namespace==s) %>%
			group_by(name_level1, name_level2) %>%
			summarise(n=sum(freq)) %>%
			arrange(desc(n)) %>%
			pull(name_level2) %>% unique()
		groups_level2 <- groups_level2[1:10]

		# dd_data$name_level1 %in% groups_level1
		# dd_plots_data[[i]] %>% distinct(name_level1)
		dd_plots_data[[i]] <- dd_data %>% 
			filter(namespace==s) %>% 
			filter(regulation==r) %>% 
			mutate(name_level2 = ifelse(name_level2 %in% groups_level2, name_level2, "other"), 
			       name_level1 = ifelse(name_level1 %in% groups_level1, name_level1, "other")) %>% 
			mutate(name_level1 = str_wrap(name_level1, width = 20), 
			       name_level2 = str_wrap(name_level2, width = 40)) 
	
		# reorganise level1 to top 5 rank
		dd_plots_data[[i]] <- dd_plots_data[[i]] %>% 
			mutate(name_level1=factor(name_level1,
			                          levels=c(str_wrap(groups_level1, width=20), "other")),
					name_level2=factor(name_level2, 
					                   levels=c(str_wrap(groups_level2, width=40), "other")) 
			)	
	
		dd_plots[[i]] <- ggplot(dd_plots_data[[i]], aes(x=name_level1, y=n)) +
			geom_bar(aes(fill=name_level2), stat="identity") +
			geom_text(aes(label=paste(n_genes, "genes")), 
			          data = dd_plots_data[[i]] %>% distinct(dd, line, n_genes), 
			          y=Inf, x=Inf, vjust="inward", hjust="inward", 
			          alpha=1/2) +
			labs(x="Category", y="# of GOs") + 
			ggtitle(s) +
			scale_y_continuous(breaks= pretty_breaks()) +
			scale_fill_manual(values=c(rev(hue_pal()(10)), "grey")) +
			theme(axis.text.x=element_text(angle = 60, hjust = 1)) +
			facet_grid(line~dd) +
			basic_style
		ggsave(paste0(s, "GO entrichment (dd ", r, ").pdf"), width=10, height=5)
		i=i+1
	}
}
