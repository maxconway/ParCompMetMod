## Make figures

# Load libraries
library(sybil)
library(sybilSBML)
library(plyr)
devtools::load_all('../GDMOr/')

# Load models
mp <- system.file(package = "sybil", "extdata")
Ec_core <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")
targetFluxes <- c('Biomass_Ecoli_core_w_GAM', 'EX_ac(e)')

## GDMO
# Prepare starting point
start <- as.list(rep(TRUE, times = length(allGenes(Ec_core))))
names(start) <- allGenes(Ec_core)
ancestors <- list(list(genotype = start))

# Prepare evaluation function
evaluate <- function(genotype){
	# avoid an error with no knockouts
	# sybil should really deal with this
	genes <- if(suppressWarnings(all(genotype))){
		NULL
	} else{
		names(genotype)[genotype==FALSE]
	}
	# call sybil to compute the fluxes
	solution <- optimizeProb(object = Ec_core, 
													 gene = genes,
													 lb = 0,
													 ub = 0,
													 retOptSol = FALSE)
	phenotype <- c(solution$fluxes[match(targetFluxes,react_id(Ec_core))], sum(genotype==FALSE))
}

# Run optimization (if necessary)
# reslist <- GDMO(200, 500, ancestors, evaluate)

# Load result
reslist <- load(file='checkpoint.RData')

# Reshape result
resdf <- ldply(.data = reslist, .progress = 'text', .fun = function(individual){
	data.frame(
		genotype = individual$genotype,
		phenotype = {names(individual$phenotype) <- c(targetFluxes, 'kos'); as.list(individual$phenotype)},
		dom = individual$front,
		crowding = individual$crowding
	)
})

strains <- list(Ec_core = list(model = Ec_core, dataset = resdf)
)

# Create graphics
l_ply(.data = strains,
			.progress = 'text',
			.fun = function(x){
				
				png(paste0('./figures/heatmap_',x$model@mod_desc,'.png'))
				heatmapify(x$dataset)
				dev.off()
				
				png(paste0('./figures/interestingByPos_',x$model@mod_desc,'.png'),height=400,width=800)
				outlyingGenes(x$dataset)
				dev.off()
			}
)
