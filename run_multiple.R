# run_multiple

# models <- llply(.data = list.files('./raw/15_final_models/', full.names=TRUE),
# 								.fun = function(filename){
# 									readSBMLmod(filename)
# 								}
# )

#save(models, file='raw/PNASMS2013-07797_DatasetS2_mutli_ecoli_55_models.RData')
load('raw/PNASMS2013-07797_DatasetS2_mutli_ecoli_55_models.RData')

targetFluxes <- c('Ec_biomass_iJO1366_core_53p95M', 'EX_ac(e)')

protoevaluate <- function(genotype,model){
	# avoid an error with no knockouts
	# sybil should really deal with this
	genes <- if(suppressWarnings(all(genotype))){
		NULL
	} else{
		names(genotype)[genotype==FALSE]
	}
	# call sybil to compute the fluxes
	solution <- optimizeProb(object = model, 
													 gene = genes,
													 lb = 0,
													 ub = 0,
													 retOptSol = FALSE)
	phenotype <- c(solution$fluxes[match(targetFluxes,react_id(model))], sum(genotype==FALSE))
}

l_ply(.data = models,
			.fun = function(model){
				start <- as.list(rep(TRUE, times = length(allGenes(model))))
				names(start) <- allGenes(model)
				ancestors <- list(list(genotype = start))
				
				evaluate <- function(x){protoevaluate(x,model)}
				
				GDMO(population=50, 
						 generations=50, 
						 startingpoint=ancestors, 
						 evaluate= evaluate, 
						 savefile=paste0('results/',model@mod_desc,'_checkpoint.RData')
				)
			}
)