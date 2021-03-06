#' eQTLsingle: Discover single-cell eQTLs from scRNA-seq data only
#'
#' A function to discover eQTLs from scRNA-seq data. The SNV-gene pair information is generated with function eQTLsingle_build_metadata.
#' @param expressionMatrix A dataframe, describes gene expressions (gene expressions by cells), rownames of this dataframe are Geneid, colnames of this dataframe are CellId
#' @param snv.gene.pair.metadata the snv-gene pair for testing, can be generated by the function eQTLsingle_build_metadata
#' @param p.adjust.method Method for adjusting P-values for multiple comparisons. Default method is bonferroni correction
#' @importFrom stringr str_split
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom maxLik maxLik
#' @importFrom pscl zeroinfl
#' @importFrom stats p.adjust pchisq plogis
#' @return A dataframe, each row describes eQTL discovering result of a snv-gene pair
#' \itemize{
#' \item SNVid: Id of SNV, represented as CHR__POS
#' \item Geneid: gene id of target gene
#' \item sample_size_Ref, sample_size_Alt: number of cells in REF and ALT group
#' \item theta_Ref, theta_Alt, mu_Ref, mu_Alt, size_Ref, size_Alt, prob_Ref, prob_Alt: estimated parameters of ZINB for REF group and ALT group
#' \item total_mean_Ref, total_mean_Alt: mean of the REF group and ALT group
#' \item foldChange: fold change of mean of REF group (total_mean_Ref) with respect to mean of ALT group (total_mean_Alt)
#' \item chi2LR1: chi square statistic for the test
#' \item pvalue, adjusted_pvalue: pvalue and adjusted pvalue of the test. If adjusted p-value is smaller than some threshold, this SNV shows significant eQTL effect on the target gene
#' \item Remark: to record abnormity
#' }
#' @examples
#' # Generate metadata for further eQTL analysis.
#' # Load snvMatrix and expressionMatrix
#' data(toy_snvMatrix)
#' data(toy_expressionMatrix)
#' # Only test the SNVs in which cells with either genotypes (REF and ALT) are at least 30 cells
#' # Only test the genes which express in > 30 cells in group REF and ALT with expression level >1
#' snv.gene.pair.metadata <- eQTLsingle_build_metadata(toy_snvMatrix,
#'                                                     toy_expressionMatrix,
#'                                                     snv.number.of.cells=30,
#'                                                     expression.min=1,
#'                                                     expression.number.of.cells=30)
#' # Discover eQTLs
#' eQTL.result <- eQTLsingle(toy_expressionMatrix, snv.gene.pair.metadata)
#' @export

eQTLsingle <- function(expressionMatrix, snv.gene.pair.metadata, p.adjust.method = "bonferroni"){
  # Invalid input control
  if(!is.matrix(expressionMatrix) & !is.data.frame(expressionMatrix) & class(expressionMatrix)[1] != "dgCMatrix")
    stop("Wrong data type of 'expressionMatrix'")
  if(sum(is.na(expressionMatrix)) > 0)
    stop("NA detected in 'expressionMatrix'");gc();
  if(sum(expressionMatrix < 0) > 0)
    stop("Negative value detected in 'expressionMatrix'");gc();
  if(all(expressionMatrix == 0))
    stop("All elements of 'expressionMatrix' are zero");gc();

  # Function of testing homogeneity of two ZINB populations
  CalleQTL <- function(i){
    # Memory management
    if(i %% 100 == 0)
      gc()

    # gene and snv extraction
    snvid <- snv.gene.pair.metadata[i, "SNVid"]
    cells_character <- snv.gene.pair.metadata[i, 'CellList']
    cells <- stringr::str_split(cells_character, ",")[[1]]
    ref_cells_character <- snv.gene.pair.metadata[i, 'Ref_cells']
    ref_cells <- stringr::str_split(ref_cells_character, ",")[[1]]
    alt_cells_character <- snv.gene.pair.metadata[i, 'Alt_cells']
    alt_cells <- stringr::str_split(alt_cells_character, ",")[[1]]
    genes <- snv.gene.pair.metadata[i, 'GeneList']
    genes <- stringr::str_split(genes, ",")[[1]]
    gene.cnt <- 0

    # for each snv, whole data.frame for result
    results_SNV <- data.frame(SNVid = character(),
                              Geneid = character(),
                              sample_size_1 = integer(),
                              sample_size_2 = integer(),
                              theta_1 = double(),
                              theta_2 = double(),
                              mu_1 = double(),
                              mu_2 = double(),
                              size_1 = double(),
                              size_2 = double(),
                              prob_1 = double(),
                              prob_2 = double(),
                              total_mean_1 = double(),
                              total_mean_2 = double(),
                              foldChange = double(),
                              chi2LR1 = double(),
                              pvalue = double(),
                              adjusted_pvalue = double(),
                              Remark = character(),
                              stringsAsFactors=FALSE)

    # for each gene
    for (gene in genes){
      gene.cnt <- gene.cnt + 1
      counts_1 <- unlist(expressionMatrix[gene, ref_cells]) # gene expression for ref group
      counts_2 <- unlist(expressionMatrix[gene, alt_cells]) # gene expression for alt group
      results_gene <- data.frame(SNVid = snvid, Geneid = gene, sample_size_1 = length(counts_1), sample_size_2 = length(counts_2), theta_1 = NA, theta_2 = NA, mu_1 = NA, mu_2 = NA, size_1 = NA, size_2 = NA,
                                 prob_1 = NA, prob_2 = NA, total_mean_1 = NA, total_mean_2 = NA, foldChange = NA, chi2LR1 = NA, pvalue = NA, adjusted_pvalue = NA, Remark = NA,
                                 stringsAsFactors=FALSE)

      # general difference
      totalMean_1 <- mean(counts_1)
      totalMean_2 <- mean(counts_2)
      foldChange <- totalMean_1/totalMean_2

      results_gene[1,"total_mean_1"] <- totalMean_1
      results_gene[1,"total_mean_2"] <- totalMean_2
      results_gene[1,"foldChange"] <- foldChange

      # Log likelihood functions
      logL <- function(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2){
        logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1, pstr0 = theta_1, log = TRUE))  # log-likelihood for count1 under parameter
        logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2, pstr0 = theta_2, log = TRUE)) # log-likelihood for count2 under parameter
        logL <- logL_1 + logL_2
        logL
      }
      logL2 <- function(param){
        theta_resL2 <- param[1]
        size_1_resL2 <- param[2]
        prob_1_resL2 <- param[3]
        size_2_resL2 <- param[4]
        prob_2_resL2 <- param[5]
        logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL2NZ <- function(param){
        theta_resL2 <- 0
        size_1_resL2 <- param[1]
        prob_1_resL2 <- param[2]
        size_2_resL2 <- param[3]
        prob_2_resL2 <- param[4]
        logL_1 <- sum(dzinegbin(counts_1, size = size_1_resL2, prob = prob_1_resL2, pstr0 = theta_resL2, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_2_resL2, prob = prob_2_resL2, pstr0 = theta_resL2, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3 <- function(param){
        theta_1_resL3 <- param[1]
        size_resL3 <- param[2]
        prob_resL3 <- param[3]
        theta_2_resL3 <- param[4]
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3NZ1 <- function(param){
        theta_1_resL3 <- 0
        size_resL3 <- param[1]
        prob_resL3 <- param[2]
        theta_2_resL3 <- param[3]
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3NZ2 <- function(param){
        theta_1_resL3 <- param[1]
        size_resL3 <- param[2]
        prob_resL3 <- param[3]
        theta_2_resL3 <- 0
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3AZ1 <- function(param){
        theta_1_resL3 <- 1
        size_resL3 <- param[1]
        prob_resL3 <- param[2]
        theta_2_resL3 <- param[3]
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3AZ2 <- function(param){
        theta_1_resL3 <- param[1]
        size_resL3 <- param[2]
        prob_resL3 <- param[3]
        theta_2_resL3 <- 1
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3NZ1AZ2 <- function(param){
        theta_1_resL3 <- 0
        size_resL3 <- param[1]
        prob_resL3 <- param[2]
        theta_2_resL3 <- 1
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      logL3NZ2AZ1 <- function(param){
        theta_1_resL3 <- 1
        size_resL3 <- param[1]
        prob_resL3 <- param[2]
        theta_2_resL3 <- 0
        logL_1 <- sum(dzinegbin(counts_1, size = size_resL3, prob = prob_resL3, pstr0 = theta_1_resL3, log = TRUE))
        logL_2 <- sum(dzinegbin(counts_2, size = size_resL3, prob = prob_resL3, pstr0 = theta_2_resL3, log = TRUE))
        logL <- logL_1 + logL_2
        logL
      }
      judgeParam <- function(param){
        if((param >= 0) & (param <= 1))
          res <- TRUE
        else
          res <- FALSE
        res
      }

      # MLE of parameters of ZINB counts_1
      if(sum(counts_1 == 0) > 0){
        if(sum(counts_1 == 0) == length(counts_1)){
          theta_1 <- 1
          mu_1 <- 0
          size_1 <- 1
          prob_1 <- size_1/(size_1 + mu_1)
        }else{
          options(show.error.messages = FALSE)
          zinb_try <- try(gamlssML(counts_1, family="ZINBI"), silent=TRUE)
          options(show.error.messages = TRUE)
          if('try-error' %in% class(zinb_try)){
            zinb_try_twice <- try(zeroinfl(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
            if('try-error' %in% class(zinb_try_twice)){
              print("MLE of ZINB failed!");
              results_gene[1,"Remark"] <- "ZINB failed!"
              return(results_gene)
            }else{
              zinb_1 <- zinb_try_twice
              theta_1 <- plogis(zinb_1$coefficients$zero);names(theta_1) <- NULL
              mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
              size_1 <- zinb_1$theta;names(size_1) <- NULL
              prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
            }
          }else{
            zinb_1 <- zinb_try
            theta_1 <- zinb_1$nu;names(theta_1) <- NULL
            mu_1 <- zinb_1$mu;names(mu_1) <- NULL
            size_1 <- 1/zinb_1$sigma;names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
          }
        }
      }else{
        op <- options(warn=2)
        nb_try <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
        options(op)
        if('try-error' %in% class(nb_try)){
          nb_try_twice <- try(fitdistr(counts_1, "Negative Binomial"), silent=TRUE)
          if('try-error' %in% class(nb_try_twice)){
            nb_try_again <- try(mle2(counts_1~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_1), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
            if('try-error' %in% class(nb_try_again)){
              nb_try_fourth <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
              if('try-error' %in% class(nb_try_fourth)){
                print("MLE of NB failed! 33");
                results_gene[1,"Remark"] <- "NB failed!"
                return(results_gene)
              }else{
                nb_1 <- nb_try_fourth
                theta_1 <- 0
                mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
                size_1 <- nb_1$theta;names(size_1) <- NULL
                prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
              }
            }else{
              nb_1 <- nb_try_again
              theta_1 <- 0
              mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
              size_1 <- 1/nb_1@coef["invk"];names(size_1) <- NULL
              prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
            }
          }else{
            nb_1 <- nb_try_twice
            theta_1 <- 0
            mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
            size_1 <- nb_1$estimate["size"];names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
          }
        }else{
          nb_1 <- nb_try
          theta_1 <- 0
          mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
          size_1 <- nb_1$theta;names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }

      # MLE of parameters of ZINB counts_2
      if(sum(counts_2 == 0) > 0){
        if(sum(counts_2 == 0) == length(counts_2)){
          theta_2 <- 1
          mu_2 <- 0
          size_2 <- 1
          prob_2 <- size_2/(size_2 + mu_2)
        }else{
          options(show.error.messages = FALSE)
          zinb_try <- try(gamlssML(counts_2, family="ZINBI"), silent=TRUE)
          options(show.error.messages = TRUE)
          if('try-error' %in% class(zinb_try)){
            zinb_try_twice <- try(zeroinfl(formula = counts_2 ~ 1 | 1, dist = "negbin"), silent=TRUE)
            if('try-error' %in% class(zinb_try_twice)){
              print("MLE of ZINB failed!");
              results_gene[1,"Remark"] <- "ZINB failed!"
              return(results_gene)
            }else{
              zinb_2 <- zinb_try_twice
              theta_2 <- plogis(zinb_2$coefficients$zero);names(theta_2) <- NULL
              mu_2 <- exp(zinb_2$coefficients$count);names(mu_2) <- NULL
              size_2 <- zinb_2$theta;names(size_2) <- NULL
              prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
            }
          }else{
            zinb_2 <- zinb_try
            theta_2 <- zinb_2$nu;names(theta_2) <- NULL
            mu_2 <- zinb_2$mu;names(mu_2) <- NULL
            size_2 <- 1/zinb_2$sigma;names(size_2) <- NULL
            prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
          }
        }
      }else{
        op <- options(warn=2)
        nb_try <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
        options(op)
        if('try-error' %in% class(nb_try)){
          nb_try_twice <- try(fitdistr(counts_2, "Negative Binomial"), silent=TRUE)
          if('try-error' %in% class(nb_try_twice)){
            nb_try_again <- try(mle2(counts_2~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_2), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
            if('try-error' %in% class(nb_try_again)){
              nb_try_fourth <- try(glm.nb(formula = counts_2 ~ 1), silent=TRUE)
              if('try-error' %in% class(nb_try_fourth)){
                print("MLE of NB failed! 11");
                results_gene[1,"Remark"] <- "NB failed!"
                return(results_gene)
              }else{
                nb_2 <- nb_try_fourth
                theta_2 <- 0
                mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
                size_2 <- nb_2$theta;names(size_2) <- NULL
                prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
              }
            }else{
              nb_2 <- nb_try_again
              theta_2 <- 0
              mu_2 <- exp(nb_2@coef["logmu"]);names(mu_2) <- NULL
              size_2 <- 1/nb_2@coef["invk"];names(size_2) <- NULL
              prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
            }
          }else{
            nb_2 <- nb_try_twice
            theta_2 <- 0
            mu_2 <- nb_2$estimate["mu"];names(mu_2) <- NULL
            size_2 <- nb_2$estimate["size"];names(size_2) <- NULL
            prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
          }
        }else{
          nb_2 <- nb_try
          theta_2 <- 0
          mu_2 <- exp(nb_2$coefficients);names(mu_2) <- NULL
          size_2 <- nb_2$theta;names(size_2) <- NULL
          prob_2 <- size_2/(size_2 + mu_2);names(prob_2) <- NULL
        }
      }

      # Restricted MLE under H0
      # to merge two counts vector together to build a uniform data
      if(sum(c(counts_1, counts_2) == 0) > 0){
        options(show.error.messages = FALSE)
        zinb_try <- try(gamlssML(c(counts_1, counts_2), family="ZINBI"), silent=TRUE)
        options(show.error.messages = TRUE)
        if('try-error' %in% class(zinb_try)){
          zinb_try_twice <- try(zeroinfl(formula = c(counts_1, counts_2) ~ 1 | 1, dist = "negbin"), silent=TRUE)
          if('try-error' %in% class(zinb_try_twice)){
            print("MLE of ZINB failed!");
            results_gene[1,"Remark"] <- "ZINB failed!"
            return(results_gene)
          }else{
            zinb_res <- zinb_try_twice
            theta_res <- plogis(zinb_res$coefficients$zero);names(theta_res) <- NULL
            mu_res <- exp(zinb_res$coefficients$count);names(mu_res) <- NULL
            size_res <- zinb_res$theta;names(size_res) <- NULL
            prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
          }
        }else{
          zinb_res <- zinb_try
          theta_res <- zinb_res$nu;names(theta_res) <- NULL
          mu_res <- zinb_res$mu;names(mu_res) <- NULL
          size_res <- 1/zinb_res$sigma;names(size_res) <- NULL
          prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
        }
      }else{
        op <- options(warn=2)
        nb_try <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
        options(op)
        if('try-error' %in% class(nb_try)){
          nb_try_twice <- try(fitdistr(c(counts_1, counts_2), "Negative Binomial"), silent=TRUE)
          if('try-error' %in% class(nb_try_twice)){
            nb_try_again <- try(mle2(c(counts_1, counts_2)~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(c(counts_1, counts_2)), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
            if('try-error' %in% class(nb_try_again)){
              nb_try_fourth <- try(glm.nb(formula = c(counts_1, counts_2) ~ 1), silent=TRUE)
              if('try-error' %in% class(nb_try_fourth)){
                print("MLE of NB failed! 22");
                results_gene[1,"Remark"] <- "NB failed!"
                return(results_gene)
              }else{
                nb_res <- nb_try_fourth
                theta_res <- 0
                mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
                size_res <- nb_res$theta;names(size_res) <- NULL
                prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
              }
            }else{
              nb_res <- nb_try_again
              theta_res <- 0
              mu_res <- exp(nb_res@coef["logmu"]);names(mu_res) <- NULL
              size_res <- 1/nb_res@coef["invk"];names(size_res) <- NULL
              prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
            }
          }else{
            nb_res <- nb_try_twice
            theta_res <- 0
            mu_res <- nb_res$estimate["mu"];names(mu_res) <- NULL
            size_res <- nb_res$estimate["size"];names(size_res) <- NULL
            prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
          }
        }else{
          nb_res <- nb_try
          theta_res <- 0
          mu_res <- exp(nb_res$coefficients);names(mu_res) <- NULL
          size_res <- nb_res$theta;names(size_res) <- NULL
          prob_res <- size_res/(size_res + mu_res);names(prob_res) <- NULL
        }
      }

      # LRT test of H0
      chi2LR1 <- 2 *(logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2) - logL(counts_1, theta_res, size_res, prob_res, counts_2, theta_res, size_res, prob_res))
      pvalue <- 1 - pchisq(chi2LR1, df = 3)

      # Format output
      results_gene[1,"theta_1"] <- theta_1
      results_gene[1,"theta_2"] <- theta_2
      results_gene[1,"mu_1"] <- mu_1
      results_gene[1,"mu_2"] <- mu_2
      results_gene[1,"size_1"] <- size_1
      results_gene[1,"size_2"] <- size_2
      results_gene[1,"prob_1"] <- prob_1
      results_gene[1,"prob_2"] <- prob_2
      results_gene[1,"chi2LR1"] <- chi2LR1
      results_gene[1,"pvalue"] <- pvalue

      # add snv-gene result to current snv result
      results_SNV <- rbind(results_SNV, results_gene)
    }
    # return SNV result
    return(results_SNV)
  }

  # Call SNV by SNV
  results <- data.frame(SNVid = character(),
                        Geneid = character(),
                        sample_size_1 = integer(),
                        sample_size_2 = integer(),
                        theta_1 = double(),
                        theta_2 = double(),
                        mu_1 = double(),
                        mu_2 = double(),
                        size_1 = double(),
                        size_2 = double(),
                        prob_1 = double(),
                        prob_2 = double(),
                        total_mean_1 = double(),
                        total_mean_2 = double(),
                        foldChange = double(),
                        chi2LR1 = double(),
                        pvalue = double(),
                        adjusted_pvalue = double(),
                        Remark = character(),
                        stringsAsFactors=FALSE)

  for(i in 1:dim(snv.gene.pair.metadata)[1]){
    message(paste0("Processing ", i, " / ", dim(snv.gene.pair.metadata)[1], " term"))
    results <- rbind(results, CalleQTL(i))
  }

  # change column names for a better understanding
  colnames(results)[colnames(results) == 'sample_size_1'] <- 'sample_size_Ref'
  colnames(results)[colnames(results) == 'sample_size_2'] <- 'sample_size_Alt'
  colnames(results)[colnames(results) == 'theta_1'] <- 'theta_Ref'
  colnames(results)[colnames(results) == 'theta_2'] <- 'theta_Alt'
  colnames(results)[colnames(results) == 'mu_1'] <- 'mu_Ref'
  colnames(results)[colnames(results) == 'mu_2'] <- 'mu_Alt'
  colnames(results)[colnames(results) == 'size_1'] <- 'size_Ref'
  colnames(results)[colnames(results) == 'size_2'] <- 'size_Alt'
  colnames(results)[colnames(results) == 'prob_1'] <- 'prob_Ref'
  colnames(results)[colnames(results) == 'prob_2'] <- 'prob_Alt'
  colnames(results)[colnames(results) == 'total_mean_1'] <- 'total_mean_Ref'
  colnames(results)[colnames(results) == 'total_mean_2'] <- 'total_mean_Alt'

  # Format output results
  results[,"adjusted_pvalue"] <- p.adjust(results[,"pvalue"], method=p.adjust.method)

  if(exists("lastFuncGrad") & exists("lastFuncParam")){
    lastFuncGrad <- NULL
    lastFuncParam <- NULL
    remove(lastFuncGrad, lastFuncParam, envir=.GlobalEnv)
  }

  return(results)
}
