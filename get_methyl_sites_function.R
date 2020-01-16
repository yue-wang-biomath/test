get_methyl_function   <- function(sample_name){
  SN_original_data    <- read.csv('E:/Workspace/Workspace/SingleNucleotideReselutionData.csv');
  sample_column       <- SN_original_data[, sample_name]
  sample_index        <- which(sample_column == replicate(length(sample_column), 1))
  sample_list         <- SN_original_data[sample_index, ]
  sample_data         <- GRanges(Rle(sample_list[, 'chromsome']),
                                 IRanges(sample_list[, 'modEnd'],
                                 width = replicate(length(sample_list[, 'modEnd']), 1)))
  strand(sample_data) <- sample_list[, 'strand']
  return(sample_data)
}

get_trans_function       <- function(trans, sample_methyl){
  sample_methyl_pos      <- start(ranges(sample_methyl))
  sample_methyl_seqnames <- data.frame(seqnames(sample_methyl))[[1]]
  trans_metadata         <- mcols(trans)

  index_trans            <- 1:length(trans)
  index_methyl           <- data.frame(trans)[, 'xHits']
  seqnames               <- sample_methyl_seqnames[data.frame(trans)[, 'xHits']]
  methyl_pos             <- sample_methyl_pos[data.frame(trans)[, 'xHits']]
  strand                 <- data.frame(trans)[, 'strand']
  trans_ID               <- data.frame(trans)[, 'transcriptsHits']

  trans_info             <- data.frame(index_trans  = index_trans,
                                       index_methyl = index_methyl,
                                       seqnames     = seqnames,
                                       methyl_pos   = methyl_pos,
                                       strand       = strand,
                                       trans_ID     = trans_ID)
  return(trans_info)
}

get_validation_site_function <- function(tx0, num_validation_site, num_transcript_site){
  cds_by_tx0_sample          <- sample(unlist(tx0), num_validation_site)
  for (i in 1:num_transcript_site){
  validation_site_pos        <- start(cds_by_tx0_sample) + runif(num_validation_site, 0, 1) * width(cds_by_tx0_sample)
  validation_site_seqnames   <- seqnames(cds_by_tx0_sample)
  validation_site_strand     <- strand(cds_by_tx0_sample)
  validation_site_i          <- GRanges(seqnames = validation_site_seqnames,
                                        IRanges(validation_site_pos,
                                        width    = replicate(num_validation_site, 1)),
                                        validation_site_strand)
  if (i == 1){
  validation_site            <- validation_site_i
  }
  validation_site            <- c(validation_site, validation_site_i)
  }
  return(validation_site)
}


width_mtr <- matrix(1, nrow(width_mtr), length(width_mtr))
get_correct_prob_function   <- function(num_bin_sum, align_mtr, width_mtr, trans_info){

  alpha                     <- matrix(1/num_bin_sum, 1, num_bin_sum)
  index_methyl              <- trans_info[, 'index_methyl']

  for (j in 1:20){
    numerator_prob          <- align_mtr * width_mtr

    for (k in 1:num_bin_sum){
      numerator_prob[, k]   <- numerator_prob[, k] * alpha[k]
    }
    row_sum_prob            <- as.vector(rowSums(numerator_prob))
    row_sum_prob            <- data.frame(group_name = as.character(index_methyl), value = row_sum_prob)
    sum_prob                <- aggregate(row_sum_prob[,'value'],
                                         by = list(group_name = factor(trans_info[,'index_methyl'], levels = unique(trans_info[,'index_methyl']))),
                                         FUN = sum)
    names(sum_prob)         <- c('index_methyl', 'value')
    denominator_prob        <- sum_prob[, 'value']
    names(denominator_prob) <- sum_prob[, 'index_methyl']
    denominator_prob        <- denominator_prob[as.character(trans_info[, 'index_methyl'])]
    denominator_prob        <- replicate(num_bin_sum, denominator_prob)
    denominator_prob[denominator_prob == 0] <- 1

    prob_mtr                <- numerator_prob / denominator_prob
    alpha_numerator         <- colSums(prob_mtr)
    alpha                   <- alpha_numerator/sum(alpha_numerator)
  }


  return(alpha)
}

  likelihood <- data.frame(index_methyl = index_methyl,
                          likelihood_i =  rowSums(prob_mtr*width_mtr*alpha))
likelihood$index_methyl<-factor(likelihood$index_methyl)
likelihood_sum <- tapply(likelihood$likelihood_i,likelihood$index_methyl,sum)
likelihood <- sum(log(likelihood_sum)[!is.infinite(log(likelihood_sum))])
alpha <- cbind(likelihood, alpha)
d_c_f = length(unique(index_methyl))*alpha/colSums(width_mtr)
d_base = length(unique(index_methyl))*alpha/sum(colSums(width_mtr))


overlap <- data.frame(index_methyl = index_methyl, overlap_i =  1)
overlap$index_methyl<-factor(overlap$index_methyl)
overlap <- tapply(overlap$overlap_i,overlap$index_methyl,sum)
overlap_mtr <- overlap[as.character(index_methyl)]
overlap_mtr <- replicate(num_bin_sum, overlap_mtr)
