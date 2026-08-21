# Feature filtering
#
# Part of the 'proteoSE' package. Grouped by theme during the 2026 refactor
# (see AUDIT.md). Function bodies are preserved verbatim from the original
# main.R / mlasse.txt; renames and cleanup follow in later phases.


#' Filter on occurence 
#' 
#' Filters proteins with too many missing values. Parameters mimic the options
#' given by Perseus.
#'
#' @param se summarized experiment
#' @param perc_na fraction of missing values that is tolerated
#' @param filter_mode "each_group": perc_na in every group; "one_group": perc_na
#' in at least one group; "total": perc_na across all samples.
#'
#' @return summarized experiment without entries with too many missing values
#' @importFrom dplyr group_by summarize filter ungroup select mutate n
#' @importFrom BiocGenerics unique
#' @export
filter_perseus<-function(se, perc_na = 0.33, filter_mode = "each_group"){
  .assert_se(se, require_coldata = "condition")

  if(filter_mode %in% c("one_group", "each_group", "total")){

      # 1. Step: create long table with groups

      data_long <- .se_long(se)

      # 2. Step: Apply filtering

      if(filter_mode=="one_group"){

        filtered_terms<-data_long%>%
          group_by(name, condition)%>%
          summarize(number_of_na=sum(is.na(intensity)), group_size=n()) %>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          ungroup()%>%
          select(name)%>%
          BiocGenerics::unique()

      }

      if(filter_mode=="each_group"){
        filtered_terms<-data_long%>%
          group_by(name, condition)%>%
          summarize(number_of_na=sum(is.na(intensity)), group_size=n())%>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          mutate(n_data=length(unique(data_long$condition)))%>%
          group_by(name)%>%dplyr::filter(n()==n_data)%>%
          select(name)%>%unique()

      }

      if(filter_mode=="total"){

        filtered_terms<-data_long%>%
          group_by(sample, name)%>%
          summarize(m=mean(intensity))%>%ungroup()%>%
          group_by(name)%>%
          summarise(number_of_na=sum(is.na(m)), group_size=n())%>%
          filter(number_of_na<=round(group_size*perc_na))%>%
          select(name)%>%unique()

      }

      # 3. Step: return summarized experiment with only remaining genes
      return(se[filtered_terms %>% unlist(),])
  }else{
    cat("Unknown filter_mode. Select one of: \"one_group\", \"each_group\", \"total\".")
  }
}

