pipeline <- function(df = data,
                     analysis_name,
                     data_name,
                     comparison,
                     const = 0,
                     exclude = ''
                    ){
    
    df <- limma_analysis(df = df,
                         analysis_name = analysis_name,
                         data_name = data_name,
                         col_name = 'unique_name',
                         col_value = 'norm_int',
                         comparison = comparison,
                         const = const,
                         exclude = exclude
                        )
    
    
    # Export normalized intensity (metabolite levels in medium over time)
    df <- csv_total(df = df,
                    analysis_name = analysis_name,
                    data_name = data_name,
                    statistic = 'limma',
                    export = TRUE
                   )
    
    df <- csv_volcano(df = df,
                      analysis_name = analysis_name,
                      data_name = data_name,
                      statistic = 'limma',
                      export = TRUE
                     )
    
    df <- csv_metaboanalyst(df = df,
                            analysis_name = analysis_name,
                            data_name = data_name,
                            value = 'norm_int',
                            export = TRUE
                           )
    
    return(df)
    
}

load_libraries <- function(){
    library(tidyverse)
    library(limma)
    #library(MetaboAnalystR) # Loaded in the function instead to improve loading times
    library(rstatix)
    
    options(dplyr.summarise.inform = FALSE) # Suppresses "`summarise()` has grouped output by 'group'. You can override using the `.groups` argument."
}

initialize_data <- function(expID, # E.g. 'ABC042'
                            analysis_type, # E.g. 'lipids'
                            working_directory,
                            data_targeted_neg,
                            data_targeted_pos,
                            data_untargeted_neg,
                            data_untargeted_pos,
                            data_BCA,
                            metadata_samples,
                            metadata_groups,
                            metadata_library_neg,
                            metadata_library_pos
                           ){
    message("Initialize data")
    
    cat('Set working directory\n')
    setwd(working_directory)    
    cat('Success: set working directory\n')
    
    cat('Load files\n')
    # Load files
    cat('data_targeted_neg\n')
    data_targeted_neg <- read.csv(data_targeted_neg)
    cat('data_targeted_pos\n')
    data_targeted_pos <- read.csv(data_targeted_pos)
    cat('data_untargeted_neg\n')
    data_untargeted_neg <- read.csv(data_untargeted_neg)
    cat('data_untargeted_pos\n')
    data_untargeted_pos <- read.csv(data_untargeted_pos)
    cat('data_BCA\n')
    data_BCA <- read.csv(data_BCA)
    cat('metadata_samples\n')
    metadata_samples <- read.csv(metadata_samples)
    cat('metadata_groups\n')
    metadata_groups <- read.csv(metadata_groups)
    cat('metadata_library_neg\n')
    metadata_library_neg <- read.csv(metadata_library_neg)
    cat('metadata_library_pos\n')
    metadata_library_pos <- read.csv(metadata_library_pos)
    cat("Success: Load files\n")
    
    
    # Prepare raw data files
    
    cat("Adding mode\n")
    # Add mode  
    data_targeted_neg <- add_mode(data_targeted_neg, 'neg')
    data_targeted_pos <- add_mode(data_targeted_pos, 'pos')
    data_untargeted_neg <- add_mode(data_untargeted_neg, 'neg')
    data_untargeted_pos <- add_mode(data_untargeted_pos, 'pos')
    cat("Success: Adding mode\n")
    
    # Combine respective data frames
    cat("Combine data frames")
    data_targeted <- rbind(data_targeted_neg, data_targeted_pos)
    data_untargeted <- rbind(data_untargeted_neg, data_untargeted_pos)
    metadata_library <- rbind (metadata_library_neg, metadata_library_pos)
    cat("Success: Combine data frames\n")
    
    # Remove " (1)" from name if present
    # Sometimes when reopening El-Maven files, compounds get numbered. This might be a bug.
    cat("Remove trailing numbers from compounds\n")
    data_targeted <- remove_number_from_name(data_targeted)
    
    cat("Success: Remove trailing numbers from compounds\n")
    
    # Add unique name in case the same metabolite is identified in both modes
    cat("Add unique name\n")
    data_targeted <- add_unique_name(data_targeted)
    data_untargeted <- add_unique_name(data_untargeted)
    cat("Success: Add unique name\n")
    
    # Prepare data_BCA
    cat("Prepare BCA data\n")
    data_BCA <- prepare_data_BCA(data_BCA, analysis_type)
    cat("Success: Preare BCA data\n")
    
    # Rename columns from metadata file
    cat("Preare samples metadata\n")
    metadata_samples <- prepare_metadata_samples(metadata_samples)
    cat("Success: Preare samples metadata\n")
    
    # Initialize list in order to add items later
    df <- list()
    
    cat("Write variables for output\n")
    df$metadata$working_directory = working_directory
    df$metadata$expID = expID
    df$metadata$analysis_type = analysis_type
    df$metadata$samples = metadata_samples
    df$metadata$groups = metadata_groups
    df$metadata$library = metadata_library
    
    df$raw_data$targeted = data_targeted
    df$raw_data$untargeted = data_untargeted
    df$raw_data$BCA = data_BCA
    cat("Success: Write variables for output\n")
    
    message("Success: Initialize data")
    
    return(df)   
    
    }

add_mode <- function(df, selected_mode){
    df <- df %>%
    mutate(mode = as.factor(selected_mode))
    return(df)
}

remove_number_from_name <- function(df){
    df <- df %>%
        mutate(compound = case_when(grepl("[[:space:]]\\([[:digit:]]\\)", compound) ~ str_sub(compound,1,str_length(compound)-4),
                                    TRUE ~ as.character(compound)
                                   ),
               compound = as.factor(compound)
              )

    return(df)
}    

add_unique_name <- function(df){
    df <- df %>%
    mutate(unique_name = paste(compound, medRt, mode, sep='@'))
    return(df)
}

prepare_data_BCA <- function(df, analysis_type){
    
    df <- df %>%
        select(all_of(paste('sample', analysis_type, sep='_')), # Column correct column containing sample names from respective analysis
               group,
               protein,
               mc_protein
              ) %>%
        rename(sample = all_of(paste('sample', analysis_type, sep='_'))) # Rename
    
    return(df)
    
}

prepare_metadata_samples <- function(df){
    
    df <- df %>%
        select(c('Sample', 'Cohort', 'Injection.Order')) %>% # Drop unnecessary columns 'Scaling' and 'Color'
        # Rename columns for consistent name scheme
        rename(sample = Sample,
               group = Cohort,
               injection_order = Injection.Order
              )
    
}

prep_scaling <- function(df = data,
                         export = TRUE){
    
    message("Preparing scaling data frame")
    
    # Summed signal
    
    # Summed signal of all identified peaks
    cat("Get summed signal\n")
    cat("ss_identified\n")
    ss_identified_neg <- get_summed_signal(filter(df$raw_data$targeted, label == 'g'), 'ss_identified', 'neg')
    ss_identified_pos <- get_summed_signal(filter(df$raw_data$targeted, label == 'g'), 'ss_identified', 'pos')
    cat("Success: ss_identified\n")
    # Summed signal of all peaks found during targeted peak finding
    cat("ss_targeted\n")
    ss_targeted_neg <- get_summed_signal(df$raw_data$targeted, 'ss_targeted', 'neg')
    ss_targeted_pos <- get_summed_signal(df$raw_data$targeted, 'ss_targeted', 'pos')
    cat("Success: ss_targeted\n")
    
    # Summed signal of all peaks found during untargeted peak finding
    cat("ss_untargeted\n")
    ss_untargeted_neg <- get_summed_signal(df$raw_data$untargeted, 'ss_untargeted', 'neg')
    ss_untargeted_pos <- get_summed_signal(df$raw_data$untargeted, 'ss_untargeted', 'pos')
    cat("Success: ss_targeted\n")
    
    # Collect all summed signals in a single data frame
    cat("Combine in single data frame\n")
    ss <- rbind(ss_identified_neg,
                ss_identified_pos,
                
                ss_targeted_neg,
                ss_targeted_pos,
                
                ss_untargeted_neg,
                ss_untargeted_pos
               )
    cat("Success: Combine in single data frame\n")
    cat("Success: Get summed signal\n")
    
    # Collect all peaks marked as internal standard
    cat("Collect IS\n")
    IS <- collect_IS(data$raw_data$targeted)
    cat("Success: Collect IS\n")
    
    # Collect summed signal and internal standard in single data frame
    cat("Combine summed signal and IS\n")
    IS <- rbind(ss, IS) %>%
        rename(IS = compound)
    cat("Success: Combine summed signal and IS\n")
    
    
    # Add protein and save
    cat("Add protein\n")
    scaling <- merge(df$raw_data$BCA, IS, by='sample', all.y=TRUE)
    cat("Success: Add protein\n")
    
    
    # Add group metadata
    cat("Add group metadata\n")
    df$scaling <- merge(scaling, df$metadata$groups, by='group', all.x=TRUE)
    cat("Success: Add group metadata\n")
    
    
    if(export == TRUE){
        cat("Export\n")
        setwd(df$metadata$working_directory)
        dir.create(paste('../', 'report', sep=''), showWarnings = FALSE)
             
        write.csv(df$scaling, paste('../report/', df$metadata$expID, '_', df$metadata$analysis_type, '_scaling.csv', sep=''), row.names=FALSE)
        
        setwd(df$metadata$working_directory)
        cat("Success: Export\n")
    }
    
    message('Scaling dataframe creation successfull.')
    
    return(df)
    
}

# Summed signal is calculated for positive and negative mode separately
get_summed_signal <- function(df,
                              IS_name,
                              selected_mode){
    
    df <- df %>%
        filter(mode == selected_mode) %>%
        select(data$metadata$samples$sample) # only select columns containing samples
    
    v <- as.data.frame(colSums(df)) %>% # Calculate summed signal for selected mode
        rename(IS_intensity = 1) %>% # rename first column
        mutate(IS_mc_int = IS_intensity/mean(IS_intensity), # Additionally calculate mean centered summed signal for plotting
               mode = selected_mode,
               compound = IS_name
              )
    v <- tibble::rownames_to_column(v, 'sample')
    
    return(v)    
}

collect_IS <- function(df){
    df <- df %>%
        filter(label == 'b') %>%
        select(compound, mode, data$metadata$samples$sample) %>%
        pivot_longer(data$metadata$samples$sample, names_to='sample', values_to='IS_intensity') %>%
        group_by(compound, mode) %>%
        mutate(IS_mc_int = IS_intensity/mean(IS_intensity)) %>%
        ungroup

    return(df)
}

normalize <- function(df = data,
                      analysis_name,
                      data_name,
                      excluded_groups,
                      excluded_samples,
                      metadata_IS_neg,
                      metadata_IS_pos,
                      control
                     ){
    
    message("Normalizing")
    
    # Prepare list of samples that are included
    cat("Write samples metadata with selected groups and samples removed\n")
    metadata_samples <- df$metadata$samples %>%
        # Exclude specified groups
        filter(!group %in% excluded_groups) %>%
        # Exclude specified samples
        filter(!sample %in% excluded_samples)
    cat("Success: Write samples metadata with selected groups and samples removed\n")
    
    # Save information in respective normalization
    cat("Write information\n")
    df$analysis[[analysis_name]]$settings$excluded_groups = excluded_groups
    df$analysis[[analysis_name]]$settings$excluded_samples = excluded_samples
    df$analysis[[analysis_name]]$settings$metadata_samples = metadata_samples
    cat("Success: Write information\n")
    
    
    # Prepare data frame that defines which compound is normalized by which standard
    # Combine data frames
    cat("Read and combine IS metadata files\n")
    metadata_IS_neg <- read.csv(metadata_IS_neg)
    metadata_IS_pos <- read.csv(metadata_IS_pos)
    
    metadata_IS <- rbind(metadata_IS_neg, metadata_IS_pos)
    
    # Save information in respective normalization
    df$analysis[[analysis_name]]$settings$metadata_IS = metadata_IS
    cat("Success: Read and combine IS metadata files\n")
    
    
    
    # Select data from El-Maven raw data excluding specified groups and samples
    cat("Get data without excluded groups and samples\n")
    intensity <- select_data(df$raw_data$targeted, metadata_samples$sample)
    cat("Success: Get data without excluded groups and samples\n")
        
    # Ontology
    cat("Add library metadata to data\n")
    intensity <- merge(intensity, df$metadata$library, by=c('compound','mode'), all.x=TRUE)
    cat("Success: Add library metadata to data\n")

    
    # Add column of internal standard to each compound
    cat("Add column specifying normalization\n")
    intensity <- merge(intensity, metadata_IS, by=c('class','mode'), all.x=TRUE)
    cat("Success: Add column specifying normalization\n")
    
    
    # Scaling data
    cat("Add scaling data\n")
    intensity <- merge(intensity, data$scaling, by=c('sample', 'mode', 'IS'), all.x=TRUE)
    cat("Success: Add scaling data\n")
    
    

    # Normalization
    
    # Compounds that are normalized with summed signal do not need protein normalization
    cat("Summed signal normalization\n")
    normalization_ss <- intensity %>%
        filter(IS %in% c('ss_identified','ss_targeted','ss_untargeted')) %>%
        mutate(norm_int = intensity/IS_intensity)
    cat("Success: Summed signal normalization\n")
    
    # Compound normalized to IS need protein normalization
    cat("IS + protein normalization\n")
    normalization_IS <- intensity %>%
        filter(!IS %in% c('ss_identified','ss_targeted','ss_untargeted')) %>%
        mutate(norm_int = intensity/IS_intensity/protein)
    cat("Success: IS + protein normalization\n")
    
    # Reunite individual dataframes
    cat("Combine summed signal and IS + protein normalization again\n")
    intensity <- rbind(normalization_ss, normalization_IS)
    cat("Success: Combine summed signal and IS + protein normalization again\n")
    
    
    
    ## Add relative normalized intensity
    # Since intensity can not be compared between compounds the absolute intensity does not matter.\
    # Therefore, data can be displayed relative to control (HPLMax).

    # Get average of control
    cat("Get intensities relative to control\n")
    
    cat("Get control mean\n")
    control_mean_norm_int <- get_control_mean_norm_int(intensity, control)
    cat("Success: Get control mean\n")

    # Add average of HPLMax to every group
    cat("Add control mean\n")
    intensity <- merge(intensity, control_mean_norm_int, by=c('unique_name'), all.x=TRUE)
    cat("Success: Add control mean\n")

    
    # Calculate intensity relative to HPLMax
    cat("Calculate intensities relative to control\n")
    intensity <- intensity %>%
        mutate(rel_norm_int = norm_int/control_mean_int) %>%
        select(!control_mean_int) %>%
        arrange(unique_name, sample, group)
    cat("Success: Calculate intensities relative to control\n")
    
    message('Success: Normalization.')
    

    # Add quality control parameters
    intensity <- add_qc(df = intensity,
                   analysis_name = analysis_name,
                   data_name = data_name
                  )
    
    
    # Save information in respective normalization
    cat("Write variables for output\n")
    df$analysis[[analysis_name]][[data_name]]$data <- intensity
    cat("Success: Write variables for output\n")
    
    
    return(df)
    
}

select_data <- function(df,
                        samples
                       ){

    df <- df %>%
        filter(label == 'g') %>%
        select(unique_name, compound, mode, all_of(samples)) %>%
        pivot_longer(all_of(samples), names_to='sample', values_to='intensity') %>%
        mutate(sample = as.factor(sample))
    
    
    return(df)
    
}

get_control_mean_norm_int <- function(df, control){
    df <- df %>%
        group_by(group, unique_name) %>%
        summarize(control_mean_int = mean(norm_int)) %>%
        ungroup() %>%
        filter(group == control) %>%
        select(!group)
    
    return(df)
}

add_qc <- function(df,
                   analysis_name,
                   data_name
                  ){

    message("Add quality control information.")
    
    cat("Calculate.\n")
    df <- df %>%
        
        # Metabolites
        # Metabolite passes the quality control if:
            # Missingness < 10%, i.e. if the metabolite was detected in 90% of the samples
        mutate(n_all_samples = n_distinct(sample)) %>%
        group_by(unique_name) %>%
        mutate(zeros_compound = sum(norm_int==0),
               perc_zeros_compound = zeros_compound/n_all_samples,
               qc_missingness = case_when(perc_zeros_compound > 0.1 ~ 'failed',
                                         perc_zeros_compound <= 0.1 ~ 'passed'
                                        )
              ) %>%
        ungroup() %>%

        # IQR
        # Metabolite passes the quality control if:
            # IQR is equal or larger than 5% of the normalized intensity (across all samples)
        # That means metabolites which are constant across all conditions fail
        group_by(unique_name) %>%
        mutate(iqr_all_samples = IQR(norm_int),
               avg_norm_int = mean(norm_int),
               rel_iqr_compound = iqr_all_samples/avg_norm_int,
               qc_iqr = case_when(rel_iqr_compound < 0.05 ~ 'failed',
                                  rel_iqr_compound >= 0.05 ~ 'passed'
                                 )
              ) %>%
        ungroup() %>%

        # Samples
        # Sample passes the quality control if:
            # Missingness < 10%, i.e. if 90% of the metabolites were detected in this sample
        mutate(n_qc_passed_missingness = n_distinct(unique_name[qc_missingness == 'passed'])) %>%
        group_by(sample) %>%
        mutate(zeros_sample = sum(norm_int==0),
               perc_zeros_sample = zeros_sample/n_qc_passed_missingness,
               qc_sample = case_when(perc_zeros_sample > 0.1 ~ 'failed',
                                     perc_zeros_sample <= 0.1 ~ 'passed'
                                    )
              ) %>%
        ungroup()
        cat("Success: Calculate.\n")


    message("Success: Add quality control information.")
    
    return(df)
    
}

filter_compounds <- function(df = data,
                             analysis_name,
                             data_name,
                             new_data_name,
                             selection, # 'sd' or 'area'
                             qc_missingness, # 'passed' when compounds that failed qc should be excluded
                             qc_iqr # 'passed' when compounds that failed qc should be excluded
                            ){
    
    
    cat("Get data\n")
    selected <- df$analysis[[analysis_name]][[data_name]]$data
    cat("Success: Get data\n")
    
    if(qc_missingness == 'passed'){
        
        message("Remove compounds that failed missingness qc.")
        
        selected <- selected %>%
            filter(qc_missingness == 'passed')
        
        message("Success: Remove compounds that failed missingness qc.")
    }
    
    if(qc_iqr == 'passed'){
        
        message("Remove compounds that failed iqr qc.")
        
        selected <- selected %>%
            filter(qc_iqr == 'passed')
        
        message("Success: Remove compounds that failed iqr qc.")
    }
    
    
    
    message("Selecting mode for compounds")
    
    # Selection based on standard deviation
    if(selection == 'sd'){
        
        cat("sd\n")        
        selected <- selected %>%
            group_by(compound, mode, group) %>%
            mutate(perc_sd_groups = sd(norm_int)/mean(norm_int)) %>% # Percentage standard deviation for each group
            ungroup() %>%

            group_by(compound, mode) %>%
            mutate(median_perc_sd_groups = median(perc_sd_groups)) %>% # Median standard deviation of all groups
            ungroup() %>%
            arrange(median_perc_sd_groups) %>%
            distinct(compound, sample, .keep_all=TRUE) %>%
            select(!c(perc_sd_groups, median_perc_sd_groups))        
        cat("Success: sd\n")
        
    }
    
    
    # Selection based on (raw/unnnormalized) peak area
    if(selection == 'area'){
        
        cat("area\n")
        selected <- selected %>%
            group_by(unique_name) %>%
            mutate(median_int = median(intensity)) %>%
            ungroup() %>%
            arrange(median_int) %>%
            distinct(compound, sample, .keep_all=TRUE) %>%
            select(!c(median_int))
        cat("Success: area\n")
        
    }
    
    cat("Write variables for output\n")
    df$analysis[[analysis_name]][[new_data_name]]$data <- selected
    cat("Success: Write variables for output\n")
    
    message("Success: Selecting mode for compounds")
    
    return(df)
    
}

limma_analysis <- function(df = data,
                           analysis_name,
                           data_name,
                           col_name,
                           col_value,
                           comparison,
                           const = 0,
                           exclude = '' # See prep_limma_df, usually not necessary because data is filtered before it is written in data_name$data
                          ){

    message("Limma analysis")

    # Save information in respective analysis
    cat("Write comparison information\n")
    df$analysis[[analysis_name]][[data_name]]$limma$comparison = comparison
    cat("Success: Write comparison information\n")

    
    # Prepare dataframe
    cat("Prepare limma data frame\n")
    data_limma <- prep_limma_df(df$analysis[[analysis_name]][[data_name]]$data, col_name, col_value, exclude)
    cat("Success: Prepare limma data frame\n")

    # Design matrix
    cat("Prepare groups\n")
    groups <- df$analysis[[analysis_name]]$settings$metadata_samples %>%
        arrange(group) # Important sorting step to make sure that groups are in a consistent order
    groups <- factor(groups$group)
    cat("Success: Prepare groups\n")
    
    cat("Prepare design matrix\n")
    design <- model.matrix(~0+groups)
    colnames(design) <- levels(groups)
    cat("Success: Prepare design matrix\n")


    # Convert to matrix and log2-transform
    cat("Convert to matrix\n")
    data_limma <- as.data.frame(data_limma)
    rownames(data_limma) <- data_limma[,1]
    data_limma[,1] <- NULL
    m <- data.matrix(data_limma)
    m <- log2(m+const)
    cat("Success: Convert to matrix\n")


    # Contrast matrix
    # a-b means "a/b":
    #    positive FC: increased in a relative to b
    #    negative FC: decreased in a relative to b
    # comparison is defined above in settings
    cat("Make contrasts\n")
    cont_m <- makeContrasts(contrasts=comparison, levels=design)
    cat("Success: Make contrasts\n")

    # FC calculation
    cat("Fit\n")
    fit <- lmFit(m, design)
    fit <- contrasts.fit(fit, cont_m)
    fit <- eBayes(fit, trend=TRUE)
    cat("Success: Fit\n")
    
    
    # ANOVA
    cat("ANOVA\n")
    ANOVA <- topTable(fit, n=Inf)

    # Select only p-value column and move row names to column for later merging.
    ANOVA <- ANOVA %>%
        select(P.Value, adj.P.Val) %>%
        rename(ANOVA_p_val = P.Value,
               ANOVA_adj_p_val = adj.P.Val
              )

    ANOVA <- tibble::rownames_to_column(ANOVA, 'unique_name')
    cat("Success: ANOVA\n")
    
    
    cat("Aggregate pairwise comparisons\n")
    pairwise <- aggregate_pairwise_results(fit, comparison)
    cat("Success: Aggregate pairwise comparisons\n")
    
    
    cat("Add transformed data to original data frame\n")
    
    transformed <- as.data.frame(m) %>%
        tibble::rownames_to_column("unique_name") %>%
        pivot_longer(!unique_name, names_to='sample', values_to='limma_int')
    
    original <- df$analysis[[analysis_name]][[data_name]]$data
    
    added<- merge(original, transformed, by=c('unique_name','sample'), all.x=TRUE)
    
    cat("Success: Add transformed data to original data frame\n")
        

    # Add data to respective analysis
    cat("Write variables for output\n")
    df$analysis[[analysis_name]][[data_name]]$data <- added
    df$analysis[[analysis_name]][[data_name]]$limma$data_limma <- data_limma
    df$analysis[[analysis_name]][[data_name]]$limma$m <- m
    df$analysis[[analysis_name]][[data_name]]$limma$fit <- fit
    df$analysis[[analysis_name]][[data_name]]$limma$plotSA <- plotSA(fit)
    df$analysis[[analysis_name]][[data_name]]$limma$ANOVA <- ANOVA
    df$analysis[[analysis_name]][[data_name]]$limma$pairwise <- pairwise
    cat("Success: Write variables for output\n")
    
    message("Success: Limma analysis")
    
    return(df)
}

prep_limma_df <- function(df,
                          col_name,
                          col_value,
                          exclude # Vector specifying all features to exclude
                                       # Create manually, e.g. by:
                                       # exclude <- data$analysis$D8Phe_prot$intensity$data %>%
                                       #     filter(qc_missingness_metabolite == 'failed') %>%
                                       #     distinct(unique_name) %>%
                                       #     pull() # Creates a vector
                         ){
    
    df <- df%>%
        filter(!unique_name %in% exclude) %>% # Remove features in 'exclude'        
        arrange(group) %>% # Important sorting step to make sure that groups are in a consistent order
        select(sample, all_of(col_name), all_of(col_value)) %>%
        pivot_wider(names_from='sample', values_from=all_of(col_value))
    
    return(df)
    
}

aggregate_pairwise_results <- function(df_fit, 
                                    comparison
                                      ){

    res <- data.frame() # Initialize variable

    for(i in comparison){
        df <- topTable(df_fit, number = Inf, coef=i) %>%
            select(1,4,5) %>%
            mutate(comparison = i,
                   neglog10p = -log10(P.Value)
                  ) %>%
            relocate(comparison)
        df <- tibble::rownames_to_column(df, 'unique_name') %>%
            separate(unique_name, c('compound',NA,NA), sep='@', remove=FALSE)
        res <- rbind(res, df)
    }
    
    return(res)
}

csv_total <- function(df = data,
                      analysis_name,
                      data_name,
                      statistic, # e.g. limma
                      export = TRUE
                     ){
    
    message("Create csv for total (= peak intensity)")
    
    cat("Get data\n")
    total <- df$analysis[[analysis_name]][[data_name]]$data
    cat("Success: Get data\n")
    
    cat("Get statistics\n")
    test <- c('limma','anova_posthoc\n') %in% statistic # Returns a string that contains TRUE when one is present
    if (TRUE %in% test){
        cat("Statistics added\n")
        total <- merge(total, df$analysis[[analysis_name]][[data_name]][[statistic]]$ANOVA, by='unique_name', all.x=TRUE)
    } else {
        cat("No statistics added\n")
    }
    cat("Success: Get statistics\n")
    
    
    cat("Write variables for output\n")
    df$analysis[[analysis_name]][[data_name]]$csv$total <- total   
    cat("Success: Write variables for output\n") 
    
    
    if(export == TRUE){
        cat("Exporting csv\n")
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv', sep=''), showWarnings = FALSE)
        
        write.csv(df$analysis[[analysis_name]][[data_name]]$csv$total, paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv/', df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_total.csv', sep=''), row.names=FALSE)
        cat("Success: Exporting csv\n")
    }

    message("Success: Create csv for total (= peak intensity)")
    
    return(df)
}

csv_volcano <- function(df = data,
                        analysis_name,
                        data_name,
                        statistic, # e.g. limma
                        export = TRUE
                       ){    
    
    message("Create csv for vulcano")
    
    cat("Get data\n")
    pairwise <- df$analysis[[analysis_name]][[data_name]][[statistic]]$pairwise %>%
            separate(unique_name, c(NA,NA,'mode'), sep='@', remove=FALSE) %>%
            mutate(mode = as.factor(mode))
    cat("Success: Get data\n")
    
    cat("Add metadata\n")
    volcano <- merge(pairwise, df$metadata$library, by=c('compound', 'mode'), all.x=TRUE)
    cat("Success: Add metadata\n")
    
       
    cat("Write variables for output\n") 
    df$analysis[[analysis_name]][[data_name]]$csv$volcano <- volcano   
    cat("Success: Write variables for output\n") 
    
    
    
    if(export == TRUE){
        cat("Exporting csv\n")
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv', sep=''), showWarnings = FALSE)
        
        write.csv(df$analysis[[analysis_name]][[data_name]]$csv$volcano, paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv/', df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_volcano.csv', sep=''), row.names=FALSE)
        cat("Success: Exporting csv\n")
    }

    message("Success: Create csv for vulcano")
        
    return(df)
}

csv_IPA <- function(df = data,
                       analysis_name,
                       data_name, # KeggID needs to be unique -> run select_mode_for_compounds first if required
                       export = TRUE
                      ){
    
    message("Create csv for IPA")
    
    cat("Get data from $csv$volcano\n")
    selected <- df$analysis[[analysis_name]][[data_name]]$csv$volcano %>%
        filter(id != '') # If metabolites have no id this field is empty but not NA
    cat("Success: Get data from $csv$volcano\n")
    
    cat("Get log2fc\n")
    IPA_log2fc <- selected %>%
        select(comparison, id, logFC) %>%
        pivot_wider(names_from='comparison', values_from='logFC', names_prefix='log2fc_')
    cat("Success: Get log2fc\n")

    cat("Get p-value\n")
    IPA_pval <- selected %>%
        select(comparison, id, P.Value) %>%
        pivot_wider(names_from='comparison', values_from='P.Value', names_prefix='P.Value_')
    cat("Success: Get p-value\n")

    cat("Get adjusted p-value\n")
    IPA_adjpval <- selected %>%
        select(comparison, id, adj.P.Val) %>%
        pivot_wider(names_from='comparison', values_from='adj.P.Val', names_prefix='adj.P.Val_')
    cat("Success: Get adjusted p-value\n")

    
    cat("Merge data\n")
    IPA <- cbind(IPA_log2fc,IPA_pval[,-1],IPA_adjpval[,-1])
    cat("Success: Merge data\n")
    
    
    if(export == TRUE){
        cat("Exporting csv\n")
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv', sep=''), showWarnings = FALSE)

        write.csv(IPA, paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv/', df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_IPA.csv', sep=''), row.names=FALSE)
        cat("Success: Exporting csv\n")
    }
    
    cat("Write variables for output\n")
    df$analysis[[analysis_name]][[data_name]]$csv$IPA <- IPA
    cat("Success: Write variables for output\n")
    
    message("Success: Create csv for IPA")
    
    return(df)
    
}

csv_metaboanalyst <- function(df = data,
                              analysis_name,
                              data_name,
                              value,
                              export = TRUE
                             ){    
    
    message("Create csv for metaboanalyst")
    
    
    cat("Get data\n")
    metaboanalyst <- df$analysis[[analysis_name]][[data_name]]$data %>%
        arrange(group) %>% # Important sorting step 1 to make sure that the group row is in the correct order
        select(sample, unique_name, all_of(value)) %>%
        pivot_wider(names_from='sample', values_from=all_of(value))
    cat("Success: Get data\n")

    groups <- df$analysis[[analysis_name]]$settings$metadata_samples %>%
        arrange(group) # Important sorting step 2
    
    cat("Add group row\n")
    metaboanalyst <- rbind(c('', as.character(groups$group)), metaboanalyst)
    cat("Success: Add group row\n")
    
       
    cat("Write variables for output\n") 
    df$analysis[[analysis_name]][[data_name]]$csv$metaboanalyst <- metaboanalyst
    cat("Success: Write variables for output\n") 
    
    if(export == TRUE){
        cat("Exporting csv\n")
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, sep=''), showWarnings = FALSE)
        dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv', sep=''), showWarnings = FALSE)
        
        write.csv(df$analysis[[analysis_name]][[data_name]]$csv$metaboanalyst, paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/csv/', df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_metaboanalyst.csv', sep=''), row.names=FALSE)
        cat("Success: Exporting csv\n")
    }

    message("Success: Create csv for metaboanalyst")
    
    return(df)
}

metaboanalyst <- function(df = data,
                          analysis_name,
                          data_name,
                          norm
                         ){

    message("Metaboanalyst")
    
    library(MetaboAnalystR)
    
    cat("Change working directory\n")        
    dir.create(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/metaboanalyst', sep=''), showWarnings = FALSE)
    setwd(paste(df$metadata$working_directory, '/', analysis_name, '/', data_name, '/metaboanalyst', sep=''))
    cat("Success: Change working directory\n")

    # General setup
    cat("Initialization\n")
    mSet<-InitDataObjects("pktable", "stat", FALSE)
    mSet<-Read.TextData(mSet, paste('../csv/', df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_metaboanalyst.csv', sep=''), "colu", "disc");
    mSet<-SanityCheckData(mSet)
    mSet<-ReplaceMin(mSet);
    mSet<-SanityCheckData(mSet)
    mSet<-FilterVariable(mSet, "none", "F", 25)
    mSet<-PreparePrenormData(mSet)
    cat("Success: Initialization\n")
    
    # Normalization
    cat("Normalization\n")
    mSet<-Normalization(mSet, "NULL", norm, "NULL", ratio=FALSE, ratioNum=20)
    #mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
    #mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
    cat("Success: Normalization\n")

    # Heatmap
    # Only works from command line! Does not work in jupyter.
    cat("Heatmap\n")
    mSet<-PlotHeatMap(mSet, paste(df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_heatmap_', sep=''), "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NULL, T, F)
    cat("Success: Heatmap\n")
    
    # PCA
    cat("PCA\n")
    mSet<-PCA.Anal(mSet)
    #mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
    #mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
    mSet<-PlotPCA2DScore(mSet, paste(df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_PCA_', sep=''), "png", 72, width=NA, 1,2,0.95,0,0)
    #mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
    #mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
    #mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
    cat("Success: PCA\n")

    # PLSDA
    cat("PLSDA\n")
    mSet<-PLSR.Anal(mSet, reg=FALSE)
    #mSet<-PlotPLSPairSummary(mSet, "pls_pair_1_", "png", 72, width=NA, 5)
    mSet<-PlotPLS2DScore(mSet, paste(df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_PLSDA_', sep=''), "png", 72, width=NA, 1,2,0.95,1,0)
    #mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_1_", "png", 72, width=NA, 1,2,3, 40)
    #mSet<-PlotPLSLoading(mSet, "pls_loading_1_", "png", 72, width=NA, 1, 2);
    #mSet<-PLSDA.CV(mSet, "L",5, "Q2")
    #mSet<-PlotPLS.Classification(mSet, "pls_cv_1_", "png", 72, width=NA)
    #mSet<-PlotPLS.Imp(mSet, "pls_imp_1_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
    cat("Success: PLSDA\n")

    # sPLSDA
    cat("sPLSDA\n")
    mSet<-SPLSR.Anal(mSet, 2, 10, "same", "Mfold")
    #mSet<-PlotSPLSPairSummary(mSet, "spls_pair_0_", "png", 72, width=NA, 5)
    mSet<-PlotSPLS2DScore(mSet, paste(df$metadata$expID, '_', df$metadata$analysis_type, '_', analysis_name, '_', data_name, '_sPLSDA_', sep=''), "png", 72, width=NA, 1,2,0.95,0,0)
    #mSet<-PlotSPLS3DScoreImg(mSet, "spls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
    #mSet<-PlotSPLSLoading(mSet, "spls_loading_0_", "png", 72, width=NA, 1,"overview");
    #mSet<-PlotSPLSDA.Classification(mSet, "spls_cv_0_", "png", 72, width=NA)
    #mSet<-PlotSPLS3DLoading(mSet, "spls_loading3d_0_", "json", 1,2,3)
    cat("Success: sPLSDA\n")
    
    #O-PLSDA
    #mSet<-OPLSR.Anal(mSet, reg=FALSE)
    #mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
    #mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
    #mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE) # Does currently not work
    #mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)

    cat("Reset working directory\n")
    setwd(df$metadata$working_directory)
    cat("Success: Reset working directory\n")
    
    #df$analysis[[analysis_name]][[data_name]]$metaboanalyst$mSet = mSet
    
    message("Success: Metaboanalyst")
    
    #return(df)
    
}

combine_csvs <- function(df = data,
                         analysis_name, # vector of analysis names
                         data_name, # vector of data names
                         csv_source_name,
                         csv_target_name,
                         export = TRUE
                        ){
    
    message("Combine csvs")
    
    cat("Loop over csvs\n")
    df_combined <- NULL
    i = 1
    while (i <= length(analysis_name)){
        analysis <- df$analysis[[analysis_name[i]]][[data_name[i]]]$csv[[csv_source_name]] %>%
            mutate(analysis = paste(analysis_name[i], data_name[i], sep='_'))
        df_combined <- rbind(df_combined, analysis)
        
        i = i+1
    }
    cat("Success: Loop over csvs\n")

    cat("Write variables for output\n") 
    df$report[[csv_target_name]] <- df_combined
    cat("Success: Write variables for output\n")

    if(export == TRUE){
        cat("Exporting csv\n")
        setwd(df$metadata$working_directory)
        dir.create(paste('../', 'report', sep=''), showWarnings = FALSE)
        
        write.csv(df_combined, paste('../report/', df$metadata$expID, '_', df$metadata$analysis_type, '_', csv_target_name, '_summary.csv', sep=''), row.names=FALSE)
    
        setwd(df$metadata$working_directory)
        cat("Success: Exporting csv\n")
    }
    
    message("Success: Combine csvs\n")
    
    return(df)
    
}

pipeline2 <- function(df = data,
                     analysis_name,
                     data_name,
                     excluded_groups,
                     excluded_samples,
                     metadata_IS_neg,
                     metadata_IS_pos,
                     control,
                     comparison
                    ){
    
    # Normalization
    
    # No normalization because data has already been normalized
    

    
    df <- limma_analysis(df = df,
                         analysis_name = analysis_name,
                         data_name = data_name,
                         col_name = 'unique_name',
                         col_value = 'norm_int',
                         comparison = comparison
                        )
    
    
    # Export normalized intensity (metabolite levels in medium over time)
    df <- csv_total(df = df,
                    analysis_name = analysis_name,
                    data_name = data_name,
                    statistic = 'limma',
                    export = TRUE
                   )
    
    df <- csv_volcano(df = df,
                      analysis_name = analysis_name,
                      data_name = data_name,
                      statistic = 'limma',
                      export = TRUE
                     )
    
    df <- csv_metaboanalyst(df = df,
                            analysis_name = analysis_name,
                            data_name = data_name,
                            value = 'norm_int',
                            export = TRUE
                           )
       
    metaboanalyst(df = df,
                  analysis_name = analysis_name,
                  data_name = data_name,
                  norm='LogNorm')
    
    return(df)
    
}
