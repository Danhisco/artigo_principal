dinamica_coalescente <- function(U, S=0, N_simul, disp_range, landscape, seed = as.numeric(Sys.time()), disp_kernel = 2, decay = 2){
    # Runs coalescent simulations for a given heterogeneous landscape
    #
    # Parameters:
    # U: speciation rate
    # S: observed richness (integer) - used to fit the value of U, or set to
    #       0 (default) if that is not desired
    # N_simul: number of simulations
    # seed: seed of the RNG (an integer)
    # disp_range: width of the dispersal kernel
    # disp_kernel: an integer corresponding to the type of dispersal kernel. One of
    #               0: uniform
    #               1: normal
    #               2: Laplacian
    # landscape: either a filename containing the landscape data, or a
    #   bidimensional R array or matrix.
    #   TODO: describe the format of the input - trinary matrix)
    #
    # Returns:
    # r: an array of dimension N_simul x landscape dimensions, that is, each
    #   r[i.,] is a bidimensional array of the same shape as the landscape.
    #   Each site is labeled according to the identity of the species occupying
    #   that site.
    # U_est: estimated speciation rate. This is returned only if input parameter S > 0
    if (is.character(landscape)){
        l <- as.matrix(read.table(landscape))
        infile <- landscape
        land_dims <- dim(l)
    } else {
        land_dims <- dim(landscape)
        infile <- tempfile()
        # input file *must* be clean: no comments, headers or anything
        write.table(landscape, infile, col.names=F, row.names=F)
    }
    outfile <- tempfile()
    repeat {
        system(paste('./dinamica_coalescente', land_dims[1], land_dims[2], U, S, N_simul,
                 seed, disp_range, disp_kernel, infile, outfile))
        if (file.exists(outfile) || S == 0)
            break
        U <- U/decay
        print(paste("Decreasing value of U to", U))
        # set some lowest boundary here so simulations don't take forever
        if (U < 1e-9){
            print("Richness value too low, giving up...")
            return(NA)
        }
    }
    r <- as.matrix(read.table(outfile))
    # transpose each grid, as output is written along lines but R reads it along columns
    # TODO: I thought I got it right, but it was wrong... please DO re-check
    #r <- aperm(r, c(1,3,2))

    # recover estimated speciation rate
    if (S > 0){
        out_con <- file(outfile)
        U_line <- strsplit(readLines(out_con, 2)[2], ' ')[[1]]
        close(out_con)
        U_est <- as.double(U_line[length(U_line)])
        return(list(r = r, U_est = U_est))
    }
    return(r)
}
rec_distribuicao_espacial <- function(r, landscape){
    # receives a vector `r` with species identities and the original landscape
    # (either a filename or a matrix), and returns a 2d matrix with the 2's
    # substituted by the species identity
    if (is.character(landscape))
        l <- as.matrix(read.table(landscape))
    else
        l <- landscape
    if (length(l[l==2]) != length(r)){
        print("Identity species vector does not correspond to landscape file! Wrong number of individuals!")
        return()
    }
    # TODO: check if this is correct or if we need to transpose the landscape
    # matrix first
    l[l==2] <- r
    return(l)
}

f_simU <- function(df_1row,land_mat=FALSE){
  # receives a data frame with 
  ## S_obs (observed species richness), 
  ## d (sd laplace dispersal kernel), 
  ## txt.file (path to the habitat spatial configuration matrix)
  ## land_mat=FALSE: o txt da simulação
  # returns a estimated U rate
  if(is.matrix(land_mat)){
    v <- with(df_1row,
              dinamica_coalescente(U = 1.25e-06, 
                                   S=S_obs, 
                                   disp_range = d, 
                                   N_simul=1,
                                   landscape = land_mat))  
  }else{
    v <- with(df_1row,
              dinamica_coalescente(U = 1.25e-06, 
                                   S=S_obs, 
                                   disp_range = d, 
                                   N_simul=1,
                                   landscape = txt.path))  
  }
  v$U_est
}

f_writeUcsv <- function(df_bySite, land_matrix = FALSE,
                        n_replicas=10, Umin = 1.25e-06, path_csv){
  l_dfU <- list()
  for(i in 1:nrow(df_bySite)){
    v_U <- replicate(n=n_replicas,f_simU(df_1row = df_bySite[i,],
                                         land_mat = land_matrix))
    l_dfU[[i]] <- cbind(select(df_bySite[i,],SiteCode,land_type,k,d),
                        matrix(v_U,nrow=1))
  }
  df_write <- rbind.fill(l_dfU)
  write_csv(df_write,file = paste0(path_csv,df_bySite$SiteCode[1],".csv"))
}

f_simSAD <- function(df_1row,n=100,land_mat=FALSE){
    # receives a data frame with 
    ## S_obs = 0
    ## d (sd laplace dispersal kernel), 
    ## txt.file (path to the habitat spatial configuration matrix)
    ## land_mat=FALSE: o txt da simulação
    # returns a estimated U rate
    if(is.matrix(land_mat)){
      with(df_1row,dinamica_coalescente(U = Umed, 
                                        S=0, 
                                        disp_range = d, 
                                        N_simul=n,
                                        landscape = land_mat))  
    }else{
      with(df_1row,dinamica_coalescente(U = Umed, 
                                        S=0, 
                                        disp_range = d, 
                                        N_simul=n,
                                        landscape = txt.path))  
    }
}
f_writeSADcsv <- function(df_bySite, land_matrix = FALSE,
                          n_replicas=100, path_csv){
  f_adply <- \(X) f_simSAD(df_1row = X, land_mat = land_matrix,n = n_replicas)
  v_names <- names(df_bySite)[!names(df_bySite)%in%c("SiteCode","k")]
  m_SADreplicas <- adply(df_bySite,1,f_adply)
  write_csv(select(m_SADreplicas,-all_of(v_names)),
            file = paste0(path_csv,df_bySite$SiteCode[1],".csv"))
}

f_simMNEE <- function(df,
                      U_rep=10,
                      SAD_rep=100,
                      Umin = 1.25e-06,
                      general_path){
  # nessa função já é feito todo o processo: U -> SADs
  f_simUeSAD <- \(df_exti,
                  land_type,
                  m_landi){
    folder_path <- paste0(general_path,
                          "/csv_SoE/taxaU/MNEE/",
                          land_type,
                          "/")
    path_df_simSAD <- paste0(folder_path,
                             df_exti$SiteCode[1],
                             ".csv")
    if(!file.exists(path_df_simSAD)){
      # taxa U
      f_writeUcsv(df_bySite = df_exti,
                  land_matrix = m_landi,
                  path_csv = folder_path,
                  n_replicas = U_rep)
    }
    # síntese taxa U
    df_simSAD <- read_csv(path_df_simSAD) |> 
      pivot_longer(-c(SiteCode:d)) |> 
      summarise(Umed = mean(value),.by = "k") |> 
      inner_join(x=df_exti,by="k")
    # SAD
    folder_path <- paste0(general_path,
                          "/csv_SoE/SADs_neutras/MNEE/",
                          land_type,
                          "/")
    f_writeSADcsv(df_bySite = df_simSAD,
                  land_matrix = m_landi, 
                  path_csv = folder_path,
                  n_replicas = SAD_rep)
  }
  f_singleext <- \(mland,dfi){
    # frag
    
    # non frag
    
    # pristine
    
  }
  f_SoE <- \(m_maxext){
    v_maxext <- max(df$SoE)
    ldf <- split(df,df$SoE)
    v_namesSoE <- abs(sort(desc(as.numeric(names(ldf)))))
    # 1o ciclo com a paisagem máxima
    teste <- f_singleext(mland = m_maxext,
                         dfi=ldf[[as.character(v_maxext)]]) 
      
    
  }
  #
  SAD_path <- list.files(path="../csv/SADs_neutras/MNEE",
                         pattern = df$SiteCode[1],
                         recursive = TRUE)
    
  if(length(SAD_path)==0){
    return(NA)
  }
  # tipo da paisagem
  m_land <- read.table(df$txt.path[1]) |> as.matrix()
  if(df$land_type[1] == "ideal"){
    m_land[m_land==0] <- 1
  }else if(df$land_type[1] == "non_frag"){
    m_land <- f_nonFragLand(m_land)
  }

  
  }

f_simMNEE2 <- \(df, 
                SAD_rep=100, 
                repo_path="../csv/SADs_neutras/MNEE_taxaU_idealizado/", 
                U_path="../csv/taxaU/df_U.csv"){
  df_U <- read_csv(U_path) %>% 
    filter(land_type=="ideal",SiteCode==df$SiteCode[1]) %>% select(-land_type)
  df_simSAD <- inner_join(df,df_U) %>% select(SiteCode,k:p,Umed)
  l_mland <- list()
  l_mland[[1]] <- read.table(df$txt.path[1]) |> as.matrix()
  l_mland[[2]] <- f_nonFragLand(l_mland[[1]])
  names(l_mland) <- c("contemp","non_frag")
  a_ply(names(l_mland),1,\(vname){
    f_writeSADcsv(df_bySite = df_simSAD,
                  land_matrix = l_mland[[vname]], 
                  path_csv = paste0(repo_path,vname,"__"),
                  n_replicas = SAD_rep)  
  })
  
}

f_simMNEE3 <- \(df_bySiteLand, 
                SAD_rep=100, 
                repo_path="../csv/SADs_neutras/MNEE_taxaU_"){
  # preparação
  ## caminho do respositorio
  land_U <- df_bySiteLand$land_type[1]
  repo_path <- paste0(repo_path,land_U,"/")
  ## paisagens contrafactuais
  l_mland <- list()
  l_mland[[1]] <- read.table(df_bySiteLand$txt.path[1]) |> as.matrix()
  l_mland[[2]] <- f_nonFragLand(l_mland[[1]])
  l_mland[[3]] <- l_mland[[1]]
  l_mland[[3]][l_mland[[3]]==0] <- 1
  names(l_mland) <- c("contemp","non_frag","ideal")
  ## vetor de indicação
  vname <- names(l_mland)[names(l_mland)!=land_U]
  # rotinas
  lapply(vname,\(i){
    f_writeSADcsv(df_bySite = df_bySiteLand,
                  land_matrix = l_mland[[i]], 
                  path_csv = paste0(repo_path,i,"__"),
                  n_replicas = SAD_rep)
  })
}

