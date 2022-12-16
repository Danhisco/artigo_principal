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
f_simU <- function(df){
  v <- dinamica_coalescente(U = 1.25e-06, 
                            S=df$S_obs, 
                            disp_range = df$d, 
                            N_simul=1,
                            landscape = df$txt.file)
  v$U_est
}
