
# Help:
# We use MLE to estimate the parameters
# you just need to input the phenotype, moderartor, GRM and the identity matrix I
# Function "ras_MLE" will maximuze the loglike function
# See the "test_model_GxE_*" function for more info.
#
# model 1:
# log-like function: "loglike_model_GxE_a_ap_e_ep"
#
# model 2:
# log-like function: "loglike_model_GxE_a_e_ep"
#
# model 3:
# log-like function: "loglike_model_GxE_a_ap_e"
#
# model 4:
# log-like function: "likeloglike_GCTA_model"

# model 5:
# log-like function: "loglike_model_GxE_h2_costant"





################################
ras_MLE <- function(loglike_fun, param0)
{
    res <- optim(param0, loglike_fun, method="L-BFGS-B", control=list(REPORT=1,trace=6), hessian = T)
    return(res)
}




#################################################################
#################################################################
# CGTA estimation model
# model 4 in paper
# input: phen, I=diag(rep(1,nind)), GRM1=1/m*ZZp
loglike_GCTA_model <- function(phen,GRM1,I)
{
    ras_ret <- function(param)
    {
        s_a <- param[1]
        s_e <- param[2]
        
        Sigma <- s_a^2 * GRM1 + s_e^2 * I
        
        log_det <- determinant(Sigma)$modulus[1]
        if(is.na(log_det))
        {
            print("--NA")
            assign("not_conerged", T, envir = .GlobalEnv)
            0
        } else {
            if (log_det == -Inf)
            {
                print("--singular")
                diag(Sigma) <- diag(Sigma)+.01
                log_det <- determinant(Sigma)$modulus[1]
            }
            return(-(-.5 * log_det - .5 * t(phen) %*% solve(Sigma) %*% (phen)))
        }
    }
}


# GRM is n*n matrix
ras_init_GCTA <- function(phen,GRM1)
{
    phen_mean_extrcat <- phen - mean(phen)
    Phi  <- phen_mean_extrcat %*% t(phen_mean_extrcat)
    phi <- Phi[upper.tri(Phi, diag = FALSE)]
    
    grm1 <- GRM1[upper.tri(GRM1, diag = FALSE)]
    fit_a2 <- lm(phi~grm1)
    est_init_a  <- 0
    if(coef(fit_a2)[2]>0)
    {
        est_init_a=sqrt(coef(fit_a2)[2])
    }
    
    
    #for e
    y <- diag(Phi) - est_init_a^2 * diag(GRM1)
    est_init_e  <- 0
    if(mean(y)>0)
    {
        est_init_e <- sqrt(mean(y))
    }
    return(c(est_init_a,est_init_e))
}





#################################################################
#################################################################
# GxE model 1
# input: phen, moderator, I=diag(rep(1,nind)), GRM1=1/m*ZZp
loglike_model_GxE_a_ap_e_ep <- function(phen,moderator,GRM1,I)
{
    ras_ret <- function(param)
    {
        s_a    <- param[1]
        s_ap   <- param[2]
        s_e    <- param[3]
        s_ep   <- param[4]
        temp   <- s_a + s_ap * moderator
        T      <- temp %*% t(temp)
        temp_e <- s_e + s_ep * moderator
        K      <- temp_e %*% t(temp_e)
        
        Sigma  <- T * GRM1 + K * I
        log_det <- determinant(Sigma)$modulus[1]
        if(is.na(log_det))
        {
            print("--NA")
            assign("not_conerged", T, envir = .GlobalEnv)
            return(0)
        } else {
            if (log_det == -Inf)
            {
                print("--singular")
                diag(Sigma) <- diag(Sigma)+.01
                log_det <- determinant(Sigma)$modulus[1]
            }
            return(-(-.5 * log_det - .5 * t(phen) %*% solve(Sigma) %*% (phen)))
        }
    }

}

# To estimate
test_model_GxE_a_ap_e_ep <- function()
{
    # define phen,moderator,GRM1,I
    loglike_fun <- loglike_model_GxE_a_ap_e_ep(phen,moderator,GRM1,I)
    # find the initial value, see the paper
    param0      <- c(0,0,0,0)
    est         <- ras_MLE(loglike_fun, param0)
    print(est)
}


#################################################################
#################################################################
# GxE model 2
# input: phen, moderator, I=diag(rep(1,nind)), GRM1=1/m*ZZp
loglike_model_GxE_a_e_ep <- function(phen,moderator,GRM1,I)
{
    ras_ret <- function(param)
    {
        s_a    <- param[1]
        s_e    <- param[2]
        s_ep   <- param[3]
        temp_e <- s_e + s_ep * moderator
        K      <- temp_e %*% t(temp_e)
        
        Sigma  <- s_a^2 * GRM1 + K * I
        
        log_det <- determinant(Sigma)$modulus[1]
        if(is.na(log_det))
        {
            print("--NA")
            assign("not_conerged", T, envir = .GlobalEnv)
            return(0)
        } else {
            if (log_det == -Inf)
            {
                print("--singular")
                diag(Sigma) <- diag(Sigma)+.01
                log_det <- determinant(Sigma)$modulus[1]
            }
            return(-(-.5 * log_det - .5 * t(phen) %*% solve(Sigma) %*% (phen)))
        }
    }
    
}


# To estimate
test_model_GxE_a_e_ep <- function()
{
    # define phen,moderator,GRM1,I
    loglike_fun <- loglike_model_GxE_a_e_ep(phen,moderator,GRM1,I)
    # find the initial value, see the paper
    param0      <- c(0,0,0)
    est         <- ras_MLE(loglike_fun, param0)
    print(est)
}


#################################################################
#################################################################
# GxE model 3
# input: phen, moderator, I=diag(rep(1,nind)), GRM1=1/m*ZZp
loglike_model_GxE_a_ap_e <- function(phen,moderator,GRM1,I)
{
    ras_ret <- function(param)
    {
        s_a   <- param[1]
        s_ap  <- param[2]
        s_e   <- param[3]
        temp  <- s_a + s_ap * moderator
        T     <- temp %*% t(temp)
        
        Sigma <- T * GRM1 + s_e^2 * I
        
        log_det <- determinant(Sigma)$modulus[1]
        if(is.na(log_det))
        {
            print("--NA")
            assign("not_conerged", T, envir = .GlobalEnv)
            return(0)
        } else {
            if (log_det == -Inf)
            {
                print("--singular")
                diag(Sigma) <- diag(Sigma)+.01
                log_det <- determinant(Sigma)$modulus[1]
            }
            return(-(-.5 * log_det - .5 * t(phen) %*% solve(Sigma) %*% (phen)))
        }
    }
    
}


# To estimate
test_model_GxE_a_ap_e <- function()
{
    # define phen,moderator,GRM1,I
    loglike_fun <- loglike_model_GxE_a_ap_e(phen,moderator,GRM1,I)
    # find the initial value, see the paper
    param0      <- c(0,0,0)
    est         <- ras_MLE(loglike_fun, param0)
    print(est)
}


#################################################################
#################################################################
# GxE model 5
# input: phen, moderator, I=diag(rep(1,nind)), GRM1=1/m*ZZp
loglike_model_GxE_h2_costant <- function(phen,moderator,GRM1,I)
{
    ras_ret <- function(param)
    {
        s_a    <- param[1]
        s_ap   <- param[2]
        s_e    <- param[3]
        s_ep   <- s_e * s_ap / s_a
        
        temp   <- s_a + s_ap * moderator
        T      <- temp %*% t(temp)
        temp_e <- s_e + s_ep * moderator
        K      <- temp_e %*% t(temp_e)
        
        Sigma  <- T * ZZp + K * I
        
        log_det <- determinant(Sigma)$modulus[1]
        if(is.na(log_det))
        {
            print("--NA")
            assign("not_conerged", T, envir = .GlobalEnv)
            return(0)
        } else {
            if (log_det == -Inf)
            {
                print("--singular")
                diag(Sigma) <- diag(Sigma)+.01
                log_det <- determinant(Sigma)$modulus[1]
            }
            return(-(-.5 * log_det - .5 * t(phen) %*% solve(Sigma) %*% (phen)))
        }
    }
    
}


# To estimate
test_model_GxE_h2_costant <- function()
{
    # define phen,moderator,GRM1,I
    loglike_fun <- loglike_model_GxE_h2_costant(phen,moderator,GRM1,I)
    # find the initial value, see the paper
    param0      <- c(0,0,0)
    est         <- ras_MLE(loglike_fun, param0)
    print(est)
}








