setwd(here::here())
source("code/startup.R")
# Draw a slice sample for given position and momentum variables

# Calculate the Hamiltonian value for position and momentum variables.
#
# @details This function currently assumes iid standard normal momentum
# variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)

# Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
  runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))

# Test whether a "U-turn" has occured in a branch of the binary tree
# created by \ref\code{.buildtree} function.
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- (theta.plus-theta.minus)
  as.numeric( theta.temp %*% r.minus >= 0 | theta.temp %*% r.plus >= 0)
}
# A recursive function that builds a leapfrog trajectory using a balanced
# binary tree.
#
# @references This is from the No-U-Turn sampler with dual averaging
# (algorithm 6) of Hoffman and Gelman (2014).
#
# @details The function repeatedly doubles (in a random direction) until
# either a U-turn occurs or the trajectory becomes unstable. This is the
# 'efficient' version that samples uniformly from the path without storing
# it. Thus the function returns a single proposed value and not the whole
# trajectory.v


# Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
  runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))

# Test whether a "U-turn" has occured in a branch of the binary tree
# created by \ref\code{.buildtree} function.
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- (theta.plus-theta.minus)
  as.numeric( theta.temp %*% r.minus >= 0 | theta.temp %*% r.plus >= 0)
}
.buildtree <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v
    eps <- v*eps
    r <- r+(eps/2)*gr(theta)
    theta <- theta+eps*r
    r <- r+(eps/2)*gr(theta)
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## Useful code for debugging. Returns entire path to global env.
    # if(!exists('theta.trajectory')){
    #     theta.trajectory <<- theta
    # } else {
    #     theta.trajectory <<- rbind(theta.trajectory, theta)
    # }
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$v <- c(info$v, v)
    info$n.calls <- info$n.calls + 5
    info$theta <- rbind(info$theta, theta)
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}


get_NUTS_trajectory <-
  function(obj, init0, r0, covar, eps=.01, max_treedepth=15) {
  chd <- t(chol(covar))               # lower triangular Cholesky decomp.
  chd.inv <- solve(chd)               # inverse
  theta.cur <- as.vector(chd.inv %*% init0)
  fn2 <- function(theta) -obj$fn(chd %*% theta)
  gr2 <- function(theta) -as.vector( t( obj$gr(as.vector(chd %*% theta)) %*% chd ))
  theta.minus <- theta.plus <- theta0 <- theta.cur
  r.cur <- r.plus <- r.minus <- r0
  ## Draw a slice variable u
  u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
  j <- 0; n <- 1; s <- 1
  info <- as.environment( list(n.calls = 0, v=0, theta=theta.cur) )
  while(s==1) {
    v <- sample(x=c(1,-1), size=1)
    if(v==1){
      ## move in right direction
      res <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                        j=j, eps=eps, theta0=theta0, r0=r0,
                        fn=fn2, gr=gr2, info=info)
      theta.plus <- res$theta.plus
      r.plus <- res$r.plus
    } else {
      ## move in left direction
      res <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                        j=j, eps=eps, theta0=theta0, r0=r0,
                        fn=fn2, gr=gr2, info=info)
      theta.minus <- res$theta.minus
      r.minus <- res$r.minus
    }
    ## test whether to accept this state
    if(is.na(res$s) | is.nan(res$s))  res$s <- 0
    # does this need to be here?
    s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
    if(s==1) {
      info.keep <- as.list(info) # this keeps the doubled trajectory out
      if(runif(n=1, min=0,max=1) <= res$n/n){
        theta.out <- res$theta.prime
        #message("updating theta at ", j, " ", theta.out)
        #theta.out[m,] <- res$theta.prime
      }
    }
    n <- n+res$n
    s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
    ## Stop trajectory if there are any problems, probably happens
    ## when jumping way too far into the tails and the model isn't
    ## defined
    if(is.na(s) | is.nan(s))  s <- 0
    j <- j+1
    ## Stop doubling if too many or it's diverged enough
    if(j>max_treedepth & s) {
       warning("j larger than max_doublings, skipping to next m")
      break
    }
  }
  return(list(traj=info.keep$theta, v=info.keep$v, theta=theta.out))
  }
# for plotting trajectories as lines need to reorder the time-reversed parts
reorder_traj <- function(traj){
  v <- traj$v[-1]
  # initial point ignored here
  traj1 <- traj$traj[-1,]
  # split the trajectory into forward and backward time pieces, as detected via 'v'
  forward <- traj1[v==1,]
  backward <- traj1[v==-1,]
  # reverse order to flip time forward
  backward <- backward[rev(seq_len(nrow(backward))),]
  traj$traj <- rbind(backward, forward)
  return(traj)
}
print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
  tmp <- par("usr")
  text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
  text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
  label <- paste0('(',label, ')')
  text(x=text.x, y=text.y, labels=label, ...)
}

## Make exploratory plots of trajectories
png('plots/NUTS_trajectories.png', width=5, height=6.5, units='in', res=300)
par(mfrow=c(3,2), oma=c(0,1.1,0,0), mar=c(1.5,1.5,.25,.25), mgp=c(1,.25,0), tck=-.01)
init0 <- c(-1.05, -1.15) # init in iid Z space
init0 <- c(-2.1112354,.1124123)
#init0 <- c(-1,-1)
nreps <- 3
init0 <- lapply(1:nreps, \(x) init0)
# init0 <- lapply(1:nreps, \(x) rnorm(2))
# z <- seq(-2,2, len=3)
# z <- expand.grid(z,z, KEEP.OUT.ATTRS = FALSE) |> as.matrix()
# init0 <- split(z, row(z))
# nreps <- length(init0)
ses <- c(1,1)
eps <- .01
td <- 12
seed <- 1124261
col2 <- rgb(1,0,0,alpha=1)
col1 <- rgb(0,0,0,alpha=1)
kk <- 1
for(mycor in c(.1, .5, .9)){
  corr <- matrix(c(1,mycor,mycor,1), 2)
  covar <- diag(ses) %*% corr %*% diag(ses)
  covar.inv <- solve(covar)
  chd <- t(chol(covar))               # lower triangular Cholesky decomp.
  chd.inv <- solve(chd)               # inverse
  # transform to correlated space
  init <- lapply(init0, \(x) as.numeric(chd %*% x))
  library(RTMB)
  func <- function(pars){ -RTMB::dmvnorm(x = pars$x, mu = c(0,0), Sigma = covar, log=TRUE)}
  obj0 <- RTMB::MakeADFun(func = func, parameters = list(x=c(0,0)), silent=TRUE)
   # test <- get_NUTS_trajectory(obj0, init0=init, r0=c(1,1), eps=eps, covar=diag(2),
   #                             max_treedepth = 20)
  # get a set of random trajectories
  traj1 <- traj2 <- list()
  for(ii in 1:nreps){
    ## without mass matrix (covar)
    set.seed(seed+ii)
    r.cur <- rnorm(2)
    set.seed(seed+ii)
    traj1[[ii]] <- get_NUTS_trajectory(obj0, init0=init[[ii]], r0=r.cur, eps=eps, covar=diag(2), max_treedepth = td)
    traj1[[ii]] <- reorder_traj(traj1[[ii]])
    set.seed(seed+ii)
    traj2[[ii]] <- get_NUTS_trajectory(obj0, init0=init[[ii]], r0=r.cur, eps=eps, covar=covar, max_treedepth = td)
    traj2[[ii]] <- reorder_traj(traj2[[ii]])
  }
  # plot them in untransformed and transformed space
  for(transformed in c(FALSE,TRUE)){
    if(!transformed){
      plot(ellipse::ellipse(covar), type='l', xlab=NA, ylab=NA)
      tmp <-  init
      a1 <- lapply(traj1, \(y) y$traj)
      a2 <- lapply(traj2, \(y) t(apply(y$traj, 1, function(x) chd %*% x)))
      p1 <- lapply(traj1, \(y) y$theta)
      p2 <- lapply(traj2, \(y) as.numeric(chd %*% y$theta))
      mtext(bquote(.(as.name('rho'))==.(mycor)), side=2, line=1.1)
    } else {
      plot(ellipse::ellipse(diag(2)), type='l', xlab=NA, ylab=NA)
      tmp <-  init0
      a1 <-  lapply(traj1, \(y) t(apply(y$traj, 1, function(x) chd.inv %*% x)))
      a2 <-  lapply(traj2, \(y) y$traj)
      p1 <- lapply(traj1, \(y) as.numeric(chd.inv %*% y$theta))
      p2 <- lapply(traj2, \(y) y$theta)
    }
    print.letter(label=letters[kk], xy=c(.075,.925), cex=1.25)
    kk <- kk +1
    if(kk==6) legend('bottomright', legend=c('Correlated', 'Decorrelated'), bty='n', lty=1, col=c(col1,col2))
    lapply(a1, \(a) {
      N <- nrow(a)
      segments(a[-N,1], a[-N,2], a[-1,1], a[-1,2], lwd=2, col=col1)
      #points(a, pch=16, cex=.15,  col=col1)
      #points(a[1,1], a[1,2], col='green', pch=16)
    })
    lapply(a2, \(a) {
      N <- nrow(a)
      segments(a[-N,1], a[-N,2], a[-1,1], a[-1,2], lwd=2, col=col2)
      #points(a, pch=16, cex=.15,  col=col2)
      # points(a[1,1], a[1,2], col='green', pch=16)
    })

    #lapply(p1, \(p) points(p[1], p[2], col=col1, pch=15, cex=1.5))
    #lapply(p2, \(p) points(p[1], p[2], col=col2, pch=15, cex=1.5))
    lapply(tmp, \(x) points(x[1], x[2], pch=15, col='blue', cex=1.5))
  }
}
dev.off()

sapply(traj1, \(x) nrow(x$traj)*eps) |> mean()
sapply(traj2, \(x) nrow(x$traj)*eps) |> mean()

