

kr.1d <- function(x, y, h, eval.points, poly.deg=1)
{
  n <- length(x)
  nsumy <- sum(y)/n
  wy <- y/nsumy
  if (missing(eval.points))
    s0 <- kde(x=x, h=h)
  else
    s0 <- kde(x=x, h=h, eval.points=eval.points)

  fhat.wy <- kde(x=x, h=h, w=wy, eval.points=s0$eval.points)
  fhat.wy$estimate <- fhat.wy$estimate * nsumy
  mhat <- fhat.wy
  if (poly.deg==0)
  {  
    mhat$estimate <- mhat$estimate/s0$estimate
  }

  if (poly.deg==1)
  {
    fhat1.wy <- kdde(x=x, h=h, w=wy, eval.points=s0$eval.points, deriv.order=1)
    fhat1.wy$estimate <- fhat1.wy$estimate * nsumy*h
    s1 <- -kdde(x=x, h=h, deriv=1, eval.points=s0$eval.points)$estimate*h
    s2 <- kdde(x=x, h=h, deriv=2, eval.points=s0$eval.points)$estimate*h^2 - s0$estimate

    mhat$estimate <- s2/(s2*s0$estimate-s1^2) * fhat.wy$estimate + s1/(s2*s0$estimate-s1^2) * fhat1.wy$estimate
  }
  return(mhat)
}


kr.1d.extrap1 <- function(x, y, h, poly.deg)
{
  dx <- tail(diff(x), n=1)
  mhat <- kr.1d(x=x, y=y, h=h, eval.points=c(x, tail(x, n=1)+dx), poly.deg=poly.deg)

  return(list(x=mhat$eval.points, y=mhat$estimate))
}


kr.1d.extrap <- function(x, y, h, ne, poly.deg)
{
  mhat1 <- kr.1d.extrap1(x=x, y=y, h=h, poly.deg=poly.deg)
  x1 <- mhat1[[1]]
  y1 <- mhat1[[2]]

  browser()
  if (ne > 1)
    for (i in 2:ne)
    {
      mhat1 <- kr.1d.extrap1(x=x1, y=y1, h=h, poly.deg=poly.deg)
      x1 <- mhat1[[1]]
      y1 <- mhat1[[2]]
    }

  return(list(x=x1,y=y1))
}
