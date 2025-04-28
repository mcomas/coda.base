library(coda.base)
dat <- matrix(c(404.0555,738.4495,236.9774,57.90058,2.616169,
                522.9290,573.0399,253.1222,87.79391,3.115526), ncol = 5, byrow = TRUE,
              dimnames=list(NULL,c("SLEEP","SB","LPA","MPA","VPA")))

sbp_basis(b1=SLEEP+SB ~ LPA+MPA, data = dat, fill = TRUE)


partition = matrix(c(1,1,-1,-1,0,0,0,
                     0,0,1,-1,0,0,0), ncol = 2)

sbp_basis(partition)

