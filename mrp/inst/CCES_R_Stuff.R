library(mrpdata) 
data(CCES.complete)

data(mrp.census)
data(mrp.regions)
data(spmap.states)


<<prepare-data>>=
CCES.complete <- within (CCES.complete, {
  education <- factor(education,exclude=NA)
  female <- factor(sex=="Female",exclude=NA)
  race <- factor(race,exclude=NA)
  f.race <- interaction(female,race)
})
@

<<prepare-census>>=
  mrp.census <- na.omit(mrp.census[mrp.census$state %in% CCES.complete$state,])
  mrp.census <- na.omit(mrp.census[mrp.census$education %in% CCES.complete$education,])

  mrp.census <- within(mrp.census,{
    age <- factor(age,exclude=NA,labels=c("18-29","30-44","45-64","65+"))

    education[education=="postgraduate"] <- "college graduate"
    education <- factor(education,exclude=NA)
    edu <- factor(education,exclude=NA,labels=c("< High School",
                                         "High School",
                                         "Some College",
                                         "Graduated College"))
    state <- factor(state,exclude=NA)
    #race[race=="Other"] <- NA
    race <- factor(race,exclude=NA)
    f.race <- interaction(sex,race)
  })
  mrp.census <- na.omit(mrp.census)
@ 


<<simple-model,echo=TRUE,results=hide>>=
mrp.simple <- mrp(ban.gaymarr ~ state+age+education, 
                  data=CCES.complete,
                  population=mrp.census,
                  use="weighted2004",cov.prior="none")
@ 

<<simple-model-table-code,eval=FALSE>>=
xtable(100*poststratify(mrp.simple, ~ education+age), digits=0)
@
  
print(spplot(mrp.simple, state~state))
  

<<eval=FALSE>>=
mr.formula= .~.+ (1|region) + (1|age.edu) + z.age + p.relig.full + p.kerry.full
@ 

  
<<fullmodelcall,eval=FALSE>>=
mrp.statelevel <- mrp(ban.gaymarr~
                      state+f.race+age+education,
                      data=CCES.complete,
                      population=mrp.census,use="weighted2008",
                      #population.formula=.~.-poll,
                      add=list(
                        Statelevel,
                        mrp.regions,
                        expression(age.edu <- interaction(age,education))
                        ),
                      mr.formula=.~.+(1|region)+ (1|age.edu)+
                       p.relig.full+p.kerry.full
                      )
@
  
spplot(mrp.statelevel, state ~ education+age,
                   subset=TRUE,
                   spmap.states, "STATE", exclude=c("AK","DC","HI"),
                   stroke=list(expression(hasmarriage2010==TRUE),
                     "CA"),
                   center=poststratify(mrp.statelevel), cuts=50,
                   sub=paste("National average:",
                     format(poststratify(mrp.statelevel),digits=2)),
                   add.settings=list(
                     regions=list(col=fBasics:::divPalette(51,"BrBG")),
                     superpose.line=list(col=c("black","#00000066"),lwd=c(.3,1.3))
                     ),
                   colorkey=list(
                     space="bottom",height=.5,width=.5,
                     labels=list(at=c(.04,.34,.64),
                       labels=c("-30%","|","+30%"), cex=.7)
                     )
                   )
