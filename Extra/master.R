# Lese inn data
data1 = read.csv(file.choose(), header=TRUE) # eller file="plassering",
data2 = read.table(file.choose(), header=TRUE, sep=",") # sep=type seperasjon

data3 = read.delim(file.choose(), header=TRUE) # innlese .txt-fil
data4 = read.table(file.choose(), header=TRUE, sep="\t") # sep ved tab

# Jobbe med data
rm(data1) # fjerne data1 fra workspace
dim(data1) # antall rader og kolonner
head(data1) # de første 6 radene
tail(data1) # de siste 6 radene
data1[c(1,2,3,4,5), ] # eller data1[1:5, ] gir rad nr. 1,2,3,4,5 og alle kolonner
data1[-(5:10), ] # gir alle rader utenom radene 5-10
names(data1) # gir navn i header
length(data1$variabel1) # lengde av vektor
data1$variabel1[11:14] # subsett
data1[11:14, ]
summary(data1)

# Jobbe med variabler
mean(data1$variabel1) # dersom ikke attach() brukes
class(data$variabel1) # gir type
levels(data1$variabel1) # gir kategoriene i variabel1 dersom kategorisk
x = c(0,1,1,1,0,0,1,0,1) # lagres som en numerisk vektor
x = as.factor(x) # lagres som en faktor -> gir frekvenser ved summary()
mean(data1$variabel1[data1$variabel2 == "noe"])
data = data1[data1$variabel1=="noe" & variabel2>13, ]

temp = Age>15
temp2 = as.numeric(temp)
temp3 = variabel1=="noe" & variabel2=="nehgee"
moreData = cbind(data1, temp3) # legge til kolonne

rm(list=ls()) # clear

# Working directory
getwd()
setwd(~/Desktop/Dropboxblabla)

# Boxplot brukes for å oppsumere fordelingen av en numerisk variabel
boxplot(data1$variabel1, main="Tittel", ylab="y-akse", ylim=c(0,16), ylas=1) # variabel1 må være numerisk
quantile(data1$variabel1, probs=c(0, 0.25, 0.5, 0.75, 1))
# Ønsker å sammenligne fordelingen av en numeriske variabl for ulike kategorigrupper 
boxplot(data1$variabel1 ~ data1$variabel2, blabla) 

AgeGroups = cut(Age, breaks=c(0,13,15,17,25), labels=c("<13", "14/15", "16/17", "18+")) # konverterer numerisk variabel til kategorisk variabel
boxplot(data1$variabel1 ~ data1$variabel2**AgeGroups, las=2, blabla) # las=2 roterer x-verdiene

# Histogrammer: brukes for å oppsumere fordelingen av en numerisk variabel
hist(data1$variabel1, freq=FALSE, ylim=c(0,0.2), breaks=7, las=1) # 8 binn
lines(density(data1$variabel1), col=2, lwd=3)

# Scatterplott: brukes for å se på det lineære forholdet mellom to numeriske variabler
plot(data1$variabel1, data1$variabel2, main="Tittel", xlab="x", ylab="y", las=1, xlim=c(0,25), cex=0.5, col=2) # cex er pnktstørrelse
abline(lm(data1$variabel1 ~ data1$variabel2)) # lin.reg-linje
cor(data1$variabel1, data1$variabel2) # Pearson's correlation
cov(data1$variabel1, data1$variabel2) # kovarians
var(data1$variabel1) # varians
sd(data1$variabel1)
mean(data1$variabel1, trim=0.10) # fjerne øverste og nederste 10% av observasjoner

# Modifisere plott
# cex: endre størrelse av plottpunkter
# cex.main: endre tittelstørrelse
# cex.axis: akseverdistørrelse 
# font.main: endre font på tittelen
# col: endre farge på punkter
# col.main: endre farge på tittel
# pch="w": endre punktform
plot(data1$variabel1, data1$variabel2, main="Tittel", 
     cex=0.5, cex.main=2, cex.lab=1.5, cex.axis=0.7, 
     font.main=3)

# lty: linjetype
# lwd: linjetykkelse
abline(lm(data1$variabel1~data1$variabel2), col=2, lty=2, lwd=6)

plot(data1$variabel1[data1$variabel4=="noe"], 
     data1$variabel2[data1$variabel4=="noe"],
     col=4, pch=16,
     xlab="x", ylab="y")

# points(): legge til flere punkter i eksisterende plott
# bty="n": ingen boks rundt legend
points(data1$variabel1[data1$variabel4=="noe2"], 
       data1$variabel2[data1$variabel4=="noe2"],
       col=2, pch=17)
lines(smooth.spline(data1$variabel1[data1$variabel4=="noe"], 
                    data1$variabel2[data1$variabel4=="noe"),
      col=4, lwd=3, lty=2)
lines(smooth.spline(data1$variabel1[data1$variabel4=="noe2"], 
                    data1$variabel2[data1$variabel4=="noe2"),
      col=2, lwd=3, lty=3)
legend(x=3.5, y=14, legend=c("Noe", "Noe2"), 
       col=c(4,2), pch=c(16,17), bty="n")
legend(x=3.5, y=14, legend=c("Noe", "Noe2"), 
       col=c(4,2), lty=c(2,3), bty="n") # lty: linjetype


# Flere plott på en skjerm
# mfrow: 1 rad og 2 kolonner
par(mfrow=c(1,2))
plot(data1$variabel1[data1$variabel4=="noe"], 
     data1$variabel2[data1$variabel4=="noe"],
     xlab="x", ylab="y",
     xlim=c(0,20), ylim=c(45,85))
plot(data1$variabel1[data1$variabel4=="noe2"], 
     data1$variabel2[data1$variabel4=="noe2"],
     xlab="x", ylab="y",
     xlim=c(0,20), ylim=c(45,85))

# For å få ett plott på skjermen igjen
# axes=FALSE fjerner akse-verdier
# side=1: ønsker å merke x-aksen
# side=4: høyre side
par(mfrow=c(1,1))
plot(data1$variabel1, data1$variabel2,
     xlab="x", ylab="y", axes=FALSE)
axis(side=1, at=c(7, 12.3, 15), labels=c("sev", "mean", "15"))
axis(side=2, at=c(65, 70, 87), labels=c("65", "70", "87"))
box()

# Legge til tekst i et plott
# x,y: indikerer tekstplassering
# adj=0: teksten starter i x,y
plot(data1$variabel1, data1$variabel2)
text(x=5, y=11, label="r=0.82", adj=0, col=4, cex=0.5, font=4)
abline(h=mean(data1$variabel1), col=2, lwd=2) # horisontal linje ved mu

# Binomialfordeling
# n = 20 forsøk, p=1/6 for suksess 
# X ~ bin(n=20, p=1/6), der X er binoimalt fordelt
# dbinom: brukes for å finne verdier for pdf av X, f(x)
dbinom(x=3, size=20, prob=1/6) # P(X=3)
dbinom(x=0:3, size=20, prob=1/6) # P(X=0) & ... & P(X=3)
sum(dbinom(x=0:3, size=20, prob=1/6)) # P(X<=3)
pbinom(q=3, size=20, prob=1/6, lower.tail=TRUE) # også P(X<=3)

# Poisson fordeling
# lambda = 7
# X ~ POISSON(lambda=7)
# dpois: brukes for å finne verdier for pdf av X, f(x)
dpois(x=4, lambda=7) # P(X=4)
dpois(x=0:4, lambda=7) # P(X=0) & ... & P(X=4)
sum(dpois(x=0:4, lambda=7)) # P(X<=4)
ppois(q=4, lambda=7, lower.tail=TRUE) # P(X<=4)
ppois(q=12, lambda=7, lower.tail=FALSE) # P(X>12)

# Normalfordeling
# mu = 75, sd = 5
# X ~ N(mu=75, var=5^2)
pnorm(q=70, mean=75, sd=5, lower.tail=TRUE) # P(X<=70)
pnorm(q=85, mean=75, sd=5, lower.tail=FALSE) # P(X>85)
pnorm(q=1, mean=0, sd=1, lower.tail=FALSE) # P(Z>=1)
qnorm(p=0.25, mean=75, sd=5, lower.tail=TRUE) # finn Q1

x= seq(from=55, to=95, by=0.25)
dens = dnorm(x, mean=75, sd=5)
plot(x, dens, type="l")
abline(x=75) # legg til en vertikal linje ved gjennomsnittet

rand = rnorm(n=40, mean=75, sd=5)
hist(rand)

# t-fordeling og t-verdier
# mu = 0, sd = 1, df = 25
# t ~ t_df=25, mu = 0, sd = 1
# t-stat=2.3, df = 25, vil finne én-sidet p-verdi -> P(t>2.3)
pt(q=2.3, df=25, lower.tail=FALSE) # P(t>2.3)
# to-sidet p-verdi: finne arealet over t=2.3 og under t=-2.3
pt(q=2.3, df=25, lower.tail=FALSE) + pt(q=-2.3, df=25, lower.tail=TRUE)
pt(q=2.3, df=25, lower.tail=FALSE)*2 # også 

# Finn t-verdi for 95% konfidens
# -> finn t-verdi med 2.5% i hver hale
# dette er t-verdien brukt som kritisk verdi for 95% konfidens
qt(p=0.025, df=25, lower.tail=TRUE) # dette er t-verdien

# En-samplet t-test og konstruere en-samplet konfidensintervall
# En-samplet t-test og konfidensintervall er parametriske metoder for å 
# undersøke en enkel numerisk variabel

# Ho: mu>=8, Ha: mu<8
# en-side 95% konfidensintervall for mu
t.test(data1$variabel1, mu=8, alternative="less", conf.level=0.95)
# to-sidet
TEST = t.test(data1$variabel1, mu=8, 
              alternative="two.sided", conf.level=0.95)
attr(TEST) # gir attributtene i variabelen
# TEST$conf.int

# 2-samplet t-test og konfidensintervall
# Dette er parametriske metoder for å undersøke forskjellen i mu for
# 2 populasjoner. 
# Disse er måter å eksaminere forhold mellom en numerisk utfall variabel (Y)
# og en kategorisk variabel (X, med 2 nivåer)

# Ho: mean lung cap of smokers = of non-smokers
# to-sidet t-test, 
# antar ikke like varianser, var(LungCap[Smoke=="yes"]) mot "no"
# paired: gruppene mellom smokers og non-smokers er uavhengig
t-test(LungCap~Smoke, mu=0, alt="two.sided", 
       conf=0.95, var.eq=FALSE, paired=FALSE)


# Korrelasjon og Kovarians
# Pearson's korrelasjonskoeff. er et parametrisk mål av det
# lineære forholdet mellom 2 numeriske variabler
# Kendall's rank korrelasjon er et ikke-parametrisk mål på forholdet
# basert på samsvar mellom x-y par
cor.test(Age, LungCap, method="perason") # person er default, "kendall"

# Konfidensintervall kan returneres og samtidig teste om cor=0
# Ha: cor > 0
cor.test(Age, LungCap, method="pearson", alt="greater", conf.level=0.99)
cov(Age, LungCap)
# produsere par-vise plot for kun numeriske variabler og ikke kategoriske
# alle rader, for kolonnene 1-3
pairs(LungCapData[,1:3])
# Korrelasjon-matrise. Kan også gjøres for cov
cor(LungCapData[,1:3], method="pearson")


# Lineær regresjon
# LungCap er utfallsvariabelen (Y), avhengig variabel
plot(Age, LungCap, main="Scatterplot")
cor(Age, LungCap)
mod = lm(LungCap ~ Age) # lineær modell
abline(mod) # legge til regresjonslinjen. Denne linjen er hat(Y)=aX + b, tilpasset y-verdi
# Residual standard error: variasjoner av observasjoner rundt linjen, sqrt(MSE)
summary(mod)
attr(mod)
confint(mod, level=0.99) # konfidensintervall for modell-koeffisientene

# Residual: e(error) = Y(observert verdi) - hat(Y)
# 1. antagelse: Y-verdiene (eller errorene, "e") er uavhengige
# 2. antagelse: Y-verdiene kan uttrykkes som en lineær funksjon av X
# 3. antagelse: Variasjonene av observasjonene rundt en regresjonslinje (residual SE) er konstant (homoscedasticity)
# 4. antagelse: For en gitt X-verdi, er Y-verdiene (eller "e") normalfordelt
# Antagelse 2-4 kan undersøkes ved å betrakte residualer eller errorer

# Diagnostisk plot, for å evaluere modellen
plot(mod)
# Første plot viser Residualer (eller error e) vs hat(Y)
# Om antagelse nr. 2 er oppfylt, må det ikke være et mønster i plottet
# Linjen burde være flat, horisontal
# Om variasjonen er konstant, burde punktene være en sky av punkter
# I andre plot (Q-Q) burde linjen være diagonal hvis error/residualer er normalfordelt


# Multippel lineær regresjon
mod = lm(LungCap ~ Age + Height)
summary(mod)
# Multiple R-squared: prosentvis variasjon i LungCap kan forklares fra modellen (Age og Height)
cor(Age, Height, method="pearson")
# Kolinearitet mellom Age og Height betyr at vi burde ikke tolke stigningstallet (f.eks. Age), som en effekt av Age på LungCap
# Høy korrelasjon betyr at disse to effektene er noe bundet sammen
confint(mod, conf.level=0.95)
# konfidensintervaller for modellkoeff.

# Konvertere en numerisk variabel til en kategorisk variabel
############################################################
# dette gjøres for kryss-tabulering for en variabel, 
# eller tilpasse en regresjonsmodell når linearitetantagelsen er ikke gyldig for en variabel

# konverterer Height (en numerisk variabel) til en kategorisk
# kategorier: A<=50, B=50-55, C=55-60, D=60-65, E=65-70, F=70+
# øvre og nedre grense er 100 og 0, hhv.
# intervallene er (a,b]: 60 faller i kategorien C=(55,60]
# right=FALSE gir [a,b)
# breaks=4 f.eks gir fire kategorier eller nivåer
catHeight = cut(Height, breaks=c(0,50,55,60,65,70,100),
                labels=c("A","B","C","D","E","F"))

# Dummy eller indikator variabler
#################################

# kategorisk variabel med k-nivåer eller kategorier, krever k-1 indikator variabler
# Smoke har 2 nivåer "no" og "yes", og krever 1 indikator variabel for å representere røykestatus
# catHeight har 6 nivåer, som vil si 5 indikator variabler; X_B, ..., X_F
# kategori A er referanse/baseline gruppe; alltid første kategori som oppstår alfabetisk eller numerisk
# hat(mu)_(LungCap gitt x) = hat(Y) = b0 + b_BX_B + ... + b_FX_F
# der b_B er endring i gjennomsnittet for kategori B relativt til kategori A
mod = lm(LungCap, catHeight)

# Smoke har to nivåer "yes" og "no" (som er refereanse)
# kan endre referansenivået ved Smoke = relevel(Smoke, ref="yes")
# dette gir nye koeffisienter og kalles for reparametrisering av en modell

# Polynomial regresjon
mod = lm(LungCap ~ Height + I(Height^2))
summary(mod)
lines(smooth.spline(Height, predict(mod), col="blue", lwd=3)





































