#!/usr/bin/env Rscript

options(java.parameters = "-Xmx32G")
workDir = tempdir()

library(SELEX)
library(SelexGLM)
library(grid)

selex.config(workingDir=workDir, maxThreadNumber=800)

R0.file = commandArgs(TRUE)[1]
data.file = commandArgs(TRUE)[2]
round = strtoi(commandArgs(TRUE)[3])
consensus = commandArgs(TRUE)[4]
kLen = nchar(consensus)
varLen = strtoi(commandArgs(TRUE)[5])
leftFixed = commandArgs(TRUE)[6]
rightFixed = commandArgs(TRUE)[7]

selex.defineSample('R0', R0.file, 'R0', 0, varLen, '', '')
selex.defineSample('RN', data.file, 'RN', round, varLen, '', '')

r0 = selex.sample(seqName = 'R0', sampleName = 'R0', round = 0)
r0.split = selex.split(r0)
r0.train = r0.split$train
r0.test = r0.split$test
dataSample = selex.sample(seqName = 'RN', sampleName = 'RN', round = round)

# MARKOV MODEL BUILT
kmax = selex.kmax(sample = r0.test, threshold = 100)
mm = selex.mm(sample = r0.train, order = NA, crossValidationSample = r0.test, Kmax = kmax, mmMethod = "TRANSITION")
mmscores = selex.mmSummary(sample = r0.train)
ido = which(mmscores$R==max(mmscores$R))
mm.order = mmscores$Order[ido]

selex.mmSummary()

# INFOGAIN USED TO CALCULATE KLEN
libLen = as.numeric(as.character(selex.getAttributes(dataSample)$VariableRegionLength))
selex.infogain(sample = dataSample, k = c(kLen), markovModel = mm)

data.kmerTable = getKmerCountAffinities(sample=dataSample, k=kLen, minCount = 0, markovModel=mm)
data.probeCounts = getProbeCounts(dataSample, markovModel = mm)

# Inputs about library are data specific 
glm_model = model("model",
             varRegLen = libLen,
             leftFixedSeq =  leftFixed, 
             rightFixedSeq = rightFixed, 
             consensusSeq = consensus,
             affinityType = "AffinitySym",
             leftFixedSeqOverlap = 10,
             minAffinity = 0.01,
             missingValueSuppression = 0.5,
             minSeedValue = .01, 
             upFootprintExtend = 5,
             confidenceLevel = 0.1, 
             rounds = list(c(round)),
             # Original paper forces symmetry
            #  rcSymmetric = TRUE,
             rcSymmetric = FALSE,
             verbose = FALSE, 
          )

summary(glm_model)

addSeedPsam(glm_model) = seedTable2psam(glm_model, data.kmerTable)
print(getValues(getN(glm_model)))

#Use this definition of data for complete analysis
data = data.probeCounts
data = topModelMatch(data, glm_model)

# Uses aligned probes to build design matrix
data = addDesignMatrix(data, glm_model)
designMatrixSummary = getDesignMatrix(glm_model, data)
print("Round summary: ")
print (designMatrixSummary$Round)
print("View/strand orientation summary: ")
print (designMatrixSummary$Intercept)
print("Mono-nucleotide summary: ")
print (designMatrixSummary$N)

regressionFormula = updatedRegressionFormula(data, glm_model)
print("Regression Formula: ")
print (regressionFormula)

fit = glm(regressionFormula, 
          data=data, 
          family = poisson(link="log"))

model = addNewBetas(glm_model, data, fit)

data = data.probeCounts
data.nrow = nrow(data)

for (i in 2:20) {
  print (paste("i =",i))

  data = topModelMatch(data, glm_model)
  data = addDesignMatrix(data, glm_model)

  designMatrixSummary = getDesignMatrix(glm_model, data)
  print("\n")
  print("Round summary: ")
  print (designMatrixSummary$Round)
  print("\n")
  print("Mono-nucleotide summary: ")
  print (designMatrixSummary$N)
  print("\n")
  print("View/strand orientation summary: ")
  print (designMatrixSummary$Intercept)

  if (data.nrow == nrow(data)) {
    print ("Stability Reached")
    break
  } else {
    data.nrow = nrow(data)
  }
          
  regressionFormula = updatedRegressionFormula(data, glm_model)
  print("\n")
  print("Regression Formula: ")
  print (regressionFormula)

  fit = glm(regressionFormula, 
            data=data, 
            family = poisson(link="log"))

  glm_model = addNewBetas(glm_model,data,fit)
}

save(glm_model, file = paste0("SelexGLM_R", round, ".RData"))

pwm = getPSAM(glm_model)
write.table(t(pwm), file=paste0("SelexGLM_k", kLen, "_R", round, ".pwm"), row.names=FALSE, col.names=c('A', 'C', 'G', 'T'), sep='\t', quote=FALSE)