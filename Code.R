
library('e1071')

set.seed(22200149)                      
file <- read.csv('Milk_MIR_Traits_data_2023.csv')     #read in dataset
N <- nrow(file)                                      #Number of rows in our file
del_obs <- sample(1:N,1)                       #sample a random integer here 310
file <- file[-del_obs,]                          #Delete that random integer row 
del_obs

#colnames(file)   #to find corresponing column name
#beta_lactoglobulin_b     # required column name
# We will now replace any missing values with NA
file_n <- file[!is.na(file$beta_lactoglobulin_b),] #identify missing rows in beta_l_b column and then remove those rows from our object file
nrow(file_n)
sum(is.na(file$beta_lactoglobulin_b))

#We will also reindex our dataframe as some observations have been removed
#row.names(file_n) <-1 :nrow(file_n)


ncol(file_n) - 531           #Total columns - Spectral columns = non spectral
spectra <- file_n[,52:582]   #Extracting columns of containing spectrum  
spectra_names <- as.numeric( gsub("X", "", colnames(spectra) )) # extract wavelength number for plotting by replacing all X, now spectrum is in corresponding wavelength (integer)

easy_df <-  data.frame(file_n$beta_lactoglobulin_b,spectra)
colnames(easy_df) <- c('beta_lac_b',spectra_names)

summary(file_n$beta_lactoglobulin_b)

temp1 <- subset(easy_df,(easy_df$beta_lac_b<3))
temp2 <- subset(easy_df,((easy_df$beta_lac_b>=3)&(easy_df$beta_lac_b<6)))
temp3 <- subset(easy_df,(easy_df$beta_lac_b>=6))

nrow(temp1)+nrow(temp2)+nrow(temp3)      #verified our observations our correctly divided

par(mfrow = c(1,3))
matplot(t(temp1[,-1]),x = spectra_names,pch=1:25, xlim = c(940,3825),ylim = c(-0.2,0.7), main = "Beta in [0,3)")
matplot(t(temp2[,-1]),x = spectra_names,pch=1:25, xlim = c(940,3825),ylim = c(-0.2,0.7), main = "Beta in [3,6)")
matplot(t(temp3[,-1]),x = spectra_names,pch=1:25, xlim = c(940,3825),ylim = c(-0.2,0.7), main = "Beta in [6,9.7)")

Plotting the graph as such we can see that observations with beta_lac_b values falling in [0,3) have much variation than others, we can see lots of non overlapping regions for examples at x near 1300,2200 and so on, wheras for other beta_lac_b values we can see the sprectrum returns similar result and huge overlapping makes the plot look dense following the same pattern. So if we are looking for explanatory variables I would go for beta with values that explains much variation.

Let us standardise our data with mean 0 and standard deviation 1

mean_easy_df <- apply(easy_df,2,mean)     #find the mean along columns
sd_easy_df <- apply(easy_df,2,sd)         #find sd along columns
std_mean_easy_df <- sweep(easy_df,2,mean_easy_df,'-')    #subtract respective means from col
std_easy_df <- sweep(std_mean_easy_df,2,sd_easy_df,'/')     #divide respective cols by their sd

Now we have standardised our data, let us remove all occurences in which abs(beta_lac_b)>3.

std_easy_df <- subset(std_easy_df,abs(std_easy_df$beta_lac_b)<=3*sd(std_easy_df$beta_lac_b))
summary(std_easy_df$beta_lac_b)
#row.names(std_easy_df) <- 1:nrow(std_easy_df)
So we have 303 observations after the preprocessing. 

df.eucl <- dist(std_easy_df,method = "euclidean")
#c1.average <- hclust(df.eucl,method = 'average')
#plot(c1.average, xlab = "Average Linkage")
#c1.single <- hclust(df.eucl,method = 'single')
#plot(c1.single, xlab = "Single Linkage")


c1.complete <- hclust(df.eucl,method = 'complete')
plot(c1.complete, xlab = "complete Linkage")

hcl = cutree(c1.complete, k=5)
table(hcl)WGSS = rep(0,10)    #placeholder fr within sum of square
n = nrow(std_easy_df[,-1])
WGSS[1] = (n-1) * sum(apply(std_easy_df[,-1], 2, var))   # WSS = sum(variances of points to their centroids, if n =1 then the mean is the mean of the data)

for(k in 2:10){
WGSS[k] = sum(kmeans(std_easy_df[,-1], centers = k)$withinss)
}

plot(1:10, WGSS, type="b", xlab="k", ylab="Within group sum of squares", main = "WSS vs clusters")

hcl = matrix(NA,nrow = 303, ncol = 9)    #make placeholders for clustering result
k_m = matrix(NA,nrow = 303, ncol = 9)
rand <- rep(NA, 9)                # initialising vector for adjusted random index
library(e1071)

  for (i in 2:10){ 
  hcl[,i-1] <- cutree(c1.complete, k = i)   #calculate cluster allocation using hierachial clustering 
  k_m [,i-1]<- kmeans(std_easy_df[,-1], centers = i)$cluster   #calculate cluster allocation using k means
  rand[i-1] <- classAgreement(table(k_m[,i-1],hcl[,i-1]))$crand   #find randIndex using cross Tab agree for each k
  }

plot(2:10,rand, main = "Adjusted rand index for clusters", ylab = "Adj Rand index", xlab = "Clusters")

set.seed(22200149)
k = 5
k_5 = kmeans(std_easy_df[,-1], center=5)
plot(std_easy_df[,c(2,70,100,178,200,349)],col = k_5$cluster, main = "Clusters in the data using 5 random spectrum")

k_5_df <- as.data.frame(k_5$cluster)

k_5_df_imp <- subset(k_5_df,k_5_df$`k_5$cluster`==3)  #we subset the rows with cluster 3

des_rows <- row.names(k_5_df_imp)
des_file <- file[des_rows,1:5]
table(des_file$Breed)

#table(des_file$Date_of_sampling)   # no conclusion here
#table(des_file$DaysInMilk) #no conclusion
table(des_file$Parity)



fit = prcomp(std_easy_df[,-1])    #303 x 531 dimension
#names(fit)
s <- (fit$sdev[1:10])^2/sum((fit$sdev)^2)     #finding variation from eigenvalues
sum = rep(0,10)
sum[1] = s[1]

for (i in 2:10) {
  sum[i] = s[i] + sum[i-1]
  
}
plot(1:10,sum, type = 'b', col = 'blue', ylab = "Cummulative prop sum", xlab = "Principal Components", main="Total Variation explained by Principal Components")
From here we can see that approximately 99% of variation in data is explained by taking first 5 principal components, these values can be referred from the cumulatve proportion table of fit or even the plot shows that after 4 Principal Components almost 98% variation is explained by the components, so 4-5 principal components can be considered. Even though there is not much difference in jumping from 4 to 5 I would like to take 4 principal components as there is not much change for going from 4 PC to 5 PC.

p <- predict(fit,std_easy_df) #Let us store the results to compare with our numerically derived solution

data_temp <- scale(std_easy_df[,-1],center = TRUE, scale = FALSE) #We need to scale our data as fit$scale = FALSE

# dim(data_temp)      check if the dimensions match
# dim(fit$rotation)

pq <- data_temp%*%fit$rotation   #using lecture notes Principle components * Variables == Scores

all(p==pq)
pairs(pq[,1:4], main = "Pairs plot of first 4 principal Components")
#install.packages('pls')
set.seed(22200149)
library('pls')
row_names_std <- row.names(std_easy_df)  # store in case
N = nrow(std_easy_df)  #total rows in our dataset
train_i =  2/3*(N)   # take 2/3 of orginal data fr training
test_i = N -train_i
train <- sample(1:N,train_i)  #randomly sampling 202 points

test <- setdiff(1:N,train)
#any(train==test)   returns False hence no duplicates

train_data <- std_easy_df[train,]  #subset training data
test_data <- std_easy_df[test,]
PCR_res <- pcr(beta_lac_b~., ncomp = 10, data = train_data, validation = "LOO") #regressing beta_lac values from our 3 principal components
summary(PCR_res)
plot(RMSEP(PCR_res), legendpos = "topright")

Our plot says that 4 components are sufficient to give us a good idea of our data, the lowest RMSEP occurs at n =4.

Let us look at our data closely with 4 components
plot(PCR_res, ncomp = 4, asp = 1, line = TRUE)

plot(PCR_res, plottype = "scores", comps = 1:4)
y_pred <-predict(PCR_res, ncomp = 4, newdata = test_data)
RMSEP(PCR_res, newdata = test_data)

milk_proteins <- file_n[,7:13] #Our data set with 7 milk proteins
row_index_milk <- row.names(milk_proteins)
#row.names(milk_proteins) <- 1:nrow(file_n)
beta0 <- milk_proteins[milk_proteins$beta_lactoglobulin_b==0,] # we can see there are 32 rows in which beta values after the missings values are removed.
names0 <- row.names(beta0)   #store the row names of 0 values in df0
milk_proteins[names0, "beta_lactoglobulin_b"] <- mean(milk_proteins$beta_lactoglobulin_b)

milk_fit <- prcomp(milk_proteins)   #Principle Components of milk prot
summary(milk_fit)
mu = colMeans(milk_proteins)
ncomp = 5  # I am taking 5 components into account for better accuracy

milk_hat <- milk_fit$x[,1:ncomp]%*%t(milk_fit$rotation[,1:ncomp])
milk_hat = scale(milk_hat,center = -mu, scale = FALSE)
head(milk_hat)
head(milk_proteins)
milk_proteins2 <- milk_proteins 
mu <- colMeans(milk_proteins)
mssold <-   mean(apply(((milk_hat - milk_proteins2)[,1:6])^2,2, mean))+ mean(((milk_hat - milk_proteins2)[-as.numeric(names0),7])^2)
thresh <- 1e-7
rel_err <- 1
iter <- 0
while (rel_err > thresh) {
# milk_proteins2[names0, "beta_lactoglobulin_b"] <- mean(milk_hat_df[,7])

  iter <- iter + 1
 # Step 2(a)
milk_hat <- prcomp(milk_proteins2)   #PCA
 # Step 2(b) 
milk_hat_df <- milk_fit$x[,1:5]%*%t(milk_fit$rotation[,1:5])  #Reconstructing original mmatrix

 # Step 2(c)
milk_hat_df = scale(milk_hat_df,center = -mu, scale = FALSE)
 

mss <- mean(apply(((milk_hat_df - milk_proteins2)[,1:6])^2,2, mean))  + mean(((milk_hat_df - milk_proteins2)[-as.numeric(names0),7])^2)  #Finding mean squared error seperately for non missing values 

 rel_err <- (mssold - mss) / mssold    #Relative error
 mssold <- mss

milk_proteins2[names0,7] <- milk_hat_df[names0,7] 
 
  cat (" Iter :", iter , " MSS :", mss ,
 " Rel . Err :", rel_err , "\n")
}

milk_proteins2[names0,7]
PCR_reg <- function(dataset){
  #dataset should have 1st column for target and rest for predictors
 N= nrow(dataset) 
 row_names <- row.names(dataset)
 train_i =  round(2/3*(N))   # take 2/3 of orginal data fr training
 test_i = N -train_i
 train <- sample(row_names,train_i)  #randomly sampling train_i points
 test <- setdiff(row_names,train)  #test data indexes
 
 dataset <- scale(dataset,center = TRUE, scale = TRUE) #scaled data
 dataset <- data.frame(dataset)
 train_data <- dataset[train,]  #subset training data
 test_data <- dataset[test,]
 
 PCR_result <- pcr(beta_lactoglobulin_b ~., ncomp = 10, data = train_data, validation = "LOO") #regressing beta_lac values from our 10 principal components
 
 y_pred <-predict(PCR_result, ncomp = 10, newdata = test_data)
 
 predplot(PCR_result, ncomp = 1:5)
 
 
 df_ <- data.frame(y_pred,test_data$beta_lactoglobulin_b)
 df_$MSE <- (df_[,1]-df_[,2])^2
 mean_MSE <- mean(df_$MSE)
 return(df_)
 }

#file_n contains non missing values
data_0 <- subset(file_n, file_n$beta_lactoglobulin_b==0)
data_non0 <- subset(file_n, file_n$beta_lactoglobulin_b!=0)
data_non0 <- data_non0[,c(13,52:582)]
non0 <- PCR_reg(data_non0)
mean(non0$MSE)
data_0$beta_lactoglobulin_b <- mean(file_n[,13])
data_0 <- data_0[,c(13,52:582)]
data_set_imp_mean <- rbind(data_0,data_non0)

mean0 <- PCR_reg(data_set_imp_mean)
mean(mean0$MSE)
#any(milk_proteins2$beta_lactoglobulin_b==0) gives FALSE
rowdata_imp_pca <- row.names(milk_proteins2)

data_set_imp_pca = data.frame()
data_set_imp_pca[rowdata_imp_pca,1:6] <- file_n[rowdata_imp_pca,1:6]
data_set_imp_pca[rowdata_imp_pca,7:13]<- milk_proteins2[rowdata_imp_pca,]
data_set_imp_pca[rowdata_imp_pca,14:582] <- file_n[rowdata_imp_pca,14:582]
data_set_pca <- data_set_imp_pca[,c(13,53:582)]
mean_pca <- PCR_reg(data_set_pca)
mean(mean_pca$MSE)